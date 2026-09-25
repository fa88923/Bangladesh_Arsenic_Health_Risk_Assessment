"""Phase 7: district x population Monte Carlo engine (reused by Phases 8 and 9).

Seeds: one stream per run,

    SeedSequence(master_seed, spawn_key=(district_index, population_index, n_index, replicate)),

where the indices are positions in ``arsenic.selected_districts``, ``simulation.populations``
and ``simulation.convergence_N``. The primary run is (N = primary_N, replicate 0).

Sampling: explicit inverse CDF. From each run's Generator, six independent uniform
vectors are drawn in the fixed order ``simulation.input_order`` (C, IR, BW, EF, ET, SA),
u = Generator.uniform(tiny, 1, N), then X = F^-1(u). Inputs are independent (approved).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from arsenic_hra import risk_equations as re_
from arsenic_hra import risk_parameters as rp
from arsenic_hra import simulation_inputs as si

U_LOW = np.finfo(float).tiny  # keeps u > 0 so Lognormal IR/SA ppf(u) > 0
SUMMARY_PROBS = (0.05, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
SUMMARY_METRICS = ("HI", "ELCR", "ELCR_ED_AT", "HQ_ing", "HQ_dermal")


@dataclass(frozen=True)
class RunSpec:
    district: str
    population: str
    N: int
    replicate: int
    district_index: int
    population_index: int
    n_index: int

    @property
    def run_id(self) -> str:
        return f"{self.district.lower()}_{self.population}_N{self.N}_r{self.replicate}"

    @property
    def spawn_key(self) -> tuple[int, int, int, int]:
        return (self.district_index, self.population_index, self.n_index, self.replicate)


def run_spec(cfg: dict, district: str, population: str, N: int, replicate: int = 0) -> RunSpec:
    sim = cfg["simulation"]
    return RunSpec(district, population, int(N), int(replicate),
                   cfg["arsenic"]["selected_districts"].index(district),
                   sim["populations"].index(population),
                   sim["convergence_N"].index(int(N)))


def generator(cfg: dict, spec: RunSpec) -> np.random.Generator:
    seq = np.random.SeedSequence(cfg["simulation"]["master_seed"], spawn_key=spec.spawn_key)
    return np.random.default_rng(seq)


@dataclass(frozen=True)
class ModelInputs:
    """Everything needed to simulate: final C and BW records plus base-paper parameter sets."""
    arsenic: dict        # district -> record
    bodyweight: dict     # population -> record
    params: dict         # population -> PopulationParameters

    def distributions(self, district: str, population: str) -> dict:
        p = self.params[population]
        a, b = self.arsenic[district], self.bodyweight[population]
        return {
            "C_mg_per_L": si.frozen(a["distribution"], a["params"]),
            "BW_kg": si.frozen(b["distribution"], b["params"]),
            **{name: d.frozen() for name, d in p.stochastic.items()},
        }


def load_model_inputs(cfg: dict) -> ModelInputs:
    bw, _ = si.bodyweight_inputs(cfg)
    return ModelInputs(si.arsenic_inputs(cfg), bw, {pop: rp.population_parameters(pop, cfg)
                                                   for pop in cfg["simulation"]["populations"]})


def draw_uniforms(rng: np.random.Generator, order, N: int) -> dict:
    return {name: rng.uniform(U_LOW, 1.0, size=N) for name in order}


def simulate(cfg: dict, model: ModelInputs, spec: RunSpec) -> pd.DataFrame:
    """One run: N iterations of sampled inputs and every risk output."""
    order = cfg["simulation"]["input_order"]
    if tuple(order) != rp.STOCHASTIC_INPUTS:
        raise ValueError(f"simulation.input_order must be {rp.STOCHASTIC_INPUTS}")
    dists = model.distributions(spec.district, spec.population)
    u = draw_uniforms(generator(cfg, spec), order, spec.N)
    x = {name: np.asarray(dists[name].ppf(u[name]), dtype=float) for name in order}
    for name, v in x.items():
        if not np.all(np.isfinite(v)):
            raise FloatingPointError(f"{spec.run_id}: non-finite {name}")
    p = model.params[spec.population]
    out = re_.evaluate_risk(x, p)
    # Base-paper convention (AT = ED * 365) kept alongside the lifetime ELCR for comparison.
    elcr_ed_at = re_.elcr(out["ADD_ing_mg_per_kg_day"], p.fixed["CSF_per_mg_per_kg_day"])
    frame = pd.DataFrame({"iteration": np.arange(spec.N, dtype=np.int64), **x, **out, "ELCR_ED_AT": elcr_ed_at})
    return frame


def wilson_interval(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return float(max(0.0, centre - half)), float(min(1.0, centre + half))


def summarize(frame: pd.DataFrame, spec: RunSpec, hi_thresholds=(1.0, 2.0), elcr_threshold: float = 1e-4) -> list[dict]:
    """One row per metric; threshold probabilities come directly from the iterations."""
    rows = []
    n = len(frame)
    for metric in SUMMARY_METRICS:
        v = frame[metric].to_numpy()
        q = np.quantile(v, SUMMARY_PROBS)
        row = {"run_id": spec.run_id, "district": spec.district, "population": spec.population, "N": n,
               "replicate": spec.replicate, "metric": metric, "mean": float(v.mean()), "sd": float(v.std(ddof=1)),
               "min": float(v.min()), "max": float(v.max()),
               **{f"P{int(round(p * 100))}": float(x) for p, x in zip(SUMMARY_PROBS, q)}}
        thresholds = hi_thresholds if metric == "HI" else (elcr_threshold,) if metric.startswith("ELCR") else ()
        for t in thresholds:
            k = int(np.sum(v > t))
            p_hat = k / n
            lo, hi = wilson_interval(k, n)
            tag = f"gt_{t:g}"
            row[f"P_{tag}"] = p_hat
            row[f"P_{tag}_mc_se"] = float(np.sqrt(p_hat * (1 - p_hat) / n))
            row[f"P_{tag}_wilson_lo"], row[f"P_{tag}_wilson_hi"] = lo, hi
        rows.append(row)
    return rows


def sampling_check(frame: pd.DataFrame, model: ModelInputs, spec: RunSpec) -> list[dict]:
    """Compare each sampled input with its target distribution (KS distance and quantiles)."""
    dists = model.distributions(spec.district, spec.population)
    n = len(frame)
    crit = 1.358 / np.sqrt(n)  # 5% two-sided KS critical value, large-n
    rows = []
    for name, d in dists.items():
        v = frame[name].to_numpy()
        ks = stats.kstest(v, d.cdf).statistic
        rows.append({"run_id": spec.run_id, "input": name, "N": n, "ks_distance": float(ks),
                     "ks_5pct_critical": float(crit), "ks_within_critical": bool(ks <= crit),
                     "target_P50": float(d.ppf(0.5)), "sample_P50": float(np.median(v)),
                     "target_P95": float(d.ppf(0.95)), "sample_P95": float(np.quantile(v, 0.95)),
                     "min": float(v.min()), "all_finite": bool(np.all(np.isfinite(v))),
                     "all_physical": bool(np.all(v >= 0) if name == "C_mg_per_L" else np.all(v > 0))})
    return rows
