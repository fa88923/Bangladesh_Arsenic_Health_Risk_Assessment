"""Phase 9: sensitivity (Spearman), parameter uncertainty (2-D MC), and model scenarios.

Sensitivity uses the exact Phase 7 primary-run iterations (loaded from Parquet), with the
tie-aware ``scipy.stats.spearmanr``. Only the six sampled inputs enter; fixed ED, Kp, CF,
RfDs and CSF have zero variance and are excluded (design doc 16.2).

Parameter uncertainty and scenarios reuse the primary run's uniforms (common random numbers),
so each comparison with the primary result differs only in the fitted parameters or model.
"""
from __future__ import annotations

import copy
from dataclasses import replace

import numpy as np
import pandas as pd
from scipy import stats

from arsenic_hra import monte_carlo as mc
from arsenic_hra import risk_equations as re_
from arsenic_hra import risk_parameters as rp
from arsenic_hra import simulation_inputs as si

INPUTS = rp.STOCHASTIC_INPUTS
# Expected sign of d(output)/d(input) from the equations: every input is in the numerator
# except BW; ELCR (ingestion only) does not depend on ET or SA.
EXPECTED_SIGN = {
    "HI": {"C_mg_per_L": 1, "IR_L_per_day": 1, "BW_kg": -1, "EF_days_per_year": 1, "ET_h_per_day": 1, "SA_m2": 1},
    "ELCR": {"C_mg_per_L": 1, "IR_L_per_day": 1, "BW_kg": -1, "EF_days_per_year": 1, "ET_h_per_day": 0, "SA_m2": 0},
}
# A sign contradicting the equations is a defect only if the correlation is distinguishable from
# zero; the dermal inputs ET and SA carry <0.3% of HI, so their true rho is below sampling noise.
SIGN_ALPHA = 1e-3
STAT_NAMES = ("mean_HI", "P50_HI", "P95_HI", "P_HI_gt_1", "P_HI_gt_2", "P95_ELCR")


# ---------------------------------------------------------------------------
# Spearman sensitivity
# ---------------------------------------------------------------------------

def spearman_rows(frame: pd.DataFrame, district: str, population: str, outputs=("HI", "ELCR")) -> list[dict]:
    rows = []
    for out in outputs:
        y = frame[out].to_numpy()
        for name in INPUTS:
            rho, p = stats.spearmanr(frame[name].to_numpy(), y)
            expected = EXPECTED_SIGN[out][name]
            rows.append({"district": district, "population": population, "output": out, "input": name,
                         "spearman_rho": float(rho), "p_value": float(p), "abs_rho": abs(float(rho)),
                         "expected_sign": expected,
                         "sign_matches": bool(expected == 0 or np.sign(rho) == expected),
                         "distinguishable_from_zero": bool(p < SIGN_ALPHA),
                         "sign_consistent": bool(expected == 0 or np.sign(rho) == expected or p >= SIGN_ALPHA)})
    return rows


def rank_inputs(table: pd.DataFrame) -> pd.DataFrame:
    table = table.copy()
    table["rank"] = table.groupby(["district", "population", "output"])["abs_rho"].rank(ascending=False,
                                                                                        method="first")
    return table


# ---------------------------------------------------------------------------
# Shared statistics
# ---------------------------------------------------------------------------

def risk_statistics(out: dict) -> dict:
    hi, elcr = np.asarray(out["HI"]), np.asarray(out["ELCR"])
    return {"mean_HI": float(hi.mean()), "P50_HI": float(np.quantile(hi, 0.5)),
            "P95_HI": float(np.quantile(hi, 0.95)), "P_HI_gt_1": float(np.mean(hi > 1)),
            "P_HI_gt_2": float(np.mean(hi > 2)), "P95_ELCR": float(np.quantile(elcr, 0.95))}


def primary_uniforms(cfg: dict, district: str, population: str) -> tuple[mc.RunSpec, dict]:
    spec = mc.run_spec(cfg, district, population, int(cfg["simulation"]["primary_N"]), 0)
    return spec, mc.draw_uniforms(mc.generator(cfg, spec), cfg["simulation"]["input_order"], spec.N)


# ---------------------------------------------------------------------------
# Two-dimensional Monte Carlo (parameter uncertainty)
# ---------------------------------------------------------------------------

def two_dimensional(cfg: dict, model: mc.ModelInputs, district: str, population: str,
                    arsenic_sets: list[dict], bw_sets: pd.DataFrame) -> pd.DataFrame:
    """One row per outer replicate b: risk statistics with (C, BW) parameters from replicate b."""
    spec, u = primary_uniforms(cfg, district, population)
    p = model.params[population]
    fixed_inputs = {name: np.asarray(d.frozen().ppf(u[name]), dtype=float) for name, d in p.stochastic.items()}
    a_family = model.arsenic[district]["distribution"]
    b_rec = model.bodyweight[population]
    bw_rows = bw_sets[bw_sets.population == population].reset_index(drop=True)
    n_outer = min(len(arsenic_sets), len(bw_rows))
    rows = []
    for b in range(n_outer):
        a = arsenic_sets[b]
        bw_params = {k: float(bw_rows.loc[b, k]) for k in b_rec["params"]}
        x = {"C_mg_per_L": np.asarray(si.frozen(a_family, a["params"]).ppf(u["C_mg_per_L"]), dtype=float),
             "BW_kg": np.asarray(si.frozen(b_rec["distribution"], bw_params).ppf(u["BW_kg"]), dtype=float),
             **fixed_inputs}
        out = re_.evaluate_risk(x, p)
        rows.append({"district": district, "population": population, "outer": b,
                     "arsenic_replicate": a["replicate"], "bw_replicate": int(bw_rows.loc[b, "replicate"]),
                     **{f"C_{k}": v for k, v in a["params"].items()},
                     **{f"BW_{k}": v for k, v in bw_params.items()}, **risk_statistics(out)})
    return pd.DataFrame(rows)


def uncertainty_summary(two_d: pd.DataFrame, primary: pd.DataFrame) -> pd.DataFrame:
    """Percentile intervals across outer replicates next to the primary (variability-only) value."""
    rows = []
    prim = primary.set_index(["district", "population"])
    for (d, pop), g in two_d.groupby(["district", "population"], sort=False):
        for stat in STAT_NAMES:
            v = g[stat].to_numpy()
            point = float(prim.loc[(d, pop), stat])
            lo, med, hi = np.quantile(v, [0.025, 0.5, 0.975])
            rows.append({"district": d, "population": pop, "statistic": stat, "primary_value": point,
                         "outer_replicates": len(v), "param_unc_median": float(med),
                         "param_unc_q025": float(lo), "param_unc_q975": float(hi),
                         "interval_width_over_primary": float((hi - lo) / point) if point else np.nan})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Model scenarios (common random numbers with the primary run)
# ---------------------------------------------------------------------------

def scenario_model(model: mc.ModelInputs, arsenic: dict | None = None, bodyweight: dict | None = None,
                   params: dict | None = None) -> mc.ModelInputs:
    return replace(model, arsenic={**model.arsenic, **(arsenic or {})},
                   bodyweight={**model.bodyweight, **(bodyweight or {})},
                   params={**model.params, **(params or {})})


def run_scenario(cfg: dict, model: mc.ModelInputs, cells: list[tuple[str, str]]) -> list[dict]:
    rows = []
    for district, pop in cells:
        spec = mc.run_spec(cfg, district, pop, int(cfg["simulation"]["primary_N"]), 0)
        frame = mc.simulate(cfg, model, spec)
        rows.append({"district": district, "population": pop, **risk_statistics(frame)})
    return rows


def arsenic_alternative(district: str, family: str) -> dict:
    fits = pd.read_csv(si.ARSENIC_FITS)
    row = fits[(fits.DISTRICT == district) & (fits.distribution == family)].iloc[0]
    params = {n: float(row[f"ls_{n}"]) for n in si.af_param_names(family)}
    return {"district": district, "distribution": family, "params": params,
            "scipy": si.scipy_spec(family, params), "source": "Phase 4 fit table (scenario)"}


def child_bw_gamma() -> dict:
    import json
    prov = {r["population"]: r for r in json.loads(si.BW_SELECTED.read_text(encoding="utf-8"))}["child"]
    return {"population": "child", "distribution": prov["distribution"], "params": prov["params"]}


def child_bw_age_extended(lower: float = 1.6) -> tuple[dict, dict]:
    """0-71-month child BW: 60/72 weight on the MICS 0-59 fit, 12/72 on an extrapolated 60-71 group.

    The 60-71-month group is the MICS 48-59-month sample shifted by the observed one-year gain
    (weighted mean of 48-59 minus 36-47 months), refitted with the same truncated-Normal LS.
    """
    child = pd.read_csv(si.CHILD_BW_CLEAN)
    x, w, age = child["BW_kg"].to_numpy(), child["survey_weight"].to_numpy(), child["CAGE_months"].to_numpy()

    def wmean(mask):
        return float(np.sum(x[mask] * w[mask]) / np.sum(w[mask]))

    older, prev = (age >= 48) & (age <= 59), (age >= 36) & (age <= 47)
    gain = wmean(older) - wmean(prev)
    ext = si.fit_truncated_normal_bw(x[older] + gain, w[older], lower)
    base = si.fit_truncated_normal_bw(x, w, lower)
    params = {"w1": 60 / 72, "mu1": base["params"]["mu"], "sigma1": base["params"]["sigma"],
              "mu2": ext["params"]["mu"], "sigma2": ext["params"]["sigma"], "lower": lower}
    info = {"gain_kg_per_year": gain, "mean_48_59_kg": wmean(older), "mean_36_47_kg": wmean(prev),
            "extrapolated_60_71_mean_kg": float(si.frozen("TruncatedNormal", ext["params"]).mean()),
            "mixture_mean_kg": float(si.frozen("TruncatedNormalMixture", params).mean())}
    return {"population": "child", "distribution": "TruncatedNormalMixture", "params": params}, info


def log_space_params(cfg: dict) -> dict:
    alt = copy.deepcopy(cfg)
    alt["risk_model"]["plus_minus_interpretation"] = "log_space"
    return {pop: rp.population_parameters(pop, alt) for pop in cfg["simulation"]["populations"]}
