"""Phase 9: fitted-parameter uncertainty (two-dimensional Monte Carlo).

Outer loop: B parameter sets for the fitted inputs.
  - Arsenic C: the Phase 4 well-level bootstrap, regenerated bit-for-bit with Member 2's
    ``arsenic_bootstrap.run_bootstrap`` and seeds, refitting the *final* family of each district
    (for Barisal the Gamma override) by the same CDF least squares.
  - Body weight: a survey-design (Rao-Wu rescaling) bootstrap. Within each stratum with
    n_h PSUs, n_h - 1 PSUs are drawn with replacement; each drawn PSU's weights are multiplied
    by n_h / (n_h - 1) times its draw count. The final family is refitted by weighted CDF LS.
Inner loop: N iterations with **common random numbers** - every replicate reuses the uniforms
of the Phase 7 primary run, so the spread across replicates is due to the parameters alone.

The base-paper inputs IR, EF, ET, SA have no data from which to estimate parameter uncertainty
and are held at their specified distributions (limitation).
"""
from __future__ import annotations

import os
from multiprocessing import get_context

import numpy as np
import pandas as pd

from arsenic_hra import arsenic_bootstrap as ab
from arsenic_hra import arsenic_fitting as af
from arsenic_hra import bodyweight_fitting as bf
from arsenic_hra import simulation_inputs as si

BW_DESIGN = {"adult": ("psu", "stratum"), "child": ("PSU", "stratum")}
BW_CLEAN = {"adult": si.RESULTS_BODYWEIGHT / "adult_bw_clean.csv", "child": si.CHILD_BW_CLEAN}
# Spawn-key tag for BW bootstrap streams: (tag, population_index, replicate); length 3 so it can
# never equal a simulation key (length 4) or a Phase 4 arsenic bootstrap key (length 2).
BW_STREAM_TAG = 9001


# ---------------------------------------------------------------------------
# Arsenic: regenerate the Phase 4 bootstrap
# ---------------------------------------------------------------------------

def regenerate_arsenic_bootstrap(cfg: dict, jobs: int | None = None) -> pd.DataFrame:
    acfg = cfg["arsenic"]
    data = pd.read_csv(si.ARSENIC_PREPROCESSED)
    districts = acfg["selected_districts"]
    return ab.run_bootstrap(data[data["DISTRICT"].isin(districts)], replicates=int(acfg["bootstrap_replicates"]),
                            master_seed=int(acfg["bootstrap_seed"]),
                            tie_tolerance=float(acfg["bootstrap_tie_tolerance"]),
                            jobs=jobs or max(1, (os.cpu_count() or 2) - 2), district_order=districts)


def compare_with_phase4(boot: pd.DataFrame) -> pd.DataFrame:
    """Recompute Member 2's summary from the regenerated table and diff it against the saved one."""
    saved = pd.read_csv(si.RESULTS_ARSENIC / "arsenic_bootstrap_stability.csv")
    mine = ab.summarise_bootstrap(boot)
    cols = [c for c in saved.columns if c not in ("DISTRICT", "distribution") and saved[c].dtype.kind == "f"]
    merged = saved.merge(mine, on=["DISTRICT", "distribution"], suffixes=("_saved", "_regen"))
    rows = []
    for c in cols:
        a, b = merged[f"{c}_saved"].to_numpy(float), merged[f"{c}_regen"].to_numpy(float)
        both = np.isfinite(a) & np.isfinite(b)
        nan_mismatch = int(np.sum(np.isfinite(a) != np.isfinite(b)))
        diff = float(np.max(np.abs(a[both] - b[both]) / np.maximum(np.abs(a[both]), 1e-300))) if both.any() else 0.0
        rows.append({"column": c, "max_rel_diff": diff, "nan_pattern_mismatches": nan_mismatch})
    return pd.DataFrame(rows)


def arsenic_parameter_sets(boot: pd.DataFrame, final: dict[str, dict]) -> dict[str, list[dict]]:
    """District -> list of (replicate, params) for the final family; failed refits are dropped."""
    out = {}
    for district, rec in final.items():
        family = rec["distribution"]
        names = af.PARAM_NAMES[family]
        part = boot[(boot.DISTRICT == district) & (boot.status == "ok")]
        sets = []
        for _, r in part.iterrows():
            params = {n: float(r[f"ls_{n}"]) for n in names}
            if np.all(np.isfinite(list(params.values()))) and ab._admissible_in_replicate(family, params):
                sets.append({"replicate": int(r.replicate), "params": params})
        out[district] = sets
    return out


# ---------------------------------------------------------------------------
# Body weight: survey-design bootstrap
# ---------------------------------------------------------------------------

def rao_wu_weights(psu: np.ndarray, stratum: np.ndarray, weight: np.ndarray,
                   rng: np.random.Generator) -> np.ndarray:
    """One Rao-Wu (n_h - 1) with-replacement PSU bootstrap of the survey weights."""
    new = np.zeros_like(weight, dtype=float)
    for h in np.unique(stratum):
        in_h = stratum == h
        psus = np.unique(psu[in_h])
        n_h = psus.size
        if n_h < 2:
            new[in_h] = weight[in_h]
            continue
        drawn, counts = np.unique(rng.choice(psus, size=n_h - 1, replace=True), return_counts=True)
        mult = dict(zip(drawn.tolist(), (counts * n_h / (n_h - 1)).tolist()))
        idx = np.flatnonzero(in_h)
        new[idx] = weight[idx] * np.array([mult.get(p, 0.0) for p in psu[idx].tolist()])
    return new


def _bw_refit(job) -> dict:
    population, family, params0, replicate, master_seed, pop_index, lower = job
    df = pd.read_csv(BW_CLEAN[population])
    psu_col, stratum_col = BW_DESIGN[population]
    rng = np.random.default_rng(np.random.SeedSequence(master_seed, spawn_key=(BW_STREAM_TAG, pop_index, replicate)))
    w = rao_wu_weights(df[psu_col].to_numpy(), df[stratum_col].to_numpy(), df["survey_weight"].to_numpy(), rng)
    keep = w > 0
    x, w = df["BW_kg"].to_numpy()[keep], w[keep]
    row = {"population": population, "replicate": replicate, "family": family, "n_kept": int(keep.sum())}
    if family == "TruncatedNormal":
        fit = si.fit_truncated_normal_bw(x, w, lower)
        row.update(fit["params"], converged=fit["ls_converged"])
    else:
        fit = bf.fit_cdf_least_squares(x, w, family)
        row.update(fit.params, converged=bool(fit.success))
    return row


def bw_parameter_sets(cfg: dict, final_bw: dict[str, dict], B: int, jobs: int | None = None) -> pd.DataFrame:
    seed = int(cfg["simulation"]["master_seed"])
    pops = cfg["simulation"]["populations"]
    jobs_list = []
    for pop in pops:
        rec = final_bw[pop]
        lower = rec["params"].get("lower")
        jobs_list += [(pop, rec["distribution"], rec["params"], b, seed, pops.index(pop), lower) for b in range(B)]
    n = jobs or max(1, (os.cpu_count() or 2) - 2)
    with get_context("fork").Pool(n) as pool:
        rows = pool.map(_bw_refit, jobs_list, chunksize=4)
    return pd.DataFrame(rows).sort_values(["population", "replicate"]).reset_index(drop=True)
