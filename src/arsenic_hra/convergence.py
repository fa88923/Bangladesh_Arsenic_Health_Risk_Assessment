"""Phase 8: convergence and Monte Carlo error analysis.

For every district x population the full simulation is repeated at each N in
``simulation.convergence_N`` with ``replicate_seeds_per_N`` independent streams.
Tracked statistics T: mean HI, P50 HI, P95 HI, P(HI>1), P(HI>2), P95 ELCR.

Per (district, population, metric, N), over the R replicates:

    T_bar_N  = mean_r T_{N,r}                 across-seed mean
    s_N      = sd_r T_{N,r}                   between-seed spread (MC error of one run)
    Delta(N) = |T_bar_N - T_bar_Nprev| / |T_bar_N| * 100%    successive relative change

For probabilities p = k/N the binomial MC standard error sqrt(p(1-p)/N) is reported
next to the observed between-seed SD as a cross-check.

Convergence criterion (Member 4 decision, reports/convergence_report.md D16):
a statistic is adequate at N when
  - continuous statistics: Delta(N) <= 5% and CV_N = s_N / |T_bar_N| <= 5%;
  - probabilities:        |T_bar_N - T_bar_Nprev| <= 0.01 and s_N <= 0.01 (1 percentage point);
and the smallest adequate N is the smallest N from which every larger N is also adequate.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

TRACKED = ("mean_HI", "P50_HI", "P95_HI", "P_HI_gt_1", "P_HI_gt_2", "P95_ELCR")
PROBABILITY_METRICS = ("P_HI_gt_1", "P_HI_gt_2")
REL_TOL_PCT = 5.0
ABS_TOL_PROB = 0.01


def run_statistics(frame: pd.DataFrame) -> dict:
    hi = frame["HI"].to_numpy()
    elcr = frame["ELCR"].to_numpy()
    return {"mean_HI": float(hi.mean()), "P50_HI": float(np.quantile(hi, 0.50)),
            "P95_HI": float(np.quantile(hi, 0.95)), "P_HI_gt_1": float(np.mean(hi > 1.0)),
            "P_HI_gt_2": float(np.mean(hi > 2.0)), "P95_ELCR": float(np.quantile(elcr, 0.95))}


def aggregate(runs: pd.DataFrame) -> pd.DataFrame:
    """Long table: one row per (district, population, metric, N) with across-seed statistics."""
    long = runs.melt(id_vars=["district", "population", "N", "replicate"], value_vars=list(TRACKED),
                     var_name="metric", value_name="value")
    g = long.groupby(["district", "population", "metric", "N"], sort=False)["value"]
    agg = g.agg(replicates="count", across_seed_mean="mean", between_seed_sd="std",
                seed_min="min", seed_max="max").reset_index()
    agg["between_seed_cv_pct"] = 100 * agg["between_seed_sd"] / agg["across_seed_mean"].abs()
    agg = agg.sort_values(["district", "population", "metric", "N"], kind="stable").reset_index(drop=True)
    key = ["district", "population", "metric"]
    prev = agg.groupby(key, sort=False)["across_seed_mean"].shift(1)
    agg["N_prev"] = agg.groupby(key, sort=False)["N"].shift(1)
    agg["abs_change"] = (agg["across_seed_mean"] - prev).abs()
    agg["successive_rel_change_pct"] = 100 * agg["abs_change"] / agg["across_seed_mean"].abs()
    is_prob = agg["metric"].isin(PROBABILITY_METRICS)
    p = agg["across_seed_mean"].clip(0, 1)
    agg["binomial_mc_se"] = np.where(is_prob, np.sqrt(p * (1 - p) / agg["N"]), np.nan)
    agg["adequate"] = adequate(agg)
    return agg


def adequate(agg: pd.DataFrame) -> pd.Series:
    is_prob = agg["metric"].isin(PROBABILITY_METRICS)
    has_prev = agg["N_prev"].notna()
    # With T_bar = 0 (for example P(HI>2) in a very low-risk cell) relative measures are undefined;
    # the absolute rule is used whenever T_bar is exactly zero.
    zero = agg["across_seed_mean"] == 0
    cont_ok = (agg["successive_rel_change_pct"] <= REL_TOL_PCT) & (agg["between_seed_cv_pct"] <= REL_TOL_PCT)
    prob_ok = (agg["abs_change"] <= ABS_TOL_PROB) & (agg["between_seed_sd"] <= ABS_TOL_PROB)
    return has_prev & np.where(is_prob | zero, prob_ok, cont_ok)


def conclusions(agg: pd.DataFrame, primary_N: int) -> pd.DataFrame:
    """Smallest adequate N per (district, population, metric), and whether the primary N passes."""
    rows = []
    for (d, pop, metric), g in agg.groupby(["district", "population", "metric"], sort=False):
        g = g.sort_values("N")
        ok = g["adequate"].to_numpy()
        Ns = g["N"].to_numpy()
        smallest = None
        for i in range(len(Ns)):
            if ok[i:].all():
                smallest = int(Ns[i])
                break
        at_primary = g[g["N"] == primary_N].iloc[0]
        rows.append({"district": d, "population": pop, "metric": metric,
                     "smallest_adequate_N": smallest, "adequate_at_primary_N": bool(at_primary["adequate"]),
                     "primary_rel_change_pct": at_primary["successive_rel_change_pct"],
                     "primary_between_seed_cv_pct": at_primary["between_seed_cv_pct"],
                     "primary_between_seed_sd": at_primary["between_seed_sd"],
                     "N20000_rel_change_pct": g["successive_rel_change_pct"].iloc[-1],
                     "N20000_between_seed_cv_pct": g["between_seed_cv_pct"].iloc[-1],
                     "stabilises_within_grid": smallest is not None})
    return pd.DataFrame(rows)
