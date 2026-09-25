"""Step D: bootstrap stability of the district arsenic fits (plan.md section 3.D).

Per district, ``bootstrap_replicates`` times: resample the wells **with
replacement**, preserving each row's exact/censored status and its censoring bound
(``< L`` rows keep carrying ``L``), rebuild the reverse-KM curve, refit all three
families both ways, and record the parameters, fit error and winner.

Seeds follow the project convention
``numpy.random.SeedSequence(master_seed, spawn_key=(district_index, replicate))``.

Two implementation notes:

- Resampling duplicates exact concentrations. Ties are harmless for the reverse-KM
  curve (lifelines handles them) but make SciPy's censored fitters slow, so
  duplicated values are separated by ~1e-9 of the local spacing. On the real
  districts this is a no-op: 27 of the 726 (district, replicate) draws contained no
  duplicate value at all, and the perturbation is 10 orders of magnitude below the
  smallest measured concentration. It touches the fitted curve, never the data.
- Replicates are independent, so the work is a deterministic map over
  ``(district, replicate)`` and is safe to spread over processes; results are
  re-sorted by index and are bit-identical to a serial run.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .arsenic_fitting import (BOOTSTRAP_LS_OPTIONS, FAMILIES, NORMAL_MAX_NONPOSITIVE_PROB, PARAM_NAMES,
                              defective_starts, fit_cdf_least_squares, fit_censored_mle, frozen,
                              objective_mask, reverse_km_cdf)

WINNER_COLUMNS = ("ls_rmse_winner", "ls_rmse_winner_admissible", "aicc_winner")

# Report quantiles for every recorded quantity.
Q_LOW, Q_MED, Q_HIGH = 2.5, 50.0, 97.5
PARAMETER_STATISTICS = (
    ("_median", 50.0),
    ("_q025", 2.5),
    ("_q975", 97.5),
    ("_rmse_spread", None),  # sqrt(mean((theta_b - theta_bar)^2))
)


@dataclass(frozen=True)
class BootstrapSpec:
    """Everything a single replicate needs; picklable for process pools."""

    district: str
    district_index: int
    values_mgL: np.ndarray
    censored: np.ndarray
    master_seed: int
    tie_tolerance: float = 1e-9


def break_ties(values_mgL, tolerance: float = 1e-9) -> np.ndarray:
    """Separate duplicated concentrations so the censored fitters stay well behaved.

    Returns ``values_mgL`` unchanged when every value is distinct (the real-data case).
    """
    x = np.asarray(values_mgL, dtype=float)
    if np.unique(x).size == x.size:
        return x
    order = np.argsort(x, kind="stable")
    sorted_x = x[order]
    unique_x, first_index = np.unique(sorted_x, return_index=True)
    spacing = np.diff(unique_x)
    step = (spacing.min() if spacing.size else unique_x[0]) * tolerance
    if not np.isfinite(step) or step <= 0:
        step = tolerance * max(abs(unique_x[0]), 1.0)
    # Rank of each element inside its own tie group.
    offsets = np.arange(sorted_x.size) - np.repeat(first_index, np.diff(np.append(first_index, sorted_x.size)))
    out = np.empty_like(sorted_x)
    out[:] = sorted_x + offsets * step
    result = np.empty_like(x)
    result[order] = out
    return result


def replicate_seed(master_seed: int, district_index: int, replicate: int) -> np.random.Generator:
    """``SeedSequence(master_seed, spawn_key=(district_index, replicate))`` -> Generator."""
    seq = np.random.SeedSequence(master_seed, spawn_key=(int(district_index), int(replicate)))
    return np.random.default_rng(seq)


def bootstrap_replicate(spec: BootstrapSpec, replicate: int) -> dict:
    """One resampled refit of every candidate family. Never raises: degenerate
    replicates are reported with ``status`` instead so the run cannot be killed by
    an unlucky draw."""
    rng = replicate_seed(spec.master_seed, spec.district_index, replicate)
    n = spec.values_mgL.size
    index = rng.integers(0, n, size=n)
    values = break_ties(spec.values_mgL[index], spec.tie_tolerance)
    censored = spec.censored[index]

    row: dict = {"DISTRICT": spec.district, "district_index": spec.district_index,
                 "replicate": int(replicate), "n": int(n),
                 "n_censored": int(censored.sum()), "status": "ok",
                 "skipped_reason": ""}

    if censored.all() or not censored.any():
        row["status"] = "skipped"
        row["skipped_reason"] = ("every resampled row is censored" if censored.all()
                                 else "no resampled row is censored")
        return _fill_empty_families(row)

    try:
        curve = reverse_km_cdf(values, censored)
    except ValueError as exc:
        row["status"] = "skipped"
        row["skipped_reason"] = f"reverse-KM failed: {exc}"
        return _fill_empty_families(row)

    row["rkm_points"] = int(curve.x.size)
    row["rkm_censoring_fraction"] = curve.censoring_fraction
    try:
        row["objective_points"] = int(objective_mask(curve).sum())
    except ValueError:
        row["status"] = "skipped"
        row["skipped_reason"] = "objective region leaves too few points"
        return _fill_empty_families(row)

    per_family: dict[str, dict] = {}
    for family in FAMILIES:
        per_family[family] = _fit_replicate_family(values, censored, curve, family)
    return _assemble_row(row, per_family)


def _fit_replicate_family(values: np.ndarray, censored: np.ndarray, curve, family: str) -> dict:
    """LS + censored-MLE refit for one family on one replicate."""
    extra = defective_starts(family, values, censored)
    try:
        ls = fit_cdf_least_squares(values, censored, family, curve=curve,
                                   options=BOOTSTRAP_LS_OPTIONS, extra_starts=extra)
    except (ValueError, FloatingPointError) as exc:
        return {"ls_ok": False, "ls_note": str(exc), "mle_ok": False, "mle_note": str(exc)}
    mle = fit_censored_mle(values, censored, family)
    return {"ls_ok": ls.success and np.all(np.isfinite(list(ls.params.values()))),
            "ls_params": ls.params, "ls_sse": ls.sse, "ls_points": ls.n_points,
            "ls_rmse": float(np.sqrt(ls.sse / ls.n_points)) if ls.n_points else np.nan,
            "ls_note": ls.message,
            "mle_ok": mle.success and np.all(np.isfinite(list(mle.params.values()))),
            "mle_params": mle.params, "aicc": mle.aicc, "mle_note": mle.message}


def _fill_empty_families(row: dict) -> dict:
    for family in FAMILIES:
        row[f"ls_rmse_{family.lower()}"] = np.nan
        row[f"aicc_{family.lower()}"] = np.nan
        for name in PARAM_NAMES[family]:
            row[f"ls_{name}"] = np.nan
            row[f"mle_{name}"] = np.nan
    for column in WINNER_COLUMNS:
        row[column] = None
    return row


def _admissible_in_replicate(family: str, params: dict) -> bool:
    """The same physical gate as the main fit, applied to a replicate."""
    if not np.all(np.isfinite([params[k] for k in PARAM_NAMES[family]])):
        return False
    if family in ("Lognormal", "Gamma"):
        lo, _ = frozen(family, params).support()
        if lo != 0.0:
            return False
    if family == "Normal" and float(frozen(family, params).cdf(0.0)) >= NORMAL_MAX_NONPOSITIVE_PROB:
        return False
    return True


def _assemble_row(row: dict, per_family: dict[str, dict]) -> dict:
    for family, fit in per_family.items():
        key = family.lower()
        row[f"ls_rmse_{key}"] = fit.get("ls_rmse", np.nan)
        row[f"aicc_{key}"] = fit.get("aicc", np.nan)
        for name in PARAM_NAMES[family]:
            row[f"ls_{name}"] = fit.get("ls_params", {}).get(name, np.nan)
            row[f"mle_{name}"] = fit.get("mle_params", {}).get(name, np.nan)

    row["ls_rmse_winner"] = _argmin_family(per_family, "ls_rmse", admissible_only=False)
    row["ls_rmse_winner_admissible"] = _argmin_family(per_family, "ls_rmse", admissible_only=True)
    row["aicc_winner"] = _argmin_family(per_family, "aicc", admissible_only=False)
    return row


def _argmin_family(per_family: dict[str, dict], metric: str, admissible_only: bool) -> str | None:
    best, best_value = None, np.inf
    for family, fit in per_family.items():
        if not fit.get("ls_ok" if metric == "ls_rmse" else "mle_ok", False):
            continue
        value = fit.get(metric, np.nan)
        if not np.isfinite(value):
            continue
        if admissible_only and not _admissible_in_replicate(family, fit.get("ls_params", {})):
            continue
        if value < best_value:
            best, best_value = family, value
    return best


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def run_bootstrap(districts: pd.DataFrame, district_column: str = "DISTRICT",
                  value_column: str = "As_model_value_mgL", censored_column: str = "As_censored",
                  replicates: int = 1000, master_seed: int = 20260923,
                  tie_tolerance: float = 1e-9, jobs: int = 1,
                  district_order: list[str] | None = None, progress=None) -> pd.DataFrame:
    """Run the full bootstrap and return one long table of replicate results.

    ``districts`` holds every selected district's rows; ``district_order`` fixes the
    ``district_index`` used in the seed derivation (defaults to sorted order, which
    must match the order the runner passes so seeds are reproducible).
    """
    order = district_order if district_order is not None else sorted(districts[district_column].unique())
    specs = []
    for district_index, district in enumerate(order):
        part = districts[districts[district_column] == district]
        if part.empty:
            raise ValueError(f"no rows for district {district!r}")
        specs.append(BootstrapSpec(
            district=district,
            district_index=district_index,
            values_mgL=part[value_column].to_numpy(dtype=float),
            censored=part[censored_column].to_numpy(dtype=bool),
            master_seed=master_seed,
            tie_tolerance=tie_tolerance,
        ))

    jobs_list = [(spec, rep) for spec in specs for rep in range(replicates)]
    if jobs is None or jobs <= 1:
        rows = [bootstrap_replicate(spec, rep) for spec, rep in jobs_list]
        if progress is not None:
            progress(len(rows), len(rows))
    else:
        rows = _run_parallel(jobs_list, jobs, progress)

    table = pd.DataFrame(rows)
    return table.sort_values(["district_index", "replicate"], kind="stable").reset_index(drop=True)


def _run_parallel(jobs_list, jobs: int, progress) -> list[dict]:
    """Deterministic process-pool map: chunks are ordered, so results re-sort cleanly."""
    import os
    from multiprocessing import get_context

    chunk = max(1, len(jobs_list) // (jobs * 8) or 1)
    ctx = get_context("fork") if "fork" in __import__("multiprocessing").get_all_start_methods() else get_context()
    done = 0
    rows: list[dict] = []
    env_jobs = min(jobs, os.cpu_count() or 1)
    with ctx.Pool(env_jobs) as pool:
        for result in pool.imap(_replicate_worker, jobs_list, chunksize=chunk):
            rows.append(result)
            done += 1
            if progress is not None and done % 50 == 0:
                progress(done, len(jobs_list))
    if progress is not None:
        progress(len(rows), len(jobs_list))
    return rows


def _replicate_worker(job) -> dict:
    spec, replicate = job
    return bootstrap_replicate(spec, replicate)


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def summarise_bootstrap(table: pd.DataFrame) -> pd.DataFrame:
    """Per district x family: parameter median/IQR, RMSE spread, and winner frequency."""
    rows = []
    total = table.groupby("DISTRICT")["replicate"].nunique()
    for (district, ), part in table.groupby(["DISTRICT"], sort=True):
        usable = part[part["status"] == "ok"]
        n_total = int(total.loc[district])
        winner_counts = {
            column: usable[column].value_counts(dropna=True) for column in WINNER_COLUMNS
        }
        for family in FAMILIES:
            key = family.lower()
            row: dict = {
                "DISTRICT": district,
                "distribution": family,
                "bootstrap_replicates_requested": n_total,
                "bootstrap_replicates_used": int(len(usable)),
                "bootstrap_replicates_skipped": int(len(part) - len(usable)),
                "ls_rmse_median": _safe_quantile(usable[f"ls_rmse_{key}"], 50.0),
                "ls_rmse_q025": _safe_quantile(usable[f"ls_rmse_{key}"], 2.5),
                "ls_rmse_q975": _safe_quantile(usable[f"ls_rmse_{key}"], 97.5),
                "aicc_median": _safe_quantile(usable[f"aicc_{key}"], 50.0),
            }
            for name in PARAM_NAMES[family]:
                values = usable[f"ls_{name}"].to_numpy(dtype=float)
                finite = values[np.isfinite(values)]
                row[f"ls_{name}_median"] = float(np.median(finite)) if finite.size else np.nan
                row[f"ls_{name}_q025"] = _safe_quantile(values, 2.5)
                row[f"ls_{name}_q975"] = _safe_quantile(values, 97.5)
                row[f"ls_{name}_rmse_spread"] = (float(np.sqrt(np.mean((finite - finite.mean()) ** 2)))
                                                 if finite.size else np.nan)
                mle_values = usable[f"mle_{name}"].to_numpy(dtype=float)
                row[f"mle_{name}_median"] = float(np.median(mle_values[np.isfinite(mle_values)])) \
                    if np.any(np.isfinite(mle_values)) else np.nan
            for column in WINNER_COLUMNS:
                count = int(winner_counts[column].get(family, 0))
                row[f"{column}_frequency"] = count / len(usable) if len(usable) else np.nan
                row[f"{column}_count"] = count
            rows.append(row)
    return pd.DataFrame(rows)


def _safe_quantile(values, q: float) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    return float(np.percentile(finite, q)) if finite.size else np.nan


def bootstrap_rationale(summary: pd.DataFrame, district: str) -> str:
    """One-line stability verdict for the selected family of one district."""
    part = summary[summary["DISTRICT"] == district]
    if part.empty:
        return ""
    pieces = []
    for row in part.itertuples():
        pieces.append(
            f"{row.distribution}: LS-RMSE winner {row.ls_rmse_winner_frequency:.1%}, "
            f"admissible LS-RMSE winner {row.ls_rmse_winner_admissible_frequency:.1%}, "
            f"AICc winner {row.aicc_winner_frequency:.1%} "
            f"({row.bootstrap_replicates_used} usable replicates)"
        )
    return " | ".join(pieces)
