"""Run the body-weight track end to end (Phases 3 and 5).

Reads the raw STEPS and MICS files (never modifies them) and writes every
output to results/bodyweight/:

    adult_bw_clean.csv, adult_bw_cleaning_audit.csv
    child_bw_clean.csv, child_bw_cleaning_audit.csv
    bodyweight_weighted_summary.csv
    bodyweight_distribution_fit_results.csv
    bodyweight_selected_distributions.csv / .json
    bodyweight_run_manifest.json
    figures/{adult,child}/*.png

Usage:  python scripts/run_bodyweight.py [--no-figures]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import matplotlib
import numpy as np
import pandas as pd
import pyreadstat
import scipy

from bodyweight import fitting, preprocess
from bodyweight.paths import FIGURES_DIR, MICS_CH_PATH, PROJECT_ROOT, RESULTS_DIR, STEPS_PATH


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rel(path: Path) -> str:
    return path.resolve().relative_to(PROJECT_ROOT).as_posix()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--no-figures", action="store_true", help="skip diagnostic figures")
    args = parser.parse_args(argv)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    input_hashes_before = {rel(p): sha256(p) for p in (STEPS_PATH, MICS_CH_PATH)}

    # Phase 3: preprocessing
    adult, adult_audit = preprocess.clean_steps_bw(preprocess.load_steps_bw())
    child, child_audit = preprocess.clean_mics_bw(preprocess.load_mics_bw())
    adult.to_csv(RESULTS_DIR / "adult_bw_clean.csv", index=False)
    adult_audit.to_csv(RESULTS_DIR / "adult_bw_cleaning_audit.csv", index=False)
    child.to_csv(RESULTS_DIR / "child_bw_clean.csv", index=False)
    child_audit.to_csv(RESULTS_DIR / "child_bw_cleaning_audit.csv", index=False)
    summary = pd.DataFrame([
        preprocess.weighted_summary(adult["BW_kg"], adult["survey_weight"], "adult"),
        preprocess.weighted_summary(child["BW_kg"], child["survey_weight"], "child"),
    ])
    summary.to_csv(RESULTS_DIR / "bodyweight_weighted_summary.csv", index=False)

    # Phase 5: fitting and provisional selection
    fit_tables, selected, figures = [], [], []
    for population, df in (("adult", adult), ("child", child)):
        table = fitting.fit_population(df["BW_kg"], df["survey_weight"], population)
        chosen, rationale = fitting.select_bw_distribution(table)
        table["selected"] = table["distribution"].eq(chosen["distribution"])
        fit_tables.append(table)
        selected.append(fitting.selected_record(chosen, rationale))
        if not args.no_figures:
            from bodyweight.plots import plot_population

            figures += plot_population(df["BW_kg"], df["survey_weight"], table, chosen["distribution"],
                                       population, FIGURES_DIR / population)

    fit_results = pd.concat(fit_tables, ignore_index=True)
    fit_results.to_csv(RESULTS_DIR / "bodyweight_distribution_fit_results.csv", index=False)
    selected_flat = pd.DataFrame([{
        "population": s["population"],
        "distribution": s["distribution"],
        **{f"param_{k}": v for k, v in s["params"].items()},
        **{f"scipy_{k}": v for k, v in s["scipy"].items()},
        "weighted_cdf_sse": s["weighted_cdf_sse"],
        "weighted_cdf_max_discrepancy": s["weighted_cdf_max_discrepancy"],
        "n": s["n"],
        "status": s["status"],
        "rationale": s["rationale"],
    } for s in selected])
    selected_flat.to_csv(RESULTS_DIR / "bodyweight_selected_distributions.csv", index=False)
    (RESULTS_DIR / "bodyweight_selected_distributions.json").write_text(json.dumps(selected, indent=2))

    input_hashes_after = {rel(p): sha256(p) for p in (STEPS_PATH, MICS_CH_PATH)}
    if input_hashes_after != input_hashes_before:
        raise RuntimeError("raw input file changed during the run")

    manifest = {
        "run_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": rel(Path(__file__)),
        "inputs_sha256": input_hashes_before,
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "pyreadstat": pyreadstat.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "populations": {
            "adult": {"source": "STEPS 2018 m12 (kg), weight wstep2", "age": "18-69 years",
                      "rows_raw": int(adult_audit["rows_after"].iloc[0]), "rows_clean": len(adult)},
            "child": {"source": "MICS7 ch.sav AN8 (kg), weight chweight", "age": "0-59 months (under-5 only)",
                      "rows_raw": int(child_audit["rows_after"].iloc[0]), "rows_clean": len(child)},
        },
        "fitting": {
            "candidates": list(fitting.FAMILIES),
            "primary": "survey-weighted CDF least squares, midpoint plotting positions, scipy.optimize.least_squares",
            "robustness": "survey-weighted pseudo-MLE, scipy.optimize.minimize",
            "normal_max_negative_prob": fitting.NORMAL_MAX_NEGATIVE_PROB,
            "triangular_max_excluded_mass": fitting.TRIANGULAR_MAX_EXCLUDED_MASS,
            "selection_rule": "lowest weighted CDF SSE among physically admissible candidates; provisional",
        },
        "random_seed": None,
        "outputs": sorted(rel(p) for p in RESULTS_DIR.glob("*.*")) + sorted(rel(p) for p in figures),
    }
    (RESULTS_DIR / "bodyweight_run_manifest.json").write_text(json.dumps(manifest, indent=2))

    for s in selected:
        print(f"{s['population']:>5}: {s['distribution']} {s['params']}  [{s['status']}]")
        print(f"       {s['rationale']}")
    print(f"outputs written to {rel(RESULTS_DIR)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
