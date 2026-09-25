"""Phases 3 and 5: body-weight preprocessing and weighted distribution fitting.

Reads the raw STEPS (bgd2018.csv) and MICS7 (ch.sav) files named in
config/run_config.json (never modified) and writes, under results/bodyweight/:

    adult_bw_clean.csv, adult_bw_cleaning_audit.csv, adult_bw_bmi_review.csv
    child_bw_clean.csv, child_bw_cleaning_audit.csv, child_bw_an8_codebook_review.csv
    bodyweight_weighted_summary.csv
    bodyweight_distribution_fit_results.csv
    bodyweight_selected_distributions.csv / .json
    bodyweight_manifest.json

and diagnostic figures results/figures/bodyweight_<subject>_<population>.png.

Exits non-zero if a raw input differs from results/provenance/raw_data.sha256,
changes during the run, or an audit fails to reconcile.

Usage:  python scripts/run_bodyweight.py [--no-figures]
"""
from __future__ import annotations

import argparse
import importlib.metadata
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd  # noqa: E402

from arsenic_hra import bodyweight_fitting as bf  # noqa: E402
from arsenic_hra import bodyweight_preprocessing as bp  # noqa: E402
from arsenic_hra.paths import (RESULTS_BODYWEIGHT, RESULTS_FIGURES, RESULTS_PROVENANCE,  # noqa: E402
                               ensure_output_dir, load_run_config, raw_path, rel)
from arsenic_hra.provenance import read_hash_manifest, sha256_file  # noqa: E402

OUT = {
    "adult_clean": RESULTS_BODYWEIGHT / "adult_bw_clean.csv",
    "adult_audit": RESULTS_BODYWEIGHT / "adult_bw_cleaning_audit.csv",
    "adult_bmi_review": RESULTS_BODYWEIGHT / "adult_bw_bmi_review.csv",
    "child_clean": RESULTS_BODYWEIGHT / "child_bw_clean.csv",
    "child_audit": RESULTS_BODYWEIGHT / "child_bw_cleaning_audit.csv",
    "child_an8_review": RESULTS_BODYWEIGHT / "child_bw_an8_codebook_review.csv",
    "summary": RESULTS_BODYWEIGHT / "bodyweight_weighted_summary.csv",
    "fit_results": RESULTS_BODYWEIGHT / "bodyweight_distribution_fit_results.csv",
    "selected_csv": RESULTS_BODYWEIGHT / "bodyweight_selected_distributions.csv",
    "selected_json": RESULTS_BODYWEIGHT / "bodyweight_selected_distributions.json",
    "manifest": RESULTS_BODYWEIGHT / "bodyweight_manifest.json",
}
PACKAGES = ["numpy", "pandas", "scipy", "matplotlib", "pyreadstat"]


def check_config(bcfg: dict) -> None:
    """The module's column choices must match the run configuration."""
    expected = {
        "adult_weight_column": "m12",
        "adult_survey_weight_column": "wstep2",
        "child_weight_column": "AN8",
        "child_survey_weight_column": "chweight",
    }
    for key, column in expected.items():
        if bcfg[key] != column:
            raise ValueError(f"config bodyweight.{key}={bcfg[key]!r}, module uses {column!r}")
    if [c.lower() for c in bcfg["candidate_distributions"]] != [f.lower() for f in bf.FAMILIES]:
        raise ValueError(f"config candidates {bcfg['candidate_distributions']} != {list(bf.FAMILIES)}")


def audit_reconciles(audit: pd.DataFrame, n_clean: int) -> bool:
    chained = (audit["rows_before"].iloc[1:].to_numpy() == audit["rows_after"].iloc[:-1].to_numpy()).all()
    steps = (audit["rows_before"] - audit["removed"] == audit["rows_after"]).all()
    total = audit["rows_after"].iloc[0] - audit["removed"].sum() == n_clean
    return bool(chained and steps and total)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--no-figures", action="store_true", help="skip diagnostic figures")
    args = parser.parse_args(argv)

    cfg = load_run_config()
    check_config(cfg["bodyweight"])
    inputs = {name: raw_path(cfg["raw_inputs"][name]) for name in ("adult_bw_steps", "child_bw_mics")}
    frozen_hashes = read_hash_manifest(RESULTS_PROVENANCE / "raw_data.sha256")
    hash_before = {name: sha256_file(p) for name, p in inputs.items()}
    mismatched = [rel(p) for name, p in inputs.items() if frozen_hashes.get(rel(p)) != hash_before[name]]
    if mismatched:
        print(f"raw inputs differ from the frozen manifest: {mismatched}")
        return 1

    # Phase 3: preprocessing
    steps_raw = bp.load_steps_bw(inputs["adult_bw_steps"])
    mics_raw = bp.load_mics_bw(inputs["child_bw_mics"])
    adult, adult_audit = bp.clean_steps_bw(steps_raw)
    child, child_audit = bp.clean_mics_bw(mics_raw)
    summary = pd.DataFrame([
        bp.weighted_summary(adult["BW_kg"], adult["survey_weight"], "adult"),
        bp.weighted_summary(child["BW_kg"], child["survey_weight"], "child"),
    ])

    ensure_output_dir(RESULTS_BODYWEIGHT)
    adult.to_csv(OUT["adult_clean"], index=False)
    adult_audit.to_csv(OUT["adult_audit"], index=False)
    bp.adult_bmi_review(adult).to_csv(OUT["adult_bmi_review"], index=False)
    child.to_csv(OUT["child_clean"], index=False)
    child_audit.to_csv(OUT["child_audit"], index=False)
    bp.child_an8_codebook_review(mics_raw).to_csv(OUT["child_an8_review"], index=False)
    summary.to_csv(OUT["summary"], index=False)

    # Phase 5: fitting and provisional selection
    fit_tables, selected, figures = [], [], []
    for population, df in (("adult", adult), ("child", child)):
        table = bf.fit_population(df["BW_kg"], df["survey_weight"], population)
        chosen, rationale = bf.select_bw_distribution(table)
        table["selected"] = table["distribution"].eq(chosen["distribution"])
        fit_tables.append(table)
        selected.append(bf.selected_record(chosen, rationale))
        if not args.no_figures:
            from arsenic_hra.bodyweight_plots import plot_population

            figures += plot_population(df["BW_kg"], df["survey_weight"], table, chosen["distribution"],
                                       population, RESULTS_FIGURES)

    pd.concat(fit_tables, ignore_index=True).to_csv(OUT["fit_results"], index=False)
    pd.DataFrame([{
        "population": s["population"],
        "distribution": s["distribution"],
        **{f"param_{k}": v for k, v in s["params"].items()},
        **{f"scipy_{k}": v for k, v in s["scipy"].items()},
        "weighted_cdf_sse": s["weighted_cdf_sse"],
        "weighted_cdf_max_discrepancy": s["weighted_cdf_max_discrepancy"],
        "n": s["n"],
        "status": s["status"],
        "rationale": s["rationale"],
    } for s in selected]).to_csv(OUT["selected_csv"], index=False)
    OUT["selected_json"].write_text(json.dumps(selected, indent=2) + "\n", encoding="utf-8")

    hash_after = {name: sha256_file(p) for name, p in inputs.items()}
    reconciled = {"adult": audit_reconciles(adult_audit, len(adult)), "child": audit_reconciles(child_audit, len(child))}
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/run_bodyweight.py",
        "config_version": cfg["config_version"],
        "raw_inputs": {rel(p): hash_before[name] for name, p in inputs.items()},
        "raw_matches_frozen_manifest": True,
        "raw_unchanged_during_run": hash_before == hash_after,
        "audits_reconcile": reconciled,
        "software": {"python": platform.python_version(),
                     **{name: importlib.metadata.version(name) for name in PACKAGES}},
        "populations": {
            "adult": {"source": "STEPS 2018 m12 (kg), survey weight wstep2", "age": "18-69 years",
                      "rows_raw": len(steps_raw), "rows_clean": len(adult)},
            "child": {"source": "MICS7 2025 ch.sav AN8 (kg), survey weight chweight",
                      "age": "0-59 months (under-5 only)", "rows_raw": len(mics_raw), "rows_clean": len(child)},
        },
        "fitting": {
            "candidates": list(bf.FAMILIES),
            "primary": "survey-weighted CDF least squares on midpoint plotting positions (scipy.optimize.least_squares)",
            "robustness": "survey-weighted pseudo-MLE (scipy.optimize.minimize)",
            "normal_max_negative_prob": bf.NORMAL_MAX_NEGATIVE_PROB,
            "triangular_max_excluded_mass": bf.TRIANGULAR_MAX_EXCLUDED_MASS,
            "selection_rule": "lowest weighted CDF SSE among physically admissible candidates; provisional",
        },
        "random_seed": None,
        "outputs": {k: rel(v) for k, v in OUT.items() if k != "manifest"},
        "figures": [rel(p) for p in figures],
    }
    OUT["manifest"].write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    for s in selected:
        print(f"{s['population']:>5}: {s['distribution']} {s['params']}  [{s['status']}]")
        print(f"       {s['rationale']}")
    if hash_before != hash_after or not all(reconciled.values()):
        print(f"FAILED: raw unchanged={hash_before == hash_after}; audits reconcile={reconciled}")
        return 1
    print(f"audits reconcile; raw inputs match the frozen manifest; outputs in {rel(RESULTS_BODYWEIGHT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
