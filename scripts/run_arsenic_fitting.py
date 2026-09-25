"""Phase 4: fit a censoring-aware arsenic distribution to each approved district.

Reads results/arsenic/arsenic_selected_districts_preprocessed.csv (Phase 2 output;
read-only) and writes, under results/arsenic/:

    arsenic_distribution_fit_results.csv   10 districts x 3 families
    arsenic_bootstrap_stability.csv        parameter CIs, RMSE spread, winner frequency
    arsenic_selected_distributions.csv     one provisional choice per district
    arsenic_selected_distributions.json    the same, machine-readable
    arsenic_fitting_manifest.json          provenance for the whole run

and diagnostic figures results/figures/arsenic_<subject>_<district>.png (4 per district).

Method (plan.md section 3): reverse Kaplan-Meier empirical CDF -> CDF least squares
(primary) + censored maximum likelihood (robustness) -> bootstrap stability ->
admissibility gate -> provisional selection. Censored wells are never substituted.

Exits non-zero if the raw BGS file differs from results/provenance/raw_data.sha256,
changes during the run, or a selection is missing.

Usage:  python scripts/run_arsenic_fitting.py [--bootstrap N] [--no-figures]
                                               [--jobs N] [--districts A B ...]
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

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from arsenic_hra import arsenic_bootstrap as ab  # noqa: E402
from arsenic_hra import arsenic_fitting as af  # noqa: E402
from arsenic_hra.paths import (RESULTS_ARSENIC, RESULTS_FIGURES, RESULTS_PROVENANCE,  # noqa: E402
                               ensure_output_dir, load_run_config, raw_path, rel)
from arsenic_hra.provenance import read_hash_manifest, sha256_file  # noqa: E402

OUT = {
    "input": RESULTS_ARSENIC / "arsenic_selected_districts_preprocessed.csv",
    "summary": RESULTS_ARSENIC / "arsenic_district_preprocessing_summary.csv",
    "fit_results": RESULTS_ARSENIC / "arsenic_distribution_fit_results.csv",
    "bootstrap": RESULTS_ARSENIC / "arsenic_bootstrap_stability.csv",
    "selected_csv": RESULTS_ARSENIC / "arsenic_selected_distributions.csv",
    "selected_json": RESULTS_ARSENIC / "arsenic_selected_distributions.json",
    "manifest": RESULTS_ARSENIC / "arsenic_fitting_manifest.json",
}
PACKAGES = ["numpy", "pandas", "scipy", "matplotlib", "lifelines"]


def check_config(acfg: dict) -> None:
    """The module's families, objective region and selection rule must match config."""
    configured = [c.lower() for c in acfg["candidate_distributions"]]
    if configured != [f.lower() for f in af.FAMILIES]:
        raise ValueError(f"config candidates {acfg['candidate_distributions']} != {list(af.FAMILIES)}")
    if acfg["objective_region"] != af.OBJECTIVE_REGION:
        raise ValueError(f"config objective_region {acfg['objective_region']!r} != {af.OBJECTIVE_REGION!r}")
    if acfg["selection_rule"] != af.SELECTION_RULE:
        raise ValueError(f"config selection_rule {acfg['selection_rule']!r} != {af.SELECTION_RULE!r}")
    if float(acfg["lognormal_gamma_location"]) != 0.0:
        raise ValueError("the plan fixes the Lognormal/Gamma location at 0")
    if acfg["censored_substitution"] != "none":
        raise ValueError("the plan forbids substituting censored values")


def _progress(done: int, total: int) -> None:
    print(f"  bootstrap {done}/{total} replicates", flush=True)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bootstrap", type=int, default=None,
                        help="override the frozen bootstrap replicate count (for quick runs)")
    parser.add_argument("--no-figures", action="store_true", help="skip diagnostic figures")
    parser.add_argument("--jobs", type=int, default=1,
                        help="parallel processes for the bootstrap (default 1, deterministic)")
    parser.add_argument("--districts", nargs="+", default=None,
                        help="restrict the run to these districts (smoke runs only)")
    args = parser.parse_args(argv)

    cfg = load_run_config()
    acfg = cfg["arsenic"]
    check_config(acfg)
    af.configure_constants(acfg)

    districts = list(acfg["selected_districts"])
    if args.districts:
        unknown = [d for d in args.districts if d not in districts]
        if unknown:
            raise SystemExit(f"--districts contains names outside the approved list: {unknown}")
        districts = [d for d in districts if d in set(args.districts)]
    replicates = acfg["bootstrap_replicates"] if args.bootstrap is None else int(args.bootstrap)
    if replicates < 1:
        raise SystemExit("--bootstrap must be at least 1")

    # Provenance: the raw BGS file must match the frozen manifest, before and after.
    raw = raw_path(cfg["raw_inputs"]["arsenic_bgs"])
    frozen_hashes = read_hash_manifest(RESULTS_PROVENANCE / "raw_data.sha256")
    hash_before = sha256_file(raw)
    if frozen_hashes.get(rel(raw)) != hash_before:
        print(f"raw input differs from the frozen manifest: {rel(raw)}")
        return 1

    if not OUT["input"].is_file():
        print(f"missing Phase 2 output {rel(OUT['input'])}; run scripts/preprocess_arsenic.py first")
        return 1
    data = pd.read_csv(OUT["input"])
    summary = pd.read_csv(OUT["summary"]).set_index("DISTRICT") if OUT["summary"].is_file() else None

    missing = [d for d in districts if not (data["DISTRICT"] == d).any()]
    if missing:
        print(f"approved districts absent from the preprocessed file: {missing}")
        return 1

    # ---- Steps A, B, C, E: fit and select every district ----
    fit_tables, selected, figures = [], [], []
    for district in districts:
        part = data[data["DISTRICT"] == district]
        values = part["As_model_value_mgL"].to_numpy(dtype=float)
        censored = part["As_censored"].to_numpy(dtype=bool)

        table = af.fit_district(district, values, censored)
        if summary is not None:
            src = summary.loc[district]
            table["high_censoring_warning"] = bool(src["high_censoring_warning"])
            table["minimum_detected_As_mgL"] = float(src["minimum_detected_As_mgL"])
        chosen, rationale = af.select_arsenic_distribution(
            table, warning_percent=float(acfg["high_censoring_warning_percent"]),
            min_detected_mgL=float(table["minimum_detected_As_mgL"].iloc[0])
            if "minimum_detected_As_mgL" in table else None,
            curve=af.reverse_km_cdf(values, censored))
        table["selected"] = table["distribution"].eq(chosen["distribution"])
        table["district_index"] = districts.index(district)
        fit_tables.append(table)
        selected.append(af.selected_record(chosen, rationale))

        if not args.no_figures:
            from arsenic_hra.arsenic_plots import plot_district

            figures += plot_district(district, values, censored, table, chosen["distribution"],
                                     summary.loc[district] if summary is not None else None, RESULTS_FIGURES)

    fit_results = pd.concat(fit_tables, ignore_index=True)

    # ---- Step D: bootstrap stability ----
    print(f"bootstrapping {len(districts)} districts x {replicates} replicates "
          f"(seed {acfg['bootstrap_seed']}, jobs {args.jobs})", flush=True)
    boot_rows = ab.run_bootstrap(
        data[data["DISTRICT"].isin(districts)],
        replicates=replicates,
        master_seed=int(acfg["bootstrap_seed"]),
        tie_tolerance=float(acfg["bootstrap_tie_tolerance"]),
        jobs=args.jobs,
        district_order=districts,
        progress=_progress if args.jobs == 1 else None,
    )
    stability = ab.summarise_bootstrap(boot_rows)

    # Fold the bootstrap verdict into each selection's rationale.
    for record in selected:
        stability_note = None
        part = stability[stability["DISTRICT"] == record["district"]]
        row = part[part["distribution"] == record["distribution"]]
        if len(row):
            r = row.iloc[0]
            stability_note = (f"bootstrap ({int(r['bootstrap_replicates_used'])} usable replicates): "
                              f"{r['distribution']} wins admissible LS-RMSE in "
                              f"{r['ls_rmse_winner_admissible_frequency']:.1%} of replicates and AICc in "
                              f"{r['aicc_winner_frequency']:.1%}; median RMSE "
                              f"{r['ls_rmse_median']:.4f} "
                              f"[{r['ls_rmse_q025']:.4f}, {r['ls_rmse_q975']:.4f}]")
            record["rationale"] = f"{record['rationale']}; {stability_note}"
        record["bootstrap_stability"] = _stability_record(stability, record["district"], record["distribution"])

    # ---- Outputs ----
    ensure_output_dir(RESULTS_ARSENIC)
    float_fmt = "%.10g"
    fit_results.to_csv(OUT["fit_results"], index=False, float_format=float_fmt)
    stability.to_csv(OUT["bootstrap"], index=False, float_format=float_fmt)
    af.selection_table(selected).to_csv(OUT["selected_csv"], index=False, float_format=float_fmt)
    OUT["selected_json"].write_text(json.dumps(selected, indent=2) + "\n", encoding="utf-8")

    hash_after = sha256_file(raw)
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/run_arsenic_fitting.py",
        "config_version": cfg["config_version"],
        "input": rel(OUT["input"]),
        "input_sha256": sha256_file(OUT["input"]),
        "raw_input": rel(raw),
        "raw_sha256": hash_before,
        "raw_matches_frozen_manifest": True,
        "raw_unchanged_during_run": hash_before == hash_after,
        "districts": districts,
        "software": {"python": platform.python_version(),
                     **{name: importlib.metadata.version(name) for name in PACKAGES}},
        "fitting": {
            "candidates": list(af.FAMILIES),
            "primary": "censoring-aware (reverse Kaplan-Meier) empirical CDF + CDF least squares",
            "robustness": "censored maximum likelihood (scipy.stats.CensoredData)",
            "censored_substitution": "none",
            "lognormal_gamma_location": 0.0,
            "unweighted": True,
            "objective_region": af.OBJECTIVE_REGION,
            "objective_region_rationale": acfg.get("objective_region_rationale", ""),
            "ls_multi_start": ["moments_detected", "censored_mle"],
            "normal_max_nonpositive_prob": af.NORMAL_MAX_NONPOSITIVE_PROB,
            "parameter_magnitude_band": {"min": af.PARAM_MIN_MAGNITUDE, "max": af.PARAM_MAX_MAGNITUDE},
            "selection_rule": af.SELECTION_RULE,
            "selection_status": af.SELECTION_STATUS,
            "high_censoring_warning_percent": acfg["high_censoring_warning_percent"],
        },
        "bootstrap": {
            "replicates_requested": replicates,
            "replicates_in_config": acfg["bootstrap_replicates"],
            "master_seed": int(acfg["bootstrap_seed"]),
            "seed_derivation": "numpy.random.SeedSequence(master_seed, spawn_key=(district_index, replicate))",
            "districts_order": districts,
            "jobs": args.jobs,
            "resampling": "rows with replacement, preserving each row's As_censored flag and value",
            "tie_tolerance": float(acfg["bootstrap_tie_tolerance"]),
            "usable_replicates_per_district": {
                d: int(v) for d, v in stability.groupby("DISTRICT")["bootstrap_replicates_used"].max().items()},
        },
        "rows": {"districts": len(districts), "families": len(af.FAMILIES),
                 "fit_rows": len(fit_results), "bootstrap_rows": len(boot_rows)},
        "outputs": {k: rel(v) for k, v in OUT.items() if k != "manifest"},
        "figures": [rel(p) for p in figures],
    }
    OUT["manifest"].write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    # ---- Report ----
    print()
    for record in selected:
        params = ", ".join(f"{k}={v:.6g}" for k, v in record["params"].items())
        print(f"{record['district']:>11}: {record['distribution']:<9} {params}")
        print(f"             RMSE {record['cdf_rmse']:.4f} on {record['objective_points']} points; "
              f"censoring {record['censoring_percent']:.2f}%")
        print(f"             {record['rationale']}")

    if hash_before != hash_after:
        print("FAILED: the raw BGS file changed during the run")
        return 1
    if len(selected) != len(districts):
        print(f"FAILED: {len(selected)} selections for {len(districts)} districts")
        return 1
    print(f"\n{len(selected)} districts fitted; raw input matches the frozen manifest and is unchanged; "
          f"outputs in {rel(RESULTS_ARSENIC)}")
    return 0


def _stability_record(stability: pd.DataFrame, district: str, family: str) -> dict:
    part = stability[(stability["DISTRICT"] == district) & (stability["distribution"] == family)]
    if part.empty:
        return {}
    row = part.iloc[0].to_dict()
    return {k: (float(v) if isinstance(v, (int, float, np.floating, np.integer)) and not pd.isna(v) else v)
            for k, v in row.items()}


if __name__ == "__main__":
    raise SystemExit(main())
