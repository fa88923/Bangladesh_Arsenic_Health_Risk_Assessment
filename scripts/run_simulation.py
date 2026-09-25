"""Phase 7: primary Monte Carlo simulation (N = primary_N, replicate 0) for every
district x population, plus the deterministic district benchmark.

Writes:

    results/simulation/simulation_input_distributions.json / .csv   final C and BW sampling models
    results/simulation/arsenic_selection_review.csv                  Lognormal vs Gamma evidence per district
    results/simulation/child_bw_selection_review.csv                 Gamma vs truncated Normal (child BW)
    results/simulation/iterations/<run_id>.parquet                   iteration-level inputs and outputs
    results/simulation/manifests/run_manifest_<run_id>.json          seed, spawn key, distributions, hashes
    results/simulation/simulation_summary.csv                        summaries of HI, ELCR, HQ per run
    results/simulation/sampling_validation.csv                       sampled inputs vs target distributions
    results/simulation/simulation_manifest.json
    results/tables/risk_probabilistic_p95.csv                        base-paper Table 5 layout
    results/tables/risk_exceedance.csv                               P(HI>1), P(HI>2), P(ELCR>1e-4) with MC SE
    results/tables/risk_deterministic_district.csv                   base-paper Table 4 layout
    results/tables/risk_deterministic_vs_probabilistic.csv
    results/figures/simulation_hi_distribution_<district>.png

Exits non-zero if a rerun with the same seed differs, a sampled input is
non-physical, or a sampled input fails the KS check at the 0.1% level.

Usage:  python scripts/run_simulation.py [--no-figures]
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

from arsenic_hra import monte_carlo as mc  # noqa: E402
from arsenic_hra import risk_equations as re_  # noqa: E402
from arsenic_hra import simulation_inputs as si  # noqa: E402
from arsenic_hra import simulation_plots as sp  # noqa: E402
from arsenic_hra.paths import (RESULTS_FIGURES, RESULTS_SIMULATION, RESULTS_TABLES, RUN_CONFIG_PATH,  # noqa: E402
                               ensure_output_dir, load_run_config, rel)
from arsenic_hra.provenance import sha256_file  # noqa: E402

ITER_DIR = RESULTS_SIMULATION / "iterations"
MANIFEST_DIR = RESULTS_SIMULATION / "manifests"
OUT = {
    "inputs_json": RESULTS_SIMULATION / "simulation_input_distributions.json",
    "inputs_csv": RESULTS_SIMULATION / "simulation_input_distributions.csv",
    "arsenic_review": RESULTS_SIMULATION / "arsenic_selection_review.csv",
    "child_bw_review": RESULTS_SIMULATION / "child_bw_selection_review.csv",
    "summary": RESULTS_SIMULATION / "simulation_summary.csv",
    "sampling": RESULTS_SIMULATION / "sampling_validation.csv",
    "p95": RESULTS_TABLES / "risk_probabilistic_p95.csv",
    "exceedance": RESULTS_TABLES / "risk_exceedance.csv",
    "deterministic": RESULTS_TABLES / "risk_deterministic_district.csv",
    "det_vs_prob": RESULTS_TABLES / "risk_deterministic_vs_probabilistic.csv",
    "manifest": RESULTS_SIMULATION / "simulation_manifest.json",
}
BW_SUMMARY = si.RESULTS_BODYWEIGHT / "bodyweight_weighted_summary.csv"
UPSTREAM = [si.ARSENIC_SELECTED, si.ARSENIC_FITS, si.ARSENIC_PREPROCESSED, si.BW_SELECTED, si.CHILD_BW_CLEAN,
            BW_SUMMARY]
PACKAGES = ["numpy", "pandas", "scipy", "pyarrow", "matplotlib"]
KS_GATE = 1.949  # 0.1% two-sided KS critical value coefficient: D <= 1.949 / sqrt(N)


def arsenic_review(model: mc.ModelInputs) -> pd.DataFrame:
    fits = pd.read_csv(si.ARSENIC_FITS)
    rows = []
    for district, rec in model.arsenic.items():
        for _, r in fits[(fits.DISTRICT == district) & (fits.distribution != "Normal")].iterrows():
            names = si.af_param_names(r.distribution)
            d = si.frozen(r.distribution, {n: float(r[f"ls_{n}"]) for n in names})
            rows.append({
                "district": district, "family": r.distribution,
                "provisional": bool(r.selected), "final": r.distribution == rec["distribution"],
                "cdf_rmse": r.cdf_rmse, "cdf_max_discrepancy": r.cdf_max_discrepancy, "aicc": r.aicc,
                "empirical_P95_mgL": r.empirical_P95_mgL, "fitted_P95_mgL": r.fitted_P95_mgL,
                "empirical_P99_mgL": r.empirical_P99_mgL, "fitted_P99_mgL": r.fitted_P99_mgL,
                "P95_ratio": r.fitted_P95_mgL / r.empirical_P95_mgL, "P99_ratio": r.fitted_P99_mgL / r.empirical_P99_mgL,
                "reverse_km_mean_mgL": rec["reverse_km_mean_mgL"], "fitted_mean_mgL": float(d.mean()),
                "mean_ratio": float(d.mean()) / rec["reverse_km_mean_mgL"], "fitted_P99_9_mgL": float(d.ppf(0.999)),
            })
    return pd.DataFrame(rows)


def deterministic_district(cfg: dict, model: mc.ModelInputs, bw_means: dict) -> pd.DataFrame:
    rows = []
    for district, rec in model.arsenic.items():
        row = {"district": district, "C_reverse_km_mean_mgL": rec["reverse_km_mean_mgL"]}
        for pop in cfg["simulation"]["populations"]:
            p = model.params[pop]
            out = re_.evaluate_risk(re_.deterministic_inputs(p, rec["reverse_km_mean_mgL"], bw_means[pop]), p)
            row[f"HI_{pop}"] = float(out["HI"])
            row[f"ELCR_{pop}"] = float(out["ELCR"])
            row[f"ELCR_ED_AT_{pop}"] = float(re_.elcr(out["ADD_ing_mg_per_kg_day"], p.fixed["CSF_per_mg_per_kg_day"]))
        rows.append(row)
    return pd.DataFrame(rows)


def write_parquet(frame: pd.DataFrame, path: Path) -> None:
    frame.to_parquet(path, engine="pyarrow", compression="zstd", index=False)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args()

    cfg = load_run_config()
    sim = cfg["simulation"]
    for d in (RESULTS_SIMULATION, ITER_DIR, MANIFEST_DIR, RESULTS_TABLES, RESULTS_FIGURES):
        ensure_output_dir(d)
    upstream_before = {rel(p): sha256_file(p) for p in UPSTREAM}

    model = mc.load_model_inputs(cfg)
    _, child_review = si.bodyweight_inputs(cfg)
    records = {"arsenic": model.arsenic, "bodyweight": model.bodyweight}
    OUT["inputs_json"].write_text(json.dumps(records, indent=2) + "\n", encoding="utf-8")
    pd.DataFrame([{"input": "C_mg_per_L", "key": k, "distribution": v["distribution"],
                   "params": json.dumps(v["params"]), "decision": v["decision"]} for k, v in model.arsenic.items()]
                 + [{"input": "BW_kg", "key": k, "distribution": v["distribution"],
                     "params": json.dumps(v["params"]), "decision": v["decision"]}
                    for k, v in model.bodyweight.items()]).to_csv(OUT["inputs_csv"], index=False)
    arsenic_review(model).to_csv(OUT["arsenic_review"], index=False)
    child_review.to_csv(OUT["child_bw_review"], index=False)

    N = int(sim["primary_N"])
    summaries, checks, runs = [], [], []
    hi_by_run = {}
    for district in cfg["arsenic"]["selected_districts"]:
        for pop in sim["populations"]:
            spec = mc.run_spec(cfg, district, pop, N, replicate=0)
            frame = mc.simulate(cfg, model, spec)
            path = ITER_DIR / f"{spec.run_id}.parquet"
            write_parquet(frame, path)
            summaries += mc.summarize(frame, spec, tuple(sim["hi_thresholds"]), sim["elcr_threshold"])
            checks += mc.sampling_check(frame, model, spec)
            hi_by_run[(district, pop)] = frame["HI"].to_numpy()
            run_manifest = {
                "run_id": spec.run_id, "district": district, "population": pop, "N": N, "replicate": 0,
                "master_seed": sim["master_seed"], "spawn_key": list(spec.spawn_key),
                "seed_derivation": sim["seed_derivation"], "input_order": sim["input_order"],
                "arsenic_distribution": model.arsenic[district], "bw_distribution": model.bodyweight[pop],
                "base_paper_inputs": {k: {"family": v.family, "params": v.params, "scipy": v.scipy}
                                      for k, v in model.params[pop].stochastic.items()},
                "fixed": model.params[pop].fixed,
                "output": rel(path), "output_sha256": sha256_file(path),
                "config_version": cfg["config_version"], "config_sha256": sha256_file(RUN_CONFIG_PATH),
            }
            (MANIFEST_DIR / f"run_manifest_{spec.run_id}.json").write_text(
                json.dumps(run_manifest, indent=2) + "\n", encoding="utf-8")
            runs.append({"run_id": spec.run_id, "spawn_key": list(spec.spawn_key), "output": rel(path),
                         "sha256": run_manifest["output_sha256"]})

    # Reproducibility: re-simulate the first run with the same seed; it must be bit-identical.
    first = mc.run_spec(cfg, cfg["arsenic"]["selected_districts"][0], sim["populations"][0], N, 0)
    rerun = mc.simulate(cfg, model, first)
    saved = pd.read_parquet(ITER_DIR / f"{first.run_id}.parquet")
    reproducible = bool(rerun.equals(saved))

    summary = pd.DataFrame(summaries)
    summary.to_csv(OUT["summary"], index=False)
    sampling = pd.DataFrame(checks)
    sampling["ks_0.1pct_critical"] = KS_GATE / np.sqrt(sampling["N"])
    sampling["ks_within_0.1pct"] = sampling["ks_distance"] <= sampling["ks_0.1pct_critical"]
    sampling.to_csv(OUT["sampling"], index=False)

    s = summary.set_index(["district", "population", "metric"])
    districts = cfg["arsenic"]["selected_districts"]
    p95 = pd.DataFrame([{"district": d,
                         "P95_HI_adult": s.loc[(d, "adult", "HI"), "P95"],
                         "P95_ELCR_adult": s.loc[(d, "adult", "ELCR"), "P95"],
                         "P95_HI_child": s.loc[(d, "child", "HI"), "P95"],
                         "P95_ELCR_child": s.loc[(d, "child", "ELCR"), "P95"],
                         "P95_ELCR_child_ED_AT": s.loc[(d, "child", "ELCR_ED_AT"), "P95"]} for d in districts])
    p95.to_csv(OUT["p95"], index=False)
    exc_rows = []
    for d in districts:
        for pop in sim["populations"]:
            h, e = s.loc[(d, pop, "HI")], s.loc[(d, pop, "ELCR")]
            exc_rows.append({"district": d, "population": pop, "N": N,
                             "P_HI_gt_1": h["P_gt_1"], "P_HI_gt_1_mc_se": h["P_gt_1_mc_se"],
                             "P_HI_gt_1_wilson_lo": h["P_gt_1_wilson_lo"], "P_HI_gt_1_wilson_hi": h["P_gt_1_wilson_hi"],
                             "P_HI_gt_2": h["P_gt_2"], "P_HI_gt_2_mc_se": h["P_gt_2_mc_se"],
                             "P_ELCR_gt_1e-4": e["P_gt_0.0001"], "P_ELCR_gt_1e-4_mc_se": e["P_gt_0.0001_mc_se"]})
    pd.DataFrame(exc_rows).to_csv(OUT["exceedance"], index=False)

    bw_means = pd.read_csv(BW_SUMMARY).set_index("population")["weighted_mean_kg"].to_dict()
    det = deterministic_district(cfg, model, bw_means)
    det.to_csv(OUT["deterministic"], index=False)
    cmp_rows = []
    for _, r in det.iterrows():
        for pop in sim["populations"]:
            for metric in ("HI", "ELCR"):
                row = s.loc[(r.district, pop, metric)]
                cmp_rows.append({"district": r.district, "population": pop, "metric": metric,
                                 "deterministic": r[f"{metric}_{pop}"], "prob_mean": row["mean"],
                                 "prob_P50": row["P50"], "prob_P95": row["P95"],
                                 "P95_over_deterministic": row["P95"] / r[f"{metric}_{pop}"]})
    pd.DataFrame(cmp_rows).to_csv(OUT["det_vs_prob"], index=False)

    figures = []
    if not args.no_figures:
        for d in districts:
            figures.append(rel(sp.plot_district_hi(d, hi_by_run[(d, "adult")], hi_by_run[(d, "child")],
                                                   RESULTS_FIGURES)))

    upstream_after = {rel(p): sha256_file(p) for p in UPSTREAM}
    physical = bool(sampling["all_finite"].all() and sampling["all_physical"].all())
    ks_ok = bool(sampling["ks_within_0.1pct"].all())
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/run_simulation.py",
        "config_version": cfg["config_version"], "config_sha256": sha256_file(RUN_CONFIG_PATH),
        "upstream_inputs": upstream_before, "upstream_unchanged_during_run": upstream_before == upstream_after,
        "software": {"python": platform.python_version(),
                     **{name: importlib.metadata.version(name) for name in PACKAGES}},
        "master_seed": sim["master_seed"], "seed_derivation": sim["seed_derivation"],
        "N": N, "replicate": 0, "input_dependence": sim["input_dependence"],
        "input_selection": sim["input_selection"],
        "risk_model": {k: cfg["risk_model"][k] for k in ("plus_minus_interpretation", "cancer_averaging_time")},
        "checks": {"same_seed_rerun_identical": reproducible, "sampled_inputs_finite_and_physical": physical,
                   "ks_within_0.1pct_all": ks_ok,
                   "ks_within_5pct": f"{int(sampling['ks_within_critical'].sum())}/{len(sampling)}",
                   "distinct_spawn_keys": len({tuple(r['spawn_key']) for r in runs}) == len(runs)},
        "runs": runs,
        "outputs": {k: rel(v) for k, v in OUT.items() if k != "manifest"},
        "figures": figures,
    }
    OUT["manifest"].write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    with pd.option_context("display.width", 200, "display.float_format", "{:.4g}".format):
        print(p95.to_string(index=False))
        print(pd.DataFrame(exc_rows)[["district", "population", "P_HI_gt_1", "P_HI_gt_2", "P_ELCR_gt_1e-4"]]
              .to_string(index=False))
    print(f"checks: {manifest['checks']}")
    if not (reproducible and physical and ks_ok and manifest["upstream_unchanged_during_run"]):
        print("FAILED")
        return 1
    print(f"{len(runs)} runs; outputs in {rel(RESULTS_SIMULATION)} and {rel(RESULTS_TABLES)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
