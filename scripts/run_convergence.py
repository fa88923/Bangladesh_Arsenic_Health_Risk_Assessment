"""Phase 8: convergence experiment for every district x population.

Runs the Phase 7 simulation at every N in simulation.convergence_N with
simulation.replicate_seeds_per_N independent seeds (400 runs), and writes:

    results/convergence/convergence_runs.csv          one row per run (design doc table 19.4)
    results/convergence/convergence_summary.csv       across-seed mean, SD, CV, range, successive change
    results/convergence/repeated_seed_stability.csv   the N = primary_N replicates side by side
    results/convergence/convergence_conclusions.csv   smallest adequate N per district/population/metric
    results/convergence/convergence_manifest.json
    results/figures/convergence_<metric>_vs_N.png, convergence_mc_error_decay.png

Iterations are not stored (decision D17): every run is reproducible from its seed,
and the primary N = 10,000 replicate-0 iterations are already saved by Phase 7.

Exits non-zero if the N = 10,000 replicate-0 run does not reproduce the Phase 7 primary run.

Usage:  python scripts/run_convergence.py [--no-figures]
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

from arsenic_hra import convergence as cv  # noqa: E402
from arsenic_hra import convergence_plots as cp  # noqa: E402
from arsenic_hra import monte_carlo as mc  # noqa: E402
from arsenic_hra.paths import (RESULTS_CONVERGENCE, RESULTS_FIGURES, RESULTS_SIMULATION, RUN_CONFIG_PATH,  # noqa: E402
                               ensure_output_dir, load_run_config, rel)
from arsenic_hra.provenance import sha256_file  # noqa: E402

OUT = {
    "runs": RESULTS_CONVERGENCE / "convergence_runs.csv",
    "summary": RESULTS_CONVERGENCE / "convergence_summary.csv",
    "seeds": RESULTS_CONVERGENCE / "repeated_seed_stability.csv",
    "conclusions": RESULTS_CONVERGENCE / "convergence_conclusions.csv",
    "manifest": RESULTS_CONVERGENCE / "convergence_manifest.json",
}
PACKAGES = ["numpy", "pandas", "scipy", "matplotlib"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args()

    cfg = load_run_config()
    sim = cfg["simulation"]
    ensure_output_dir(RESULTS_CONVERGENCE)
    ensure_output_dir(RESULTS_FIGURES)
    model = mc.load_model_inputs(cfg)
    districts = cfg["arsenic"]["selected_districts"]
    primary_N = int(sim["primary_N"])

    rows = []
    primary_matches = []
    for district in districts:
        for pop in sim["populations"]:
            for N in sim["convergence_N"]:
                for r in range(int(sim["replicate_seeds_per_N"])):
                    spec = mc.run_spec(cfg, district, pop, N, r)
                    frame = mc.simulate(cfg, model, spec)
                    rows.append({"district": district, "population": pop, "N": N, "replicate": r,
                                 "run_id": spec.run_id, "spawn_key": str(spec.spawn_key),
                                 **cv.run_statistics(frame)})
                    if N == primary_N and r == 0:
                        saved = pd.read_parquet(RESULTS_SIMULATION / "iterations" / f"{spec.run_id}.parquet")
                        primary_matches.append(bool(frame.equals(saved)))
    runs = pd.DataFrame(rows)
    runs.to_csv(OUT["runs"], index=False)

    agg = cv.aggregate(runs)
    agg.to_csv(OUT["summary"], index=False)
    seeds = runs[runs.N == primary_N].pivot_table(index=["district", "population"], columns="replicate",
                                                  values=list(cv.TRACKED))
    seeds.columns = [f"{m}_r{r}" for m, r in seeds.columns]
    seeds = seeds.reset_index()
    for m in cv.TRACKED:
        cols = [f"{m}_r{r}" for r in range(int(sim["replicate_seeds_per_N"]))]
        seeds[f"{m}_cv_pct"] = 100 * seeds[cols].std(axis=1, ddof=1) / seeds[cols].mean(axis=1).abs()
    seeds.to_csv(OUT["seeds"], index=False)
    concl = cv.conclusions(agg, primary_N)
    concl.to_csv(OUT["conclusions"], index=False)

    # Cross-check: observed between-seed SD of probabilities vs binomial MC SE.
    probs = agg[agg.metric.isin(cv.PROBABILITY_METRICS) & (agg.binomial_mc_se > 0)]
    ratio = (probs.between_seed_sd / probs.binomial_mc_se).describe()

    figures = []
    if not args.no_figures:
        figures = [rel(cp.plot_metric_vs_N(agg, m, districts, RESULTS_FIGURES)) for m in cv.TRACKED]
        figures.append(rel(cp.plot_error_decay(agg, RESULTS_FIGURES)))

    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/run_convergence.py",
        "config_version": cfg["config_version"], "config_sha256": sha256_file(RUN_CONFIG_PATH),
        "phase7_manifest_sha256": sha256_file(RESULTS_SIMULATION / "simulation_manifest.json"),
        "software": {"python": platform.python_version(),
                     **{name: importlib.metadata.version(name) for name in PACKAGES}},
        "master_seed": sim["master_seed"], "seed_derivation": sim["seed_derivation"],
        "N_grid": sim["convergence_N"], "replicates_per_N": sim["replicate_seeds_per_N"],
        "runs": len(runs), "iterations_total": int(runs.N.sum()),
        "criterion": {"continuous": f"successive change <= {cv.REL_TOL_PCT}% and between-seed CV <= {cv.REL_TOL_PCT}%",
                      "probabilities": f"absolute change <= {cv.ABS_TOL_PROB} and between-seed SD <= {cv.ABS_TOL_PROB}"},
        "checks": {"N10000_r0_reproduces_phase7_primary": all(primary_matches) and len(primary_matches) == 20,
                   "between_seed_sd_over_binomial_se": {k: float(v) for k, v in ratio.items()}},
        "adequate_at_primary_N": f"{int(concl.adequate_at_primary_N.sum())}/{len(concl)}",
        "outputs": {k: rel(v) for k, v in OUT.items() if k != "manifest"},
        "figures": figures,
    }
    OUT["manifest"].write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(concl.pivot_table(index=["district", "population"], columns="metric", values="smallest_adequate_N",
                            aggfunc="first").to_string())
    print(f"adequate at N={primary_N}: {manifest['adequate_at_primary_N']}; "
          f"between-seed SD / binomial SE: median {ratio['50%']:.2f}")
    if not manifest["checks"]["N10000_r0_reproduces_phase7_primary"]:
        print("FAILED: N=10000 r0 does not reproduce the Phase 7 primary run")
        return 1
    print(f"{len(runs)} runs; outputs in {rel(RESULTS_CONVERGENCE)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
