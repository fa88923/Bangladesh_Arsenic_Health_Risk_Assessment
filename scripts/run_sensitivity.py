"""Phase 9: sensitivity and uncertainty analysis.

1. Spearman rank sensitivity of HI (and ELCR) to the six sampled inputs, computed on the exact
   Phase 7 primary-run iterations; tornado plots per district.
2. Fitted-parameter uncertainty: 2-D Monte Carlo with the Phase 4 arsenic bootstrap (regenerated
   bit-for-bit) and a Rao-Wu survey-PSU bootstrap of BW; common random numbers with the primary run.
3. Model scenarios (common random numbers): alternative arsenic families, child-BW Gamma,
   the 0-71-month child-BW age extension, and the old log-space reading of IR/SA.

Writes results/sensitivity/* and results/figures/sensitivity_*.png.
Exits non-zero if a Spearman sign contradicts the equations, the regenerated arsenic bootstrap
does not reproduce Member 2's saved summary, or an outer replicate fails.

Usage:  python scripts/run_sensitivity.py [--reuse-bootstrap] [--no-figures] [--jobs J]
"""
from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
import platform
import sys
from datetime import datetime, timezone
from multiprocessing import get_context
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from arsenic_hra import monte_carlo as mc  # noqa: E402
from arsenic_hra import parameter_uncertainty as pu  # noqa: E402
from arsenic_hra import sensitivity as sn  # noqa: E402
from arsenic_hra import sensitivity_plots as spl  # noqa: E402
from arsenic_hra.paths import (RESULTS_FIGURES, RESULTS_SENSITIVITY, RESULTS_SIMULATION, RUN_CONFIG_PATH,  # noqa: E402
                               ensure_output_dir, load_run_config, rel)
from arsenic_hra.provenance import sha256_file  # noqa: E402

OUT = {
    "spearman": RESULTS_SENSITIVITY / "sensitivity_spearman.csv",
    "spearman_table": RESULTS_SENSITIVITY / "sensitivity_spearman_HI_table.csv",
    "arsenic_boot": RESULTS_SENSITIVITY / "arsenic_bootstrap_replicates.csv",
    "arsenic_boot_check": RESULTS_SENSITIVITY / "arsenic_bootstrap_reproduction_check.csv",
    "bw_boot": RESULTS_SENSITIVITY / "bw_survey_bootstrap_replicates.csv",
    "two_d": RESULTS_SENSITIVITY / "parameter_uncertainty_replicates.csv",
    "two_d_summary": RESULTS_SENSITIVITY / "parameter_uncertainty_summary.csv",
    "scenarios": RESULTS_SENSITIVITY / "scenario_results.csv",
    "manifest": RESULTS_SENSITIVITY / "sensitivity_manifest.json",
}
PACKAGES = ["numpy", "pandas", "scipy", "matplotlib", "lifelines"]
# Member 2 ran Phase 4 under numpy 2.5.3 / scipy 1.18.1 (arsenic_fitting_manifest.json); this
# environment pins numpy 1.26.4 / scipy 1.17.1. Optimizer end points then differ at ~1e-6 relative,
# so continuous summaries are compared at 1e-5 and the discrete winner frequencies exactly.
BOOT_REPRO_RTOL = 1e-5
WINNER_COLUMNS_RTOL = 1e-9

_STATE: dict = {}


def _two_d_worker(cell):
    d, pop = cell
    return sn.two_dimensional(_STATE["cfg"], _STATE["model"], d, pop, _STATE["arsenic_sets"][d], _STATE["bw_sets"])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reuse-bootstrap", action="store_true", help="load saved bootstrap tables if present")
    ap.add_argument("--no-figures", action="store_true")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    args = ap.parse_args()

    cfg = load_run_config()
    sim = cfg["simulation"]
    ensure_output_dir(RESULTS_SENSITIVITY)
    ensure_output_dir(RESULTS_FIGURES)
    model = mc.load_model_inputs(cfg)
    districts = cfg["arsenic"]["selected_districts"]
    pops = sim["populations"]
    cells = [(d, p) for d in districts for p in pops]
    iter_dir = RESULTS_SIMULATION / "iterations"
    figures = []

    # 1. Spearman sensitivity on the exact primary iterations ------------------------------
    rows, primary_stats = [], []
    for d, pop in cells:
        frame = pd.read_parquet(iter_dir / f"{d.lower()}_{pop}_N{sim['primary_N']}_r0.parquet")
        rows += sn.spearman_rows(frame, d, pop)
        primary_stats.append({"district": d, "population": pop, **sn.risk_statistics(frame)})
    spear = sn.rank_inputs(pd.DataFrame(rows))
    spear.to_csv(OUT["spearman"], index=False)
    wide = spear[spear.output == "HI"].pivot_table(index=["district", "population"], columns="input",
                                                   values="spearman_rho").reindex(columns=list(sn.INPUTS))
    wide.columns = [f"rho_{c}" for c in wide.columns]
    wide.reset_index().to_csv(OUT["spearman_table"], index=False)
    primary = pd.DataFrame(primary_stats)
    signs_ok = bool(spear.sign_consistent.all())
    if not args.no_figures:
        figures += [rel(spl.plot_tornado(spear, d, RESULTS_FIGURES)) for d in districts]

    # 2. Parameter uncertainty -----------------------------------------------------------------
    if args.reuse_bootstrap and OUT["arsenic_boot"].exists():
        boot = pd.read_csv(OUT["arsenic_boot"])
    else:
        boot = pu.regenerate_arsenic_bootstrap(cfg, jobs=args.jobs)
        boot.to_csv(OUT["arsenic_boot"], index=False)
    check = pu.compare_with_phase4(boot)
    check.to_csv(OUT["arsenic_boot_check"], index=False)
    winners = check.column.str.contains("winner")
    boot_reproduced = bool((check.max_rel_diff <= BOOT_REPRO_RTOL).all()
                           and (check.loc[winners, "max_rel_diff"] <= WINNER_COLUMNS_RTOL).all()
                           and (check.nan_pattern_mismatches == 0).all())
    arsenic_sets = pu.arsenic_parameter_sets(boot, model.arsenic)

    B = int(cfg["arsenic"]["bootstrap_replicates"])
    if args.reuse_bootstrap and OUT["bw_boot"].exists():
        bw_sets = pd.read_csv(OUT["bw_boot"])
    else:
        bw_sets = pu.bw_parameter_sets(cfg, model.bodyweight, B, jobs=args.jobs)
        bw_sets.to_csv(OUT["bw_boot"], index=False)
    bw_converged = bool(bw_sets.converged.all())

    _STATE.update(cfg=cfg, model=model, arsenic_sets=arsenic_sets, bw_sets=bw_sets)
    with get_context("fork").Pool(min(args.jobs, len(cells))) as pool:
        two_d = pd.concat(pool.map(_two_d_worker, cells), ignore_index=True)
    two_d.to_csv(OUT["two_d"], index=False)
    unc = sn.uncertainty_summary(two_d, primary)
    unc.to_csv(OUT["two_d_summary"], index=False)
    if not args.no_figures:
        figures.append(rel(spl.plot_parameter_uncertainty(unc, "P95_HI", "P95 HI", districts, RESULTS_FIGURES)))
        figures.append(rel(spl.plot_parameter_uncertainty(unc, "P_HI_gt_1", "P(HI > 1)", districts,
                                                          RESULTS_FIGURES, log=False)))
        figures.append(rel(spl.plot_parameter_uncertainty(unc, "P95_ELCR", "P95 ELCR (lifetime AT)", districts,
                                                          RESULTS_FIGURES)))

    # 3. Scenarios -------------------------------------------------------------------------------
    age_rec, age_info = sn.child_bw_age_extended(model.bodyweight["child"]["params"]["lower"])
    scenarios = {
        "S1_barisal_lognormal": (sn.scenario_model(model, arsenic={"Barisal": sn.arsenic_alternative("Barisal", "Lognormal")}),
                                 [("Barisal", p) for p in pops]),
        "S2_bogra_gamma": (sn.scenario_model(model, arsenic={"Bogra": sn.arsenic_alternative("Bogra", "Gamma")}),
                           [("Bogra", p) for p in pops]),
        "S3_rangpur_lognormal": (sn.scenario_model(model, arsenic={"Rangpur": sn.arsenic_alternative("Rangpur", "Lognormal")}),
                                 [("Rangpur", p) for p in pops]),
        "S4_child_bw_gamma": (sn.scenario_model(model, bodyweight={"child": sn.child_bw_gamma()}),
                              [(d, "child") for d in districts]),
        "S5_child_bw_0_71_months": (sn.scenario_model(model, bodyweight={"child": age_rec}),
                                    [(d, "child") for d in districts]),
        "S6_ir_sa_log_space": (sn.scenario_model(model, params=sn.log_space_params(cfg)), cells),
    }
    prim = primary.set_index(["district", "population"])
    scen_rows = []
    for name, (smodel, scells) in scenarios.items():
        for r in sn.run_scenario(cfg, smodel, scells):
            base = prim.loc[(r["district"], r["population"])]
            scen_rows.append({"scenario": name, **r,
                              **{f"ratio_{k}": (r[k] / base[k]) if base[k] else np.nan for k in sn.STAT_NAMES},
                              **{f"primary_{k}": base[k] for k in sn.STAT_NAMES}})
    scen = pd.DataFrame(scen_rows)
    scen.to_csv(OUT["scenarios"], index=False)

    outer_ok = bool(two_d.groupby(["district", "population"]).size().min() >= 0.95 * B)
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/run_sensitivity.py",
        "config_version": cfg["config_version"], "config_sha256": sha256_file(RUN_CONFIG_PATH),
        "phase7_manifest_sha256": sha256_file(RESULTS_SIMULATION / "simulation_manifest.json"),
        "software": {"python": platform.python_version(),
                     **{name: importlib.metadata.version(name) for name in PACKAGES}},
        "master_seed": sim["master_seed"],
        "spearman": {"iterations": "Phase 7 primary runs (N = 10,000, replicate 0)", "method": "scipy.stats.spearmanr (tie-aware)",
                     "inputs": list(sn.INPUTS), "signs_consistent_with_equations": signs_ok},
        "parameter_uncertainty": {
            "arsenic": "Phase 4 bootstrap regenerated with Member 2's code and seeds; final family refitted per replicate",
            "arsenic_bootstrap_reproduces_phase4_summary": boot_reproduced,
            "bw": "Rao-Wu (n_h - 1) PSU-within-stratum bootstrap of survey weights; final family refitted",
            "bw_stream": f"SeedSequence(master_seed, spawn_key=({pu.BW_STREAM_TAG}, population_index, replicate))",
            "bw_replicates": B, "bw_all_converged": bw_converged,
            "outer_replicates_per_cell": two_d.groupby(["district", "population"]).size().to_dict().__repr__(),
            "inner": "common random numbers: the primary run's uniforms reused for every outer replicate",
            "held_fixed": ["IR_L_per_day", "EF_days_per_year", "ET_h_per_day", "SA_m2"],
        },
        "scenarios": list(scenarios), "age_extension": age_info,
        "checks": {"spearman_signs": signs_ok, "arsenic_bootstrap_reproduced": boot_reproduced,
                   "bw_bootstrap_converged": bw_converged, "outer_replicates_sufficient": outer_ok},
        "outputs": {k: rel(v) for k, v in OUT.items() if k != "manifest"},
        "figures": figures,
    }
    OUT["manifest"].write_text(json.dumps(manifest, indent=2, default=float) + "\n", encoding="utf-8")

    with pd.option_context("display.width", 220, "display.float_format", "{:.3f}".format):
        print(wide.round(3).to_string())
        print(unc[unc.statistic.isin(["P95_HI", "P_HI_gt_1"])]
              [["district", "population", "statistic", "primary_value", "param_unc_q025", "param_unc_q975"]].to_string(index=False))
        print(scen[["scenario", "district", "population", "P95_HI", "ratio_P95_HI", "P_HI_gt_1", "ratio_P_HI_gt_1"]]
              .to_string(index=False))
    print(f"age extension: {age_info}")
    print(f"checks: {manifest['checks']}")
    if not all(manifest["checks"].values()):
        print("FAILED")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
