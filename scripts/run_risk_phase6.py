"""Phase 6: risk-equation specification, parameter table, and deterministic benchmarks.

Writes, under results/tables/:

    risk_equation_specification.csv     equations (1)-(6) with input and output units
    risk_parameter_table.csv            every fixed and sampled input, both populations
    risk_paper_table4_reproduction.csv  Yadav & Kalkal (2024) Table 4 recomputed from its own inputs
    risk_deterministic_benchmark.csv    Bangladesh point inputs at reference concentrations
    risk_phase6_manifest.json

No random numbers are drawn. Exits non-zero if the Table 4 reproduction deviates
from the published values by more than 0.5%, other than the three documented cells.

Usage:  python scripts/run_risk_phase6.py
"""
from __future__ import annotations

import copy
import importlib.metadata
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd  # noqa: E402

from arsenic_hra import risk_equations as re_  # noqa: E402
from arsenic_hra import risk_parameters as rp  # noqa: E402
from arsenic_hra.paths import (RESULTS_BODYWEIGHT, RESULTS_TABLES, RUN_CONFIG_PATH,  # noqa: E402
                               ensure_output_dir, load_run_config, rel)
from arsenic_hra.provenance import sha256_file  # noqa: E402

OUT = {
    "equations": RESULTS_TABLES / "risk_equation_specification.csv",
    "parameters": RESULTS_TABLES / "risk_parameter_table.csv",
    "paper_table4": RESULTS_TABLES / "risk_paper_table4_reproduction.csv",
    "benchmark": RESULTS_TABLES / "risk_deterministic_benchmark.csv",
    "manifest": RESULTS_TABLES / "risk_phase6_manifest.json",
}
BW_SUMMARY = RESULTS_BODYWEIGHT / "bodyweight_weighted_summary.csv"
PACKAGES = ["numpy", "pandas", "scipy"]

# Yadav & Kalkal (2024) Table 1 means (mg/L) and Table 4 results.
PAPER_TABLE4 = {
    "Moga": (0.0043, 0.41, 1.72, 1.84e-4, 7.74e-4),
    "Faridkot": (0.017, 1.62, 6.83, 7.29e-4, 3.06e-3),
    "Fazilka": (0.037, 3.54, 14.87, 1.59e-3, 6.66e-3),
    "Patiala": (0.021, 2.01, 8.44, 9.01e-4, 3.78e-3),
    "Ferozepur": (0.027, 2.58, 10.85, 1.16e-3, 4.80e-3),
    "Rupnagar": (0.035, 3.35, 14.07, 1.50e-3, 6.30e-3),
    "Amritsar": (0.083, 7.94, 33.26, 3.50e-3, 1.40e-2),
}
PAPER_BW_KG = {"adult": 70.0, "child": 15.0}
KNOWN_DISCREPANCIES = {("Ferozepur", "child", "ELCR"), ("Amritsar", "adult", "ELCR"), ("Amritsar", "child", "ELCR")}
REPRODUCTION_RTOL = 0.005

# Reference concentrations for the Bangladesh benchmark (mg/L).
REFERENCE_C = {"WHO guideline": 0.01, "Bangladesh standard": 0.05}


def paper_table4_reproduction(cfg: dict) -> pd.DataFrame:
    paper_cfg = copy.deepcopy(cfg)
    paper_cfg["risk_model"]["cancer_averaging_time"] = "exposure_duration"
    rows = []
    for district, (C, hi_a, hi_c, elcr_a, elcr_c) in PAPER_TABLE4.items():
        refs = {"adult": {"HI": hi_a, "ELCR": elcr_a}, "child": {"HI": hi_c, "ELCR": elcr_c}}
        for pop in rp.POPULATIONS:
            p = rp.population_parameters(pop, paper_cfg)
            out = re_.evaluate_risk(re_.deterministic_inputs(p, C, PAPER_BW_KG[pop]), p)
            for metric in ("HI", "ELCR"):
                ours, paper = float(out[metric]), refs[pop][metric]
                rel_dev = (ours - paper) / paper
                known = (district, pop, metric) in KNOWN_DISCREPANCIES
                rows.append({"district": district, "population": pop, "metric": metric, "C_mg_per_L": C,
                             "BW_kg": PAPER_BW_KG[pop], "recomputed": ours, "paper_value": paper,
                             "relative_deviation": rel_dev, "known_paper_discrepancy": known,
                             "within_rounding": abs(rel_dev) <= REPRODUCTION_RTOL})
    return pd.DataFrame(rows)


def bangladesh_benchmark(cfg: dict, bw_means: dict) -> pd.DataFrame:
    rows = []
    for label, C in REFERENCE_C.items():
        for pop in rp.POPULATIONS:
            p = rp.population_parameters(pop, cfg)
            x = re_.deterministic_inputs(p, C, bw_means[pop])
            out = re_.evaluate_risk(x, p)
            rows.append({"reference": label, "population": pop, **x, **p.fixed,
                         **{k: float(v) for k, v in out.items()}})
    return pd.DataFrame(rows)


def main() -> int:
    cfg = load_run_config()
    ensure_output_dir(RESULTS_TABLES)

    pd.DataFrame(re_.EQUATIONS).to_csv(OUT["equations"], index=False)
    pd.DataFrame(rp.parameter_table(cfg)).to_csv(OUT["parameters"], index=False)

    repro = paper_table4_reproduction(cfg)
    repro.to_csv(OUT["paper_table4"], index=False)
    unexplained = repro[~repro.within_rounding & ~repro.known_paper_discrepancy]
    stale_known = repro[repro.within_rounding & repro.known_paper_discrepancy]

    bw = pd.read_csv(BW_SUMMARY).set_index("population")["weighted_mean_kg"].to_dict()
    bench = bangladesh_benchmark(cfg, bw)
    bench.to_csv(OUT["benchmark"], index=False)

    rcfg = cfg["risk_model"]
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/run_risk_phase6.py",
        "config_version": cfg["config_version"],
        "config_sha256": sha256_file(RUN_CONFIG_PATH),
        "inputs": {rel(BW_SUMMARY): sha256_file(BW_SUMMARY)},
        "software": {"python": platform.python_version(),
                     **{name: importlib.metadata.version(name) for name in PACKAGES}},
        "risk_model": {
            "plus_minus_interpretation": rcfg["plus_minus_interpretation"],
            "cancer_averaging_time": rcfg["cancer_averaging_time"],
            "base_paper": rcfg["base_paper"],
        },
        "deterministic_BW_kg": {"source": "Phase 5 survey-weighted mean", **bw},
        "paper_table4": {"cells": len(repro), "within_rounding": int(repro.within_rounding.sum()),
                         "known_paper_discrepancies": sorted("/".join(k) for k in KNOWN_DISCREPANCIES),
                         "unexplained_deviations": len(unexplained), "rtol": REPRODUCTION_RTOL},
        "random_seed": None,
        "outputs": {k: rel(v) for k, v in OUT.items() if k != "manifest"},
    }
    OUT["manifest"].write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    cols = ["reference", "population", "C_mg_per_L", "BW_kg", "HQ_ing", "HQ_dermal", "HI", "ELCR"]
    print(bench[cols].to_string(index=False))
    print(f"paper Table 4: {int(repro.within_rounding.sum())}/{len(repro)} cells within "
          f"{REPRODUCTION_RTOL:.1%}; {len(KNOWN_DISCREPANCIES)} documented paper discrepancies")
    if len(unexplained) or len(stale_known):
        print("FAILED: unexplained Table 4 deviations:\n" + unexplained.to_string(index=False)
              + ("\nknown discrepancies that now reproduce:\n" + stale_known.to_string(index=False)))
        return 1
    print(f"outputs in {rel(RESULTS_TABLES)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
