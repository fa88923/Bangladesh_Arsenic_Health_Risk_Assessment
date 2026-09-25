"""Phase 2: preprocess the BGS arsenic data for the approved districts.

Reads data/raw/NationalSurveyData.csv (never modified) and writes, under
results/arsenic/:

    arsenic_all_districts_classified.csv      every sample row, classified (QC basis)
    arsenic_selected_districts_preprocessed.csv
    arsenic_district_preprocessing_summary.csv
    arsenic_duplicate_review.csv
    arsenic_malformed_value_review.csv
    arsenic_reconciliation_audit.csv
    arsenic_preprocessing_manifest.json

Exits non-zero if any reconciliation check fails.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from arsenic_hra import arsenic_preprocessing as ap  # noqa: E402
from arsenic_hra.paths import RESULTS_ARSENIC, ensure_output_dir, load_run_config, raw_path, rel  # noqa: E402
from arsenic_hra.provenance import sha256_file  # noqa: E402

OUT = {
    "all": RESULTS_ARSENIC / "arsenic_all_districts_classified.csv",
    "selected": RESULTS_ARSENIC / "arsenic_selected_districts_preprocessed.csv",
    "summary": RESULTS_ARSENIC / "arsenic_district_preprocessing_summary.csv",
    "duplicates": RESULTS_ARSENIC / "arsenic_duplicate_review.csv",
    "malformed": RESULTS_ARSENIC / "arsenic_malformed_value_review.csv",
    "audit": RESULTS_ARSENIC / "arsenic_reconciliation_audit.csv",
    "manifest": RESULTS_ARSENIC / "arsenic_preprocessing_manifest.json",
}


def main() -> int:
    cfg = load_run_config()
    acfg = cfg["arsenic"]
    districts = acfg["selected_districts"]
    src = raw_path(cfg["raw_inputs"]["arsenic_bgs"])
    hash_before = sha256_file(src)

    raw, units = ap.read_bgs_raw(src)
    if units.get("As", "").lower() != acfg["source_unit"].lower():
        raise ValueError(f"As unit is {units.get('As')!r}, config expects {acfg['source_unit']!r}")

    classified = ap.add_typed_metadata(ap.classify_arsenic(raw))
    classified = classified[ap.OUTPUT_COLUMNS]
    selected = ap.filter_districts(classified, districts)

    summary = ap.district_summary(selected, districts, acfg["high_censoring_warning_percent"])
    duplicates = ap.duplicate_review(classified)
    duplicates["in_selected_district"] = duplicates["DISTRICT"].isin(districts)
    malformed = ap.malformed_review(classified)
    malformed["in_selected_district"] = malformed["DISTRICT"].isin(districts)
    audit = ap.reconciliation_audit(raw, classified, selected, districts, ap.count_lines(src))

    ensure_output_dir(RESULTS_ARSENIC)
    float_fmt = "%.10g"
    classified.to_csv(OUT["all"], index=False, float_format=float_fmt)
    selected.to_csv(OUT["selected"], index=False, float_format=float_fmt)
    summary.to_csv(OUT["summary"], index=False, float_format=float_fmt)
    duplicates.to_csv(OUT["duplicates"], index=False, float_format=float_fmt)
    malformed.to_csv(OUT["malformed"], index=False, float_format=float_fmt)
    audit.to_csv(OUT["audit"], index=False)

    hash_after = sha256_file(src)
    failed = audit.loc[audit["step"].str.startswith("check_") & (audit["count"] != 1), "step"].tolist()
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/preprocess_arsenic.py",
        "config_version": cfg["config_version"],
        "raw_input": rel(src),
        "raw_sha256": hash_before,
        "raw_unchanged_during_run": hash_before == hash_after,
        "source_units": units,
        "conversion": "mg/L = ug/L / 1000 (detected values and censoring limits)",
        "selected_districts": districts,
        "rows": {"raw_samples": len(raw), "selected": len(selected)},
        "failed_checks": failed,
        "outputs": {k: rel(v) for k, v in OUT.items() if k != "manifest"},
    }
    OUT["manifest"].write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(summary[["DISTRICT", "total_n", "detected_n", "censored_n", "censoring_percent",
                   "high_censoring_warning"]].to_string(index=False))
    print(f"\n{len(raw)} raw sample rows -> {len(selected)} selected rows")
    print(f"duplicate-review rows: {len(duplicates)}; malformed-review rows: {len(malformed)}")
    if failed or hash_before != hash_after:
        print(f"RECONCILIATION FAILED: {failed}; raw unchanged: {hash_before == hash_after}")
        return 1
    print("all reconciliation checks passed; raw file unchanged")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
