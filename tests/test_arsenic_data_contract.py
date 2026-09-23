"""Data-contract tests for the Phase 1-2 outputs.

Run `python scripts/build_provenance.py` and `python scripts/preprocess_arsenic.py`
first; these tests check the saved files against the raw source.
"""
import json

import numpy as np
import pandas as pd
import pytest

from arsenic_hra import arsenic_preprocessing as ap
from arsenic_hra.paths import PROJECT_ROOT, RESULTS_ARSENIC, RESULTS_PROVENANCE, load_run_config, raw_path
from arsenic_hra.provenance import read_hash_manifest, verify_raw_hashes

CFG = load_run_config()
DISTRICTS = CFG["arsenic"]["selected_districts"]
SELECTED_CSV = RESULTS_ARSENIC / "arsenic_selected_districts_preprocessed.csv"
SUMMARY_CSV = RESULTS_ARSENIC / "arsenic_district_preprocessing_summary.csv"
AUDIT_CSV = RESULTS_ARSENIC / "arsenic_reconciliation_audit.csv"
HASHES = RESULTS_PROVENANCE / "raw_data.sha256"

pytestmark = pytest.mark.skipif(not SELECTED_CSV.exists(), reason="run scripts/preprocess_arsenic.py first")


@pytest.fixture(scope="module")
def selected():
    return pd.read_csv(SELECTED_CSV, dtype={"As_raw": str, "SAMPLE_ID": str, "DISTRICT": str,
                                            "GEOCODE": str}, keep_default_na=False,
                       na_values={c: [""] for c in ["As_bound_ugL", "As_detected_ugL", "As_model_value_ugL",
                                                     "As_bound_mgL", "As_detected_mgL", "As_model_value_mgL",
                                                     "LAT_DEG", "LONG_DEG", "WELL_DEPTH"]})


@pytest.fixture(scope="module")
def raw():
    return ap.read_bgs_raw(raw_path(CFG["raw_inputs"]["arsenic_bgs"]))[0]


# ---- raw integrity -------------------------------------------------------
def test_raw_files_unchanged():
    assert HASHES.exists(), "run scripts/build_provenance.py first"
    assert verify_raw_hashes(HASHES) == []


def test_every_analysis_input_is_frozen():
    frozen = set(read_hash_manifest(HASHES))
    for relative in CFG["raw_inputs"].values():
        assert relative in frozen, relative


# ---- schema --------------------------------------------------------------
def test_required_columns(selected):
    assert list(selected.columns) == ap.OUTPUT_COLUMNS


def test_sample_ids_unique(selected):
    assert selected["SAMPLE_ID"].is_unique


# ---- reconciliation ------------------------------------------------------
def test_audit_checks_all_pass():
    audit = pd.read_csv(AUDIT_CSV)
    checks = audit[audit["step"].str.startswith("check_")]
    assert len(checks) >= 8
    assert (checks["count"] == 1).all(), checks[checks["count"] != 1]


def test_row_counts_match_raw(selected, raw):
    expected = raw["DISTRICT"].isin(DISTRICTS).sum()
    assert len(selected) == expected
    counts = selected["DISTRICT"].value_counts()
    raw_counts = raw["DISTRICT"].value_counts()
    for d in DISTRICTS:
        assert counts[d] == raw_counts[d], d


def test_summary_reconciles_with_rows(selected):
    summary = pd.read_csv(SUMMARY_CSV).set_index("DISTRICT")
    assert list(summary.index) == DISTRICTS
    for d, grp in selected.groupby("DISTRICT"):
        s = summary.loc[d]
        assert s["total_n"] == len(grp)
        assert s["detected_n"] + s["censored_n"] + s["missing_n"] + s["malformed_n"] == s["total_n"]
        assert s["censored_n"] == (grp["As_status"] == ap.STATUS_CENSORED).sum()


def test_rajshahi_carries_high_censoring_warning():
    summary = pd.read_csv(SUMMARY_CSV).set_index("DISTRICT")
    assert bool(summary.loc["Rajshahi", "high_censoring_warning"])


# ---- value contracts -----------------------------------------------------
def test_as_raw_matches_source_string(selected, raw):
    source = raw.set_index("source_line")
    assert (selected["As_raw"].to_numpy() == source.loc[selected["source_line"], "As"].to_numpy()).all()
    assert (selected["SAMPLE_ID"].to_numpy() == source.loc[selected["source_line"], "SAMPLE_ID"].to_numpy()).all()


def test_district_names_exactly_match_source(selected):
    assert set(selected["DISTRICT"]) == set(DISTRICTS)


def test_no_censored_value_in_detected_field(selected):
    cens = selected["As_status"] == ap.STATUS_CENSORED
    assert selected.loc[cens, "As_detected_ugL"].isna().all()
    assert selected.loc[cens, "As_detected_mgL"].isna().all()
    assert selected.loc[cens, "As_raw"].str.strip().str.startswith("<").all()
    assert selected.loc[~cens, "As_bound_ugL"].isna().all()
    assert not selected.loc[~cens, "As_raw"].str.contains("<").any()


def test_censor_flag_consistent_with_status(selected):
    flag = selected["As_censored"].astype(str).str.lower()
    assert ((flag == "true") == (selected["As_status"] == ap.STATUS_CENSORED)).all()


def test_mgL_is_ugL_over_1000(selected):
    for stem in ("bound", "detected", "model_value"):
        ug, mg = selected[f"As_{stem}_ugL"], selected[f"As_{stem}_mgL"]
        assert (ug.isna() == mg.isna()).all()
        np.testing.assert_allclose(mg.dropna(), ug.dropna() / 1000.0, rtol=1e-12)


def test_model_values_positive_and_finite(selected):
    usable = selected["As_status"].isin([ap.STATUS_DETECTED, ap.STATUS_CENSORED])
    vals = selected.loc[usable, "As_model_value_mgL"]
    assert np.isfinite(vals).all() and (vals > 0).all()
    # Physical sanity for mg/L: the survey maximum is 1.66 mg/L (1660 ug/L).
    assert vals.max() < 2.0


def test_manifest_records_unchanged_raw():
    manifest = json.loads((RESULTS_ARSENIC / "arsenic_preprocessing_manifest.json").read_text())
    assert manifest["raw_unchanged_during_run"] is True
    assert manifest["failed_checks"] == []
    assert manifest["selected_districts"] == DISTRICTS
    assert (PROJECT_ROOT / manifest["raw_input"]).exists()
