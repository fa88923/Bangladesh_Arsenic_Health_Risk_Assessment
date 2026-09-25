"""Data-contract tests on the saved body-weight outputs (run scripts/run_bodyweight.py first)."""

import json

import numpy as np
import pandas as pd
import pytest

from arsenic_hra import bodyweight_fitting as ft
from arsenic_hra import bodyweight_preprocessing as bp
from arsenic_hra.paths import PROJECT_ROOT, RESULTS_BODYWEIGHT, RESULTS_PROVENANCE, load_run_config, raw_path
from arsenic_hra.provenance import read_hash_manifest, sha256_file

MANIFEST = RESULTS_BODYWEIGHT / "bodyweight_manifest.json"
CFG = load_run_config()
POPULATIONS = ("adult", "child")

pytestmark = pytest.mark.skipif(not MANIFEST.exists(), reason="run scripts/run_bodyweight.py first")


def load(name):
    return pd.read_csv(RESULTS_BODYWEIGHT / name)


@pytest.fixture(scope="module")
def manifest():
    return json.loads(MANIFEST.read_text(encoding="utf-8"))


@pytest.mark.parametrize("pop", POPULATIONS)
def test_audit_reconciles_with_clean_file(pop):
    clean = load(f"{pop}_bw_clean.csv")
    audit = load(f"{pop}_bw_cleaning_audit.csv")
    assert (audit["rows_before"] - audit["removed"] == audit["rows_after"]).all()
    assert (audit["rows_before"].iloc[1:].to_numpy() == audit["rows_after"].iloc[:-1].to_numpy()).all()
    assert audit["rows_after"].iloc[0] - audit["removed"].sum() == len(clean) == audit["rows_after"].iloc[-1]


@pytest.mark.parametrize("pop", POPULATIONS)
def test_clean_values_are_physical(pop):
    clean = load(f"{pop}_bw_clean.csv")
    bw, sw = clean["BW_kg"], clean["survey_weight"]
    assert np.isfinite(bw).all() and (bw > 0).all()
    assert np.isfinite(sw).all() and (sw > 0).all()
    codes = [666, 888] if pop == "adult" else [99.3, 99.4, 99.5, 99.6]
    assert not np.isclose(bw.to_numpy()[:, None], codes).any()


def test_adult_population_and_source_column():
    clean = load("adult_bw_clean.csv")
    assert clean["age"].between(18, 69).all()
    raw = pd.read_csv(raw_path(CFG["raw_inputs"]["adult_bw_steps"]), usecols=["pid", "m12", "wstep2"]).set_index("pid")
    joined = clean.set_index("pid").join(raw)
    assert np.allclose(joined["BW_kg"], joined["m12"])
    assert np.allclose(joined["survey_weight"], joined["wstep2"])


def test_child_population_and_flags():
    clean = load("child_bw_clean.csv")
    assert clean["CAGE_months"].between(0, 59).all()
    assert (clean["WAZFLAG"] == 0).all()
    assert (clean["BW_kg"] < 90).all()
    assert clean["height_cm"].dropna().lt(bp.MICS_AN11_SPECIAL_MIN).all()


def test_review_files_match_audit_counts():
    bmi = load("adult_bw_bmi_review.csv")
    adult_audit = load("adult_bw_cleaning_audit.csv").set_index("rule")
    assert len(bmi) == adult_audit.loc["bmi_outside_review_bounds", "raw_rows_matching"]
    assert bmi["retained_in_fit"].all()
    assert set(bmi["pid"]) <= set(load("adult_bw_clean.csv")["pid"])
    an8 = load("child_bw_an8_codebook_review.csv")
    child_audit = load("child_bw_cleaning_audit.csv").set_index("rule")
    assert len(an8) == child_audit.loc["AN8_unexpected_ge_90", "removed"]


def test_saved_selection_reproduces_saved_fit_metrics():
    selected = json.loads((RESULTS_BODYWEIGHT / "bodyweight_selected_distributions.json").read_text(encoding="utf-8"))
    fits = load("bodyweight_distribution_fit_results.csv")
    assert sorted(s["population"] for s in selected) == list(POPULATIONS)
    for s in selected:
        clean = load(f"{s['population']}_bw_clean.csv")
        ecdf = ft.weighted_ecdf_midpoints(clean["BW_kg"], clean["survey_weight"])
        dist = ft.frozen(s["distribution"], s["params"])
        assert ft.weighted_sse(dist, ecdf) == pytest.approx(s["weighted_cdf_sse"], rel=1e-9)
        assert ft.weighted_max_discrepancy(dist, ecdf) == pytest.approx(s["weighted_cdf_max_discrepancy"], rel=1e-9)
        row = fits[(fits["population"] == s["population"]) & fits["selected"]]
        assert len(row) == 1 and row["distribution"].iloc[0] == s["distribution"] and bool(row["admissible"].iloc[0])
        assert s["status"].startswith("provisional")
        draws = ft.sample_bw(s, 10_000, np.random.default_rng(0))
        assert np.isfinite(draws).all() and (draws > 0).all()


def test_every_fit_reports_convergence_state():
    fits = load("bodyweight_distribution_fit_results.csv")
    assert len(fits) == 8
    assert fits["ls_converged"].dtype == bool and fits["pml_converged"].dtype == bool
    assert fits["ls_message"].notna().all() and fits["pml_message"].notna().all()


def test_manifest_contract(manifest):
    assert manifest["config_version"] == CFG["config_version"]
    assert manifest["raw_unchanged_during_run"] and manifest["raw_matches_frozen_manifest"]
    assert all(manifest["audits_reconcile"].values())
    for path in list(manifest["outputs"].values()) + manifest["figures"]:
        assert (PROJECT_ROOT / path).is_file(), path
    assert len(manifest["figures"]) == 14
    assert all(p.startswith("results/figures/bodyweight_") for p in manifest["figures"])


def test_raw_inputs_match_frozen_hashes(manifest):
    frozen = read_hash_manifest(RESULTS_PROVENANCE / "raw_data.sha256")
    for rel_path, digest in manifest["raw_inputs"].items():
        assert frozen[rel_path] == digest, rel_path
        assert sha256_file(PROJECT_ROOT / rel_path) == digest, rel_path
