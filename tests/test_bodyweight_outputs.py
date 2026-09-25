"""Data-contract tests on the saved results/bodyweight outputs (run scripts/run_bodyweight.py first)."""

import hashlib
import json

import numpy as np
import pandas as pd
import pytest

from bodyweight import fitting as ft
from bodyweight.paths import PROJECT_ROOT, RESULTS_DIR, STEPS_PATH

pytestmark = pytest.mark.skipif(not (RESULTS_DIR / "bodyweight_run_manifest.json").exists(),
                                reason="body-weight outputs not generated")

POPULATIONS = ("adult", "child")


def load(name):
    return pd.read_csv(RESULTS_DIR / name)


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
    raw = pd.read_csv(STEPS_PATH, usecols=["pid", "m12", "wstep2"]).set_index("pid")
    joined = clean.set_index("pid").join(raw)
    assert np.allclose(joined["BW_kg"], joined["m12"])
    assert np.allclose(joined["survey_weight"], joined["wstep2"])


def test_child_population_and_flags():
    clean = load("child_bw_clean.csv")
    assert clean["CAGE_months"].between(0, 59).all()
    assert (clean["WAZFLAG"] == 0).all()
    assert (clean["BW_kg"] < 90).all()


def test_saved_selection_reproduces_saved_fit_metrics():
    selected = json.loads((RESULTS_DIR / "bodyweight_selected_distributions.json").read_text())
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


def test_raw_inputs_unchanged_since_run():
    manifest = json.loads((RESULTS_DIR / "bodyweight_run_manifest.json").read_text())
    for rel_path, digest in manifest["inputs_sha256"].items():
        h = hashlib.sha256((PROJECT_ROOT / rel_path).read_bytes()).hexdigest()
        assert h == digest, rel_path
