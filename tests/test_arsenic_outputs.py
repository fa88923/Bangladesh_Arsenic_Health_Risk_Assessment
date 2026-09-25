"""Data-contract tests on the saved arsenic fitting outputs.

Run ``python scripts/run_arsenic_fitting.py`` first. The whole module is skipped when
the manifest is absent, mirroring ``test_bodyweight_outputs.py``.
"""

import json

import numpy as np
import pandas as pd
import pytest

from arsenic_hra import arsenic_bootstrap as ab
from arsenic_hra import arsenic_fitting as af
from arsenic_hra.paths import PROJECT_ROOT, RESULTS_ARSENIC, RESULTS_PROVENANCE, load_run_config, raw_path
from arsenic_hra.provenance import read_hash_manifest, sha256_file

MANIFEST = RESULTS_ARSENIC / "arsenic_fitting_manifest.json"
CFG = load_run_config()
ACFG = CFG["arsenic"]
DISTRICTS = list(ACFG["selected_districts"])

pytestmark = pytest.mark.skipif(not MANIFEST.exists(),
                                reason="run scripts/run_arsenic_fitting.py first")


def load(name):
    return pd.read_csv(RESULTS_ARSENIC / name)


@pytest.fixture(scope="module")
def manifest():
    return json.loads(MANIFEST.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def selections():
    return json.loads((RESULTS_ARSENIC / "arsenic_selected_distributions.json").read_text(encoding="utf-8"))


def _district_rows(district):
    data = load("arsenic_selected_districts_preprocessed.csv")
    part = data[data["DISTRICT"] == district]
    return part["As_model_value_mgL"].to_numpy(dtype=float), part["As_censored"].to_numpy(dtype=bool)


# ---------------------------------------------------------------------------
# 1-2. Table shapes
# ---------------------------------------------------------------------------

def test_fit_table_has_thirty_rows_covering_every_district_and_family():
    """Plan section 6 test 1: 10 districts x 3 families."""
    fits = load("arsenic_distribution_fit_results.csv")
    assert len(fits) == len(DISTRICTS) * len(af.FAMILIES) == 30
    assert set(fits["DISTRICT"]) == set(DISTRICTS)
    assert set(fits["distribution"]) == set(af.FAMILIES)
    assert not fits.duplicated(["DISTRICT", "distribution"]).any()
    combos = fits.groupby("DISTRICT")["distribution"].apply(set)
    assert all(c == set(af.FAMILIES) for c in combos)


def test_exactly_ten_selected_rows_one_per_approved_district():
    """Plan section 6 test 2."""
    selected = load("arsenic_selected_distributions.csv")
    assert len(selected) == len(DISTRICTS) == 10
    assert set(selected["DISTRICT"]) == set(DISTRICTS)
    assert not selected["DISTRICT"].duplicated().any()
    assert set(selected["distribution"]) <= set(af.FAMILIES)


# ---------------------------------------------------------------------------
# 3-4. Admissibility and round-trip
# ---------------------------------------------------------------------------

def test_every_selected_family_passes_admissibility():
    """Plan section 6 test 3."""
    fits = load("arsenic_distribution_fit_results.csv")
    selected = load("arsenic_selected_distributions.csv")
    for row in selected.itertuples():
        match = fits[(fits["DISTRICT"] == row.DISTRICT) & (fits["distribution"] == row.distribution)]
        assert len(match) == 1
        assert bool(match["admissible"].iloc[0]) is True
        assert match["inadmissible_reason"].iloc[0] == "" or pd.isna(match["inadmissible_reason"].iloc[0])


def test_normal_is_rejected_in_every_district_by_the_gate():
    """The verified consequence in plan.md section 3.E: Normal wins raw SSE but is
    physically impossible, so the real contest is Lognormal vs Gamma."""
    fits = load("arsenic_distribution_fit_results.csv")
    normal = fits[fits["distribution"] == "Normal"]
    assert (~normal["admissible"]).all()
    assert (normal["prob_arsenic_nonpositive"] >= af.NORMAL_MAX_NONPOSITIVE_PROB).all()
    assert normal["inadmissible_reason"].str.contains("P\\(C<=0\\)").all()


def test_selected_parameters_reproduce_the_stored_metrics():
    """Plan section 6 test 4: the saved parameters regenerate the saved curve metrics.

    The CSVs are written with the project's ``%.10g`` float format, so the round-trip
    tolerance is 1e-9 relative (the plan's figure) with a 1e-9 absolute floor for the
    smallest quantities; both are far tighter than any modelling uncertainty.
    """
    selected = load("arsenic_selected_distributions.csv")
    fits = load("arsenic_distribution_fit_results.csv").set_index(["DISTRICT", "distribution"])
    for row in selected.itertuples():
        x, censored = _district_rows(row.DISTRICT)
        curve = af.reverse_km_cdf(x, censored)
        mask = af.objective_mask(curve, row.objective_region)
        params = {name: getattr(row, f"param_{name}") for name in af.PARAM_NAMES[row.distribution]}
        dist = af.frozen(row.distribution, params)
        sse = af.cdf_sse(dist, curve.x[mask], curve.F[mask])
        assert sse == pytest.approx(row.cdf_sse, rel=1e-9), row.DISTRICT
        assert af.cdf_rmse_from_sse(sse, int(mask.sum())) == pytest.approx(row.cdf_rmse, rel=1e-9)
        assert af.cdf_max_discrepancy(dist, curve.x[mask], curve.F[mask]) == pytest.approx(
            row.cdf_max_discrepancy, rel=1e-8)
        # The same parameters must also be recorded in the fit table.
        fit_row = fits.loc[(row.DISTRICT, row.distribution)]
        for name in af.PARAM_NAMES[row.distribution]:
            assert fit_row[f"ls_{name}"] == pytest.approx(params[name], rel=1e-9)
        assert int(row.objective_points) == int(mask.sum())


def test_selected_rows_are_the_lowest_rmse_admissible_family():
    """The predeclared selection rule must actually have been applied."""
    fits = load("arsenic_distribution_fit_results.csv")
    selected = load("arsenic_selected_distributions.csv")
    for row in selected.itertuples():
        part = fits[fits["DISTRICT"] == row.DISTRICT]
        admissible = part[part["admissible"]].sort_values("cdf_rmse")
        assert row.distribution == admissible.iloc[0]["distribution"], row.DISTRICT


def test_lognormal_and_gamma_selections_use_zero_location():
    """Plan section 6: Lognormal/Gamma must be fitted with loc = 0."""
    selected = load("arsenic_selected_distributions.csv")
    for row in selected.itertuples():
        if row.distribution in af.ZERO_LOCATION_FAMILIES:
            dist = af.frozen(row.distribution, {n: getattr(row, f"param_{n}")
                                                for n in af.PARAM_NAMES[row.distribution]})
            lo, _ = dist.support()
            assert lo == 0.0, row.DISTRICT
        fits = load("arsenic_distribution_fit_results.csv")
        for family in af.ZERO_LOCATION_FAMILIES:
            part = fits[(fits["DISTRICT"] == row.DISTRICT) & (fits["distribution"] == family)]
            assert bool(part["zero_location_enforced"].iloc[0]) is True


# ---------------------------------------------------------------------------
# 5. Sampling
# ---------------------------------------------------------------------------

def test_sampling_from_every_selection_yields_finite_non_negative_draws(selections):
    """Plan section 6 test 5."""
    assert len(selections) == len(DISTRICTS)
    for record in selections:
        draws = af.sample_arsenic(record, 10_000, np.random.default_rng(0))
        assert draws.shape == (10_000,)
        assert np.isfinite(draws).all(), record["district"]
        assert (draws >= 0).all(), record["district"]
        assert draws.std() > 0
        # The SciPy spec must reconstruct the same distribution.
        spec = dict(record["scipy"])
        name = spec.pop("scipy_dist")
        from scipy import stats

        rebuilt = getattr(stats, name)(**spec)
        u = np.linspace(1e-6, 1 - 1e-6, 501)
        assert np.allclose(rebuilt.ppf(u), af.frozen(record["distribution"],
                                                     record["params"]).ppf(u), rtol=1e-9)


# ---------------------------------------------------------------------------
# 6-7. Convergence reporting and high-censoring labelling
# ---------------------------------------------------------------------------

def test_every_fit_row_reports_an_explicit_convergence_state():
    """Plan section 6 test 6."""
    fits = load("arsenic_distribution_fit_results.csv")
    assert fits["ls_converged"].dtype == bool
    assert fits["mle_converged"].dtype == bool
    assert fits["ls_message"].notna().all() and (fits["ls_message"].astype(str).str.len() > 0).all()
    assert fits["mle_message"].notna().all() and (fits["mle_message"].astype(str).str.len() > 0).all()
    assert fits["ls_n_starts"].min() >= 2, "multi-start is required"
    assert (fits["objective_points"] >= 2).all()
    assert (fits["objective_points"] <= fits["objective_points_total"]).all()


def test_high_censoring_districts_are_labelled_in_the_rationale():
    """Plan section 6 test 7: Rajshahi >= 50%, plus Dhaka and Mymensingh near 48%."""
    selected = load("arsenic_selected_distributions.csv").set_index("DISTRICT")
    summary = load("arsenic_district_preprocessing_summary.csv").set_index("DISTRICT")

    rajshahi = selected.loc["Rajshahi"]
    assert rajshahi["censoring_percent"] >= ACFG["high_censoring_warning_percent"]
    assert "HIGH-CENSORING DISTRICT" in rajshahi["rationale"]
    assert "0.711" in rajshahi["rationale"]
    assert summary.loc["Rajshahi", "high_censoring_warning"]

    for district in ("Dhaka", "Mymensingh"):
        row = selected.loc[district]
        assert row["censoring_percent"] >= 45.0
        assert "censoring" in row["rationale"]
        assert "decision 9.1 V2" in row["rationale"]

    # Every high-censoring district carries the caveat; low-censoring ones need not.
    for district, row in selected.iterrows():
        if row["censoring_percent"] >= 30.0:
            assert "reverse-KM curve is flat" in row["rationale"], district


def test_selection_status_is_provisional_and_rationale_is_substantive():
    """Selection is a hold point: nothing may be presented as frozen."""
    selected = load("arsenic_selected_distributions.csv")
    assert selected["status"].str.startswith("provisional").all()
    for text in selected["rationale"]:
        assert len(text) > 120
    assert selected["rationale"].str.contains("lowest CDF RMSE among admissible").all()


def test_selection_uses_more_than_one_statistic():
    """Plan section 8: selection used more than one statistic."""
    selected = load("arsenic_selected_distributions.csv")
    for text in selected["rationale"]:
        assert "AICc" in text
        assert "tail quantile agreement" in text
        assert "max CDF discrepancy" in text


def test_ls_and_aicc_agreement_or_documented_dissent():
    """Plan section 7 step 6: agreement or an explicit record of the disagreement.

    The AICc comparison is made among *admissible* families, because Normal is
    rejected by the physical gate in every district and can never be the benchmark
    the section has in mind.
    """
    fits = load("arsenic_distribution_fit_results.csv")
    selected = load("arsenic_selected_distributions.csv").set_index("DISTRICT")
    for district, part in fits.groupby("DISTRICT"):
        admissible = part[part["admissible"]]
        best_aicc = admissible.sort_values("aicc").iloc[0]["distribution"]
        chosen = selected.loc[district, "distribution"]
        rationale = selected.loc[district, "rationale"]
        if best_aicc == chosen:
            assert "AICc corroborates" in rationale, district
        else:
            assert "DISSENT" in rationale, district
            assert best_aicc in rationale, district
            assert "recorded, not used as a tie-break" in rationale, district


# ---------------------------------------------------------------------------
# 8-9. Manifest, figures, bootstrap
# ---------------------------------------------------------------------------

def test_manifest_records_the_inputs_method_and_outputs(manifest):
    assert manifest["config_version"] == CFG["config_version"]
    assert manifest["raw_unchanged_during_run"] and manifest["raw_matches_frozen_manifest"]
    assert manifest["districts"] == DISTRICTS
    assert manifest["fitting"]["objective_region"] == ACFG["objective_region"]
    assert manifest["fitting"]["selection_rule"] == ACFG["selection_rule"]
    assert manifest["fitting"]["selection_status"].startswith("provisional")
    assert manifest["fitting"]["censored_substitution"] == "none"
    assert manifest["fitting"]["lognormal_gamma_location"] == 0.0
    assert manifest["fitting"]["normal_max_nonpositive_prob"] == ACFG["normal_max_nonpositive_prob"]
    assert manifest["fitting"]["ls_multi_start"] == ["moments_detected", "censored_mle"]
    assert manifest["bootstrap"]["master_seed"] == ACFG["bootstrap_seed"]
    assert manifest["rows"]["fit_rows"] == 30
    assert set(manifest["software"]) >= {"python", "numpy", "scipy", "pandas", "matplotlib", "lifelines"}
    assert all(manifest["software"][k] for k in ("python", "numpy", "scipy", "lifelines"))


def test_manifest_outputs_and_figures_all_exist(manifest):
    """Plan section 6 test 8."""
    for key, path in manifest["outputs"].items():
        assert (PROJECT_ROOT / path).is_file(), (key, path)
    assert len(manifest["figures"]) == len(DISTRICTS) * 4 == 40
    for path in manifest["figures"]:
        assert (PROJECT_ROOT / path).is_file(), path
        assert path.startswith("results/figures/arsenic_")
        assert (PROJECT_ROOT / path).stat().st_size > 5_000, path
    for subject in ("ecdf_cdf", "histogram_pdf", "qq", "tails"):
        for district in DISTRICTS:
            assert f"results/figures/arsenic_{subject}_{district}.png" in manifest["figures"]


def test_raw_input_matches_the_frozen_hash_manifest(manifest):
    """Plan section 6 test 8: the raw BGS file is unchanged."""
    frozen = read_hash_manifest(RESULTS_PROVENANCE / "raw_data.sha256")
    raw = raw_path(CFG["raw_inputs"]["arsenic_bgs"])
    assert frozen[manifest["raw_input"]] == manifest["raw_sha256"]
    assert sha256_file(raw) == manifest["raw_sha256"]
    assert manifest["raw_input"] == "data/raw/NationalSurveyData.csv"


def test_bootstrap_table_covers_every_district_and_family(manifest):
    """Plan section 6 test 9."""
    boot = load("arsenic_bootstrap_stability.csv")
    assert len(boot) == len(DISTRICTS) * len(af.FAMILIES) == 30
    assert set(boot["DISTRICT"]) == set(DISTRICTS)
    assert set(boot["distribution"]) == set(af.FAMILIES)
    requested = manifest["bootstrap"]["replicates_requested"]
    assert requested == ACFG["bootstrap_replicates"]
    for district, part in boot.groupby("DISTRICT"):
        assert (part["bootstrap_replicates_requested"] == requested).all()
        assert (part["bootstrap_replicates_used"] > 0).all()
        assert (part["bootstrap_replicates_used"] <= requested).all()
        assert manifest["bootstrap"]["usable_replicates_per_district"][district] == int(
            part["bootstrap_replicates_used"].max())
    # Winner frequencies across families must sum to 1 per district.
    for district, part in boot.groupby("DISTRICT"):
        assert part["ls_rmse_winner_frequency"].sum() == pytest.approx(1.0, abs=1e-9)
        assert part["aicc_winner_frequency"].sum() == pytest.approx(1.0, abs=1e-9)
        assert part["ls_rmse_winner_admissible_frequency"].sum() == pytest.approx(1.0, abs=1e-9)


def test_bootstrap_reports_parameter_intervals_and_rmse_spread():
    boot = load("arsenic_bootstrap_stability.csv")
    assert (boot["bootstrap_replicates_used"] >= 100).all()
    for row in boot.itertuples():
        for name in af.PARAM_NAMES[row.distribution]:
            q025, med, q975 = (getattr(row, f"ls_{name}_q025"), getattr(row, f"ls_{name}_median"),
                               getattr(row, f"ls_{name}_q975"))
            assert np.isfinite([q025, med, q975]).all(), (row.DISTRICT, row.distribution, name)
            assert q025 <= med <= q975, (row.DISTRICT, row.distribution, name)
            assert np.isfinite(getattr(row, f"ls_{name}_rmse_spread"))
        assert np.isfinite(row.ls_rmse_median) and row.ls_rmse_q025 <= row.ls_rmse_q975


def test_bootstrap_stability_of_each_selection_is_recorded(selections):
    boot = load("arsenic_bootstrap_stability.csv").set_index(["DISTRICT", "distribution"])
    for record in selections:
        summary = record.get("bootstrap_stability") or {}
        assert summary, record["district"]
        direct = boot.loc[(record["district"], record["distribution"])]
        assert summary["bootstrap_replicates_used"] == direct["bootstrap_replicates_used"]
        assert summary["aicc_winner_frequency"] == pytest.approx(direct["aicc_winner_frequency"])
        assert "bootstrap" in record["rationale"]
        assert "wins admissible LS-RMSE" in record["rationale"]


def test_rajshahi_bootstrap_is_reported_and_dominant():
    """Plan section 7 step 5: inspect Rajshahi's bootstrap winner frequency."""
    boot = load("arsenic_bootstrap_stability.csv")
    part = boot[boot["DISTRICT"] == "Rajshahi"].set_index("distribution")
    assert part["bootstrap_replicates_used"].iloc[0] >= 100
    gamma = part.loc["Gamma"]
    assert gamma["ls_rmse_winner_admissible_frequency"] >= 0.9
    # Rajshahi's censoring burden is the highest of the ten districts.
    summary = load("arsenic_district_preprocessing_summary.csv").set_index("DISTRICT")
    assert summary["censoring_percent"].idxmax() == "Rajshahi"


def test_fit_metrics_are_physically_sane():
    """CDF RMSE and the KS-style distance are bounded, and every row is populated."""
    fits = load("arsenic_distribution_fit_results.csv")
    for column in ("cdf_sse", "cdf_rmse", "cdf_max_discrepancy", "aic", "aicc",
                   "mle_loglikelihood", "prob_arsenic_nonpositive"):
        assert np.isfinite(fits[column]).all(), column
    assert (fits["cdf_sse"] >= 0).all()
    assert (fits["cdf_rmse"] >= 0).all()
    assert fits["cdf_max_discrepancy"].between(0, 1).all()
    assert (fits["cdf_rmse"] <= fits["cdf_max_discrepancy"] + 1e-12).all()
    assert (fits["aicc"] >= fits["aic"]).all()
    # Lognormal/Gamma always start their support at exactly zero.
    zero_loc = fits[fits["distribution"].isin(af.ZERO_LOCATION_FAMILIES)]
    assert (zero_loc["support_lower"] == 0.0).all()
    assert (zero_loc["prob_arsenic_nonpositive"] == 0.0).all()


def test_censored_counts_in_the_fit_table_match_the_preprocessing_summary():
    """The censoring burden driving every censoring-aware step must reconcile."""
    fits = load("arsenic_distribution_fit_results.csv")
    summary = load("arsenic_district_preprocessing_summary.csv").set_index("DISTRICT")
    for district, part in fits.groupby("DISTRICT"):
        row = part.iloc[0]
        assert int(row["n"]) == int(summary.loc[district, "total_n"])
        assert int(row["n_censored"]) == int(summary.loc[district, "censored_n"])
        assert int(row["n_detected"]) == int(summary.loc[district, "detected_n"])
        assert row["censoring_percent"] == pytest.approx(summary.loc[district, "censoring_percent"], abs=0.01)
        # F_hat_0 never exceeds the raw censored share: detections below the largest
        # bound push the plateau *below* it, and bounds above every detection (a
        # case that breaks the lifelines estimator) are flagged rather than absorbed.
        assert 0.0 <= row["rkm_censoring_fraction"] <= row["n_censored"] / row["n"] + 1e-9, district
        assert bool(row["censoring_above_detection"]) is False, district
        if row["rkm_censoring_fraction"] > 0:
            assert row["rkm_censoring_fraction"] == pytest.approx(
                row["n_censored"] / row["n"], abs=0.031), district


def test_saved_curves_reproduce_the_reverse_km_plateau(manifest):
    """The reverse-KM curve a fresh recomputation produces must match what was saved,
    so the saved fits are reproducible from the input file alone."""
    fits = load("arsenic_distribution_fit_results.csv")
    for district, part in fits.groupby("DISTRICT"):
        x, censored = _district_rows(district)
        curve = af.reverse_km_cdf(x, censored)
        row = part.iloc[0]
        assert int(row["rkm_points"]) == int(curve.x.size)
        assert row["rkm_censoring_fraction"] == pytest.approx(curve.censoring_fraction, abs=1e-9)
        assert bool(row["censoring_above_detection"]) == curve.censoring_above_detection
        assert not curve.censoring_above_detection, district


def test_plateau_notes_describe_the_unresolved_region_honestly():
    """Districts whose low tail is resolved must not be given a plateau caveat, and
    a plateau below the raw censored share must say why."""
    selected = load("arsenic_selected_distributions.csv").set_index("DISTRICT")
    summary = load("arsenic_district_preprocessing_summary.csv").set_index("DISTRICT")
    for district, row in selected.iterrows():
        x, censored = _district_rows(district)
        curve = af.reverse_km_cdf(x, censored)
        if curve.censoring_fraction == 0.0:
            assert "fully resolved" in af.censoring_plateau_note(curve)
            # 11.36% censoring is below the 30% caveat threshold, so no high-censoring
            # note is attached at all for such a district.
            assert row["censoring_percent"] < 30.0
        else:
            note = af.censoring_plateau_note(curve)
            assert "unresolved" in note
            if curve.censoring_fraction + 0.005 < curve.n_censored / curve.n:
                assert "below the raw censored share" in note
        assert row["rkm_censoring_fraction"] == pytest.approx(curve.censoring_fraction, abs=1e-9)
        assert row["censoring_percent"] == pytest.approx(summary.loc[district, "censoring_percent"], abs=0.01)
