"""Unit tests for the censoring-aware arsenic fitting core (plan.md section 6).

Every test here is synthetic or hand-computed: it needs no real dataset, so the tests
pin the mathematics independently of whatever the BGD districts happen to show.
"""

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose

from arsenic_hra import arsenic_bootstrap as ab
from arsenic_hra import arsenic_fitting as af

TRUE = {
    "Normal": {"mu": 0.05, "sigma": 0.02},
    "Lognormal": {"mu_log": np.log(0.05), "sigma_log": 0.6},
    "Gamma": {"shape_k": 3.0, "scale_theta": 0.02},
}


# ---------------------------------------------------------------------------
# 1-3. Reverse Kaplan-Meier
# ---------------------------------------------------------------------------

def test_reverse_km_matches_hand_computed_step_heights():
    """A hand-worked fixture: 4 detected values and 2 censored bounds.

    The bounds sit below every detection, which is the pattern the reverse-KM
    estimator is valid for (and the pattern the BGS data shows). With the two
    censored rows tied at a single low bound, the curve opens flat at the censoring
    fraction and rises over the four detections.
    """
    values = np.array([0.0005, 0.0005, 0.001, 0.002, 0.003, 0.004])
    censored = np.array([True, True, False, False, False, False])
    curve = af.reverse_km_cdf(values, censored)

    assert curve.n == 6 and curve.n_detected == 4 and curve.n_censored == 2
    assert_allclose(curve.F[0], 2 / 6)        # 2/6 censored -> initial plateau
    assert_allclose(curve.F[-1], 1.0)         # the curve must reach 1
    assert np.all(np.diff(curve.F) >= -1e-12)
    # The plateau spans everything below the smallest detection.
    assert curve.x[curve.F <= curve.censoring_fraction].max() <= curve.smallest_detected
    assert curve.smallest_detected == 0.001
    assert curve.censoring_above_detection is False
    # Four detections produce four further steps above the plateau.
    assert_allclose(curve.F[curve.F > curve.censoring_fraction],
                    [3 / 6, 4 / 6, 5 / 6, 1.0])


def test_reverse_km_flags_censoring_above_the_largest_detection():
    """A bound above the largest detection breaks the lifelines estimator; the
    condition must be reported, never silently absorbed."""
    values = np.array([0.001, 0.002, 0.003, 0.006, 0.006])
    censored = np.array([False, False, False, True, True])
    curve = af.reverse_km_cdf(values, censored)
    assert curve.censoring_above_detection is True
    assert curve.n_censored == 2 and curve.n_detected == 3


def test_reverse_km_is_monotone_and_bounded():
    rng = np.random.default_rng(0)
    bounds = np.full(40, 0.0005)
    detected = rng.uniform(0.001, 0.5, 60)
    x = np.concatenate([detected, bounds])
    censored = np.concatenate([np.zeros(60, dtype=bool), np.ones(40, dtype=bool)])
    curve = af.reverse_km_cdf(x, censored)
    assert np.all(np.diff(curve.F) >= -1e-12)
    assert curve.F.min() >= 0.0 and curve.F.max() <= 1.0
    assert np.all(curve.x >= 0.0)
    assert 0.0 < curve.censoring_fraction < 1.0


def test_censored_rows_are_never_treated_as_detections():
    """The censoring mass must survive: the plateau must reflect the censored share
    and no bound may enter the curve as if it were a measurement."""
    values = np.array([0.0005, 0.0005, 0.0005, 0.002, 0.004, 0.008])
    censored = np.array([True, True, True, False, False, False])
    curve = af.reverse_km_cdf(values, censored)
    assert curve.n_censored == 3
    assert_allclose(curve.censoring_fraction, 3 / 6)
    assert curve.smallest_detected == 0.002
    # Every plateau point is below the smallest detection, i.e. unresolved.
    assert curve.x[curve.F <= curve.censoring_fraction].max() < curve.smallest_detected


@pytest.mark.parametrize("values,censored,match", [
    ([0.01, 0.02], [True, True], "every value is censored"),
    ([0.01, 0.02], [False, False], "no censored values"),
    ([0.01, np.nan], [False, True], "finite"),
    ([0.01, -0.02], [False, True], "non-negative"),
    ([0.01, 0.02], [False], "equal shape"),
])
def test_reverse_km_guards(values, censored, match):
    with pytest.raises(ValueError, match=match):
        af.reverse_km_cdf(values, censored)


def test_all_censored_district_raises_clear_error():
    """Plan section 6 test 10."""
    with pytest.raises(ValueError, match="every value is censored"):
        af.fit_district("Nowhere", np.full(10, 0.006), np.ones(10, dtype=bool))


# ---------------------------------------------------------------------------
# Objective region (frozen decision 9.1)
# ---------------------------------------------------------------------------

def test_objective_mask_variants_agree_with_the_plan_table():
    """The V1/V2/V3 masks must reproduce the point counts quoted in plan.md section 9.1."""
    values = np.array([0.001, 0.0036, 0.0059, 0.0006, 0.0006, 0.0006, 0.0006])
    censored = np.array([False, False, False, True, True, True, True])
    curve = af.reverse_km_cdf(values, censored)
    assert curve.censoring_above_detection is False
    mask_plateau = curve.F <= curve.censoring_fraction
    assert mask_plateau.sum() == 2            # x=0 placeholder plus the bound itself
    assert af.objective_mask(curve, "all").sum() == curve.x.size
    assert af.objective_mask(curve, "above_initial_mass").sum() == curve.x.size - mask_plateau.sum()
    assert af.objective_mask(curve, "above_max_limit").sum() == 3
    with pytest.raises(ValueError, match="unknown objective region"):
        af.objective_mask(curve, "nonsense")


def test_objective_mask_drops_the_x_zero_placeholder():
    values = np.array([0.0005, 0.0005, 0.002, 0.004])
    censored = np.array([True, True, False, False])
    curve = af.reverse_km_cdf(values, censored)
    assert curve.x[0] == 0.0
    assert not af.objective_mask(curve, "above_initial_mass")[0]


# ---------------------------------------------------------------------------
# 4-5, 7. Fits recover known parameters
# ---------------------------------------------------------------------------

def _censored_sample(family, n=6000, censor_prob=0.25, seed=0):
    """A censored sample drawn from a known family.

    Detection limits are drawn *below* the detected values so the reverse-KM
    estimator stays inside its validity region (bounds above the largest detection
    are flagged by ``censoring_above_detection`` rather than silently absorbed).
    """
    rng = np.random.default_rng(seed)
    x = np.asarray(af.frozen(family, TRUE[family]).ppf(rng.uniform(size=n)), dtype=float)
    limit = float(np.quantile(x, censor_prob))
    censored = x <= limit
    # Move the censored rows onto a single bound just below the surviving minimum.
    bound = 0.99 * float(x[~censored].min())
    x = np.where(censored, bound, x)
    return x, censored


@pytest.mark.parametrize("family", af.FAMILIES)
def test_ls_recovers_known_parameters_from_a_large_sample(family):
    """Plan section 6 test 4: LS on a large near-uncensored sample recovers truth."""
    x, censored = _censored_sample(family, seed=1, censor_prob=0.05)
    fit = af.fit_cdf_least_squares(x, censored, family)
    assert fit.success
    for name, expected in TRUE[family].items():
        assert fit.params[name] == pytest.approx(expected, abs=0.15 * abs(expected)), name


@pytest.mark.parametrize("family", af.FAMILIES)
def test_censored_mle_recovers_known_parameters(family):
    """Plan section 6 test 5."""
    x, censored = _censored_sample(family, seed=2, censor_prob=0.05)
    fit = af.fit_censored_mle(x, censored, family)
    assert fit.success and np.isfinite(fit.loglikelihood)
    assert np.isfinite(fit.aic) and np.isfinite(fit.aicc)
    assert fit.aicc >= fit.aic
    for name, expected in TRUE[family].items():
        assert fit.params[name] == pytest.approx(expected, abs=0.15 * abs(expected)), name


def test_censored_loglikelihood_uses_logcdf_for_censored_rows():
    """The censored likelihood is log f(x) for detections and log F(L) for bounds."""
    x = np.array([0.0, 0.0])
    x = np.array([0.03, 0.04, 0.01, 0.02])
    censored = np.array([False, False, True, True])
    params = {"mu": 0.03, "sigma": 0.01}
    dist = af.frozen("Normal", params)
    expected = (np.log(dist.pdf(x[~censored])).sum() + np.log(dist.cdf(x[censored])).sum())
    assert af.censored_loglikelihood("Normal", params, x, censored) == pytest.approx(expected)


def test_censored_mle_rejects_a_sample_with_nothing_censored():
    with pytest.raises(ValueError, match="no censored values"):
        af.fit_censored_mle([0.01, 0.02, 0.03], [False, False, False], "Normal")


@pytest.mark.parametrize("family", af.ZERO_LOCATION_FAMILIES)
def test_lognormal_and_gamma_ls_fits_use_zero_location(family):
    """Plan section 6 test 7."""
    x, censored = _censored_sample(family, seed=3)
    fit = af.fit_cdf_least_squares(x, censored, family)
    lo, _ = fit.dist().support()
    assert lo == 0.0
    spec = af.scipy_spec(family, fit.params)
    assert spec["loc"] == 0.0
    assert af.frozen(family, fit.params).cdf(0.0) == pytest.approx(0.0, abs=1e-12)


def test_multi_start_is_recorded_and_beats_single_starts():
    """Plan section 3.B: two starts per family, lowest SSE kept."""
    x, censored = _censored_sample("Gamma", seed=4)
    fit = af.fit_cdf_least_squares(x, censored, "Gamma")
    assert fit.n_starts >= 2
    assert set(fit.start_labels) >= {"moments_detected", "censored_mle"}
    for start in ("moments_detected", "censored_mle"):
        other = af.fit_cdf_least_squares(x, censored, "Gamma",
                                         extra_starts=[],
                                         mle_start=af.moment_start("Gamma", x[~censored])
                                         if start == "moments_detected" else None)
        assert fit.sse <= other.sse + 1e-12


# ---------------------------------------------------------------------------
# 6. Round-trip
# ---------------------------------------------------------------------------

def test_saved_parameters_reproduce_the_saved_objective():
    """Plan section 6 test 6: the stored parameters regenerate the stored SSE/RMSE."""
    x, censored = _censored_sample("Lognormal", seed=5)
    region = af.OBJECTIVE_REGION
    curve = af.reverse_km_cdf(x, censored)
    mask = af.objective_mask(curve, region)
    fit = af.fit_cdf_least_squares(x, censored, "Lognormal", curve=curve)
    dist = af.frozen("Lognormal", fit.params)
    recomputed = af.cdf_sse(dist, curve.x[mask], curve.F[mask])
    assert recomputed == pytest.approx(fit.sse, rel=1e-9)
    assert af.cdf_rmse_from_sse(recomputed, int(mask.sum())) == pytest.approx(
        af.cdf_rmse_from_sse(fit.sse, fit.n_points), rel=1e-9)


def test_evaluate_district_fit_reports_every_required_column():
    x, censored = _censored_sample("Gamma", seed=6)
    row = af.evaluate_district_fit("Synthetic", x, censored, "Gamma")
    for column in ("distribution", "cdf_sse", "cdf_rmse", "cdf_max_discrepancy", "aic", "aicc",
                   "admissible", "inadmissible_reason", "ls_converged", "ls_message",
                   "mle_converged", "mle_message", "objective_region", "objective_points"):
        assert column in row, column
    assert row["admissible"] is True
    assert row["inadmissible_reason"] == ""
    assert row["objective_region"] == af.OBJECTIVE_REGION


# ---------------------------------------------------------------------------
# 8. Admissibility gate
# ---------------------------------------------------------------------------

def test_admissibility_gate_rejects_a_normal_with_mass_at_or_below_zero():
    """Plan section 6 test 8, and the verified consequence in section 3.E."""
    x, censored = _censored_sample("Normal", seed=7)
    row = af.evaluate_district_fit("Synthetic", x, censored, "Normal")
    assert row["prob_arsenic_nonpositive"] >= af.NORMAL_MAX_NONPOSITIVE_PROB
    assert row["admissible"] is False
    assert "P(C<=0)" in row["inadmissible_reason"]


def test_admissibility_gate_accepts_a_strictly_positive_normal():
    x, censored = _censored_sample("Normal", seed=8)
    params = {"mu": 0.30, "sigma": 0.01}       # P(C<=0) is far below the threshold
    fit = af.LSFit("Normal", params, True, "ok", 0.0, 10, 1)
    assert af.admissibility_reasons("Normal", fit, None) == []


def test_admissibility_gate_rejects_non_convergence_and_bad_magnitudes():
    params = {"mu": 0.05, "sigma": 0.02}
    not_converged = af.LSFit("Normal", params, False, "max iterations reached", 1.0, 10, 1)
    assert any("did not report convergence" in r for r in af.admissibility_reasons("Normal", not_converged, None))

    huge = af.LSFit("Normal", {"mu": 0.05, "sigma": 1e12}, True, "ok", 0.0, 10, 1)
    assert any("magnitude band" in r for r in af.admissibility_reasons("Normal", huge, None))

    nan_params = af.LSFit("Gamma", {"shape_k": np.nan, "scale_theta": 0.02}, True, "ok", 0.0, 10, 1)
    assert any("not finite" in r for r in af.admissibility_reasons("Gamma", nan_params, None))


def test_gamma_shape_below_one_is_admissible():
    """k < 1 is legal and is what heavily censored districts produce; it must not be
    rejected by the magnitude band (which covers scale-like parameters only)."""
    fit = af.LSFit("Gamma", {"shape_k": 0.068, "scale_theta": 0.152}, True, "ok", 0.0, 10, 2)
    assert af.admissibility_reasons("Gamma", fit, None) == []


# ---------------------------------------------------------------------------
# Selection
# ---------------------------------------------------------------------------

def _synthetic_district(family, seed=11, censor_prob=0.3, n=800):
    """A synthetic district whose censoring bounds sit below every detection."""
    rng = np.random.default_rng(seed)
    x = np.asarray(af.frozen(family, TRUE[family]).ppf(rng.uniform(size=n)), dtype=float)
    raw_limit = float(np.quantile(x, censor_prob))
    censored = x <= raw_limit
    bound = 0.99 * float(x[~censored].min())
    x = np.where(censored, bound, x)
    return x, censored


def test_selection_picks_the_lowest_rmse_admissible_family():
    x, censored = _synthetic_district("Gamma", seed=12)
    table = af.fit_district("Synthetic", x, censored)
    chosen, rationale = af.select_arsenic_distribution(table)
    admissible = table[table["admissible"]].sort_values("cdf_rmse")
    assert chosen["distribution"] == admissible.iloc[0]["distribution"]
    assert "lowest CDF RMSE among admissible" in rationale


def test_selection_never_returns_an_inadmissible_family():
    x, censored = _synthetic_district("Normal", seed=13, censor_prob=0.45)
    table = af.fit_district("Synthetic", x, censored)
    chosen, rationale = af.select_arsenic_distribution(table)
    assert bool(chosen["admissible"])
    if not table.set_index("distribution").loc["Normal", "admissible"]:
        assert chosen["distribution"] != "Normal"
        assert "Normal was rejected by the admissibility gate" in rationale

def test_selection_raises_when_nothing_is_admissible():
    table = pd.DataFrame([
        {"DISTRICT": "X", "distribution": "Normal", "admissible": False, "cdf_rmse": 0.01,
         "cdf_max_discrepancy": 0.01, "tail_max_rel_quantile_error": 0.1, "aicc": 1.0,
         "prob_arsenic_nonpositive": 0.5, "inadmissible_reason": "P(C<=0) too large",
         "objective_points": 10, "objective_points_total": 12, "objective_region": "above_initial_mass",
         "ls_mle_max_rel_param_diff": np.nan},
    ])
    with pytest.raises(RuntimeError, match="no admissible arsenic candidate"):
        af.select_arsenic_distribution(table)


def test_selection_documents_aicc_dissent_instead_of_hiding_it():
    table = pd.DataFrame([
        {"DISTRICT": "X", "distribution": "Gamma", "admissible": True, "cdf_rmse": 0.010,
         "cdf_max_discrepancy": 0.02, "tail_max_rel_quantile_error": 0.05, "aicc": 50.0,
         "prob_arsenic_nonpositive": 0.0, "inadmissible_reason": "", "objective_points": 20,
         "objective_points_total": 22, "objective_region": "above_initial_mass",
         "ls_mle_max_rel_param_diff": 0.1, "censoring_percent": 10.0, "rkm_censoring_fraction": 0.1},
        {"DISTRICT": "X", "distribution": "Lognormal", "admissible": True, "cdf_rmse": 0.020,
         "cdf_max_discrepancy": 0.01, "tail_max_rel_quantile_error": 0.04, "aicc": 10.0,
         "prob_arsenic_nonpositive": 0.0, "inadmissible_reason": "", "objective_points": 20,
         "objective_points_total": 22, "objective_region": "above_initial_mass",
         "ls_mle_max_rel_param_diff": 0.2, "censoring_percent": 10.0, "rkm_censoring_fraction": 0.1},
    ])
    chosen, rationale = af.select_arsenic_distribution(table)
    assert chosen["distribution"] == "Gamma"
    assert "DISSENT" in rationale and "Lognormal" in rationale


def test_high_censoring_note_is_attached_above_thirty_percent():
    low = {"censoring_percent": 10.0, "rkm_censoring_fraction": 0.1}
    high = {"censoring_percent": 71.79, "rkm_censoring_fraction": 0.711}
    assert af.high_censoring_note(low) == ""
    note = af.high_censoring_note(high, warning_percent=50.0, min_detected_mgL=0.0006)
    assert "HIGH-CENSORING DISTRICT" in note and "0.711" in note and "decision 9.1 V2" in note
    mid = af.high_censoring_note({"censoring_percent": 45.0, "rkm_censoring_fraction": 0.45},
                                 min_detected_mgL=0.001)
    assert mid.startswith("notable censoring")


def test_plateau_note_distinguishes_resolved_from_unresolved_low_tails():
    """A plateau at 0 arises only when *every* bound sits above every detection;
    a positive plateau below the raw censored share must say why."""
    # Every bound above every detection: F_hat_0 collapses to 0 (the flagged case).
    resolved = af.reverse_km_cdf(np.array([0.001, 0.002, 0.006, 0.006, 0.006]),
                                 np.array([False, False, True, True, True]))
    assert resolved.censoring_above_detection is True
    assert resolved.censoring_fraction == 0.0
    assert "fully resolved" in af.censoring_plateau_note(resolved)

    # Khulna/Sylhet pattern: bounds at the largest limit, some detections below it.
    # Those low detections constrain the tail, so F_hat_0 drops below the raw share.
    partial = af.reverse_km_cdf(np.array([0.0005, 0.0005, 0.006, 0.006, 0.0008, 0.001, 0.002, 0.003, 0.02, 0.05]),
                                np.array([True, True, True, True, False, False, False, False, False, False]))
    assert partial.censoring_above_detection is False
    assert 0.0 < partial.censoring_fraction < partial.n_censored / partial.n
    note = af.censoring_plateau_note(partial)
    assert "unresolved" in note
    assert "below the raw censored share" in note


# ---------------------------------------------------------------------------
# 9. Sampling interface
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("family", af.FAMILIES)
def test_sampling_is_reproducible_finite_and_non_negative(family):
    """Plan section 6 test 9."""
    x, censored = _synthetic_district(family, seed=14)
    table = af.fit_district("Synthetic", x, censored)
    chosen, rationale = af.select_arsenic_distribution(table)
    # Sample from whichever admissible candidates exist, plus the chosen one.
    records = [af.selected_record(chosen, rationale)]
    for row in table.itertuples():
        if row.admissible and row.distribution != chosen["distribution"]:
            params = {name: getattr(row, f"ls_{name}") for name in af.PARAM_NAMES[row.distribution]}
            records.append({"district": "Synthetic", "distribution": row.distribution, "params": params})
    assert records
    for record in records:
        seed = np.random.SeedSequence(7)
        first = af.sample_arsenic(record, 5_000, np.random.default_rng(seed))
        second = af.sample_arsenic(record, 5_000, np.random.default_rng(seed))
        assert_allclose(first, second)
        assert np.isfinite(first).all()
        assert (first >= 0).all()
        assert first.size == 5_000 and first.std() > 0


def test_sampling_enforces_the_same_marginal_as_the_frozen_distribution():
    record = {"district": "X", "distribution": "Gamma", "params": {"shape_k": 3.0, "scale_theta": 0.02}}
    u = np.linspace(1e-6, 1 - 1e-6, 2001)
    assert_allclose(af.sample_arsenic(record, u=u), af.frozen("Gamma", record["params"]).ppf(u))


def test_sampling_clips_the_normal_branch_to_non_negative():
    record = {"district": "X", "distribution": "Normal",
              "params": {"mu": 0.0, "sigma": 0.02}}   # half the mass is below zero
    u = np.linspace(0.001, 0.999, 999)
    draws = af.sample_arsenic(record, u=u)
    assert (draws >= 0).all() and np.isfinite(draws).all()
    assert (draws == 0).any()


def test_sampling_rejects_bad_arguments():
    record = {"district": "X", "distribution": "Gamma", "params": {"shape_k": 3.0, "scale_theta": 0.02}}
    with pytest.raises(ValueError, match="provide u, or both N and rng"):
        af.sample_arsenic(record, N=10)
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        af.sample_arsenic(record, u=np.array([0.5, 1.5]))
    with pytest.raises(ValueError, match="finite"):
        af.sample_arsenic(record, u=np.array([0.5, np.nan]))


# ---------------------------------------------------------------------------
# Bootstrap helpers
# ---------------------------------------------------------------------------

def test_bootstrap_seed_follows_the_project_convention():
    expected = np.random.default_rng(np.random.SeedSequence(20260923, spawn_key=(3, 17)))
    got = ab.replicate_seed(20260923, 3, 17)
    assert_allclose(got.uniform(size=5), expected.uniform(size=5))
    assert not np.allclose(got.uniform(size=5), ab.replicate_seed(20260923, 4, 17).uniform(size=5))


def test_break_ties_is_identity_when_values_are_distinct():
    x = np.array([0.001, 0.002, 0.003, 0.005])
    assert_allclose(ab.break_ties(x), x)


def test_break_ties_separates_duplicates_and_stays_close():
    x = np.array([0.006, 0.001, 0.006, 0.002, 0.006])
    separated = ab.break_ties(x, tolerance=1e-9)
    assert np.unique(separated).size == x.size
    assert np.max(np.abs(separated - x)) < 1e-9
    assert_allclose(np.sort(separated), np.sort(separated))  # order preserved within ties


def test_bootstrap_replicate_reports_a_usable_row():
    x, censored = _synthetic_district("Gamma", seed=15, censor_prob=0.5)
    spec = ab.BootstrapSpec("Synthetic", 0, x, censored, master_seed=20260923)
    row = ab.bootstrap_replicate(spec, 0)
    assert row["status"] == "ok"
    assert row["DISTRICT"] == "Synthetic" and row["replicate"] == 0
    assert row["ls_rmse_winner_admissible"] in (None, *af.FAMILIES)
    for family in af.FAMILIES:
        assert np.isfinite(row[f"ls_rmse_{family.lower()}"])


def test_bootstrap_marks_degenerate_replicates_as_skipped_not_crashed():
    x = np.full(6, 0.006)
    censored = np.ones(6, dtype=bool)
    spec = ab.BootstrapSpec("AllCensored", 0, x, censored, master_seed=20260923)
    row = ab.bootstrap_replicate(spec, 0)
    assert row["status"] == "skipped"
    assert row["skipped_reason"] == "every resampled row is censored"
    assert row["ls_rmse_winner"] is None


def test_bootstrap_driver_covers_every_district_and_replicate():
    frames = []
    for i, district in enumerate(("Alpha", "Beta")):
        x, censored = _synthetic_district("Gamma", seed=20 + i, censor_prob=0.4, n=300)
        frames.append(pd.DataFrame({"DISTRICT": district, "As_model_value_mgL": x, "As_censored": censored}))
    data = pd.concat(frames, ignore_index=True)
    table = ab.run_bootstrap(data, replicates=6, master_seed=20260923, jobs=1,
                             district_order=["Alpha", "Beta"])
    assert len(table) == 2 * 6
    assert set(table["DISTRICT"]) == {"Alpha", "Beta"}
    assert sorted(table["replicate"].unique()) == list(range(6))
    stability = ab.summarise_bootstrap(table)
    assert len(stability) == 2 * len(af.FAMILIES)
    assert set(stability["DISTRICT"]) == {"Alpha", "Beta"}
    for district in ("Alpha", "Beta"):
        part = stability[stability["DISTRICT"] == district]
        assert_allclose(part["ls_rmse_winner_frequency"].sum(), 1.0, atol=1e-9)
    assert "usable replicates" in ab.bootstrap_rationale(stability, "Alpha")


def test_bootstrap_driver_is_seed_reproducible():
    x, censored = _synthetic_district("Gamma", seed=30, censor_prob=0.4, n=300)
    data = pd.DataFrame({"DISTRICT": "Alpha", "As_model_value_mgL": x, "As_censored": censored})
    a = ab.run_bootstrap(data, replicates=5, master_seed=20260923, jobs=1, district_order=["Alpha"])
    b = ab.run_bootstrap(data, replicates=5, master_seed=20260923, jobs=1, district_order=["Alpha"])
    pd.testing.assert_frame_equal(a, b)
    c = ab.run_bootstrap(data, replicates=5, master_seed=999, jobs=1, district_order=["Alpha"])
    assert not np.allclose(a["ls_rmse_gamma"], c["ls_rmse_gamma"])


# ---------------------------------------------------------------------------
# Config pinning
# ---------------------------------------------------------------------------

def test_configure_constants_pins_the_frozen_values_and_restores_them():
    from arsenic_hra.paths import load_run_config

    acfg = load_run_config()["arsenic"]
    before = (af.NORMAL_MAX_NONPOSITIVE_PROB, af.OBJECTIVE_REGION, af.SELECTION_RULE)
    try:
        af.configure_constants({**acfg, "normal_max_nonpositive_prob": 1e-6,
                                "objective_region": "all",
                                "parameter_magnitude_band": {"min": 1e-9, "max": 1e9}})
        assert af.NORMAL_MAX_NONPOSITIVE_PROB == 1e-6
        assert af.OBJECTIVE_REGION == "all"
        assert af.CONFIGURED["objective_region"] == "all"
    finally:
        af.configure_constants(acfg)
    assert (af.NORMAL_MAX_NONPOSITIVE_PROB, af.OBJECTIVE_REGION, af.SELECTION_RULE) == before
