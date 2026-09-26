"""Phase 7 tests: seeds, inverse-CDF sampling, input selection, summaries, and saved outputs."""
import json

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from arsenic_hra import monte_carlo as mc
from arsenic_hra import paths
from arsenic_hra import risk_equations as re_
from arsenic_hra import simulation_inputs as si

CFG = paths.load_run_config()
SIM = CFG["simulation"]
DISTRICTS = CFG["arsenic"]["selected_districts"]


@pytest.fixture(scope="module")
def model():
    return mc.load_model_inputs(CFG)


# ---------------------------------------------------------------------------
# Seeds and streams
# ---------------------------------------------------------------------------

def test_same_seed_reproduces_exactly(model):
    spec = mc.run_spec(CFG, "Rajshahi", "child", 1000, replicate=3)
    a, b = mc.simulate(CFG, model, spec), mc.simulate(CFG, model, spec)
    pd.testing.assert_frame_equal(a, b, check_exact=True)


def test_every_run_in_the_design_has_a_distinct_stream():
    keys, first_draws = set(), set()
    for d in DISTRICTS:
        for p in SIM["populations"]:
            for N in SIM["convergence_N"]:
                for r in range(SIM["replicate_seeds_per_N"]):
                    spec = mc.run_spec(CFG, d, p, N, r)
                    keys.add(spec.spawn_key)
                    first_draws.add(float(mc.generator(CFG, spec).uniform()))
    n_runs = len(DISTRICTS) * len(SIM["populations"]) * len(SIM["convergence_N"]) * SIM["replicate_seeds_per_N"]
    assert len(keys) == n_runs and len(first_draws) == n_runs


def test_streams_differ_from_phase4_bootstrap_streams():
    from arsenic_hra.arsenic_bootstrap import replicate_seed
    boot = replicate_seed(SIM["master_seed"], 0, 0).uniform(size=5)
    sim = mc.generator(CFG, mc.run_spec(CFG, DISTRICTS[0], "adult", 1000, 0)).uniform(size=5)
    assert not np.allclose(boot, sim)


def test_inputs_are_explicit_inverse_cdf_of_the_stream(model):
    spec = mc.run_spec(CFG, "Barisal", "adult", 1000, 0)
    frame = mc.simulate(CFG, model, spec)
    u = mc.draw_uniforms(mc.generator(CFG, spec), SIM["input_order"], spec.N)
    dists = model.distributions("Barisal", "adult")
    for name in SIM["input_order"]:
        assert_array_equal(frame[name].to_numpy(), dists[name].ppf(u[name]))
        assert np.all(u[name] > 0) and np.all(u[name] < 1)


def test_saved_outputs_recompute_from_saved_inputs(model):
    frame = pd.read_parquet(paths.RESULTS_SIMULATION / "iterations" / "comilla_child_N10000_r0.parquet")
    out = re_.evaluate_risk({k: frame[k].to_numpy() for k in SIM["input_order"]}, model.params["child"])
    for k in re_.OUTPUTS:
        assert_allclose(frame[k].to_numpy(), out[k], rtol=1e-15)


# ---------------------------------------------------------------------------
# Input selection decisions
# ---------------------------------------------------------------------------

def test_arsenic_override_uses_phase4_least_squares_candidate(model):
    fits = pd.read_csv(si.ARSENIC_FITS)
    provisional = {r["district"]: r for r in json.loads(si.ARSENIC_SELECTED.read_text())}
    for district, rec in model.arsenic.items():
        family = SIM["input_selection"]["arsenic_family_overrides"].get(district, provisional[district]["distribution"])
        assert rec["distribution"] == family
        row = fits[(fits.DISTRICT == district) & (fits.distribution == family)].iloc[0]
        for name, value in rec["params"].items():
            assert value == pytest.approx(row[f"ls_{name}"], rel=1e-9)
    assert model.arsenic["Barisal"]["distribution"] == "Gamma"


def test_truncated_normal_fit_recovers_parameters():
    rng = np.random.default_rng(0)
    true = {"mu": 11.0, "sigma": 3.0, "lower": 2.0}
    x = si.frozen("TruncatedNormal", true).ppf(rng.uniform(size=20_000))
    w = rng.uniform(0.2, 3.0, size=x.size)
    fit = si.fit_truncated_normal_bw(x, w, 2.0)
    assert fit["ls_converged"] and fit["pml_converged"]
    assert fit["params"]["mu"] == pytest.approx(11.0, rel=0.01)
    assert fit["params"]["sigma"] == pytest.approx(3.0, rel=0.02)
    scaled = si.fit_truncated_normal_bw(x, 7.5 * w, 2.0)  # survey-weight scale invariance
    assert scaled["params"]["mu"] == pytest.approx(fit["params"]["mu"], rel=1e-6)
    with pytest.raises(ValueError):
        si.fit_truncated_normal_bw(np.r_[x, 1.0], np.r_[w, 1.0], 2.0)


def test_child_bw_selection_beats_gamma_on_plan_criteria(model):
    review = pd.read_csv(paths.RESULTS_SIMULATION / "child_bw_selection_review.csv").set_index("role")
    sel, gam = review.loc["selected"], review.loc["Phase 5 provisional"]
    assert sel["weighted_cdf_sse"] < gam["weighted_cdf_sse"]
    assert sel["weighted_cdf_max_discrepancy"] < gam["weighted_cdf_max_discrepancy"]
    assert sel["tail_max_rel_quantile_error"] < gam["tail_max_rel_quantile_error"]
    assert sel["ls_pml_max_rel_param_diff"] < 0.02
    bw = model.bodyweight["child"]
    assert bw["distribution"] == "TruncatedNormal" and bw["params"]["lower"] == 1.6
    draws = si.frozen(bw["distribution"], bw["params"]).ppf(np.linspace(1e-12, 1 - 1e-12, 1001))
    assert draws.min() >= 1.6


def test_reverse_km_mean_matches_hand_calculation():
    # Data: <0.5 (censored), 1, 2, 4 mg/L. The censored row lies below the smallest detection, so the
    # reverse-KM plateau is F0 = 1/4 and the curve is the ordinary ECDF with that mass at 0.
    from arsenic_hra import arsenic_fitting as af
    v, c = np.array([0.5, 1.0, 2.0, 4.0]), np.array([True, False, False, False])
    curve = af.reverse_km_cdf(v, c)
    xs = np.r_[0.0, curve.x[curve.x > 0]]
    Fs = np.r_[curve.censoring_fraction, curve.F[curve.x > 0]]
    assert si.reverse_km_mean(v, c) == pytest.approx(np.sum((1 - Fs[:-1]) * np.diff(xs)))
    # Hand value: (0 + 1 + 2 + 4) / 4 = 1.75 (censored mass placed at 0, the lower-bound convention).
    assert si.reverse_km_mean(v, c) == pytest.approx(1.75)


# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------

def test_threshold_probabilities_are_direct_counts():
    spec = mc.RunSpec("X", "adult", 5, 0, 0, 0, 0)
    frame = pd.DataFrame({m: [0.5, 1.0, 1.5, 3.0, 0.2] for m in mc.SUMMARY_METRICS})
    frame["ELCR"] = frame["ELCR_ED_AT"] = [1e-5, 1e-4, 2e-4, 5e-4, 1e-6]
    rows = {r["metric"]: r for r in mc.summarize(frame, spec)}
    assert rows["HI"]["P_gt_1"] == pytest.approx(2 / 5)      # 1.0 is not > 1
    assert rows["HI"]["P_gt_2"] == pytest.approx(1 / 5)
    assert rows["ELCR"]["P_gt_0.0001"] == pytest.approx(2 / 5)
    assert rows["HI"]["P_gt_1_mc_se"] == pytest.approx(np.sqrt(0.4 * 0.6 / 5))


def test_wilson_interval_matches_statsmodels():
    from statsmodels.stats.proportion import proportion_confint
    for k, n in ((0, 100), (37, 1000), (9999, 10000)):
        assert_allclose(mc.wilson_interval(k, n), proportion_confint(k, n, method="wilson"), rtol=1e-10, atol=1e-15)


# ---------------------------------------------------------------------------
# Saved-output contract
# ---------------------------------------------------------------------------

def test_manifest_checks_and_runs():
    m = json.loads((paths.RESULTS_SIMULATION / "simulation_manifest.json").read_text())
    assert all(v is True for k, v in m["checks"].items() if k != "ks_within_5pct")
    assert len(m["runs"]) == len(DISTRICTS) * len(SIM["populations"])
    assert m["N"] == SIM["primary_N"] and m["master_seed"] == SIM["master_seed"]
    from arsenic_hra.provenance import sha256_file
    for run in m["runs"]:
        assert sha256_file(paths.PROJECT_ROOT / run["output"]) == run["sha256"]


def test_reported_tables_reconcile_with_iterations():
    p95 = pd.read_csv(paths.RESULTS_TABLES / "risk_probabilistic_p95.csv").set_index("district")
    exc = pd.read_csv(paths.RESULTS_TABLES / "risk_exceedance.csv").set_index(["district", "population"])
    for d in ("Dhaka", "Rajshahi"):
        for pop in SIM["populations"]:
            f = pd.read_parquet(paths.RESULTS_SIMULATION / "iterations" / f"{d.lower()}_{pop}_N10000_r0.parquet")
            assert p95.loc[d, f"P95_HI_{pop}"] == pytest.approx(np.quantile(f["HI"], 0.95), rel=1e-12)
            assert p95.loc[d, f"P95_ELCR_{pop}"] == pytest.approx(np.quantile(f["ELCR"], 0.95), rel=1e-12)
            assert exc.loc[(d, pop), "P_HI_gt_1"] == pytest.approx(np.mean(f["HI"] > 1))
            assert exc.loc[(d, pop), "P_HI_gt_2"] == pytest.approx(np.mean(f["HI"] > 2))
