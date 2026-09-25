import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose

from arsenic_hra import bodyweight_fitting as ft

TRUE = {
    "Normal": {"mu": 55.0, "sigma": 10.0},
    "Lognormal": {"mu_log": 4.0, "sigma_log": 0.2},
    "Gamma": {"shape_k": 11.0, "scale_theta": 1.0},
    "Triangular": {"a_min": 2.0, "m_mode": 11.0, "b_max": 20.0},
}


def synthetic(family, n=20_000, seed=0):
    rng = np.random.default_rng(seed)
    x = ft.frozen(family, TRUE[family]).ppf(rng.uniform(size=n))
    w = rng.uniform(0.2, 3.0, size=n)  # survey weights independent of x
    return x, w


def test_ecdf_midpoints_aggregate_ties():
    e = ft.weighted_ecdf_midpoints([2.0, 1.0, 1.0], [2.0, 1.0, 1.0])
    assert_allclose(e.x, [1.0, 2.0])
    assert_allclose(e.q, [0.5, 0.5])
    assert_allclose(e.p, [0.25, 0.75])
    assert_allclose(e.F, [0.5, 1.0])


@pytest.mark.parametrize("w", [[1.0, 0.0], [1.0, -1.0], [1.0, np.nan]])
def test_ecdf_rejects_bad_weights(w):
    with pytest.raises(ValueError):
        ft.weighted_ecdf_midpoints([1.0, 2.0], w)


@pytest.mark.parametrize("family", ft.FAMILIES)
def test_ls_recovers_parameters(family):
    x, w = synthetic(family)
    fit = ft.fit_cdf_least_squares(x, w, family)
    assert fit.success
    for k, v in TRUE[family].items():
        assert fit.params[k] == pytest.approx(v, rel=0.03), k


@pytest.mark.parametrize("family", ft.FAMILIES)
def test_pml_recovers_parameters(family):
    x, w = synthetic(family)
    fit = ft.fit_weighted_pseudo_mle(x, w, family)
    assert np.isfinite(fit.objective)
    for k, v in TRUE[family].items():
        assert fit.params[k] == pytest.approx(v, rel=0.03), k


def test_pml_matches_closed_form_weighted_mle():
    x, w = synthetic("Lognormal", n=5000, seed=3)
    q = w / w.sum()
    normal = ft.fit_weighted_pseudo_mle(x, w, "Normal")
    mu = np.sum(q * x)
    assert normal.params["mu"] == pytest.approx(mu, rel=1e-6)
    assert normal.params["sigma"] == pytest.approx(np.sqrt(np.sum(q * (x - mu) ** 2)), rel=1e-5)
    lognormal = ft.fit_weighted_pseudo_mle(x, w, "Lognormal")
    mlog = np.sum(q * np.log(x))
    assert lognormal.params["mu_log"] == pytest.approx(mlog, rel=1e-6)
    assert lognormal.params["sigma_log"] == pytest.approx(np.sqrt(np.sum(q * (np.log(x) - mlog) ** 2)), rel=1e-5)


@pytest.mark.parametrize("family", ft.FAMILIES)
def test_weight_scale_invariance(family):
    x, w = synthetic(family, n=3000, seed=1)
    for fitter in (ft.fit_cdf_least_squares, ft.fit_weighted_pseudo_mle):
        base = fitter(x, w, family).params
        scaled = fitter(x, 1234.5 * w, family).params
        for k in base:
            assert scaled[k] == pytest.approx(base[k], rel=1e-5), (fitter.__name__, k)


@pytest.mark.parametrize("family", ft.FAMILIES)
def test_saved_parameters_reproduce_objective(family):
    x, w = synthetic(family, n=3000, seed=2)
    fit = ft.fit_cdf_least_squares(x, w, family)
    ecdf = ft.weighted_ecdf_midpoints(x, w)
    assert ft.weighted_sse(ft.frozen(family, fit.params), ecdf) == pytest.approx(fit.objective, rel=1e-10)


def test_triangular_transforms_enforce_ordering_and_support():
    rng = np.random.default_rng(5)
    x = np.array([3.0, 7.0, 15.0])
    for t in rng.normal(scale=3.0, size=(200, 3)):
        p = ft._to_natural_ls("Triangular", t)
        assert 0 < p["a_min"] < p["m_mode"] < p["b_max"]
        p = ft._to_natural_pml("Triangular", t, x)
        assert 0 < p["a_min"] < x.min() and p["b_max"] > x.max()
        assert p["a_min"] <= p["m_mode"] <= p["b_max"]


def test_scipy_spec_matches_frozen():
    grid = np.linspace(1, 30, 50)
    from scipy import stats
    for family, params in TRUE.items():
        spec = dict(ft.scipy_spec(family, params))
        dist = getattr(stats, spec.pop("scipy_dist"))(**spec)
        assert_allclose(dist.cdf(grid), ft.frozen(family, params).cdf(grid))


def test_sample_bw_inverse_cdf():
    sel = {"distribution": "Gamma", "params": TRUE["Gamma"]}
    a = ft.sample_bw(sel, 10_000, np.random.default_rng(42))
    b = ft.sample_bw(sel, 10_000, np.random.default_rng(42))
    assert_allclose(a, b)
    assert np.all(np.isfinite(a)) and np.all(a > 0)
    u = np.array([0.05, 0.5, 0.95])
    assert_allclose(ft.sample_bw(sel, u=u), ft.frozen("Gamma", TRUE["Gamma"]).ppf(u))
    big = ft.sample_bw(sel, 200_000, np.random.default_rng(7))
    assert_allclose(np.quantile(big, [0.05, 0.5, 0.95]), ft.frozen("Gamma", TRUE["Gamma"]).ppf([0.05, 0.5, 0.95]),
                    rtol=0.01)
    with pytest.raises(ValueError):
        ft.sample_bw(sel)


def test_normal_with_negative_mass_is_inadmissible():
    x, w = synthetic("Normal", n=3000)
    ls = ft.fit_cdf_least_squares(x, w, "Normal")
    ls.params = {"mu": 10.0, "sigma": 3.0}  # P(BW<0) ~ 4e-4
    row = ft.evaluate_distribution_fit(x, w, ls, ft.fit_weighted_pseudo_mle(x, w, "Normal"))
    assert not row["admissible"] and "P(BW<=0)" in row["inadmissible_reason"]


def test_selection_skips_inadmissible_best_and_raises_when_none():
    x, w = synthetic("Gamma", n=4000, seed=9)
    table = ft.fit_population(x, w, "test")
    table.loc[table["distribution"] == "Gamma", "weighted_cdf_sse"] = 0.0
    table.loc[table["distribution"] == "Gamma", ["admissible", "inadmissible_reason"]] = [False, "forced"]
    chosen, rationale = ft.select_bw_distribution(table)
    assert chosen["distribution"] != "Gamma"
    assert "Gamma had lower SSE but was rejected: forced" in rationale
    table["admissible"] = False
    with pytest.raises(RuntimeError):
        ft.select_bw_distribution(table)


@pytest.mark.parametrize("family", ["Lognormal", "Gamma", "Triangular"])
def test_fit_population_selects_true_family(family):
    x, w = synthetic(family, n=20_000, seed=11)
    table = ft.fit_population(x, w, "t")
    assert table["ls_converged"].all()
    chosen, _ = ft.select_bw_distribution(table)
    assert chosen["distribution"] == family
