"""Phase 6 tests: risk equations, parameter sets, and hand-calculated benchmarks."""
import copy

import numpy as np
import pytest
from numpy.testing import assert_allclose

from arsenic_hra import paths
from arsenic_hra import risk_equations as re_
from arsenic_hra import risk_parameters as rp

CFG = paths.load_run_config()
RTOL = 1e-9


def cfg_with(**risk_overrides):
    cfg = copy.deepcopy(CFG)
    cfg["risk_model"].update(risk_overrides)
    return cfg


def inputs(C, IR, BW, EF, ET, SA):
    return {"C_mg_per_L": C, "IR_L_per_day": IR, "BW_kg": BW,
            "EF_days_per_year": EF, "ET_h_per_day": ET, "SA_m2": SA}


# ---------------------------------------------------------------------------
# Hand-calculated benchmark cases (arithmetic shown; not computed by the module)
# ---------------------------------------------------------------------------

def test_benchmark_adult_round_numbers():
    # C=0.05, IR=2, EF=365, ED=70, BW=60, AT=70*365=25550, ET=0.2, SA=1.8
    # ADD_ing    = 0.05*2*365*70 / (60*25550) = 0.05*2/60           = 1.6666667e-3
    # ADD_dermal = 0.05*0.001*365*70*0.2*1.8*10 / (60*25550)
    #            = 0.05*0.001*3.6/60                                = 3.0e-6
    # HQ_ing = 1.6666667e-3/3e-4 = 5.5555556;  HQ_dermal = 3e-6/2.85e-4 = 0.010526316
    # HI = 5.5660819;  ELCR = 1.6666667e-3*1.5 = 2.5e-3 (adult: lifetime AT = ED*365)
    out = re_.evaluate_risk(inputs(0.05, 2.0, 60.0, 365.0, 0.2, 1.8), rp.population_parameters("adult", CFG))
    assert_allclose(out["ADD_ing_mg_per_kg_day"], 1.6666666667e-3, rtol=RTOL)
    assert_allclose(out["ADD_dermal_mg_per_kg_day"], 3.0e-6, rtol=RTOL)
    assert_allclose(out["HQ_ing"], 5.5555555556, rtol=RTOL)
    assert_allclose(out["HQ_dermal"], 0.0105263158, rtol=1e-8)
    assert_allclose(out["HI"], 5.5660818713, rtol=1e-8)
    assert_allclose(out["ELCR"], 2.5e-3, rtol=RTOL)


def test_benchmark_child_round_numbers_lifetime_and_paper_AT():
    # C=0.05, IR=1, EF=365, ED=6, BW=15, AT_noncancer=6*365=2190, ET=0.2, SA=0.6
    # ADD_ing    = 0.05*1/15                   = 3.3333333e-3
    # ADD_dermal = 0.05*0.001*0.2*0.6*10/15    = 4.0e-6
    # HQ_ing = 11.111111;  HQ_dermal = 4e-6/2.85e-4 = 0.014035088;  HI = 11.125146
    # lifetime:  ADD_ing_cancer = 3.3333333e-3 * 6/70 = 2.8571429e-4 -> ELCR = 4.2857143e-4
    # base-paper AT = ED*365:  ELCR = 3.3333333e-3 * 1.5 = 5.0e-3
    x = inputs(0.05, 1.0, 15.0, 365.0, 0.2, 0.6)
    out = re_.evaluate_risk(x, rp.population_parameters("child", CFG))
    assert_allclose(out["ADD_ing_mg_per_kg_day"], 3.3333333333e-3, rtol=RTOL)
    assert_allclose(out["ADD_dermal_mg_per_kg_day"], 4.0e-6, rtol=RTOL)
    assert_allclose(out["HQ_ing"], 11.111111111, rtol=RTOL)
    assert_allclose(out["HQ_dermal"], 0.0140350877, rtol=1e-8)
    assert_allclose(out["HI"], 11.125146199, rtol=1e-8)
    assert_allclose(out["ADD_ing_cancer_mg_per_kg_day"], 2.8571428571e-4, rtol=RTOL)
    assert_allclose(out["ELCR"], 4.2857142857e-4, rtol=RTOL)

    paper_at = rp.population_parameters("child", cfg_with(cancer_averaging_time="exposure_duration"))
    assert_allclose(re_.evaluate_risk(x, paper_at)["ELCR"], 5.0e-3, rtol=RTOL)


def test_benchmark_adult_non_cancelling_exposure_frequency():
    # EF=300 != 365 so EF*ED does not cancel AT.
    # C=0.01, IR=1.5, EF=300, ED=70, BW=50, AT=25550, ET=0.25, SA=2.0
    # ADD_ing    = 0.01*1.5*300*70 / (50*25550) = 315/1,277,500           = 2.4657534e-4
    # ADD_dermal = 0.01*0.001*300*70*0.25*2*10 / 1,277,500 = 1.05/1,277,500 = 8.2191781e-7
    # HQ_ing = 0.82191781;  HQ_dermal = 2.8839221e-3;  HI = 0.82480173;  ELCR = 3.6986301e-4
    out = re_.evaluate_risk(inputs(0.01, 1.5, 50.0, 300.0, 0.25, 2.0), rp.population_parameters("adult", CFG))
    assert_allclose(out["ADD_ing_mg_per_kg_day"], 315 / 1_277_500, rtol=RTOL)
    assert_allclose(out["ADD_dermal_mg_per_kg_day"], 1.05 / 1_277_500, rtol=RTOL)
    assert_allclose(out["HQ_ing"], 0.82191780822, rtol=1e-9)
    assert_allclose(out["HQ_dermal"], 2.8839221341e-3, rtol=1e-8)
    assert_allclose(out["HI"], 0.82480173035, rtol=1e-9)
    assert_allclose(out["ELCR"], 3.6986301370e-4, rtol=1e-9)


# Yadav & Kalkal (2024) Table 4, from district means in Table 1 and deterministic
# Table 2 inputs (IR 2/1.8 L/day, EF 365, ET 0.58 h/day, BW 70/15 kg, SA 1.8/0.6 m^2,
# AT = ED*365 for both HI and ELCR).
PAPER_TABLE4 = {
    # district: (C mg/L, HI adult, HI child, ELCR adult, ELCR child)
    "Moga": (0.0043, 0.41, 1.72, 1.84e-4, 7.74e-4),
    "Faridkot": (0.017, 1.62, 6.83, 7.29e-4, 3.06e-3),
    "Fazilka": (0.037, 3.54, 14.87, 1.59e-3, 6.66e-3),
    "Patiala": (0.021, 2.01, 8.44, 9.01e-4, 3.78e-3),
    "Ferozepur": (0.027, 2.58, 10.85, 1.16e-3, 4.80e-3),
    "Rupnagar": (0.035, 3.35, 14.07, 1.50e-3, 6.30e-3),
    "Amritsar": (0.083, 7.94, 33.26, 3.50e-3, 1.40e-2),
}


# Published cells that disagree with the paper's own Table 1/Table 2 inputs by more
# than rounding. Checked against the hand arithmetic C * IR / BW * CSF instead.
PAPER_TABLE4_DISCREPANCIES = {
    ("Ferozepur", "child", "ELCR"): 0.027 * 1.8 / 15 * 1.5,   # 4.86e-3; paper prints 4.80e-3 (-1.2%)
    ("Amritsar", "adult", "ELCR"): 0.083 * 2.0 / 70 * 1.5,    # 3.557e-3; paper prints 3.50e-3 (-1.6%)
    ("Amritsar", "child", "ELCR"): 0.083 * 1.8 / 15 * 1.5,    # 1.494e-2; paper prints 1.40e-2 (-6.3%)
}


@pytest.mark.parametrize("district", PAPER_TABLE4)
def test_reproduces_base_paper_table4(district):
    """The equations reproduce the published deterministic table within its printed rounding (0.5%)."""
    C, hi_a, hi_c, elcr_a, elcr_c = PAPER_TABLE4[district]
    cfg = cfg_with(cancer_averaging_time="exposure_duration")
    for pop, bw, hi_ref, elcr_ref in (("adult", 70.0, hi_a, elcr_a), ("child", 15.0, hi_c, elcr_c)):
        p = rp.population_parameters(pop, cfg)
        out = re_.evaluate_risk(re_.deterministic_inputs(p, C, bw), p)
        assert_allclose(out["HI"], hi_ref, rtol=0.005, err_msg=f"{district} {pop} HI")
        known = PAPER_TABLE4_DISCREPANCIES.get((district, pop, "ELCR"))
        if known is None:
            assert_allclose(out["ELCR"], elcr_ref, rtol=0.005, err_msg=f"{district} {pop} ELCR")
        else:
            assert_allclose(out["ELCR"], known, rtol=1e-12)
            assert abs(out["ELCR"] - elcr_ref) / elcr_ref > 0.005


# ---------------------------------------------------------------------------
# Dimensional consistency
# ---------------------------------------------------------------------------

def test_dermal_dose_matches_si_first_principles():
    """CF = 10 L/(cm*m^2) is exactly the unit conversion, checked in SI units."""
    C, Kp, EF, ED, ET, SA, BW = 0.05, 0.001, 345.0, 70.0, 0.2, 1.8, 60.0
    AT = ED * 365
    C_si = C * 1e-6 / 1e-3               # mg/L -> kg/m^3
    Kp_si = Kp * 1e-2 / 3600.0           # cm/h -> m/s
    ET_si = ET * 3600.0                  # h/day -> s/day
    absorbed_kg_per_day = C_si * Kp_si * ET_si * SA          # kg/m^3 * m/s * s/day * m^2
    add_si = absorbed_kg_per_day * EF * ED / (BW * AT)       # kg/(kg*day)
    expected_mg_per_kg_day = add_si * 1e6
    assert_allclose(re_.add_dermal(C, Kp, EF, ED, ET, SA, 10.0, BW, AT), expected_mg_per_kg_day, rtol=1e-12)


def test_ingestion_dose_matches_si_first_principles():
    C, IR, EF, ED, BW = 0.05, 2.0, 345.0, 70.0, 60.0
    AT = ED * 365
    intake_kg_per_day = (C * 1e-6 / 1e-3) * (IR * 1e-3)      # kg/m^3 * m^3/day
    expected = intake_kg_per_day * EF * ED / (BW * AT) * 1e6
    assert_allclose(re_.add_ingestion(C, IR, EF, ED, BW, AT), expected, rtol=1e-12)


@pytest.mark.parametrize("name,exponent_ing,exponent_derm", [
    ("C_mg_per_L", 1, 1), ("IR_L_per_day", 1, 0), ("BW_kg", -1, -1),
    ("EF_days_per_year", 1, 1), ("ET_h_per_day", 0, 1), ("SA_m2", 0, 1),
])
def test_each_pathway_factor_appears_exactly_once(name, exponent_ing, exponent_derm):
    """Doubling an input scales each dose by 2**exponent: every factor enters once, in the right place."""
    p = rp.population_parameters("adult", CFG)
    base = inputs(0.05, 2.0, 60.0, 150.0, 0.25, 1.8)  # EF=150 so doubling stays <= 366
    doubled = dict(base, **{name: base[name] * 2})
    b, d = re_.evaluate_risk(base, p), re_.evaluate_risk(doubled, p)
    assert_allclose(d["ADD_ing_mg_per_kg_day"] / b["ADD_ing_mg_per_kg_day"], 2.0 ** exponent_ing, rtol=1e-12)
    assert_allclose(d["ADD_dermal_mg_per_kg_day"] / b["ADD_dermal_mg_per_kg_day"], 2.0 ** exponent_derm,
                    rtol=1e-12)


def test_fixed_constants_enter_once():
    C, IR, EF, ED, ET, SA, BW, AT = 0.05, 2.0, 300.0, 70.0, 0.25, 1.8, 60.0, 25550.0
    base = re_.add_dermal(C, 0.001, EF, ED, ET, SA, 10.0, BW, AT)
    assert_allclose(re_.add_dermal(C, 0.002, EF, ED, ET, SA, 10.0, BW, AT), 2 * base, rtol=1e-12)
    assert_allclose(re_.add_dermal(C, 0.001, EF, ED, ET, SA, 20.0, BW, AT), 2 * base, rtol=1e-12)
    assert_allclose(re_.add_dermal(C, 0.001, EF, ED, ET, SA, 10.0, BW, 2 * AT), base / 2, rtol=1e-12)
    assert_allclose(re_.hazard_quotient(1.0, 0.0003), 1 / 0.0003, rtol=1e-12)


# ---------------------------------------------------------------------------
# Identities and vectorization
# ---------------------------------------------------------------------------

def test_identities_hold_elementwise_on_arrays():
    rng = np.random.default_rng(0)
    n = 1000
    x = inputs(rng.uniform(0, 0.5, n), rng.uniform(0.5, 4, n), rng.uniform(5, 100, n),
               rng.uniform(180, 365, n), rng.uniform(0.13, 0.33, n), rng.uniform(0.3, 2.5, n))
    for pop in rp.POPULATIONS:
        p = rp.population_parameters(pop, CFG)
        out = re_.evaluate_risk(x, p)
        assert all(out[k].shape == (n,) for k in re_.OUTPUTS)
        assert_allclose(out["HI"], out["HQ_ing"] + out["HQ_dermal"], rtol=1e-15)
        assert_allclose(out["ELCR"], out["ADD_ing_cancer_mg_per_kg_day"] * 1.5, rtol=1e-15)
        # Cancer and non-cancer doses differ only by the averaging-time ratio.
        ratio = p.fixed["AT_noncancer_days"] / p.fixed["AT_cancer_days"]
        assert_allclose(out["ADD_ing_cancer_mg_per_kg_day"], out["ADD_ing_mg_per_kg_day"] * ratio, rtol=1e-12)
        scalar = re_.evaluate_risk({k: v[7] for k, v in x.items()}, p)
        assert_allclose(scalar["HI"], out["HI"][7], rtol=1e-15)


def test_zero_concentration_gives_zero_risk():
    out = re_.evaluate_risk(inputs(0.0, 2.0, 60.0, 365.0, 0.2, 1.8), rp.population_parameters("adult", CFG))
    assert all(float(out[k]) == 0.0 for k in re_.OUTPUTS)


@pytest.mark.parametrize("name,bad", [
    ("C_mg_per_L", -0.001), ("BW_kg", 0.0), ("IR_L_per_day", -1.0), ("SA_m2", np.nan),
    ("EF_days_per_year", 400.0), ("ET_h_per_day", 25.0), ("BW_kg", np.inf),
])
def test_rejects_physically_impossible_inputs(name, bad):
    x = dict(inputs(0.05, 2.0, 60.0, 365.0, 0.2, 1.8), **{name: bad})
    with pytest.raises(ValueError):
        re_.evaluate_risk(x, rp.population_parameters("adult", CFG))


def test_missing_input_is_an_error():
    x = inputs(0.05, 2.0, 60.0, 365.0, 0.2, 1.8)
    del x["ET_h_per_day"]
    with pytest.raises(KeyError):
        re_.evaluate_risk(x, rp.population_parameters("adult", CFG))


# ---------------------------------------------------------------------------
# Parameter sets
# ---------------------------------------------------------------------------

def test_fixed_constants_and_averaging_times():
    adult = rp.population_parameters("adult", CFG)
    child = rp.population_parameters("child", CFG)
    assert adult.fixed["ED_years"] == 70 and child.fixed["ED_years"] == 6
    assert adult.fixed["AT_noncancer_days"] == 25_550 and child.fixed["AT_noncancer_days"] == 2_190
    assert CFG["risk_model"]["cancer_averaging_time"] == "lifetime"
    assert adult.fixed["AT_cancer_days"] == 25_550 and child.fixed["AT_cancer_days"] == 25_550
    for p in (adult, child):
        assert p.fixed["Kp_cm_per_h"] == 0.001 and p.fixed["CF_L_per_cm_m2"] == 10.0
        assert p.fixed["RfD_ing_mg_per_kg_day"] == 0.0003 and p.fixed["RfD_dermal_mg_per_kg_day"] == 0.000285
        assert p.fixed["CSF_per_mg_per_kg_day"] == 1.5
        assert set(p.fixed) == set(rp.FIXED_INPUTS)
        assert set(p.stochastic) == set(rp.BASE_PAPER_STOCHASTIC)
        assert set(rp.FIXED_INPUTS).isdisjoint(rp.STOCHASTIC_INPUTS)


def test_adult_elcr_identical_under_both_cancer_AT_conventions():
    x = inputs(0.05, 2.0, 60.0, 345.0, 0.2, 1.8)
    life = re_.evaluate_risk(x, rp.population_parameters("adult", CFG))
    paper = re_.evaluate_risk(x, rp.population_parameters("adult", cfg_with(cancer_averaging_time="exposure_duration")))
    assert_allclose(life["ELCR"], paper["ELCR"], rtol=1e-15)


def test_lognormal_from_arithmetic_moments():
    mu, s = rp.lognormal_from_arithmetic(1.26, 0.66)
    # sigma^2 = ln(1 + (0.66/1.26)^2) = ln(1.2743764) = 0.2424570 ; mu = ln(1.26) - sigma^2/2
    #         = 0.2311117 - 0.1212285 = 0.1098832
    assert_allclose(s ** 2, 0.2424570, rtol=1e-6)
    assert_allclose(mu, 0.1098832, rtol=1e-5)
    with pytest.raises(ValueError):
        rp.lognormal_from_arithmetic(0.0, 1.0)


@pytest.mark.parametrize("pop", rp.POPULATIONS)
def test_sampled_input_moments_match_reported_values(pop):
    """Catches arithmetic-vs-log-space parameterization errors (design doc 21.2)."""
    p = rp.population_parameters(pop, CFG)
    for name, dist in p.stochastic.items():
        d = dist.frozen()
        if dist.family == "Lognormal":
            assert_allclose(d.mean(), dist.reported["reported_mean"], rtol=1e-12, err_msg=name)
            assert_allclose(d.std(), dist.reported["reported_sd"], rtol=1e-12, err_msg=name)
        else:
            a, m, b = dist.params["min"], dist.params["mode"], dist.params["max"]
            assert_allclose(d.mean(), (a + m + b) / 3, rtol=1e-12, err_msg=name)
            assert d.support() == pytest.approx((a, b))
        # Inverse-CDF sampling is finite and positive.
        u = np.random.default_rng(1).uniform(size=20_000)
        draws = dist.ppf(u)
        assert np.all(np.isfinite(draws)) and np.all(draws > 0)


def test_adult_skin_area_is_physically_plausible():
    """Exposed skin cannot exceed total body surface area (US EPA EFH 2011 Table 7-1: adult male P95 about 2.5 m^2)."""
    sa = rp.population_parameters("adult", CFG).stochastic["SA_m2"].frozen()
    assert sa.ppf(0.95) < 2.5
    assert rp.population_parameters("child", CFG).stochastic["SA_m2"].params["max"] <= 0.95


def test_log_space_switch_restores_previous_parameterization():
    p = rp.population_parameters("adult", cfg_with(plus_minus_interpretation="log_space"))
    ir = p.stochastic["IR_L_per_day"]
    assert ir.scipy["s"] == 0.66 and ir.scipy["scale"] == pytest.approx(np.exp(1.26))
    with pytest.raises(ValueError):
        rp.population_parameters("adult", cfg_with(plus_minus_interpretation="bogus"))


def test_triangular_parameterization():
    ef = rp.population_parameters("adult", CFG).stochastic["EF_days_per_year"]
    assert ef.scipy == {"scipy_dist": "triang", "c": pytest.approx(165 / 185), "loc": 180.0, "scale": 185.0}
    et = rp.population_parameters("child", CFG).stochastic["ET_h_per_day"]
    assert et.scipy["c"] == pytest.approx(0.07 / 0.20) and et.scipy["loc"] == 0.13
    with pytest.raises(ValueError):
        rp.input_distribution("x", {"family": "Triangular", "min": 1, "mode": 0.5, "max": 2}, "arithmetic_mean_sd")


def test_parameter_table_covers_every_input():
    rows = rp.parameter_table(CFG)
    for pop in rp.POPULATIONS:
        names = {r["parameter"] for r in rows if r["population"] == pop}
        assert names == set(rp.FIXED_INPUTS) | set(rp.STOCHASTIC_INPUTS)
        assert all(r["unit"] for r in rows)
