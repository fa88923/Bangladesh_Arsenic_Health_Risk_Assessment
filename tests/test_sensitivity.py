"""Phase 9 tests: Spearman sensitivity, survey bootstrap, 2-D MC, scenarios, and saved outputs."""
import json

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose

from arsenic_hra import monte_carlo as mc
from arsenic_hra import parameter_uncertainty as pu
from arsenic_hra import paths
from arsenic_hra import sensitivity as sn
from arsenic_hra import simulation_inputs as si

CFG = paths.load_run_config()
SIM = CFG["simulation"]
OUT = paths.RESULTS_SENSITIVITY


@pytest.fixture(scope="module")
def model():
    return mc.load_model_inputs(CFG)


# ---------------------------------------------------------------------------
# Spearman
# ---------------------------------------------------------------------------

def test_spearman_is_tie_aware_rank_correlation():
    frame = pd.DataFrame({name: np.arange(10.0) for name in sn.INPUTS})
    frame["BW_kg"] = -np.arange(10.0)
    frame["ET_h_per_day"] = [1, 1, 1, 2, 2, 2, 3, 3, 3, 3]  # ties
    frame["HI"] = frame["ELCR"] = np.arange(10.0) ** 3
    rows = {r["input"]: r for r in sn.spearman_rows(frame, "X", "adult", outputs=("HI",))}
    assert rows["C_mg_per_L"]["spearman_rho"] == pytest.approx(1.0)
    assert rows["BW_kg"]["spearman_rho"] == pytest.approx(-1.0)
    from scipy.stats import pearsonr, rankdata
    expected = pearsonr(rankdata(frame["ET_h_per_day"]), rankdata(frame["HI"]))[0]
    assert rows["ET_h_per_day"]["spearman_rho"] == pytest.approx(expected)


def test_saved_sensitivity_uses_exact_primary_iterations_and_signs_match_equations():
    spear = pd.read_csv(OUT / "sensitivity_spearman.csv")
    assert spear.sign_consistent.all()
    # every input whose correlation is distinguishable from zero has the sign the equations imply
    sig = spear[spear.distinguishable_from_zero & (spear.expected_sign != 0)]
    assert (np.sign(sig.spearman_rho) == sig.expected_sign).all()
    # C and IR are distinguishable from zero for HI in every cell (BW and EF are weak in a few
    # districts where the arsenic spread dominates)
    core = spear[(spear.output == "HI") & spear.input.isin(["C_mg_per_L", "IR_L_per_day"])]
    assert core.distinguishable_from_zero.all()
    frame = pd.read_parquet(paths.RESULTS_SIMULATION / "iterations" / "dhaka_child_N10000_r0.parquet")
    recomputed = {r["input"]: r["spearman_rho"] for r in sn.spearman_rows(frame, "Dhaka", "child", ("HI",))}
    saved = spear[(spear.district == "Dhaka") & (spear.population == "child") & (spear.output == "HI")]
    for _, r in saved.iterrows():
        assert r.spearman_rho == pytest.approx(recomputed[r.input], rel=1e-12)
    # fixed inputs never appear as sensitivity variables
    assert set(spear.input) == set(sn.INPUTS)


def test_arsenic_dominates_hi_sensitivity_everywhere():
    spear = pd.read_csv(OUT / "sensitivity_spearman.csv")
    top = spear[(spear.output == "HI") & (spear["rank"] == 1)]
    assert (top.input == "C_mg_per_L").all()


# ---------------------------------------------------------------------------
# Survey bootstrap
# ---------------------------------------------------------------------------

def test_rao_wu_weights_by_construction():
    rng = np.random.default_rng(0)
    psu = np.repeat(np.arange(6), 5)
    stratum = np.repeat([0, 0, 0, 1, 1, 1], 5)
    w = np.ones(psu.size)
    new = pu.rao_wu_weights(psu, stratum, w, rng)
    for h in (0, 1):
        in_h = stratum == h
        per_psu = pd.Series(new[in_h]).groupby(psu[in_h]).first()
        # n_h - 1 = 2 draws, each worth n_h/(n_h-1) = 1.5, so PSU multipliers are 0, 1.5 or 3 and sum to 3
        assert set(per_psu.round(12)) <= {0.0, 1.5, 3.0}
        assert per_psu.sum() == pytest.approx(3.0)
        # whole PSUs move together (clusters preserved)
        assert pd.Series(new[in_h]).groupby(psu[in_h]).nunique().max() == 1


def test_bw_bootstrap_centres_on_point_estimate(model):
    bw = pd.read_csv(OUT / "bw_survey_bootstrap_replicates.csv")
    assert bw.converged.all()
    assert bw.groupby("population").size().min() == CFG["arsenic"]["bootstrap_replicates"]
    for pop, rec in model.bodyweight.items():
        g = bw[bw.population == pop]
        for k, v in rec["params"].items():
            if k == "lower":
                assert (g[k] == v).all()
                continue
            assert g[k].median() == pytest.approx(v, rel=0.01)


# ---------------------------------------------------------------------------
# Arsenic bootstrap and 2-D MC
# ---------------------------------------------------------------------------

def test_regenerated_arsenic_bootstrap_reproduces_phase4():
    check = pd.read_csv(OUT / "arsenic_bootstrap_reproduction_check.csv")
    assert (check.max_rel_diff <= 1e-5).all() and (check.nan_pattern_mismatches == 0).all()
    assert (check[check.column.str.contains("winner")].max_rel_diff <= 1e-9).all()


def test_two_d_uses_common_random_numbers(model):
    """With the point-estimate parameters as the only outer set, the 2-D run equals the primary run."""
    d, pop = "Khulna", "adult"
    point_sets = [{"replicate": 0, "params": model.arsenic[d]["params"]}]
    bw = pd.DataFrame([{"population": pop, "replicate": 0, **model.bodyweight[pop]["params"]}])
    two_d = sn.two_dimensional(CFG, model, d, pop, point_sets, bw)
    frame = pd.read_parquet(paths.RESULTS_SIMULATION / "iterations" / "khulna_adult_N10000_r0.parquet")
    expected = sn.risk_statistics(frame)
    for k in sn.STAT_NAMES:
        assert two_d.loc[0, k] == pytest.approx(expected[k], rel=1e-12)


def test_uncertainty_intervals_bracket_primary_and_are_complete():
    unc = pd.read_csv(OUT / "parameter_uncertainty_summary.csv")
    assert len(unc) == len(CFG["arsenic"]["selected_districts"]) * len(SIM["populations"]) * len(sn.STAT_NAMES)
    assert (unc.outer_replicates >= 0.95 * CFG["arsenic"]["bootstrap_replicates"]).all()
    key = unc[unc.statistic.isin(["P95_HI", "P_HI_gt_1"])]
    inside = (key.primary_value >= key.param_unc_q025) & (key.primary_value <= key.param_unc_q975)
    assert inside.mean() >= 0.9


# ---------------------------------------------------------------------------
# Scenarios
# ---------------------------------------------------------------------------

def test_mixture_distribution_is_a_valid_cdf():
    rec, info = sn.child_bw_age_extended()
    d = si.frozen(rec["distribution"], rec["params"])
    u = np.linspace(1e-6, 1 - 1e-6, 1001)
    x = d.ppf(u)
    assert np.all(np.diff(x) >= 0) and x.min() >= 1.6
    assert_allclose(d.cdf(x), u, atol=1e-5)
    assert info["extrapolated_60_71_mean_kg"] > info["mean_48_59_kg"] > info["mean_36_47_kg"]
    assert rec["params"]["w1"] == pytest.approx(60 / 72)


def test_scenarios_change_only_their_cells():
    scen = pd.read_csv(OUT / "scenario_results.csv")
    s1 = scen[scen.scenario == "S1_barisal_lognormal"]
    assert set(s1.district) == {"Barisal"}
    assert (s1.ratio_P95_HI > 1).all()   # Lognormal tail is heavier than the selected Gamma
    s4 = scen[scen.scenario.isin(["S4_child_bw_gamma", "S5_child_bw_0_71_months"])]
    assert set(s4.population) == {"child"}
    s5 = scen[scen.scenario == "S5_child_bw_0_71_months"]
    assert (s5.ratio_P95_HI < 1).all()   # heavier 60-71-month children lower the dose per kg
    s6 = scen[scen.scenario == "S6_ir_sa_log_space"]
    assert (s6.ratio_P95_HI > 1).all()   # log-space reading inflates IR


def test_manifest_checks_pass():
    m = json.loads((OUT / "sensitivity_manifest.json").read_text())
    assert all(m["checks"].values())
