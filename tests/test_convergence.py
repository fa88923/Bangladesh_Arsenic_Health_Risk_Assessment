"""Phase 8 tests: convergence statistics, criterion, and saved outputs."""
import json

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose

from arsenic_hra import convergence as cv
from arsenic_hra import paths

CFG = paths.load_run_config()
SIM = CFG["simulation"]


def synthetic_runs(values_by_N: dict, metric="P95_HI"):
    rows = []
    for N, vals in values_by_N.items():
        for r, v in enumerate(vals):
            row = {"district": "X", "population": "adult", "N": N, "replicate": r}
            row.update({m: 1.0 for m in cv.TRACKED})
            row[metric] = v
            rows.append(row)
    return pd.DataFrame(rows)


def test_successive_relative_change_and_cv_by_hand():
    runs = synthetic_runs({1000: [9.0, 11.0], 5000: [10.0, 10.2], 10000: [10.1, 10.1]})
    agg = cv.aggregate(runs)
    g = agg[agg.metric == "P95_HI"].set_index("N")
    # T_bar: 10.0, 10.1, 10.1.  Delta(5000) = |10.1 - 10.0| / 10.1 = 0.990%;  Delta(10000) = 0
    assert_allclose(g.loc[5000, "successive_rel_change_pct"], 100 * 0.1 / 10.1)
    assert g.loc[10000, "successive_rel_change_pct"] == pytest.approx(0.0)
    assert np.isnan(g.loc[1000, "successive_rel_change_pct"])
    # CV(1000) = sd(9, 11) / 10 = 1.41421 / 10 = 14.14% -> first N can never be adequate
    assert_allclose(g.loc[1000, "between_seed_cv_pct"], 100 * np.sqrt(2) / 10)
    assert not g.loc[1000, "adequate"] and g.loc[5000, "adequate"] and g.loc[10000, "adequate"]


def test_probability_criterion_is_absolute():
    runs = synthetic_runs({1000: [0.02, 0.03], 5000: [0.021, 0.024], 10000: [0.022, 0.023]}, "P_HI_gt_1")
    g = cv.aggregate(runs).query("metric == 'P_HI_gt_1'").set_index("N")
    # relative change would be large (~4%), but the absolute change (< 0.01) and SD (< 0.01) pass
    assert g.loc[5000, "adequate"] and g.loc[10000, "adequate"]
    assert_allclose(g.loc[10000, "binomial_mc_se"], np.sqrt(0.0225 * 0.9775 / 10000))


def test_smallest_adequate_requires_all_larger_N_to_pass():
    runs = synthetic_runs({1000: [10, 10], 5000: [10, 10.01], 10000: [12, 12.1], 20000: [12, 12.05]})
    concl = cv.conclusions(cv.aggregate(runs), 10000).set_index("metric")
    # 5000 passes but 10000 jumps by 17%, so the smallest adequate N is 20000
    assert concl.loc["P95_HI", "smallest_adequate_N"] == 20000
    assert not concl.loc["P95_HI", "adequate_at_primary_N"]


def test_run_statistics_are_direct():
    frame = pd.DataFrame({"HI": [0.5, 1.5, 2.5, 3.0], "ELCR": [1e-5, 2e-5, 3e-5, 4e-5]})
    s = cv.run_statistics(frame)
    assert s["P_HI_gt_1"] == 0.75 and s["P_HI_gt_2"] == 0.5 and s["mean_HI"] == pytest.approx(1.875)


# ---------------------------------------------------------------------------
# Saved outputs
# ---------------------------------------------------------------------------

def test_design_is_complete_and_primary_reproduced():
    runs = pd.read_csv(paths.RESULTS_CONVERGENCE / "convergence_runs.csv")
    n_expected = (len(CFG["arsenic"]["selected_districts"]) * len(SIM["populations"])
                  * len(SIM["convergence_N"]) * SIM["replicate_seeds_per_N"])
    assert len(runs) == n_expected and runs.run_id.is_unique
    m = json.loads((paths.RESULTS_CONVERGENCE / "convergence_manifest.json").read_text())
    assert m["checks"]["N10000_r0_reproduces_phase7_primary"] is True
    # replicate-0 at the primary N equals the Phase 7 summary
    summ = pd.read_csv(paths.RESULTS_SIMULATION / "simulation_summary.csv")
    hi = summ[summ.metric == "HI"].set_index("run_id")
    r0 = runs[(runs.N == SIM["primary_N"]) & (runs.replicate == 0)].set_index("run_id")
    assert_allclose(r0.loc[hi.index, "P95_HI"], hi["P95"], rtol=1e-12)
    assert_allclose(r0.loc[hi.index, "P_HI_gt_1"], hi["P_gt_1"], rtol=1e-12)


def test_mc_error_follows_inverse_sqrt_N():
    """Median between-seed CV of P95 HI falls roughly as N^-1/2 from 1k to 20k (factor sqrt(20) = 4.47)."""
    agg = pd.read_csv(paths.RESULTS_CONVERGENCE / "convergence_summary.csv")
    med = agg[agg.metric == "P95_HI"].groupby("N")["between_seed_cv_pct"].median()
    assert 2.5 < med.loc[1000] / med.loc[20000] < 8.0


def test_headline_probability_adequate_everywhere_at_primary_N():
    concl = pd.read_csv(paths.RESULTS_CONVERGENCE / "convergence_conclusions.csv")
    assert concl[concl.metric == "P_HI_gt_1"].adequate_at_primary_N.all()
