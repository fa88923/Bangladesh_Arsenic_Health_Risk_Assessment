import numpy as np
import pandas as pd
import pytest

from arsenic_hra import bodyweight_preprocessing as pp


def assert_audit_reconciles(audit: pd.DataFrame, n_clean: int) -> None:
    assert (audit["rows_before"] - audit["removed"] == audit["rows_after"]).all()
    assert (audit["rows_before"].iloc[1:].to_numpy() == audit["rows_after"].iloc[:-1].to_numpy()).all()
    assert audit["rows_after"].iloc[0] - audit["removed"].sum() == n_clean
    assert audit["rows_after"].iloc[-1] == n_clean


@pytest.fixture
def steps_raw():
    rows = [
        # pid, age, m12, wstep2, m11, sex
        ("ok1", 30, 55.0, 1000.0, 160.0, 0),
        ("ok2", 45, 162.0, 250.0, 888.0, 1),   # legitimate extreme kept; height code masked
        ("young", 17, 50.0, 1000.0, 150.0, 0),
        ("miss", 40, np.nan, 1000.0, np.nan, 1),
        ("ref", 40, 888.0, 1000.0, 150.0, 1),
        ("big", 40, 666.0, 1000.0, 150.0, 0),
        ("zero", 40, 0.0, 1000.0, 150.0, 0),
        ("now", 40, 60.0, 0.0, 150.0, 0),
    ]
    df = pd.DataFrame(rows, columns=["pid", "age", "m12", "wstep2", "m11", "sex"])
    df["psu"], df["stratum"], df["division"], df["urbanrural"] = 1, 1, 1, "Urban"
    return df


def test_steps_cleaning_counts_and_codes(steps_raw):
    clean, audit = pp.clean_steps_bw(steps_raw)
    assert list(clean["pid"]) == ["ok1", "ok2"]
    removed = audit.set_index("rule")["removed"]
    assert removed["age_outside_18_69"] == 1
    assert removed["missing_m12"] == 1
    assert removed["m12_code_888"] == 1
    assert removed["m12_code_666"] == 1
    assert removed["m12_nonpositive"] == 1
    assert removed["wstep2_missing_or_nonpositive"] == 1
    assert_audit_reconciles(audit, len(clean))


def test_steps_bw_comes_from_m12_not_wstep2(steps_raw):
    clean, _ = pp.clean_steps_bw(steps_raw)
    assert clean["BW_kg"].tolist() == [55.0, 162.0]
    assert clean["survey_weight"].tolist() == [1000.0, 250.0]
    assert np.isnan(clean.loc[1, "height_cm"])  # 888 height is not a height
    assert clean.loc[0, "bmi_check"] == pytest.approx(55.0 / 1.6 ** 2)
    assert list(clean.columns) == pp.ADULT_OUTPUT_COLUMNS


@pytest.fixture
def mics_raw():
    rows = [
        # AN8, CAGE, chweight, WAZFLAG, HAZFLAG, AN11
        (10.5, 24, 1.2, 0, 1, 85.0),     # height flag retained
        (3.2, 0, 0.8, 0, 0, 999.5),      # height refusal code masked
        (np.nan, np.nan, 0.0, np.nan, np.nan, np.nan),
        (99.3, 12, 1.0, 1, 1, 80.0),
        (99.4, 12, 1.0, 1, 1, 80.0),
        (99.5, 12, 1.0, 1, 1, 80.0),
        (99.6, 12, 1.0, 1, 1, 80.0),
        (95.0, 12, 1.0, 0, 0, 80.0),     # undocumented >= 90 value
        (12.0, 60, 1.0, 0, 0, 90.0),
        (12.0, 30, 0.0, 0, 0, 90.0),
        (30.0, 6, 1.0, 1, 0, 70.0),      # WHO weight-for-age error
    ]
    df = pd.DataFrame(rows, columns=["AN8", "CAGE", "chweight", "WAZFLAG", "HAZFLAG", "AN11"])
    df["HH1"] = df["PSU"] = 1
    df["HH2"] = 1
    df["LN"] = range(len(df))
    df["stratum"], df["HH7A"], df["HL4"], df["WHZFLAG"] = 10, 26, 1, 0
    return df


def test_mics_cleaning_counts_and_codes(mics_raw):
    clean, audit = pp.clean_mics_bw(mics_raw)
    assert clean["BW_kg"].tolist() == [10.5, 3.2]
    removed = audit.set_index("rule")["removed"]
    assert removed["missing_AN8"] == 1
    for code in (99.3, 99.4, 99.5, 99.6):
        assert removed[f"AN8_code_{code}"] == 1
    assert removed["AN8_unexpected_ge_90"] == 1
    assert removed["CAGE_outside_0_59"] == 1
    assert removed["chweight_missing_or_nonpositive"] == 1
    assert removed["WAZFLAG_1"] == 1
    assert audit.set_index("rule").loc["HAZFLAG_1_retained", "raw_rows_matching"] == 1
    assert_audit_reconciles(audit, len(clean))


def test_mics_bw_comes_from_an8_not_chweight(mics_raw):
    clean, _ = pp.clean_mics_bw(mics_raw)
    assert clean["survey_weight"].tolist() == [1.2, 0.8]
    assert clean["height_cm"].iloc[0] == 85.0
    assert np.isnan(clean["height_cm"].iloc[1])
    assert (clean["WAZFLAG"] == 0).all()
    assert list(clean.columns) == pp.CHILD_OUTPUT_COLUMNS


def test_weighted_quantile_and_summary():
    x = np.array([1.0, 2.0, 3.0, 4.0])
    assert pp.weighted_quantile(x, np.ones(4), 0.5) == pytest.approx(2.5)
    # Doubling the weight of 4 moves the median up; scaling all weights changes nothing.
    w = np.array([1.0, 1.0, 1.0, 3.0])
    assert pp.weighted_quantile(x, w, 0.5) > 2.5
    assert pp.weighted_quantile(x, 7 * w, 0.5) == pytest.approx(pp.weighted_quantile(x, w, 0.5))
    s = pp.weighted_summary(x, w, "t")
    assert s["weighted_mean_kg"] == pytest.approx((1 + 2 + 3 + 12) / 6)
    assert s["kish_effective_n"] == pytest.approx(36 / 12)
