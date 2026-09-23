"""Unit tests: arsenic-string classification, limit extraction and unit conversion."""
import math

import numpy as np
import pandas as pd
import pytest

from arsenic_hra import arsenic_preprocessing as ap


@pytest.mark.parametrize("raw, status, detected, bound", [
    ("42", ap.STATUS_DETECTED, 42.0, None),
    (" 0.5 ", ap.STATUS_DETECTED, 0.5, None),
    ("1660", ap.STATUS_DETECTED, 1660.0, None),
    (".7", ap.STATUS_DETECTED, 0.7, None),
    ("< 6", ap.STATUS_CENSORED, None, 6.0),
    ("<6", ap.STATUS_CENSORED, None, 6.0),
    ("< 0.5", ap.STATUS_CENSORED, None, 0.5),
    ("", ap.STATUS_MISSING, None, None),
    ("   ", ap.STATUS_MISSING, None, None),
    (None, ap.STATUS_MISSING, None, None),
    (float("nan"), ap.STATUS_MISSING, None, None),
    ("<", ap.STATUS_MALFORMED, None, None),
    ("< 0", ap.STATUS_MALFORMED, None, None),
    ("-3", ap.STATUS_MALFORMED, None, None),
    ("> 500", ap.STATUS_MALFORMED, None, None),
    ("n.d.", ap.STATUS_MALFORMED, None, None),
    ("12 ug/l", ap.STATUS_MALFORMED, None, None),
])
def test_parse_arsenic_value(raw, status, detected, bound):
    got_status, got_det, got_bound = ap.parse_arsenic_value(raw)
    assert got_status == status
    assert (math.isnan(got_det) if detected is None else got_det == detected)
    assert (math.isnan(got_bound) if bound is None else got_bound == bound)


def test_unit_conversion_fixed_examples():
    np.testing.assert_allclose(ap.ugL_to_mgL(np.array([6.0, 0.5, 10.0, 50.0, 1660.0])),
                               [0.006, 0.0005, 0.01, 0.05, 1.66], rtol=0, atol=1e-15)


def _frame(values):
    return pd.DataFrame({"As": values})


def test_censored_row_keeps_bound_only():
    row = ap.classify_arsenic(_frame(["< 6"])).iloc[0]
    assert row["As_raw"] == "< 6"
    assert bool(row["As_censored"]) is True
    assert row["As_bound_ugL"] == 6.0
    assert np.isnan(row["As_detected_ugL"]) and np.isnan(row["As_detected_mgL"])
    assert row["As_model_value_ugL"] == 6.0
    assert row["As_bound_mgL"] == pytest.approx(0.006)
    assert row["As_model_value_mgL"] == pytest.approx(0.006)


def test_detected_row():
    row = ap.classify_arsenic(_frame(["42"])).iloc[0]
    assert bool(row["As_censored"]) is False
    assert row["As_detected_ugL"] == 42.0
    assert np.isnan(row["As_bound_ugL"])
    assert row["As_model_value_mgL"] == pytest.approx(0.042)


def test_missing_and_malformed_have_no_values_and_na_censor_flag():
    out = ap.classify_arsenic(_frame(["", "abc"]))
    assert out["As_censored"].isna().all()
    assert out[["As_detected_ugL", "As_bound_ugL", "As_model_value_mgL"]].isna().all().all()


def test_as_raw_is_verbatim_copy():
    values = ["< 6", " 13", "0.5", "", "x"]
    out = ap.classify_arsenic(_frame(values))
    assert out["As_raw"].tolist() == values
    assert out["As"].tolist() == values


def test_classify_does_not_mutate_input():
    df = _frame(["< 6", "13"])
    before = df.copy()
    ap.classify_arsenic(df)
    pd.testing.assert_frame_equal(df, before)


def test_filter_districts_requires_exact_names():
    df = pd.DataFrame({"DISTRICT": ["Dhaka", "Bogra", "Feni"]})
    assert ap.filter_districts(df, ["Dhaka", "Bogra"])["DISTRICT"].tolist() == ["Dhaka", "Bogra"]
    with pytest.raises(ValueError):
        ap.filter_districts(df, ["dhaka"])
    with pytest.raises(ValueError):
        ap.filter_districts(df, ["Bogura"])


def _meta_frame():
    df = pd.DataFrame({
        "source_line": [7, 8, 9, 10],
        "SAMPLE_ID": ["A", "B", "C", "D"],
        "SAMPLE_FIELD_ID": ["a", "b", "c", "d"],
        "SAMPLE_DATE": ["01/02/1998", "01/02/1998", "03/02/1998", "04/02/1998"],
        "LAT_DEG": ["23.1", "23.1", "24.0", ""],
        "LONG_DEG": ["90.1", "90.1", "91.0", ""],
        "YEAR_CONSTRUCTION": ["1990"] * 4,
        "WELL_TYPE": ["TW"] * 4,
        "WELL_DEPTH": ["10", "50", "20", "30"],
        "DIVISION": ["Dhaka"] * 4,
        "DISTRICT": ["Dhaka"] * 4,
        "THANA": ["t"] * 4, "UNION": ["u"] * 4, "MOUZA": ["m"] * 4, "GEOCODE": ["1"] * 4,
        "As": ["< 6", "12", "0.5", "7"],
    })
    return ap.add_typed_metadata(ap.classify_arsenic(df))


def test_repeated_coordinates_are_reviewed_not_removed():
    df = _meta_frame()
    review = ap.duplicate_review(df)
    xy = review[review["review_type"] == "repeated_coordinates"]
    assert sorted(xy["SAMPLE_ID"]) == ["A", "B"]
    assert (xy["group_n_depths"] == 2).all()
    assert len(df) == 4  # nothing dropped


def test_malformed_review_flags_missing_coordinates():
    review = ap.malformed_review(_meta_frame())
    assert review.loc[review["issue"] == "coordinates_missing", "SAMPLE_ID"].tolist() == ["D"]


def test_sample_date_parsed_day_first():
    assert _meta_frame()["SAMPLE_DATE_iso"].tolist()[:3] == ["1998-02-01", "1998-02-01", "1998-02-03"]


def test_district_summary_counts():
    df = _meta_frame()
    s = ap.district_summary(df, ["Dhaka"], high_censoring_pct=50.0).iloc[0]
    assert (s["total_n"], s["detected_n"], s["censored_n"]) == (4, 3, 1)
    assert s["censoring_percent"] == 25.0
    assert s["censoring_limits_mgL"] == "0.006"
    assert s["detected_below_max_limit_n"] == 1  # the 0.5 ug/L detection is below the 6 ug/L limit
    assert not s["high_censoring_warning"]
