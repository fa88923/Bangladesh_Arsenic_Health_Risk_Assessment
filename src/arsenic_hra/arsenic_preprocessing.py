"""BGS arsenic preprocessing (Phase 2).

Implements docs/arsenic_preprocessing_and_distribution_plan.md sections 2-7:
read the raw file without changing it, preserve `As` as `As_raw`, classify
each observation, convert detected values and censoring limits from ug/L to
mg/L, review duplicates, then filter the approved districts.

Censored values are never substituted: for "< 6" the 6 is stored only as the
upper censoring bound, and As_detected_* stays NaN.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

UG_PER_MG = 1000.0

STATUS_DETECTED = "detected"
STATUS_CENSORED = "left_censored"
STATUS_MISSING = "missing"
STATUS_MALFORMED = "malformed"
STATUSES = (STATUS_DETECTED, STATUS_CENSORED, STATUS_MISSING, STATUS_MALFORMED)

BGS_HEADER_LINE = 5  # 1-based file line holding the column names
BGS_FIRST_DATA_LINE = 7  # line 6 is the units row

RETAINED_COLUMNS = [
    "SAMPLE_ID", "SAMPLE_FIELD_ID", "SAMPLE_DATE", "LAT_DEG", "LONG_DEG",
    "YEAR_CONSTRUCTION", "WELL_TYPE", "WELL_DEPTH", "DIVISION", "DISTRICT",
    "THANA", "UNION", "MOUZA", "GEOCODE",
]

DERIVED_COLUMNS = [
    "As_raw", "As_status", "As_censored",
    "As_bound_ugL", "As_detected_ugL", "As_model_value_ugL",
    "As_bound_mgL", "As_detected_mgL", "As_model_value_mgL",
]

OUTPUT_COLUMNS = ["source_line"] + RETAINED_COLUMNS + ["SAMPLE_DATE_iso"] + DERIVED_COLUMNS

_NUMBER = r"[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?"
_DETECTED_RE = re.compile(rf"^({_NUMBER})$")
_CENSORED_RE = re.compile(rf"^<\s*({_NUMBER})$")


# --------------------------------------------------------------------------
# Reading
# --------------------------------------------------------------------------
def read_bgs_raw(path: Path) -> tuple[pd.DataFrame, dict[str, str]]:
    """Read the BGS CSV as strings, dropping only the verified units row.

    Returns the sample rows (with `source_line`, the 1-based line in the raw
    file) and the units mapping taken from that row.
    """
    df = pd.read_csv(path, header=BGS_HEADER_LINE - 1, dtype=str, keep_default_na=False)
    units_row = df.iloc[0]
    if units_row.get("As", "").strip().lower() != "ug/l" or units_row.get("SAMPLE_ID", "") != "":
        raise ValueError("BGS units row not found where expected; refusing to guess the layout.")
    units = {col: val for col, val in units_row.items() if val}
    samples = df.iloc[1:].reset_index(drop=True)
    samples.insert(0, "source_line", np.arange(BGS_FIRST_DATA_LINE, BGS_FIRST_DATA_LINE + len(samples)))
    return samples, units


# --------------------------------------------------------------------------
# Arsenic classification and unit conversion
# --------------------------------------------------------------------------
def parse_arsenic_value(raw: object) -> tuple[str, float, float]:
    """Classify one raw arsenic string.

    Returns (status, detected_ugL, bound_ugL); the unused value is NaN.
    """
    if raw is None or (isinstance(raw, float) and np.isnan(raw)):
        return STATUS_MISSING, np.nan, np.nan
    text = str(raw).strip()
    if text == "":
        return STATUS_MISSING, np.nan, np.nan
    if m := _CENSORED_RE.match(text):
        bound = float(m.group(1))
        if bound <= 0:
            return STATUS_MALFORMED, np.nan, np.nan
        return STATUS_CENSORED, np.nan, bound
    if m := _DETECTED_RE.match(text):
        return STATUS_DETECTED, float(m.group(1)), np.nan
    # Anything else (negative numbers, text, ">", "<" with no number) needs review.
    return STATUS_MALFORMED, np.nan, np.nan


def ugL_to_mgL(values):
    """Convert micrograms per litre to milligrams per litre."""
    return values / UG_PER_MG


def classify_arsenic(df: pd.DataFrame, column: str = "As") -> pd.DataFrame:
    """Add the As_* derived columns. The original column is copied verbatim to As_raw."""
    out = df.copy()
    out["As_raw"] = out[column]
    parsed = [parse_arsenic_value(v) for v in out[column]]
    out["As_status"] = [p[0] for p in parsed]
    out["As_detected_ugL"] = np.array([p[1] for p in parsed], dtype=float)
    out["As_bound_ugL"] = np.array([p[2] for p in parsed], dtype=float)

    status = out["As_status"]
    out["As_censored"] = pd.array(
        np.where(status.isin([STATUS_DETECTED, STATUS_CENSORED]), status == STATUS_CENSORED, pd.NA),
        dtype="boolean",
    )
    out["As_model_value_ugL"] = out["As_detected_ugL"].where(status == STATUS_DETECTED, out["As_bound_ugL"])

    for stem in ("bound", "detected", "model_value"):
        out[f"As_{stem}_mgL"] = ugL_to_mgL(out[f"As_{stem}_ugL"])
    return out


def add_typed_metadata(df: pd.DataFrame) -> pd.DataFrame:
    """Coerce coordinates/depth to numbers and add an ISO date; raw strings are kept."""
    out = df.copy()
    for col in ("LAT_DEG", "LONG_DEG", "WELL_DEPTH", "YEAR_CONSTRUCTION"):
        out[col] = pd.to_numeric(out[col].replace("", np.nan), errors="coerce")
    dates = pd.to_datetime(out["SAMPLE_DATE"].replace("", np.nan), format="%d/%m/%Y", errors="coerce")
    out["SAMPLE_DATE_iso"] = dates.dt.strftime("%Y-%m-%d")
    return out


# --------------------------------------------------------------------------
# Quality control
# --------------------------------------------------------------------------
def duplicate_review(df: pd.DataFrame) -> pd.DataFrame:
    """Rows sharing a SAMPLE_ID, a full record, or a coordinate pair.

    Nothing is deleted: every row is listed for human review with a group key.
    """
    review_cols = ["source_line", "SAMPLE_ID", "SAMPLE_FIELD_ID", "SAMPLE_DATE", "LAT_DEG", "LONG_DEG",
                   "WELL_DEPTH", "WELL_TYPE", "DISTRICT", "THANA", "UNION", "As_raw"]
    frames = []

    dup_id = df["SAMPLE_ID"].duplicated(keep=False) & (df["SAMPLE_ID"] != "")
    if dup_id.any():
        part = df.loc[dup_id, review_cols].copy()
        part.insert(0, "group_key", part["SAMPLE_ID"])
        part.insert(0, "review_type", "duplicate_sample_id")
        frames.append(part)

    record_cols = [c for c in df.columns if c != "source_line"]
    dup_row = df.duplicated(subset=record_cols, keep=False)
    if dup_row.any():
        part = df.loc[dup_row, review_cols].copy()
        part.insert(0, "group_key", df.loc[dup_row, record_cols].astype(str).agg("|".join, axis=1))
        part.insert(0, "review_type", "duplicate_full_row")
        frames.append(part)

    has_xy = df["LAT_DEG"].notna() & df["LONG_DEG"].notna()
    dup_xy = has_xy & df[["LAT_DEG", "LONG_DEG"]].duplicated(keep=False)
    if dup_xy.any():
        part = df.loc[dup_xy, review_cols].copy()
        part.insert(0, "group_key", part["LAT_DEG"].map("{:.4f}".format) + "," + part["LONG_DEG"].map("{:.4f}".format))
        grp = part.groupby("group_key")
        part["group_size"] = grp["SAMPLE_ID"].transform("size")
        part["group_n_districts"] = grp["DISTRICT"].transform("nunique")
        part["group_n_depths"] = grp["WELL_DEPTH"].transform("nunique")
        part["group_n_dates"] = grp["SAMPLE_DATE"].transform("nunique")
        part["group_n_as_values"] = grp["As_raw"].transform("nunique")
        part["reviewer_note"] = np.where(
            (part["group_n_depths"] > 1) | (part["group_n_dates"] > 1),
            "distinct wells or visits at one location (different depth/date); retain",
            "same location, depth and date; retained pending reviewer decision",
        )
        part.insert(0, "review_type", "repeated_coordinates")
        frames.append(part)

    columns = ["review_type", "group_key"] + review_cols + [
        "group_size", "group_n_districts", "group_n_depths", "group_n_dates", "group_n_as_values", "reviewer_note"]
    if not frames:
        return pd.DataFrame(columns=columns)
    return pd.concat(frames, ignore_index=True).reindex(columns=columns).sort_values(
        ["review_type", "group_key", "source_line"], kind="stable").reset_index(drop=True)


def malformed_review(df: pd.DataFrame) -> pd.DataFrame:
    """Rows whose arsenic or key metadata need a human decision. Nothing is removed."""
    issues = []

    def add(mask, issue):
        if mask.any():
            part = df.loc[mask, ["source_line", "SAMPLE_ID", "DISTRICT", "As_raw", "As_status",
                                 "LAT_DEG", "LONG_DEG", "SAMPLE_DATE"]].copy()
            part.insert(0, "issue", issue)
            issues.append(part)

    add(df["As_status"] == STATUS_MALFORMED, "arsenic_string_unparseable")
    add(df["As_status"] == STATUS_MISSING, "arsenic_missing")
    add(df["As_detected_ugL"] < 0, "arsenic_negative")
    add(df["DISTRICT"].str.strip() == "", "district_missing")
    add(~np.isfinite(df["LAT_DEG"]) | ~np.isfinite(df["LONG_DEG"]), "coordinates_missing")
    add(df["SAMPLE_DATE_iso"].isna(), "sample_date_unparseable")
    add(df["DISTRICT"] != df["DISTRICT"].str.strip(), "district_whitespace")

    columns = ["issue", "source_line", "SAMPLE_ID", "DISTRICT", "As_raw", "As_status", "LAT_DEG", "LONG_DEG", "SAMPLE_DATE"]
    if not issues:
        return pd.DataFrame(columns=columns)
    return pd.concat(issues, ignore_index=True)[columns]


# --------------------------------------------------------------------------
# District filtering and summaries
# --------------------------------------------------------------------------
def filter_districts(df: pd.DataFrame, districts: list[str]) -> pd.DataFrame:
    """Exact string match on DISTRICT; every requested name must exist in the source."""
    present = set(df["DISTRICT"])
    absent = [d for d in districts if d not in present]
    if absent:
        raise ValueError(f"Selected districts not found verbatim in DISTRICT: {absent}")
    return df[df["DISTRICT"].isin(districts)].copy()


def district_summary(df: pd.DataFrame, districts: list[str], high_censoring_pct: float) -> pd.DataFrame:
    rows = []
    for district in districts:
        d = df[df["DISTRICT"] == district]
        det = d.loc[d["As_status"] == STATUS_DETECTED, "As_detected_mgL"]
        cen = d.loc[d["As_status"] == STATUS_CENSORED, "As_bound_mgL"]
        usable = len(det) + len(cen)
        cens_pct = 100.0 * len(cen) / usable if usable else np.nan
        limits = sorted(cen.unique())
        rows.append({
            "DISTRICT": district,
            "DIVISION": ";".join(sorted(d["DIVISION"].unique())),
            "total_n": len(d),
            "detected_n": len(det),
            "censored_n": len(cen),
            "missing_n": int((d["As_status"] == STATUS_MISSING).sum()),
            "malformed_n": int((d["As_status"] == STATUS_MALFORMED).sum()),
            "censoring_percent": round(cens_pct, 2),
            "number_of_distinct_censoring_limits": len(limits),
            "censoring_limits_mgL": ";".join(f"{x:g}" for x in limits),
            "censored_at_0.0005_mgL_n": int((cen == 0.0005).sum()),
            "censored_at_0.006_mgL_n": int((cen == 0.006).sum()),
            "detected_below_max_limit_n": int((det < cen.max()).sum()) if len(cen) else 0,
            "minimum_detected_As_mgL": det.min() if len(det) else np.nan,
            "median_detected_As_mgL": det.median() if len(det) else np.nan,
            "maximum_detected_As_mgL": det.max() if len(det) else np.nan,
            "detected_gt_0.01_mgL_n": int((det > 0.01).sum()),
            "detected_gt_0.05_mgL_n": int((det > 0.05).sum()),
            "high_censoring_warning": bool(cens_pct >= high_censoring_pct),
        })
    return pd.DataFrame(rows)


def reconciliation_audit(raw: pd.DataFrame, classified: pd.DataFrame, selected: pd.DataFrame,
                         districts: list[str], raw_line_count: int) -> pd.DataFrame:
    """Countable flow from raw file lines to the selected-district file. Every check must pass."""
    status_counts = classified["As_status"].value_counts()
    sel_counts = selected["As_status"].value_counts()
    per_district_raw = raw["DISTRICT"].value_counts()
    rows = [
        ("raw_file_lines", raw_line_count, "physical lines in the raw CSV"),
        ("preamble_lines", BGS_HEADER_LINE - 1, "title/release lines above the header"),
        ("header_lines", 1, "column-name line"),
        ("units_rows", 1, "units row removed (not a sample)"),
        ("raw_sample_rows", len(raw), "rows after removing only preamble, header and units row"),
    ]
    for s in STATUSES:
        rows.append((f"all_{s}_n", int(status_counts.get(s, 0)), f"all districts, As_status == {s}"))
    rows.append(("all_status_total", int(status_counts.sum()), "must equal raw_sample_rows"))
    rows.append(("selected_districts_n", len(districts), "approved district list"))
    rows.append(("selected_rows_expected", int(sum(per_district_raw.get(d, 0) for d in districts)),
                 "raw DISTRICT counts summed over the approved list"))
    rows.append(("selected_rows_output", len(selected), "rows written to the preprocessed file"))
    for s in STATUSES:
        rows.append((f"selected_{s}_n", int(sel_counts.get(s, 0)), f"selected districts, As_status == {s}"))
    rows.append(("excluded_other_district_rows", len(raw) - len(selected), "rows outside the approved districts"))
    audit = pd.DataFrame(rows, columns=["step", "count", "description"])

    lookup = dict(zip(audit["step"], audit["count"]))
    checks = {
        "lines_reconcile": lookup["raw_file_lines"] == lookup["preamble_lines"] + lookup["header_lines"]
        + lookup["units_rows"] + lookup["raw_sample_rows"],
        "status_reconcile": lookup["all_status_total"] == lookup["raw_sample_rows"],
        "selected_reconcile": lookup["selected_rows_output"] == lookup["selected_rows_expected"],
        "selected_status_reconcile": sum(lookup[f"selected_{s}_n"] for s in STATUSES) == lookup["selected_rows_output"],
        "exclusion_reconcile": lookup["selected_rows_output"] + lookup["excluded_other_district_rows"]
        == lookup["raw_sample_rows"],
        "as_raw_preserved": bool((selected["As_raw"].to_numpy() == raw.set_index("source_line").loc[
            selected["source_line"], "As"].to_numpy()).all()),
        "no_censored_in_detected": bool(selected.loc[selected["As_censored"].fillna(False).astype(bool),
                                                     "As_detected_ugL"].isna().all()),
        "district_names_verbatim": set(selected["DISTRICT"]) == set(districts),
    }
    check_rows = pd.DataFrame(
        [(f"check_{k}", int(v), "1 = pass, 0 = FAIL") for k, v in checks.items()],
        columns=["step", "count", "description"],
    )
    return pd.concat([audit, check_rows], ignore_index=True)


def count_lines(path: Path) -> int:
    with path.open("rb") as fh:
        data = fh.read()
    return data.count(b"\n") + (0 if data.endswith(b"\n") or not data else 1)
