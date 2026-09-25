"""Phase 3: adult (STEPS 2018) and child (MICS7 ch.sav) body-weight preprocessing.

Every exclusion is applied as a separate, countable step so the cleaning audit
reconciles exactly: ``rows_before - removed == rows_after`` at each step and the
final row count equals the starting count minus all removals. Raw files are
only read, never written.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from .paths import MICS_CH_PATH, STEPS_PATH

# ---------------------------------------------------------------------------
# Source coding
# ---------------------------------------------------------------------------

STEPS_COLUMNS = ["pid", "psu", "stratum", "division", "urbanrural", "wstep2", "age", "sex", "m11", "m12"]
STEPS_SPECIAL_CODES = {888: "refused weight measurement", 666: "too large for scale"}
ADULT_AGE_RANGE = (18, 69)

MICS_COLUMNS = [
    "HH1", "HH2", "LN", "PSU", "stratum", "HH7A", "HL4", "CAGE",
    "AN8", "AN11", "WAZFLAG", "HAZFLAG", "WHZFLAG", "chweight",
]
MICS_AN8_SPECIAL_CODES = {
    99.3: "CHILD NOT PRESENT AFTER REVISITS",
    99.4: "CHILD REFUSED",
    99.5: "RESPONDENT REFUSED",
    99.6: "OTHER",
}
MICS_AN11_SPECIAL_MIN = 999.0  # AN11 codes 999.4/999.5/999.6 are refusals/other, not heights
MICS_UNEXPECTED_AN8_MIN = 90.0
CHILD_AGE_RANGE = (0, 59)

# Diagnostic-only BMI bounds; records outside are listed for review, never removed.
BMI_REVIEW_BOUNDS = (12.0, 60.0)

ADULT_OUTPUT_COLUMNS = [
    "pid", "age", "sex", "BW_kg", "survey_weight", "psu", "stratum",
    "division", "urbanrural", "height_cm", "bmi_check",
]
CHILD_OUTPUT_COLUMNS = [
    "HH1", "HH2", "LN", "CAGE_months", "sex", "BW_kg", "survey_weight", "district",
    "PSU", "stratum", "height_cm", "WAZFLAG", "HAZFLAG", "WHZFLAG",
]


# ---------------------------------------------------------------------------
# Sequential exclusion bookkeeping
# ---------------------------------------------------------------------------

class CleaningAudit:
    """Applies exclusion steps in order and records a reconcilable count table."""

    def __init__(self, df: pd.DataFrame, population: str, source: str):
        self.df = df
        self.population = population
        self.rows: list[dict] = [{
            "population": population,
            "step": 0,
            "rule": "start",
            "description": f"rows read from {source}",
            "rows_before": len(df),
            "removed": 0,
            "rows_after": len(df),
            "raw_rows_matching": len(df),
        }]
        self._raw = df

    def exclude(self, rule: str, description: str, mask_fn: Callable[[pd.DataFrame], pd.Series]) -> None:
        """Remove rows where ``mask_fn(df)`` is True, counting against current and raw data."""
        mask = mask_fn(self.df).fillna(False).astype(bool)
        raw_matching = int(mask_fn(self._raw).fillna(False).astype(bool).sum())
        before = len(self.df)
        self.df = self.df.loc[~mask]
        self.rows.append({
            "population": self.population,
            "step": len(self.rows),
            "rule": rule,
            "description": description,
            "rows_before": before,
            "removed": int(mask.sum()),
            "rows_after": len(self.df),
            "raw_rows_matching": raw_matching,
        })

    def note(self, rule: str, description: str, count: int) -> None:
        """Record a diagnostic count that does not remove rows."""
        n = len(self.df)
        self.rows.append({
            "population": self.population,
            "step": len(self.rows),
            "rule": rule,
            "description": description + " (diagnostic only; not removed)",
            "rows_before": n,
            "removed": 0,
            "rows_after": n,
            "raw_rows_matching": int(count),
        })

    def table(self) -> pd.DataFrame:
        return pd.DataFrame(self.rows)


# ---------------------------------------------------------------------------
# Adult: STEPS 2018
# ---------------------------------------------------------------------------

def load_steps_bw(path: Path = STEPS_PATH) -> pd.DataFrame:
    """Load the exact STEPS columns needed for adult BW."""
    return pd.read_csv(path, usecols=STEPS_COLUMNS)


def clean_steps_bw(raw: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Apply age, missing-value, survey-weight, and STEPS special-code rules.

    ``BW_kg`` is taken from measured weight ``m12``; ``wstep2`` is only the survey weight.
    Returns ``(clean, audit)``.
    """
    lo, hi = ADULT_AGE_RANGE
    audit = CleaningAudit(raw.copy(), "adult", "bgd2018.csv")
    audit.exclude("age_outside_18_69", f"age missing or outside {lo}-{hi}", lambda d: ~d["age"].between(lo, hi))
    audit.exclude("missing_m12", "measured weight m12 missing", lambda d: d["m12"].isna())
    for code, label in STEPS_SPECIAL_CODES.items():
        audit.exclude(f"m12_code_{code}", f"m12 == {code} ({label})", lambda d, c=code: d["m12"].eq(c))
    audit.exclude("m12_nonpositive", "m12 <= 0", lambda d: d["m12"] <= 0)
    audit.exclude("missing_sex", "sex missing", lambda d: d["sex"].isna())
    audit.exclude("wstep2_missing_or_nonpositive", "survey weight wstep2 missing or <= 0",
                  lambda d: d["wstep2"].isna() | (d["wstep2"] <= 0))

    df = audit.df.copy()
    audit.note("duplicate_pid", "rows sharing a pid", int(df["pid"].duplicated(keep=False).sum()))

    height = df["m11"].where(~df["m11"].isin(list(STEPS_SPECIAL_CODES)) & (df["m11"] > 0))
    df["height_cm"] = height
    df["bmi_check"] = df["m12"] / (height / 100.0) ** 2
    lo_bmi, hi_bmi = BMI_REVIEW_BOUNDS
    audit.note("bmi_outside_review_bounds", f"BMI outside {lo_bmi}-{hi_bmi} kg/m^2",
               int(((df["bmi_check"] < lo_bmi) | (df["bmi_check"] > hi_bmi)).sum()))
    audit.note("height_unavailable", "m11 missing or special code (BMI not computed)", int(height.isna().sum()))

    df = df.rename(columns={"m12": "BW_kg", "wstep2": "survey_weight"})
    clean = df[ADULT_OUTPUT_COLUMNS].reset_index(drop=True)
    return clean, audit.table()


# ---------------------------------------------------------------------------
# Child: MICS7 ch.sav
# ---------------------------------------------------------------------------

def load_mics_bw(path: Path = MICS_CH_PATH) -> pd.DataFrame:
    """Load ch.sav (not bh.sav) with raw numeric codes and keep child anthropometry/design columns."""
    import pyreadstat

    df, _ = pyreadstat.read_sav(str(path), apply_value_formats=False, usecols=MICS_COLUMNS)
    return df


def clean_mics_bw(raw: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Apply AN8 special-code, CAGE, chweight, and WAZFLAG rules.

    ``BW_kg`` is taken from measured weight ``AN8``; ``chweight`` is only the survey weight.
    Height flags are retained as diagnostics and never used to exclude a weight.
    Returns ``(clean, audit)``.
    """
    lo, hi = CHILD_AGE_RANGE
    audit = CleaningAudit(raw.copy(), "child", "ch.sav")
    audit.exclude("missing_AN8", "measured weight AN8 missing", lambda d: d["AN8"].isna())
    for code, label in MICS_AN8_SPECIAL_CODES.items():
        audit.exclude(f"AN8_code_{code}", f"AN8 == {code} ({label})",
                      lambda d, c=code: d["AN8"].sub(c).abs().lt(1e-6))
    audit.exclude("AN8_unexpected_ge_90", "AN8 >= 90 not a documented code (held for codebook review)",
                  lambda d: d["AN8"] >= MICS_UNEXPECTED_AN8_MIN)
    audit.exclude("AN8_nonpositive", "AN8 <= 0", lambda d: d["AN8"] <= 0)
    audit.exclude("missing_CAGE", "age in months CAGE missing", lambda d: d["CAGE"].isna())
    audit.exclude("CAGE_outside_0_59", f"CAGE outside {lo}-{hi} months", lambda d: ~d["CAGE"].between(lo, hi))
    audit.exclude("missing_HL4", "sex HL4 missing", lambda d: d["HL4"].isna())
    audit.exclude("chweight_missing_or_nonpositive", "survey weight chweight missing or <= 0",
                  lambda d: d["chweight"].isna() | (d["chweight"] <= 0))
    audit.exclude("WAZFLAG_missing", "WHO weight-for-age flag missing", lambda d: d["WAZFLAG"].isna())
    audit.exclude("WAZFLAG_1", "WAZFLAG == 1 (WHO weight-for-age error flag)", lambda d: d["WAZFLAG"].eq(1))

    df = audit.df.copy()
    audit.note("duplicate_child_id", "rows sharing (HH1, HH2, LN)",
               int(df.duplicated(["HH1", "HH2", "LN"], keep=False).sum()))
    audit.note("HAZFLAG_1_retained", "HAZFLAG == 1 retained for weight-only fit", int(df["HAZFLAG"].eq(1).sum()))
    audit.note("WHZFLAG_1_retained", "WHZFLAG == 1 retained for weight-only fit", int(df["WHZFLAG"].eq(1).sum()))

    df["height_cm"] = df["AN11"].where(df["AN11"] < MICS_AN11_SPECIAL_MIN)
    df = df.rename(columns={
        "CAGE": "CAGE_months",
        "HL4": "sex",
        "AN8": "BW_kg",
        "chweight": "survey_weight",
        "HH7A": "district",
    })
    clean = df[CHILD_OUTPUT_COLUMNS].reset_index(drop=True)
    return clean, audit.table()


# ---------------------------------------------------------------------------
# Weighted descriptive summary
# ---------------------------------------------------------------------------

def weighted_quantile(x: np.ndarray, w: np.ndarray, probs) -> np.ndarray:
    """Weighted quantiles by interpolating the weighted midpoint plotting positions."""
    x = np.asarray(x, dtype=float)
    w = np.asarray(w, dtype=float)
    order = np.argsort(x, kind="mergesort")
    xs, ws = x[order], w[order]
    q = ws / ws.sum()
    p = np.cumsum(q) - 0.5 * q
    return np.interp(probs, p, xs)


def weighted_summary(x, w, population: str) -> dict:
    x = np.asarray(x, dtype=float)
    w = np.asarray(w, dtype=float)
    q = w / w.sum()
    mean = float(np.sum(q * x))
    sd = float(np.sqrt(np.sum(q * (x - mean) ** 2)))
    p1, p5, p25, p50, p75, p95, p99 = weighted_quantile(x, w, [0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99])
    kish_neff = float(w.sum() ** 2 / np.sum(w ** 2))
    return {
        "population": population,
        "n": int(len(x)),
        "sum_survey_weight": float(w.sum()),
        "kish_effective_n": kish_neff,
        "weighted_mean_kg": mean,
        "weighted_sd_kg": sd,
        "weighted_P1_kg": p1,
        "weighted_P5_kg": p5,
        "weighted_P25_kg": p25,
        "weighted_median_kg": p50,
        "weighted_P75_kg": p75,
        "weighted_P95_kg": p95,
        "weighted_P99_kg": p99,
        "min_kg": float(x.min()),
        "max_kg": float(x.max()),
        "unweighted_mean_kg": float(x.mean()),
        "unweighted_sd_kg": float(x.std(ddof=1)),
    }
