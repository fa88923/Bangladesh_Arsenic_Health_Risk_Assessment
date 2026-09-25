"""Project-root-aware paths for the body-weight track."""

from __future__ import annotations

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

RAW_DIR = PROJECT_ROOT / "data" / "raw"
STEPS_PATH = RAW_DIR / "bgd2018.csv"
MICS_CH_PATH = (
    RAW_DIR
    / "BGD_2025_MICS7_v01_M"
    / "BGD_2025_MICS7_v01_M"
    / "BGD_2025_MICS7_Datasets"
    / "Bangladesh MICS7 Datasets"
    / "Bangladesh MICS7 Datasets"
    / "Bangladesh MICS7 SPSS Datasets"
    / "ch.sav"
)

RESULTS_DIR = PROJECT_ROOT / "results" / "bodyweight"
FIGURES_DIR = RESULTS_DIR / "figures"
