"""Project-root-aware path convention.

Every script resolves files relative to PROJECT_ROOT, never the current
working directory. Raw inputs live under data/raw/ and are read-only; all
derived files go under results/<workstream>/.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]

CONFIG_DIR = PROJECT_ROOT / "config"
RUN_CONFIG_PATH = CONFIG_DIR / "run_config.json"

DATA_RAW = PROJECT_ROOT / "data" / "raw"
RESULTS = PROJECT_ROOT / "results"
RESULTS_PROVENANCE = RESULTS / "provenance"
RESULTS_ARSENIC = RESULTS / "arsenic"
RESULTS_BODYWEIGHT = RESULTS / "bodyweight"
RESULTS_SIMULATION = RESULTS / "simulation"
RESULTS_CONVERGENCE = RESULTS / "convergence"
RESULTS_SENSITIVITY = RESULTS / "sensitivity"
RESULTS_TABLES = RESULTS / "tables"
RESULTS_FIGURES = RESULTS / "figures"

OUTPUT_DIRS = (
    RESULTS_PROVENANCE,
    RESULTS_ARSENIC,
    RESULTS_BODYWEIGHT,
    RESULTS_SIMULATION,
    RESULTS_CONVERGENCE,
    RESULTS_SENSITIVITY,
    RESULTS_TABLES,
    RESULTS_FIGURES,
)


def load_run_config() -> dict[str, Any]:
    with RUN_CONFIG_PATH.open(encoding="utf-8") as fh:
        return json.load(fh)


def raw_path(relative: str) -> Path:
    """Resolve a config-relative raw path (e.g. 'data/raw/NationalSurveyData.csv')."""
    path = (PROJECT_ROOT / relative).resolve()
    if DATA_RAW.resolve() not in path.parents:
        raise ValueError(f"{relative!r} is not under data/raw/")
    return path


def ensure_output_dir(path: Path) -> Path:
    """Create an output directory, refusing anything inside data/raw/."""
    resolved = path.resolve()
    if resolved == DATA_RAW.resolve() or DATA_RAW.resolve() in resolved.parents:
        raise ValueError(f"Refusing to write inside data/raw/: {path}")
    resolved.mkdir(parents=True, exist_ok=True)
    return resolved


def rel(path: Path) -> str:
    """Project-relative POSIX path for manifests and reports."""
    return Path(path).resolve().relative_to(PROJECT_ROOT).as_posix()
