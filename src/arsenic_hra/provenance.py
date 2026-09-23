"""Raw-file inventory, hashing and environment capture (Phase 1)."""
from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from .paths import DATA_RAW, PROJECT_ROOT, rel

INVENTORY_COLUMNS = [
    "relative_path",
    "file_name",
    "size_bytes",
    "sha256",
    "source",
    "release",
    "role",
    "used_in_analysis",
    "notes",
]

# Known provenance for the raw files. Anything under data/raw/ that is not
# listed here is still hashed and reported with role "unregistered".
_BGS = "DPHE/BGS/DFID National Hydrochemical Survey of Bangladesh (groundwater arsenic)"
_STEPS = "WHO STEPwise survey, Bangladesh 2018 (STEPS) microdata"
_MICS = "UNICEF Bangladesh MICS7 (BGD_2025_MICS7_v01_M), SPSS release"
_MICS_DIR = "BGD_2025_MICS7_v01_M/BGD_2025_MICS7_v01_M/BGD_2025_MICS7_Datasets"
_MICS_SPSS = f"{_MICS_DIR}/Bangladesh MICS7 Datasets/Bangladesh MICS7 Datasets/Bangladesh MICS7 SPSS Datasets"

KNOWN_FILES: dict[str, dict[str, str]] = {
    "NationalSurveyData.csv": dict(
        source=_BGS, release="Release date: 25 May 2000 (file preamble)",
        role="arsenic_input", used_in_analysis="yes",
        notes="5-line preamble; header on line 5; units row on line 6 (As in ug/l).",
    ),
    "bgd2018.csv": dict(
        source=_STEPS, release="STEPS 2018",
        role="adult_bw_input", used_in_analysis="yes",
        notes="Adult measured weight m12; physical-measurement weight wstep2.",
    ),
    f"{_MICS_SPSS}/ch.sav": dict(
        source=_MICS, release="Bangladesh MICS7 2025, v01_M",
        role="child_bw_input", used_in_analysis="yes",
        notes="Under-5 questionnaire; child weight AN8; survey weight chweight.",
    ),
    f"{_MICS_SPSS}/bh.sav": dict(source=_MICS, release="Bangladesh MICS7 2025, v01_M", role="reference_only",
                                 used_in_analysis="no", notes="Birth history; not used for body weight."),
    f"{_MICS_SPSS}/fs.sav": dict(source=_MICS, release="Bangladesh MICS7 2025, v01_M", role="reference_only",
                                 used_in_analysis="no", notes="Children 5-17 questionnaire; no anthropometry used."),
    f"{_MICS_SPSS}/hh.sav": dict(source=_MICS, release="Bangladesh MICS7 2025, v01_M", role="reference_only",
                                 used_in_analysis="no", notes="Household file."),
    f"{_MICS_SPSS}/hl.sav": dict(source=_MICS, release="Bangladesh MICS7 2025, v01_M", role="reference_only",
                                 used_in_analysis="no", notes="Household listing."),
    f"{_MICS_SPSS}/wm.sav": dict(source=_MICS, release="Bangladesh MICS7 2025, v01_M", role="reference_only",
                                 used_in_analysis="no", notes="Women's file."),
    f"{_MICS_SPSS}/Data_Extractor.py": dict(
        source=_MICS, release="Bangladesh MICS7 2025, v01_M", role="third_party_script",
        used_in_analysis="no", notes="Shipped inside the raw folder; hashed only, never executed.",
    ),
    f"{_MICS_DIR}/Bangladesh MICS7 Datasets/Bangladesh MICS7 Datasets/Readme_Bangladesh_MICS7.rtf": dict(
        source=_MICS, release="Bangladesh MICS7 2025, v01_M", role="documentation",
        used_in_analysis="no", notes="Distributor readme.",
    ),
    f"{_MICS_DIR}/Bangladesh MICS7 Datasets.zip": dict(
        source=_MICS, release="Bangladesh MICS7 2025, v01_M", role="archive",
        used_in_analysis="no", notes="Original distributor archive of the extracted datasets.",
    ),
}


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        while block := fh.read(chunk):
            digest.update(block)
    return digest.hexdigest()


def build_inventory(raw_dir: Path = DATA_RAW) -> list[dict[str, str]]:
    rows = []
    for path in sorted(p for p in raw_dir.rglob("*") if p.is_file()):
        key = path.relative_to(raw_dir).as_posix()
        meta = KNOWN_FILES.get(key, dict(source="", release="", role="unregistered",
                                         used_in_analysis="no", notes="Not in provenance register."))
        rows.append(dict(
            relative_path=rel(path),
            file_name=path.name,
            size_bytes=str(path.stat().st_size),
            sha256=sha256_file(path),
            **meta,
        ))
    return rows


def write_inventory(rows: list[dict[str, str]], csv_path: Path, sha_path: Path) -> None:
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=INVENTORY_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    # sha256sum-compatible manifest: `sha256sum -c` from the project root verifies it.
    with sha_path.open("w", newline="\n", encoding="utf-8") as fh:
        for row in rows:
            fh.write(f"{row['sha256']}  {row['relative_path']}\n")


def read_hash_manifest(sha_path: Path) -> dict[str, str]:
    manifest = {}
    for line in sha_path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            digest, relative = line.split("  ", 1)
            manifest[relative] = digest
    return manifest


def verify_raw_hashes(sha_path: Path) -> list[str]:
    """Return a list of problems; empty means every frozen raw file is unchanged."""
    problems = []
    for relative, digest in read_hash_manifest(sha_path).items():
        path = PROJECT_ROOT / relative
        if not path.exists():
            problems.append(f"missing: {relative}")
        elif sha256_file(path) != digest:
            problems.append(f"changed: {relative}")
    return problems


PLANNED_PACKAGES = ["numpy", "pandas", "scipy", "matplotlib", "pyarrow", "pytest",
                    "lifelines", "pyreadstat", "statsmodels"]


def environment_record() -> dict[str, object]:
    packages = {}
    for name in PLANNED_PACKAGES:
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, capture_output=True,
                                text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = None
    return {
        "recorded_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python": sys.version.split()[0],
        "implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "git_commit": commit,
        "packages": packages,
        "missing_packages": [k for k, v in packages.items() if v is None],
    }
