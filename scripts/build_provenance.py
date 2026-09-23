"""Phase 1: freeze the raw-data inventory and record the environment.

Usage (from anywhere):
    python scripts/build_provenance.py            # write inventory + hashes + environment
    python scripts/build_provenance.py --verify   # check raw files against the frozen hashes
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from arsenic_hra.paths import OUTPUT_DIRS, RESULTS_PROVENANCE, ensure_output_dir, rel  # noqa: E402
from arsenic_hra.provenance import build_inventory, environment_record, verify_raw_hashes, write_inventory  # noqa: E402

INVENTORY_CSV = RESULTS_PROVENANCE / "data_inventory.csv"
HASH_MANIFEST = RESULTS_PROVENANCE / "raw_data.sha256"
ENVIRONMENT_JSON = RESULTS_PROVENANCE / "environment.json"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true", help="only verify raw files against the manifest")
    args = parser.parse_args()

    if args.verify:
        problems = verify_raw_hashes(HASH_MANIFEST)
        for p in problems:
            print(p)
        print("raw data unchanged" if not problems else f"{len(problems)} problem(s)")
        return 1 if problems else 0

    if HASH_MANIFEST.exists():
        problems = verify_raw_hashes(HASH_MANIFEST)
        if problems:
            print("Refusing to re-freeze: raw data differ from the existing manifest:")
            for p in problems:
                print("  " + p)
            return 1

    for d in OUTPUT_DIRS:
        ensure_output_dir(d)
    rows = build_inventory()
    write_inventory(rows, INVENTORY_CSV, HASH_MANIFEST)
    env = environment_record()
    ENVIRONMENT_JSON.write_text(json.dumps(env, indent=2) + "\n", encoding="utf-8")

    print(f"{len(rows)} raw files inventoried -> {rel(INVENTORY_CSV)}")
    print(f"hash manifest -> {rel(HASH_MANIFEST)}")
    print(f"environment -> {rel(ENVIRONMENT_JSON)}; missing packages: {env['missing_packages'] or 'none'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
