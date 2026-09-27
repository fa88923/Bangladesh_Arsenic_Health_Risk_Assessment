# Fast Pipeline Run

Run from the project root:

```powershell
cd "E:\CSE 402\Project"
```

Use Python 3.11 on this machine:

```powershell
py -3.11 --version
```

## Fastest Full Pipeline

This runs the required stages in dependency order and skips figure generation where supported.

```powershell
py -3.11 scripts/build_provenance.py
py -3.11 scripts/preprocess_arsenic.py
py -3.11 scripts/run_arsenic_fitting.py --no-figures
py -3.11 scripts/run_bodyweight.py --no-figures
py -3.11 scripts/run_risk_phase6.py
py -3.11 scripts/run_simulation.py --no-figures
py -3.11 scripts/run_convergence.py --no-figures
py -3.11 scripts/run_sensitivity.py --reuse-bootstrap --no-figures --jobs 8
```

## Minimum-Time Smoke Run

Use this only when you need to check that the pipeline still executes. It is not the full final arsenic bootstrap because it uses one bootstrap replicate instead of the configured 1000.

```powershell
py -3.11 scripts/build_provenance.py
py -3.11 scripts/preprocess_arsenic.py
py -3.11 scripts/run_arsenic_fitting.py --bootstrap 1 --no-figures
py -3.11 scripts/run_bodyweight.py --no-figures
py -3.11 scripts/run_risk_phase6.py
py -3.11 scripts/run_simulation.py --no-figures
py -3.11 scripts/run_convergence.py --no-figures
py -3.11 scripts/run_sensitivity.py --reuse-bootstrap --no-figures --jobs 8
```

## Notes

- `--no-figures` saves time and does not change the CSV/JSON result tables.
- `run_arsenic_fitting.py` has no true skip-bootstrap option. The fastest valid shortcut is `--bootstrap 1`, but that changes the bootstrap evidence.
- `run_sensitivity.py --reuse-bootstrap` reuses existing bootstrap tables in `results/sensitivity/` when present.
- Increase or reduce `--jobs 8` depending on available CPU cores.
- If raw data have already been frozen, `build_provenance.py` verifies that the raw files still match the existing manifest before proceeding.
