# Phase 1: Reproducibility and Data Inventory

Owner: Member 1. Reviewer: Member 2 (cross-review) and Member 5 (QA).

This document covers the Phase 1 deliverables in `docs/project_implementation_plan.md` section 6.

| Deliverable | Location |
|---|---|
| Data inventory and provenance table | `results/provenance/data_inventory.csv` |
| Frozen raw-file hashes | `results/provenance/raw_data.sha256` |
| Environment record | `results/provenance/environment.json` |
| Parameter decision register | `config/parameter_register.csv` (summarised in section 4) |
| Run configuration | `config/run_config.json` (specified in section 3) |
| Output naming convention | section 5 |
| Reproducibility checklist | section 6 |

Regenerate the inventory with `python scripts/build_provenance.py`. Check the raw files with `python scripts/build_provenance.py --verify`, or run `sha256sum -c results/provenance/raw_data.sha256` from the project root. If a raw file has changed, the script refuses to re-freeze the manifest.

## 1. Raw-file inventory

All 11 files under `data/raw/` are hashed with SHA-256. Three of them are analysis inputs:

| Role | File | Source and release |
|---|---|---|
| Arsenic | `data/raw/NationalSurveyData.csv` | DPHE/BGS/DFID National Hydrochemical Survey. The file preamble says "Release date: 25 May 2000". |
| Adult BW | `data/raw/bgd2018.csv` | WHO STEPS Bangladesh 2018 |
| Child BW | `data/raw/BGD_2025_MICS7_v01_M/BGD_2025_MICS7_v01_M/BGD_2025_MICS7_Datasets/Bangladesh MICS7 Datasets/Bangladesh MICS7 Datasets/Bangladesh MICS7 SPSS Datasets/ch.sav` | UNICEF Bangladesh MICS7 2025, release `BGD_2025_MICS7_v01_M`, SPSS datasets |

The other MICS files (`bh`, `fs`, `hh`, `hl`, `wm`, the readme, and the distributor zip) are hashed as reference-only files. `Data_Extractor.py` came inside the raw MICS folder. It is hashed but never executed.

**BGS layout, verified:** lines 1–4 are the title and release preamble. Line 5 is the header, line 6 is the units row (`As` = `ug/l`; `WELL_DEPTH` = `m`; the other ions are in `mg/L`), and lines 7–3540 hold 3,534 sample rows. The earlier audit counted the units row as a sample with a missing district. The preprocessing removes that row explicitly.

## 2. Path convention

- All code resolves paths through `src/arsenic_hra/paths.py` (`PROJECT_ROOT` = the repository root), so scripts give the same result from any working directory.
- `data/raw/` is read-only. `paths.ensure_output_dir` refuses to create anything under it, and `paths.raw_path` only resolves inside it.
- Derived files go only under `results/<workstream>/`: `provenance`, `arsenic`, `bodyweight`, `simulation`, `convergence`, `sensitivity`, `tables`, `figures`. Raw data are never copied there. The preprocessed files are new derived tables that keep `source_line` back-references into the raw file.
- Manifests record paths relative to the project root, in POSIX form.

## 3. Run-configuration specification (`config/run_config.json`)

| Key | Meaning | Consumer |
|---|---|---|
| `config_version` | Bump this on any change; every output manifest records it | all |
| `raw_inputs.*` | Project-relative raw paths; each must appear in the hash manifest | all |
| `arsenic.bgs_header_row`, `bgs_units_row_present` | Verified BGS layout (0-based header row 4) | Member 1 |
| `arsenic.source_unit`, `model_unit`, `ugL_to_mgL_divisor` | ug/L → mg/L, divide by 1000 | Members 1, 2 |
| `arsenic.selected_districts` | Frozen district list, exact `DISTRICT` strings | Members 1, 2, 4 |
| `arsenic.high_censoring_warning_percent` | Warning threshold (50%) used in the district summary | Members 1, 2 |
| `arsenic.candidate_distributions`, `primary_fit`, `robustness_fit`, `censored_substitution`, `lognormal_gamma_location`, `bootstrap_replicates` | Phase 4 fitting policy | Member 2 |
| `bodyweight.*` | Populations, weight/survey-weight columns, BW candidates | Member 3 |
| `simulation.master_seed`, `seed_derivation` | Master seed; child streams via `numpy.random.SeedSequence(master_seed, spawn_key=(district_idx, population_idx, n_idx, replicate))` | Member 4 |
| `simulation.primary_N`, `convergence_N`, `replicate_seeds_per_N` | 10,000; {1k, 5k, 10k, 20k}; 5 replicates | Member 4 |
| `simulation.sampling`, `input_dependence`, `iteration_storage` | Explicit inverse CDF; independent inputs; Parquet | Member 4 |
| `fixed_policies.*` | Raw data read-only, no automatic trimming, survey weights required | all |

The district index for seed derivation is the position in `selected_districts`. The population index is the position in `simulation.populations`. Changing either list order changes the random streams, so reordering requires a config version bump.

## 4. Parameter decision register

`config/parameter_register.csv` lists every deterministic and stochastic input with its unit, source, type (measured, fitted, specified, fixed, or assumption), status, owner, and open question. The main items:

**Frozen:** district set; censoring policy; IR and adult SA Lognormal log-space parameterization; independent sampling; simulation sizes.

**Open and unit-critical.** These must be signed off before Phase 6 and 7 outputs are reported:

1. **Child SA Triangular triplet.** Owner: Member 4. It is not defined, and it must not be derived from "0.6800 ± 0.600".
2. **Averaging time for ELCR.** Owner: Member 4. `AT = ED × 365` is currently applied to both HQ and ELCR. A decision is needed on whether cancer risk uses a 70-year lifetime AT, which matters especially for children.
3. **CF = 10 dimensional check.** Owner: Member 4. It needs to be checked against `Kp` in cm/h, `SA` in m², and `C` in mg/L.
4. **EF and ET definitions.** Owner: Member 4. The values are recorded, but they are listed as a hold point in master plan section 5.6.
5. **Treatment of highly censored districts.** Owner: Member 2. Rajshahi is 71.8% censored. Dhaka (48.9%) and Mymensingh (48.1%) are close to half.

**Documentation corrections to raise with the document owners.** Member 1 did not edit these files:

- `docs/Bangladesh_Arsenic_Risk_Model_Parameters_Full_Design.md` section 4 names "Bangladesh MICS 2019", but the raw data are **MICS7 2025** (`BGD_2025_MICS7_v01_M`).
- In the same document, the section 7 parameter table has a corrupted `K_p` row ("Derma7. Parameter definitions, values, and unitsl permeability coefficient").

## 5. Output naming convention

- Pattern: `results/<workstream>/<subject>_<content>[_<qualifier>].<ext>`, in lower snake case. Example: `results/arsenic/arsenic_district_preprocessing_summary.csv`.
- Units go in column names, not file names: `_ugL`, `_mgL`, `_kg`. A concentration column without a unit suffix is not allowed.
- Machine-readable tables are CSV (UTF-8, header row, no index). Iteration-level simulation outputs are Parquet. Figures are `results/figures/<workstream>_<subject>_<district|population>.png`.
- Every script that writes results also writes `<subject>_manifest.json` (or `run_manifest_<id>.json` for simulations). The manifest records the UTC time, script, `config_version`, input SHA-256, seeds where relevant, row counts, and the list of output files.
- Simulation run IDs: `<district>_<population>_N<size>_r<replicate>`, for example `rajshahi_child_N10000_r0`.
- Review files, which list rows for a human decision without removing anything, end in `_review.csv`. Audit or reconciliation files end in `_audit.csv`.

## 6. Reproducibility checklist

- [x] Raw inventory frozen with SHA-256 (`results/provenance/raw_data.sha256`); verified by `tests/test_arsenic_data_contract.py::test_raw_files_unchanged`.
- [x] Exact MICS `ch.sav` path and release (`BGD_2025_MICS7_v01_M`) recorded in `config/run_config.json` and the inventory.
- [x] Project-root-aware paths; writes into `data/raw/` are refused (`tests/test_run_config.py`).
- [x] Environment recorded (`results/provenance/environment.json`): Python 3.12.1, numpy 1.26.4, pandas 2.3.3, scipy 1.17.1, matplotlib 3.11.1, pyarrow 25.0.0, pytest 9.1.1.
- [ ] **Install missing planned packages before Phases 3–4:** `lifelines` (reverse-KM, Member 2), `pyreadstat` (MICS `.sav`, Member 3), and `statsmodels` (Wilson interval, Member 4). They are listed in `requirements.txt`. Re-run `scripts/build_provenance.py` afterwards to update the environment record.
- [x] Run configuration covers districts, candidates, seeds, simulation sizes, and policies.
- [x] Parameter register created with an owner and status for every input.
- [ ] Open register items (section 4) signed off with reviewer and date.
- [ ] Master seed confirmed by Member 4 before the first recorded simulation.
- [ ] Final run repeated from a clean checkout; outputs compared by hash.
