# Phase 2: Arsenic Preprocessing Handoff

Owner: Member 1. Handoff to Member 2 (Phase 4 fitting), who is also the cross-reviewer. Member 5 does QA.

## How to reproduce

```bash
python scripts/build_provenance.py --verify   # raw data unchanged
python scripts/preprocess_arsenic.py          # writes results/arsenic/*, exits 1 on any failed check
python -m pytest -q                           # 47 unit, config and data-contract tests
```

## Outputs (`results/arsenic/`)

| File | Contents |
|---|---|
| `arsenic_selected_districts_preprocessed.csv` | **Fitting input.** 810 rows for the 10 approved districts |
| `arsenic_all_districts_classified.csv` | All 3,534 samples, with the same classification; basis for QC |
| `arsenic_district_preprocessing_summary.csv` | Per-district counts, censoring, limits, and detected range (mg/L) |
| `arsenic_reconciliation_audit.csv` | Countable flow from raw lines to output rows, plus 8 pass/fail checks |
| `arsenic_duplicate_review.csv` | Duplicate IDs, full rows, and repeated coordinates, for review; nothing deleted |
| `arsenic_malformed_value_review.csv` | Unparseable or missing arsenic, coordinates, or dates (currently empty) |
| `arsenic_preprocessing_manifest.json` | Raw SHA-256, config version, units, row counts, and the output list |

### Column contract

- `source_line` is the 1-based line number in the raw CSV. Use it to trace any row back to its raw value.
- `As_raw` is a verbatim copy of the raw `As` string, for example `< 6`.
- `As_status` is one of `detected`, `left_censored`, `missing`, or `malformed`. `As_censored` is True or False, and blank when the status is missing or malformed.
- `As_detected_*` holds exact measurements only. It is **NaN for every censored row**.
- `As_bound_*` holds the censoring limit L (for a row reported as `C < L`). It is NaN for detected rows.
- `As_model_value_*` is the detected value, or L for censored rows. **This is the `durations` input for `KaplanMeierFitter.fit_left_censoring` and for `CensoredData.left_censored`, together with `As_censored`. It is not a concentration to fit directly.**
- The `_mgL` columns equal the matching `_ugL` columns divided by 1000. Fit in mg/L.
- `SAMPLE_DATE` is the raw day/month/year string. `SAMPLE_DATE_iso` is the parsed date.

## Reconciliation

| Step | Count |
|---|---:|
| Raw file lines | 3,540 |
| Preamble, header, and units row removed | 4 + 1 + 1 |
| Raw sample rows | 3,534 |
| Detected / left-censored / missing / malformed (all districts) | 2,429 / 1,105 / 0 / 0 |
| Rows in the approved districts, expected from raw `DISTRICT` counts | 810 |
| Rows written | 810 (524 detected, 286 censored) |
| Rows excluded because they fall outside the approved districts | 2,724 |

All 8 checks pass. Row flows reconcile at every step; `As_raw` equals the raw string for every row; no censored value appears in a detected field; district names match the source exactly; and the raw SHA-256 was the same before and after the run.

## District summary (mg/L)

| District | n | Detected | Censored | Censored % | Limits (mg/L) | Detected range | Warning |
|---|---:|---:|---:|---:|---|---|---|
| Dhaka | 45 | 23 | 22 | 48.9 | 0.0005; 0.006 | 0.001–0.262 | |
| Chittagong | 44 | 39 | 5 | 11.4 | 0.006 | 0.0005–0.344 | |
| **Rajshahi** | 78 | 22 | 56 | **71.8** | 0.0005; 0.006 | 0.0006–0.0918 | **high censoring** |
| Khulna | 76 | 44 | 32 | 42.1 | 0.0005; 0.006 | 0.0005–0.538 | |
| Barisal | 92 | 76 | 16 | 17.4 | 0.0005; 0.006 | 0.0006–0.862 | |
| Sylhet | 77 | 52 | 25 | 32.5 | 0.0005; 0.006 | 0.0015–0.157 | |
| Rangpur | 86 | 47 | 39 | 45.3 | 0.0005 | 0.0005–0.298 | |
| Mymensingh | 108 | 56 | 52 | 48.1 | 0.0005 | 0.0005–0.200 | |
| Comilla | 110 | 101 | 9 | 8.2 | 0.0005; 0.006 | 0.0008–0.698 | |
| Bogra | 94 | 64 | 30 | 31.9 | 0.0005 | 0.0005–0.632 | |

Rajshahi is 71.8% censored. Only 22 detections remain for fitting there, so every fit for Rajshahi must carry the high-uncertainty label. Dhaka and Mymensingh are both close to 50% censored.

## Findings for Member 2 (Phase 4)

1. **Two detection limits appear together.** The survey uses both `< 0.5` µg/L (994 rows overall) and `< 6` µg/L (111 rows overall). In 7 of the 10 districts, some detected values fall **below** a `< 6` limit in the same district. For example, Chittagong has 19 detections below 6 µg/L but 5 rows reported as `< 6`, and Barisal has 34. A `< 6` row therefore says only that C is below 0.006 mg/L; it does not mean C is near 0.006. The reverse Kaplan–Meier and censored-likelihood methods handle this correctly. A simple "fraction below the detection limit" plotting position would not.
2. **Some detected values equal the lower limit.** The minimum detected value is 0.5 µg/L, which equals the `< 0.5` limit. Under the left-censoring Kaplan–Meier convention, a tie between a detection and a censored limit is ordered with the detection first. Check how `lifelines` resolves these ties, and document the result.
3. There are no zero, negative, missing, or malformed arsenic values, so the Lognormal and Gamma support (x > 0) is satisfied by the data.
4. Large values such as 0.862 mg/L in Barisal and 1.66 mg/L nationally are kept. No trimming was applied.

## Duplicate review

- Duplicate `SAMPLE_ID`: 0. Duplicate full rows: 0.
- Repeated coordinates: 227 rows in 114 location groups across the whole survey. **62 of those rows are in the approved districts.** Every repeated-coordinate group in the approved districts differs in well depth or sampling date, so each is treated as a distinct well or visit and retained.
- One group in Sunamganj (`S99_05026`/`S99_05027`, not an approved district) has the same location, depth, and date. The two records differ in union and arsenic value (1.3 and 6.3 µg/L), so both are retained. This does not affect the analysis.

## Assumptions

- The only non-sample row is the units row. It was verified by content (`As` = `ug/l`, empty `SAMPLE_ID`), and the script raises an error if that layout changes.
- Censored limits are stored as bounds and are never substituted as values (no LOD, LOD/2, or zero substitution).
- The BGS data are treated as a sample of wells within each district, not a population-weighted census.
- `SAMPLE_DATE` is parsed as day/month/year. All 3,534 dates parse (the malformed review is empty).

## Validation gate (master plan section 7)

| Gate | Status |
|---|---|
| Raw row counts reconcile with classified row counts | pass (audit checks) |
| Every selected row keeps its original arsenic string | pass (`test_as_raw_matches_source_string`) |
| Unit conversion verified with fixed examples | pass (`test_unit_conversion_fixed_examples`, `test_mgL_is_ugL_over_1000`) |
| No censored value in the detected field | pass (`test_no_censored_value_in_detected_field`) |
| Duplicate coordinates reviewed rather than deleted | pass (`arsenic_duplicate_review.csv`) |
| District names exactly match the raw source | pass (`test_district_names_exactly_match_source`) |

Reviewer sign-off (Member 2): ________  Date: ________
