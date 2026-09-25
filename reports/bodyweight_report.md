# Phases 3 and 5: Body-Weight Preprocessing and Fitting Handoff

Owner: Member 3. The cross-reviewer is Member 4, who also uses these outputs in Phase 7 (simulation). Member 5 does QA.

## How to reproduce

```bash
python scripts/build_provenance.py --verify   # raw data unchanged
python scripts/run_bodyweight.py              # writes results/bodyweight/* and results/figures/bodyweight_*; exits 1 on any failed check
python -m pytest -q                           # all tests; 45 of them cover body weight
```

The script reads raw paths from `config/run_config.json`. It stops before writing anything if a raw input differs from `results/provenance/raw_data.sha256`.

## Outputs (`results/bodyweight/`)

| File | Contents |
|---|---|
| `adult_bw_clean.csv` | **Adult fitting input.** 8,013 STEPS rows: `BW_kg` from `m12` and `survey_weight` from `wstep2`, plus the design variables `psu` and `stratum` |
| `child_bw_clean.csv` | **Child fitting input.** 22,576 MICS rows: `BW_kg` from `AN8` and `survey_weight` from `chweight`, plus `PSU`, `stratum`, and the WHO flags |
| `adult_bw_cleaning_audit.csv`, `child_bw_cleaning_audit.csv` | Row counts for each exclusion step, applied in order. `raw_rows_matching` counts the same rule against the full raw file |
| `adult_bw_bmi_review.csv` | 12 retained adults with a diagnostic BMI outside 12–60 kg/m². Nothing is removed |
| `child_bw_an8_codebook_review.csv` | Undocumented `AN8 >= 90` values held for codebook review. It is currently empty |
| `bodyweight_weighted_summary.csv` | Weighted mean, SD, and percentiles, plus the Kish effective n |
| `bodyweight_distribution_fit_results.csv` | 8 rows (2 populations × 4 families): LS and PML parameters, fit errors, tail quantiles, support checks, and convergence messages |
| `bodyweight_selected_distributions.json` / `.csv` | **The provisional choice for each population.** Includes natural parameters, explicit SciPy arguments, and the rationale |
| `bodyweight_manifest.json` | Config version, raw SHA-256, software versions, row counts, thresholds, and the output list |

The figures are in `results/figures/bodyweight_<subject>_<population>.png`. Each population has a histogram with PDFs, the ECDF with fitted CDFs, a Q-Q plot for each family, and lower- and upper-tail close-ups.

## Reconciliation

| Adults (STEPS 2018) | Rows |
|---|---:|
| Read from `bgd2018.csv` | 8,185 |
| Age outside 18–69 | 0 |
| `m12` missing | −167 |
| `m12 = 888` (refused) / `666` (too large) | −5 / 0 |
| `m12 <= 0`, sex missing, `wstep2` missing or <= 0 | 0 |
| **Clean** | **8,013** |

| Children (MICS7 `ch.sav`, not `bh.sav`) | Rows |
|---|---:|
| Read from `ch.sav` | 24,680 |
| `AN8` missing | −1,323 |
| `AN8 = 99.3 / 99.4 / 99.5 / 99.6` | −182 / −87 / −306 / −149 |
| Undocumented `AN8 >= 90`, `AN8 <= 0`, `CAGE` missing or outside 0–59, sex missing, `chweight` <= 0 | 0 |
| `WAZFLAG = 1` | −57 |
| **Clean** | **22,576** |

In the raw file, 781 records have `WAZFLAG = 1`, but 724 of them are the special-code rows already removed. The remaining 57 are real weight-for-age errors. 286 `HAZFLAG` and 288 `WHZFLAG` errors are retained, because they concern height and not the weight being fitted. Height codes (`m11 = 888`, `AN11 = 999.4–999.6`) are blanked in `height_cm` and are not treated as heights. There are no duplicate `pid` or `(HH1, HH2, LN)` keys.

## Weighted summary (kg)

| Population | n | Kish effective n | Mean | SD | P5 | P50 | P95 | Range |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| Adult 18–69 | 8,013 | 3,238 | 55.89 | 11.62 | 39.0 | 55.0 | 76.0 | 28–162 |
| Child 0–59 months | 22,576 | 16,600 | 10.87 | 3.22 | 5.6 | 10.9 | 15.9 | 1.0–36.0 |

## Fit results (weighted CDF least squares)

| Population | Family | LS parameters | SSE_w | D_w | Tail error* | Admissible |
|---|---|---|---:|---:|---:|---|
| Adult | Normal | μ=55.08, σ=10.91 | 2.28e-4 | 0.030 | 12.7% | yes |
| Adult | **Lognormal** | μ_log=4.0024, σ_log=0.2002 | **1.96e-5** | **0.019** | **1.0%** | yes |
| Adult | Gamma | k=25.27, θ=2.199 | 4.33e-5 | 0.020 | 3.6% | yes |
| Adult | Triangular | a=32.3, m=51.4, b=82.8 | 6.9e-5 | 0.027 | 10.0% | no: excludes 2.65% of survey mass |
| Child | Normal | μ=10.87, σ=3.22 | 2.8e-5 | 0.009 | 5.2% | no: P(BW≤0)=3.7e-4 |
| Child | Lognormal | μ_log=2.3706, σ_log=0.2993 | 5.96e-4 | 0.043 | 49.7% | yes |
| Child | **Gamma** | k=11.392, θ=0.9725 | **3.15e-4** | **0.032** | 37.2% | yes |
| Child | Triangular | a=3.17, m=11.23, b=18.09 | 3.1e-5 | 0.015 | 19.6% | no: excludes 1.91% of survey mass |

\*Tail error is the largest relative error at the weighted P1, P5, P95, and P99.

Every LS and PML optimizer reports convergence. The admissibility thresholds were set before fitting:
- A Normal fit is rejected if P(BW ≤ 0) ≥ 1/10,000, which means at least one expected negative draw in the N = 10,000 run.
- A Triangular fit is rejected if it excludes more than 0.1% of survey mass.

## Provisional selections (require team approval; master plan section 18)

- **Adult: Lognormal(μ_log = 4.002396, σ_log = 0.200174).** In SciPy: `lognorm(s=0.200174, loc=0, scale=exp(4.002396))`. It is best on SSE, D_w, and tail error, and within 1.4% of the PML fit. This selection is well supported.
- **Child: Gamma(k = 11.391855, θ = 0.972526).** In SciPy: `gamma(a=11.391855, loc=0, scale=0.972526)`. **This one needs a team decision.** The child data are close to symmetric. Normal fits about 11 times better (SSE 2.8e-5), but it was rejected under the predeclared rule. Gamma under-fits the light-weight tail: its fitted P1 is 4.89 kg against an empirical 3.56 kg. The lightest children receive the highest dose per kg, so this matters for the upper tail of the child HI. The plan allows a truncated Normal only if one is explicitly introduced and justified. The options are:
  1. accept Gamma and report the low-tail under-fit as a limitation;
  2. approve a truncated Normal (on BW > 0) as a fifth candidate for children.

Member 4 can draw samples with `arsenic_hra.bodyweight_fitting.sample_bw(record, N, rng)`, where `record` is one entry of the JSON file. The function applies U ~ Uniform(0,1) and then the selected `.ppf`.

## Findings for review

1. **Two adult weights look like typed-in heights.** `E13-37-1` has 144 kg with a height of 144 cm, and `E31-117-2` has 151 kg with a height of 151 cm. Both are listed in `adult_bw_bmi_review.csv` and retained, following the no-trimming rule. Together they carry 0.06% of the adult survey mass. They affect only the Triangular support check, not the selected Lognormal fit. The team should decide whether to treat them as source coding errors.
2. **The adult Kish effective n is 3,238**, compared with 8,013 rows. The `wstep2` weights are highly variable. They are used as supplied.
3. **Under-5 coverage only.** The child distribution covers 0–59 months. The 6-year child ED mismatch remains the approved limitation, to be carried into the Phase 9 sensitivity analysis.
4. **Documentation.** Design doc section 4 names "Bangladesh MICS 2019"; the data are MICS7 2025 (`BGD_2025_MICS7_v01_M`). Design doc section 9 selects models with AIC/KS/AD, but the body-weight plan (section 9.6) forbids ordinary AIC on weighted pseudo-likelihoods. This work follows the body-weight plan.

## Assumptions

- Survey weights are normalized: q = w / Σw for the ECDF, and w̃ = n·w / Σw for the PML fit. Weights are otherwise used as supplied, with no trimming.
- Tied weights are aggregated. The plotting positions are p_j = Σ_{r<j} q_r + q_j/2.
- Lognormal and Gamma use loc = 0. Triangular has 0 < a < m < b, and its PML support must contain every observation.
- The weighted pseudo-log-likelihood is reported but not converted to AIC.
- The optional survey-design-aware parameter bootstrap (plan section 13) is not done. `psu` and `stratum` are kept in both clean files so it can be added.

## Validation gate (master plan sections 8 and 10)

| Gate | Status |
|---|---|
| Starting, excluded, and final counts reconcile exactly | pass (`test_audit_reconciles_with_clean_file`) |
| Adult `BW_kg` comes from `m12`, never `wstep2` | pass (`test_adult_population_and_source_column`) |
| Child `BW_kg` comes from `AN8`, never `chweight` | pass (`test_mics_bw_comes_from_an8_not_chweight`) |
| No documented special code survives as body weight | pass (`test_clean_values_are_physical`) |
| Survey weights are present, finite, and positive | pass (`test_clean_values_are_physical`) |
| Clean outputs keep the design variables | pass (`psu`/`stratum`, `PSU`/`stratum` columns) |
| Scaling the weights by a constant does not change the estimates | pass (`test_weight_scale_invariance`) |
| Every optimizer reports convergence or failure | pass (`test_every_fit_reports_convergence_state`) |
| Saved parameters reproduce SSE_w and D_w | pass (`test_saved_selection_reproduces_saved_fit_metrics`) |
| Selected distributions pass the physical-support checks | pass (admissibility rule; `test_normal_with_negative_mass_is_inadmissible`) |
| Inverse-CDF samples are finite and positive | pass (`test_sample_bw_inverse_cdf`, output test) |
| Adult and child selections are justified independently | pass (separate fits and rationales); **child choice awaits team review** |

Reviewer sign-off (Member 4): ________  Date: ________
