# Phase 7: Monte Carlo Simulation Handoff

Owner: Member 4. Cross-reviewer: Member 3. Member 5 does QA. Phases 8 (convergence) and 9 (sensitivity) reuse the engine and the saved iterations.

## How to reproduce

```bash
python scripts/run_risk_phase6.py      # Phase 6 tables (equations, parameters, benchmarks)
python scripts/run_simulation.py       # 20 primary runs, tables, figures; exits 1 on any failed check
python -m pytest -q                    # 216 tests; tests/test_simulation.py covers this phase
```

The run takes about 7 seconds. The same seed reproduces the output bit for bit.

## Code

| Module | Contents |
|---|---|
| `src/arsenic_hra/simulation_inputs.py` | Final C and BW sampling models: the Phase 4 and 5 provisional records plus the decisions below; the truncated-Normal child-BW refit (weighted CDF LS + weighted PML); the reverse-KM district mean |
| `src/arsenic_hra/monte_carlo.py` | `run_spec` and `generator` (seed streams), `simulate` (explicit inverse-CDF sampling → Phase 6 equations), `summarize` (percentiles and direct-count exceedance with MC SE and Wilson CI), `sampling_check` (KS test against the target distributions) |
| `src/arsenic_hra/simulation_plots.py` | District HI figures (design doc section 14.2) |

## Outputs

| File | Contents |
|---|---|
| `results/simulation/simulation_input_distributions.json` / `.csv` | The C model for each district and the BW model for each population actually sampled, with parameters, SciPy arguments, source, and decision |
| `results/simulation/arsenic_selection_review.csv` | Lognormal vs Gamma for every district: CDF RMSE, D, AICc, fitted/empirical P95 and P99, and fitted mean / reverse-KM mean |
| `results/simulation/child_bw_selection_review.csv` | Gamma vs truncated Normal for child BW against the Phase 5 criteria |
| `results/simulation/iterations/<run_id>.parquet` | 20 files × 10,000 iterations: the 6 sampled inputs, ADD_ing, ADD_dermal, HQ_ing, HQ_dermal, HI, ADD_ing_cancer, ELCR, and ELCR_ED_AT (float64, zstd) |
| `results/simulation/manifests/run_manifest_<run_id>.json` | For each run: seed, spawn key, input order, every distribution, output SHA-256, config SHA-256 |
| `results/simulation/simulation_summary.csv` | For each run: mean, SD, min, max, P5–P99, and threshold probabilities for HI, ELCR, ELCR_ED_AT, HQ_ing, HQ_dermal |
| `results/simulation/sampling_validation.csv` | Sampled inputs against their target distributions (KS distance, P50, P95, support) |
| `results/tables/risk_probabilistic_p95.csv` | **Same layout as base-paper Table 5**: P95 HI and ELCR for adults and children |
| `results/tables/risk_exceedance.csv` | P(HI>1) with MC SE and Wilson 95% CI, P(HI>2), P(ELCR>1e-4) |
| `results/tables/risk_deterministic_district.csv` | **Same layout as base-paper Table 4**: point-input HI and ELCR |
| `results/tables/risk_deterministic_vs_probabilistic.csv` | Deterministic vs probabilistic mean, P50, P95 |
| `results/figures/simulation_hi_distribution_<district>.png` | Two panels (adult and child) showing the HI distribution, with HI = 1 and P95 marked |

## Decision log (2026-09-25, Member 4)

The user delegated these choices to Member 4: "recommend the best distribution … and use that; log your reasoning". Each decision applies the criteria of the arsenic plan (section 11), the body-weight plan (section 9.6) and the master plan (section 9) to the evidence in the Phase 4 and 5 outputs.

### D6. Accept 9 of the 10 provisional arsenic selections

The arsenic plan ranks by CDF least squares (RMSE and CDF overlay), then censored-MLE AICc. When those disagree, it looks at the tails and at bootstrap stability. The master plan also requires quantile and tail review, and says selection "is not based on a single statistic alone".

| District | Family | Evidence |
|---|---|---|
| Dhaka, Rajshahi, Khulna, Sylhet, Mymensingh, Comilla | Gamma | Gamma wins RMSE, D and AICc, and tracks the upper tail far better than the Lognormal alternative. For example, Khulna's fitted P99 is 1.9× empirical under Gamma and 50× under Lognormal. |
| Rangpur | Gamma | RMSE favours Gamma 2.26 : 1 and AICc favours Lognormal by 3.5. The tail breaks the tie: Gamma P95 is 1.18× empirical, against 2.26× for Lognormal. Gamma's extreme tail is thin (P99.9 = 0.11 mg/L, while the largest well measured 0.298 mg/L), so it is conservative only up to about P99. |
| Chittagong | Lognormal | Lognormal wins RMSE, AICc (ΔAICc 18) and bootstrap (98.8%). Both families under-predict P95 (0.58× and 0.25×); Lognormal is the less biased. |
| Bogra | Lognormal | Lognormal wins RMSE (ratio 1.36), D, AICc (Δ21) and bootstrap (99.5%). Its tail runs high: P99 is 3.9× empirical, and the fitted mean is 3.0× the reverse-KM mean. Gamma would under-predict (P99 0.73×). The primary criteria agree strongly, so Lognormal is kept, and its heavy tail is flagged (child HI maximum 15,865). |

### D7. Barisal: Lognormal → Gamma (override)

| Criterion | Lognormal (provisional) | Gamma | Verdict |
|---|---:|---:|---|
| CDF RMSE (primary) | 0.0415 | 0.0447 | near tie (ratio 1.08) |
| Max CDF discrepancy | 0.088 | 0.149 | favours Lognormal |
| Censored-MLE AICc | −223.5 | −221.0 | weakly favours Lognormal (Δ2.5) |
| Fitted P95 / empirical | 3.9× | 1.19× | **strongly favours Gamma** |
| Fitted P99 / empirical | 26× (19.6 mg/L) | 1.6× (1.20 mg/L) | **strongly favours Gamma** |
| Fitted mean / reverse-KM mean | 29× (2.71 vs 0.092 mg/L) | 1.09× | **strongly favours Gamma** |
| P(C > measured maximum 0.862 mg/L) | 7.5% | 2.2% | favours Gamma |
| Bootstrap LS-RMSE win rate | 84% | 16% | favours Lognormal |

The CDF-fit advantage of Lognormal is small, and it comes from the lower and central parts of the curve. Its upper tail is not supported by any data: 7.5% of draws exceed every well ever measured in Barisal, and its P99.9 is 291 mg/L. The reported risk quantities (P95 HI, P(HI>1), mean HI) are driven by exactly that tail, so Lognormal would inflate Barisal's P95 HI about fourfold as a fitting artefact. The master plan forbids hiding "unsupported tail extrapolation". Gamma is the admissible candidate with the strongest *overall* support. Its parameters are Member 2's saved CDF least-squares fit (k = 0.1724, θ = 0.5818), and nothing was refitted. **Phase 9 will run Barisal-Lognormal as a model-uncertainty scenario.**

### D8. Child body weight: Gamma → Normal truncated at 1.6 kg

The body-weight plan (section 9.6) ranks by SSE_w, then D_w, then the Q-Q and tail behaviour, then admissibility, then PML agreement. The design doc (section 9) allows a Normal "only if a justified truncated-normal model is explicitly used". The Normal was rejected in Phase 5 only because P(BW ≤ 0) = 3.7×10⁻⁴. Truncating it removes that problem.

| Criterion | Gamma (provisional) | Truncated Normal (L = 1.6 kg) |
|---|---:|---:|
| Weighted CDF SSE | 3.15×10⁻⁴ | **2.94×10⁻⁵** (11× better) |
| Max weighted discrepancy D_w | 0.032 | **0.0096** |
| Worst tail-quantile error (P1–P99) | 36% | **1.7%** |
| Fitted P1 / P99 (empirical 3.60 / 18.6 kg) | 4.89 / 20.1 | **3.58 / 18.4** |
| E[1/BW], which drives dose (empirical 0.1029 kg⁻¹) | 0.0989 (−3.8%) | **0.1030 (+0.1%)** |
| LS vs weighted-PML parameter difference | – | 0.8% |
| Physical support | BW > 0 | BW ≥ 1.6 kg |

**Choice of bound.** L = 1.6 kg is the lightest weight among the 22,576 cleaned MICS children, so the model never extrapolates below the observed data. Truncating at 0 fits the ECDF almost as well (SSE 2.86×10⁻⁵), but it places 7 in 10,000 draws below 1 kg, and since dose scales as 1/BW this pushes E[1/BW] up 17%. Parent parameters: μ = 10.8570, σ = 3.2281 (the parent Normal's, not the truncated moments). Gamma under-represented the lightest children, who receive the highest dose per kg, so it biased the child upper-tail HI low.

### D9. Adult body weight: Lognormal accepted

It is best on SSE_w, D_w and tail error, and within 1.4% of the PML fit. Every criterion agrees.

### D10. Highly censored districts (Rajshahi 71.8%) are retained with a warning

Dropping Rajshahi would bias the district set toward better-measured areas. Substituting values for non-detects is forbidden. Its Gamma fit is corroborated by AICc and bootstrap (98% win rate), and Rajshahi carries the high-uncertainty label in every table.

**Caveat for all districts.** Arsenic below the smallest detection (0.0005–0.0015 mg/L) is unresolved, and there the fitted distribution is pure extrapolation. With Gamma k < 1, the model puts a large share of draws at near-zero concentration. This affects **P50 and lower percentiles**; Rajshahi's median HI rounds to 0. It barely affects **P95 and the threshold probabilities**. With the other inputs at their medians, HI > 1 needs C ≈ 0.018 mg/L for adults and C ≈ 0.0036 mg/L for children, both inside the measured range. Measured directly, iterations with C below the district's smallest detection contribute **0.00 percentage points** to any adult P(HI>1), **at most 0.02 points** to child P(HI>1) in 9 districts, and **0.32 points in Sylhet** (smallest detection 0.0015 mg/L). Report P50 with this caveat.

### D11. Seeds and sampling

- **Master seed.** 20260923 is frozen before the first recorded run. It was not changed, since changing a seed after seeing results would be seed-shopping.
- **Streams.** Each run's stream is `SeedSequence(20260923, spawn_key=(district, population, N index, replicate))`, as specified in the config. All 400 runs in the Phase 8 design have distinct keys and distinct first draws. They are also distinct from the Phase 4 bootstrap streams, which use the key form `(district, replicate)`.
- **Uniforms.** Within a run, six uniform vectors are drawn in the fixed order C, IR, BW, EF, ET, SA, with u in [2.2×10⁻³⁰⁸, 1). Keeping u > 0 means the Lognormal IR and SA quantiles can never be 0, which the equations reject.
- **Dependence.** Inputs are independent, as approved.

### D12. Storage

The primary runs keep full iteration-level Parquet files: float64, all inputs and outputs, 26 MB in total. Phase 9 sensitivity analysis needs exactly these iterations. Storage for the Phase 8 convergence runs is decided in Phase 8.

### D13. Deterministic district benchmark

C is the **censoring-aware reverse-KM mean**, the empirical mean with no substitution and no dependence on the fitted model. The alternatives were rejected: the fitted-distribution mean is model-dependent (Barisal's provisional Lognormal would have given 2.7 mg/L), and a detected-only mean would ignore the non-detects. BW is the Bangladesh survey-weighted mean, and the other inputs are the base paper's deterministic values (Phase 6, D5). The reverse-KM mean puts censored mass at 0; placing it at the smallest detection instead changes C by less than 0.0005 mg/L.

### D14. Reported quantities

- **ELCR.** Lifetime AT is the primary result (Phase 6, D2). The paper-convention **ELCR_ED_AT** (AT = ED × 365) is stored in every iteration and reported as `P95_ELCR_child_ED_AT` for comparison with the paper. The two are identical for adults.
- **Thresholds.** P(HI > 1) with MC SE and Wilson interval; P(HI > 2) (design doc section 14.3); and P(ELCR > 10⁻⁴), the WHO and EPA upper tolerable risk.
- **Validation gate.** Sampled inputs must pass KS at α = 0.1% for every run and input. At α = 5% about 6 of 120 would fail by chance, which is not a defect. The observed result is 112/120 passing at 5% and 120/120 at 0.1%.

### D15. Figures

HI spans 10⁻¹⁷ to 10⁴, so the figures show the share of iterations in log-spaced bins on a log axis. A plain density on this range is dominated by the near-zero draws below the detection limit. Everything at or below HI = 0.001 is collected into the first bin, and its share is printed on the panel. The figures are visual aids only; no reported number is read from a plot.

## Results (primary run, N = 10,000)

### Probabilistic P95 (same layout as base-paper Table 5)

| District | P95 HI adult | P95 ELCR adult | P95 HI child | P95 ELCR child | P95 ELCR child (paper AT) |
|---|---:|---:|---:|---:|---:|
| Dhaka | 21.6 | 9.68×10⁻³ | 121.1 | 4.67×10⁻³ | 5.45×10⁻² |
| Chittagong | 6.66 | 2.99×10⁻³ | 34.7 | 1.34×10⁻³ | 1.56×10⁻² |
| Rajshahi ⚠ | 3.17 | 1.42×10⁻³ | 18.2 | 6.99×10⁻⁴ | 8.16×10⁻³ |
| Khulna | 15.7 | 7.03×10⁻³ | 84.3 | 3.25×10⁻³ | 3.79×10⁻² |
| Barisal | 32.1 | 1.44×10⁻² | 175.1 | 6.75×10⁻³ | 7.87×10⁻² |
| Sylhet | 6.61 | 2.96×10⁻³ | 35.8 | 1.38×10⁻³ | 1.61×10⁻² |
| Rangpur | 1.84 | 8.27×10⁻⁴ | 9.62 | 3.71×10⁻⁴ | 4.33×10⁻³ |
| Mymensingh | 7.19 | 3.23×10⁻³ | 37.9 | 1.46×10⁻³ | 1.70×10⁻² |
| Comilla | 51.3 | 2.30×10⁻² | 259.6 | 1.00×10⁻² | 1.17×10⁻¹ |
| Bogra | 6.57 | 2.95×10⁻³ | 33.3 | 1.28×10⁻³ | 1.50×10⁻² |

⚠ = high-censoring district (71.8% non-detects).

### Threshold exceedance

MC SE for P(HI>1) is at most 0.005 (0.5 percentage points) for every run.

| District | P(HI>1) adult | P(HI>1) child | P(HI>2) adult | P(HI>2) child | P(ELCR>10⁻⁴) adult | P(ELCR>10⁻⁴) child |
|---|---:|---:|---:|---:|---:|---:|
| Dhaka | 33.5% | 47.7% | 26.8% | 42.2% | 46.0% | 39.9% |
| Chittagong | 23.2% | 53.8% | 14.3% | 40.3% | 50.6% | 35.3% |
| Rajshahi ⚠ | 11.2% | 20.3% | 7.1% | 16.4% | 19.2% | 15.0% |
| Khulna | 30.0% | 44.9% | 23.1% | 39.2% | 43.4% | 36.7% |
| Barisal | 41.2% | 55.8% | 33.8% | 49.5% | 54.7% | 47.4% |
| Sylhet | 32.8% | 60.0% | 20.6% | 50.4% | 57.5% | 46.1% |
| Rangpur | 10.0% | 29.9% | 4.4% | 20.9% | 27.6% | 18.0% |
| Mymensingh | 19.1% | 32.0% | 13.8% | 26.7% | 30.5% | 25.0% |
| Comilla | 72.3% | 85.1% | 62.5% | 80.1% | 85.6% | 78.1% |
| Bogra | 15.5% | 32.4% | 10.6% | 24.6% | 30.5% | 21.9% |

### Main findings

1. **Comilla has the highest risk** (reverse-KM mean 0.142 mg/L; 72% of adults and 85% of children above HI = 1). It is followed by Barisal, Dhaka and Khulna. Rangpur and Rajshahi have the lowest risk.
2. **Children's risk is higher everywhere.** P95 HI is 5.1–5.7× the adult value, because under-5 children weigh about one fifth as much and drink comparable volumes. Child ELCR under lifetime AT is *lower* than adult ELCR, because 6 exposure years are averaged over 70. Under the paper's AT convention it is 5.1–5.7× higher than adult ELCR.
3. **Probabilistic P95 is 1.8–5.4× the deterministic point estimate.** Arsenic distributions are heavily right-skewed, so a mean-based deterministic table understates upper-tail risk. It also overstates the typical case: the median is well below the deterministic value in every district.
4. **The dermal pathway is negligible**: HQ_dermal is 0.1–0.3% of mean HI.

## Validation gate (master plan section 12)

| Gate | Status |
|---|---|
| Repeating a run with the same seed reproduces results exactly | pass (script check, bit-identical Parquet; `test_same_seed_reproduces_exactly`) |
| Different district/population streams do not reuse sequences | pass (400 distinct keys and first draws; separate from bootstrap streams) |
| Sample summaries agree with their target distributions | pass (KS at 0.1% for 120/120 run-inputs; P50 and P95 within sampling error) |
| No NaN, infinite or physically impossible value reaches the equations | pass (`sampling_validation.csv`; the equations also refuse such values) |
| Threshold probabilities are computed directly from iterations | pass (`test_reported_tables_reconcile_with_iterations`) |

## Limitations to carry into the final report

- The inputs are sampled independently. The sources give no joint data, and this is approved.
- The child BW covers **0–59 months only**, while ED is 6 years. The 60–71-month gap goes into the Phase 9 sensitivity analysis.
- Arsenic below the detection limit is extrapolated (D10). Report P50 and lower percentiles with this caveat.
- **Bogra's Lognormal tail runs high** (fitted mean 3× the empirical mean; child HI maximum 15,865), and **Rangpur's Gamma extreme tail is thin**. Both are selection-uncertainty items for Phase 9.
- The BGS data describe sampled wells, not a population-weighted census of wells.
- The Phase 4 and 5 selections were provisional. D6–D9 record Member 4's final choice, and the team should review them before the report is frozen.

Reviewer sign-off (Member 3): ________  Date: ________
