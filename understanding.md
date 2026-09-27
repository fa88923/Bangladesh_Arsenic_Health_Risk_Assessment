# Understanding the Bangladesh Arsenic Health-Risk Assessment Project

This is a self-contained study guide for the repository at commit `7dedfaf` on branch `abony-phase6` (reviewed 2026-09-26). It explains the reference paper, the Bangladesh adaptation, all ten planned phases, the data, mathematics, code, results, interpretation, and the important unresolved issues.

The shortest accurate description is:

> The project takes measured groundwater arsenic from selected Bangladesh districts, measured survey-weighted Bangladeshi adult and under-5 body weights, and exposure assumptions from Yadav and Kalkal's Punjab paper. It fits probability distributions to arsenic and body weight, repeatedly samples plausible combinations, evaluates ingestion and dermal risk equations, and reports the distribution and uncertainty of non-cancer and cancer risk.

The most defensible headline result is not a very large simulated maximum or mean. It is the estimated probability that the hazard index exceeds 1. In the primary model, that probability ranges from **10.0% to 72.3% for adults** and **20.3% to 85.1% for children**, depending on district. Arsenic concentration dominates sensitivity everywhere. However, the fitted arsenic distributions are uncertain because each district has only 44–110 sampled wells and many measurements are left-censored.

---

## 1. What is complete, and what is not

The master plan defines Phases 1–10.

| Phase | Purpose | Current status | Main evidence |
|---|---|---|---|
| 1 | Provenance and reproducibility | Implemented | `results/provenance/`, `config/` |
| 2 | Arsenic preprocessing | Implemented and tested | `reports/arsenic_preprocessing_report.md` |
| 3 | Body-weight preprocessing | Implemented and tested | `reports/bodyweight_report.md` |
| 4 | Arsenic fitting | Implemented, but no dedicated handoff report | `results/arsenic/`, fit code and tests |
| 5 | Body-weight fitting | Implemented; original selections were provisional | `reports/bodyweight_report.md` |
| 6 | Risk equations | Implemented and benchmarked | `reports/risk_equations_report.md` |
| 7 | Main Monte Carlo simulation | Implemented for 10 districts × 2 populations | `reports/simulation_report.md` |
| 8 | Convergence/numerical error | Implemented with 400 runs | `reports/convergence_report.md` |
| 9 | Sensitivity and uncertainty | Implemented | `reports/sensitivity_report.md` |
| 10 | Integrated final report and formal review | **Not complete** | No Phase 10 report; reviewer sign-offs remain blank |

Current verification:

- frozen raw-file hashes: **pass**;
- full test suite: **234 passed** in the repository virtual environment;
- Git working tree after verification: clean;
- phases 6–9 are in commit `7dedfaf`, one commit after `main` at `8de9e2a`.

Therefore, say **“the computational analysis through Phase 9 is implemented and tested”**, not **“the whole project is finalized and signed off.”**

---

## 2. The reference paper in plain language

The supplied `402project.pdf` and the checked-in `docs/references/Yadav_Kalkal_2024_JWH_arsenic_Punjab.pdf` have the same extracted textual content. Their binary hashes differ because their PDF metadata/modification records differ.

The paper is:

> Sangeeta Yadav and Sunil Kalkal (2024), “Health risk assessment using Monte-Carlo simulations due to arsenic contamination in groundwater in Punjab,” *Journal of Water and Health*, 22(12), 2304–2319.

### 2.1 What question did the paper ask?

The paper estimated non-cancer and cancer risks from arsenic-contaminated groundwater in seven Punjab districts. It compared:

1. a **deterministic calculation**, where every input is one fixed value; and
2. a **probabilistic calculation**, where uncertain or variable inputs have distributions and 10,000 Monte Carlo combinations are simulated.

The paper found that children generally had higher hazard indices than adults, most districts exceeded the non-cancer threshold, all districts exceeded the upper conventional cancer-risk benchmark, and arsenic concentration was the strongest sensitivity driver.

### 2.2 Why use Monte Carlo simulation?

A deterministic calculation answers something like:

> “What risk do we obtain for one average concentration, one body weight, and one intake rate?”

A Monte Carlo calculation answers:

> “If concentration, body weight, intake, exposure frequency, exposure time, and skin area vary, what distribution of possible risk values results?”

The algorithm draws one value from each input distribution, evaluates the equations, saves the output, and repeats this many times. This produces medians, upper percentiles, and threshold-exceedance probabilities rather than one number.

### 2.3 What this repository preserves from the paper

The repository preserves:

- the ingestion and dermal Average Daily Dose equations;
- ingestion and dermal Hazard Quotients;
- total Hazard Index;
- Excess Lifetime Carcinogenic Risk;
- the fixed toxicological constants;
- the paper's distributions for water intake, exposure frequency, exposure time, and skin area where usable;
- independent sampling of model inputs;
- a primary run size of 10,000;
- a Spearman rank-correlation sensitivity analysis;
- deterministic-versus-probabilistic result tables.

### 2.4 What this repository changes or improves

The Bangladesh project is **not a direct replication of Punjab results**. It adapts and extends the paper:

| Paper | Bangladesh repository |
|---|---|
| Seven Punjab district mean concentrations | Individual well observations from ten Bangladesh districts |
| Paper-specified adult/child body-weight distributions | Bangladesh STEPS adult and MICS7 under-5 survey data |
| Limited detail on non-detect handling | Explicit left-censoring-aware preprocessing and fitting |
| 10,000 simulations | 10,000 primary simulations plus convergence at 1k, 5k, 10k, and 20k with five seeds |
| Sensitivity of sampled inputs | Sensitivity plus fitted-parameter uncertainty and model scenarios |
| Child cancer averaging time equals the 6-year exposure duration | Primary child ELCR uses a 70-year lifetime averaging time; paper convention retained separately |
| Ambiguous “mean ± SD” lognormal notation | Treated as arithmetic mean and SD, then converted to log-space |
| Ambiguous child skin-area triangular distribution | EPA-informed `(0.29, 0.60, 0.95) m²` triplet |

Two changes are particularly important:

1. **Lognormal interpretation.** Reading `1.26 ± 0.66` as log-space parameters would imply a median intake of 3.5 L/day and a P95 of 10.4 L/day. Reading adult skin area `1.42 ± 0.31` that way gives a median of 4.14 m², larger than a plausible total adult body area. It also cannot reproduce the paper's own P95 result. The code therefore interprets the numbers as arithmetic mean and SD and converts them to log-space.
2. **Cancer averaging time.** For a lifetime cancer risk, the primary implementation averages exposure over 70 years. For children, exposure lasts six years, so their primary ELCR is reduced by `6/70` relative to the paper's convention. Both values are saved so the methodological difference is transparent.

---

## 3. The end-to-end pipeline

```text
Raw arsenic wells ──> classify detections/non-detects ──> censor-aware district fits ──┐
                                                                                       │
Raw adult survey ──> clean + survey weights ──> adult BW fit ───────────────────────────┤
                                                                                       ├─> Monte Carlo
Raw child survey ──> clean + survey weights ──> child BW fit ───────────────────────────┤   equations
                                                                                       │
Paper assumptions ──> IR, EF, ET, SA, ED, toxicological constants ─────────────────────┘
                                                                                       │
                         ┌──────────────────────────┬─────────────────────┬──────────────┘
                         v                          v                     v
                 risk summaries             convergence checks     sensitivity +
                 and exceedance              across N and seeds     parameter uncertainty
```

The unit of analysis in a simulation is a hypothetical exposure realization, not a known individual. A row combines a sampled well concentration with independently sampled exposure/body characteristics.

---

## 4. Repository map and execution order

### 4.1 Important directories

| Path | Role |
|---|---|
| `data/raw/` | Immutable source files. Never overwrite these. |
| `config/` | Frozen districts, model assumptions, seeds, candidate families, and parameter decisions. |
| `src/arsenic_hra/` | Reusable implementation modules. |
| `scripts/` | Phase entry points that load data, call modules, and write outputs/manifests. |
| `tests/` | Unit, data-contract, numerical, and end-to-end checks. |
| `results/` | Machine-readable derived data, fitted models, simulations, tables, figures, and manifests. |
| `reports/` | Human-readable handoff reports for completed phases. |
| `docs/` | Plans and detailed methodology written mostly before implementation. |

### 4.2 Dependency order

The actual dependency chain is:

1. `scripts/build_provenance.py`
2. `scripts/preprocess_arsenic.py`
3. `scripts/run_arsenic_fitting.py`
4. `scripts/run_bodyweight.py` (body-weight preprocessing and fitting together)
5. `scripts/run_risk_phase6.py`
6. `scripts/run_simulation.py`
7. `scripts/run_convergence.py`
8. `scripts/run_sensitivity.py`

Arsenic and body-weight work can run independently until the simulation stage.

### 4.3 Code modules by responsibility

- `paths.py`: root-aware paths and refusal to write under `data/raw/`.
- `provenance.py`: SHA-256 hashes and environment capture.
- `arsenic_preprocessing.py`: BGS layout, parsing, censoring classification, unit conversion, duplicate review.
- `arsenic_fitting.py`: reverse Kaplan–Meier CDF, CDF least squares, censored MLE, admissibility, selection.
- `arsenic_bootstrap.py`: well-level bootstrap and stability summaries.
- `bodyweight_preprocessing.py`: sequential adult/child cleaning and weighted summaries.
- `bodyweight_fitting.py`: weighted ECDF, weighted CDF least squares, pseudo-MLE, selection and sampling.
- `risk_parameters.py`: converts configuration to fixed and stochastic population parameters.
- `risk_equations.py`: equations (1)–(6), validation, deterministic inputs.
- `simulation_inputs.py`: final arsenic/BW choices, child truncated Normal, reverse-KM mean.
- `monte_carlo.py`: seeds, inverse-CDF sampling, equations, summaries, Wilson intervals, sampling checks.
- `convergence.py`: across-seed means, error, successive change, adequacy decisions.
- `sensitivity.py`: Spearman sensitivity, two-dimensional MC, model scenarios.
- `parameter_uncertainty.py`: arsenic bootstrap reuse and Rao–Wu survey bootstrap for body weight.
- `*_plots.py`: all generated diagnostic and presentation figures.

The scripts mostly orchestrate those modules, write CSV/JSON/Parquet/PNG artifacts, and record hashes and software versions.

---

## 5. Data: exactly what the project uses

Only three raw files are analytical inputs. Other MICS files are frozen as reference material but are not analyzed.

| Purpose | Raw source | Raw rows | Final analytical rows |
|---|---|---:|---:|
| Groundwater arsenic | DPHE/BGS/DFID National Hydrochemical Survey, released 2000 | 3,534 sample rows | 810 rows in ten selected districts |
| Adult body weight | WHO STEPS Bangladesh 2018 | 8,185 | 8,013 |
| Child body weight | UNICEF Bangladesh MICS7 2025 `ch.sav` | 24,680 | 22,576 |

### 5.1 Groundwater arsenic

The BGS CSV has four preamble lines, a header on line 5, a units row on line 6, and 3,534 sample rows from line 7 onward. Arsenic is stored in micrograms per litre (`µg/L`) and converted to milligrams per litre (`mg/L`) by dividing by 1,000.

The project purposively selects these exact district strings:

`Dhaka`, `Chittagong`, `Rajshahi`, `Khulna`, `Barisal`, `Sylhet`, `Rangpur`, `Mymensingh`, `Comilla`, and `Bogra`.

This is not a nationally representative district sample. It was chosen to balance geographic coverage, importance, sample count, and censoring burden.

#### Detected versus left-censored values

- A value such as `25` means arsenic was detected at 25 µg/L.
- A value such as `< 6` means the true value is somewhere below 6 µg/L. It does **not** mean 6, 3, or 0.
- The code stores the original string, status, exact detection if available, and censoring bound if not.
- No LOD, LOD/2, or zero substitution is used in the primary fit.

Of all 3,534 wells, 2,429 are detected and 1,105 are left-censored. In the ten districts, 524 are detected and 286 censored.

| District | n | Censored | Reverse-KM mean C (mg/L) | Final family | Main caution |
|---|---:|---:|---:|---|---|
| Dhaka | 45 | 48.9% | 0.0408 | Gamma | Nearly half censored; fitted P95 is about 2.1× empirical P95 |
| Chittagong | 44 | 11.4% | 0.0314 | Lognormal | Final fit under-predicts empirical P95 |
| Rajshahi | 78 | **71.8%** | 0.00725 | Gamma | Only 22 detections; strongest censoring warning |
| Khulna | 76 | 42.1% | 0.0349 | Gamma | Fitted P95 about 2.0× empirical P95 |
| Barisal | 92 | 17.4% | 0.0923 | **Gamma override** | Provisional Lognormal had an unsupported extreme tail |
| Sylhet | 77 | 32.5% | 0.0222 | Gamma | Low tail below detection is model-driven |
| Rangpur | 86 | 45.3% | 0.00807 | Gamma | Extreme Gamma tail may be too thin |
| Mymensingh | 108 | 48.1% | 0.0158 | Gamma | Nearly half censored |
| Comilla | 110 | 8.2% | **0.1416** | Gamma | Highest concentration and best-determined high-risk district |
| Bogra | 94 | 31.9% | 0.0176 | Lognormal | Heavy fitted tail; fitted mean about 3× reverse-KM mean |

There are no duplicate sample IDs or full duplicate rows. Repeated coordinates are retained because almost all differ in depth or date and can represent distinct wells or visits.

### 5.2 Adult body weight

Source fields:

- `m12`: measured weight in kg — this becomes `BW_kg`;
- `wstep2`: survey weight — this weights observations but is not a body weight;
- age 18–69, sex, participant ID, PSU, stratum, division, residence, and diagnostic height.

Cleaning removes 167 missing weights and five `m12 = 888` refusal codes. It does not automatically trim valid extremes. Twelve unusual-BMI records are listed for review rather than deleted.

Survey-weighted adult summary:

- n = 8,013; Kish effective n ≈ 3,238;
- mean 55.89 kg; SD 11.62 kg;
- P5 39 kg; median 55 kg; P95 76 kg;
- observed range 28–162 kg.

The selected adult distribution is:

\[
BW_{adult} \sim \operatorname{Lognormal}(\mu_{\log}=4.002396,\ \sigma_{\log}=0.200174).
\]

It won weighted CDF SSE, maximum discrepancy, and tail agreement.

### 5.3 Child body weight

Source fields:

- `AN8`: measured weight in kg;
- `chweight`: survey weight;
- `CAGE`: age 0–59 months;
- sex, district, PSU, stratum, and WHO anthropometric flags.

Cleaning separately removes:

- 1,323 missing weights;
- special codes `99.3`, `99.4`, `99.5`, `99.6`;
- 57 remaining `WAZFLAG = 1` weight-for-age errors.

Height-related error flags are retained because height is not the fitted outcome.

Survey-weighted child summary:

- n = 22,576; Kish effective n ≈ 16,600;
- mean 10.87 kg; SD 3.22 kg;
- P5 5.6 kg; median 10.9 kg; P95 15.9 kg;
- machine-readable observed range **1.6–31.5 kg**.

Phase 5 provisionally selected Gamma because the untruncated Normal assigned `3.7×10⁻⁴` probability to non-positive weights. Phase 7 explicitly introduced and refitted a truncated Normal because it fits the empirical distribution much better and represents the light-weight tail that drives high dose per kg:

\[
BW_{child} \sim \operatorname{Normal}(\mu=10.8570,\sigma=3.2281)\quad\text{conditioned on }BW\ge1.6\text{ kg}.
\]

Its weighted CDF SSE is about 11 times smaller than Gamma's, its worst P1–P99 quantile error is 1.7% instead of 36%, and it reproduces empirical `E[1/BW]` within 0.1%.

Important scope: these weights represent **children aged 0–59 months**, while the retained child exposure duration is six years. Phase 9 estimates that extending the model to 60–71 months would reduce child upper-tail risk by about 5.5–8%, so the under-5-only model is slightly conservative.

---

## 6. How the distributions are fitted

### 6.1 Why fit distributions at all?

Monte Carlo simulation requires a way to generate new plausible values. A fitted distribution smooths the observed sample and supplies an inverse CDF `F⁻¹(u)` for sampling. This also extrapolates beyond observed data, which is why upper-tail family choice is a major uncertainty.

### 6.2 Arsenic: reverse Kaplan–Meier plus CDF least squares

Ordinary empirical CDFs require exact values. A non-detect `<L` supplies only `C < L`. The reverse Kaplan–Meier estimator constructs a non-parametric CDF while respecting those upper bounds.

For each district and candidate family, the primary fit minimizes:

\[
SSE(\theta)=\sum_j \left[F(x_j;\theta)-\widehat F_{RKM}(x_j)\right]^2.
\]

The flat initial plateau below the smallest resolved detection is excluded because it contains no shape information and includes a software placeholder at zero. Candidates are:

- Normal `N(µ, σ)`;
- Lognormal with log-space `(µ_log, σ_log)` and `loc = 0`;
- Gamma with shape `k`, scale `θ`, and `loc = 0`.

A censored maximum-likelihood fit is used as a robustness comparison:

\[
\log L(\theta)=\sum_{i\in detected}\log f(x_i;\theta)
+\sum_{i\in censored}\log F(L_i;\theta).
\]

For a censored value, the likelihood contribution is the probability of being below its bound, not a density at the bound.

Selection uses the lowest CDF RMSE among physically admissible candidates, then checks maximum CDF discrepancy, censored-MLE AICc, quantiles, tails, diagnostics, and 1,000 well-level bootstrap refits. Normal is rejected in every district because it gives too much probability to negative concentration.

Final simulation families are Gamma for eight districts and Lognormal for Chittagong and Bogra. Barisal is deliberately changed from provisional Lognormal to Gamma because the Lognormal predicted:

- P99 = 19.6 mg/L while the observed maximum was 0.862 mg/L;
- fitted mean 2.71 mg/L versus reverse-KM mean 0.092 mg/L;
- P99.9 = 291 mg/L;
- 7.5% of draws above every observed Barisal well.

That tail would be mathematically allowed but scientifically unsupported.

### 6.3 Body weight: survey-weighted empirical CDF

Survey weights matter because rows do not represent equal numbers of people. The normalized mass is:

\[
q_i=\frac{w_i}{\sum_r w_r}.
\]

After aggregating tied weights, the midpoint plotting position for distinct value `x_j` is:

\[
p_j=\sum_{r<j}q_r+\frac{q_j}{2}.
\]

The primary weighted CDF objective is:

\[
SSE_w(\theta)=\sum_j q_j\left[F(x_j;\theta)-p_j\right]^2.
\]

Normal, Lognormal, Gamma, and Triangular candidates are fitted separately for adults and children. A survey-weighted pseudo-likelihood is a robustness check, not ordinary iid likelihood; therefore the project correctly does not report ordinary AIC for these fits.

Physical admissibility rules reject:

- Normal when `P(BW ≤ 0) ≥ 10⁻⁴`;
- Triangular when its support excludes more than 0.1% of survey mass;
- any failed/non-finite fit.

### 6.4 Lognormal conversion used for paper inputs

If a Lognormal variable has arithmetic mean `m` and standard deviation `s`, its log-space parameters are:

\[
\sigma_{\log}^2=\ln\left(1+\frac{s^2}{m^2}\right),
\qquad
\mu_{\log}=\ln(m)-\frac{\sigma_{\log}^2}{2}.
\]

That conversion is applied to paper-reported water intake and adult skin area.

---

## 7. Model inputs used in the risk equations

### 7.1 Stochastic inputs sampled each iteration

| Symbol | Meaning | Primary source/distribution |
|---|---|---|
| `C` | arsenic concentration | district-specific fitted Gamma or Lognormal, mg/L |
| `IR` | drinking-water ingestion rate | Lognormal with arithmetic mean 1.26, SD 0.66 L/day |
| `BW` | body weight | adult Lognormal or child truncated Normal, kg |
| `EF` | exposure frequency | Triangular(180, 345, 365) days/year |
| `ET` | dermal exposure time | Triangular(0.13, 0.20, 0.33) h/day |
| `SA` | exposed skin area | adult Lognormal(mean 1.42, SD 0.31); child Triangular(0.29, 0.60, 0.95) m² |

They are sampled independently because no joint dataset links a person's water source, intake, weight, and behavior. This is a modeling assumption, not a discovered fact.

### 7.2 Fixed inputs

| Symbol | Meaning | Adult | Child |
|---|---|---:|---:|
| `ED` | exposure duration | 70 years | 6 years |
| `AT_nc` | non-cancer averaging time | `70×365` days | `6×365` days |
| `AT_cancer` | primary cancer averaging time | `70×365` days | `70×365` days |
| `Kp` | dermal permeability | 0.001 cm/h | 0.001 cm/h |
| `CF` | dermal unit conversion | 10 L/(cm·m²) | same |
| `RfD_ing` | ingestion reference dose | 0.0003 mg/(kg·day) | same |
| `RfD_dermal` | dermal reference dose | 0.000285 mg/(kg·day) | same |
| `CSF` | cancer slope factor | 1.5 (mg/(kg·day))⁻¹ | same |

The conversion factor is dimensionally valid because `1 cm × 1 m² = 0.01 m³ = 10 L`.

Deterministic-only paper inputs are IR 2.0/1.8 L/day (adult/child), EF 365, ET 0.58 h/day, and SA 1.8/0.6 m². The paper's deterministic ET is oddly outside its own probabilistic 0.13–0.33 h/day range; the repository retains this published inconsistency transparently.

---

## 8. The health-risk mathematics

### 8.1 Average Daily Dose by ingestion

\[
ADD_{ing}=\frac{C\times IR\times EF\times ED}{BW\times AT}.
\]

Interpretation: arsenic concentration times water consumed and exposure frequency/duration, divided by body weight and averaging time.

Increasing `C`, `IR`, `EF`, or `ED` increases dose. Increasing `BW` or `AT` decreases dose.

For non-cancer risk, `AT = ED×365`; therefore ED algebraically cancels. The code intentionally retains the unsimplified equation so every paper factor remains visible and testable.

### 8.2 Average Daily Dose by dermal contact

\[
ADD_{dermal}=\frac{C\times K_p\times EF\times ED\times ET\times SA\times CF}
{BW\times AT}.
\]

This is exposure through skin contact with water. In this project it contributes only about 0.1–0.3% of mean HI, so ingestion overwhelmingly dominates.

### 8.3 Hazard Quotients and Hazard Index

\[
HQ_{ing}=\frac{ADD_{ing}}{RfD_{ing}},
\qquad
HQ_{dermal}=\frac{ADD_{dermal}}{RfD_{dermal}},
\]

\[
HI=HQ_{ing}+HQ_{dermal}.
\]

Interpretation:

- `HI < 1`: modeled exposure is below the combined reference-dose benchmark;
- `HI > 1`: non-cancer effects cannot be ruled out and concern increases;
- HI is **not** a probability of illness and `HI = 2` does not mean twice the disease incidence.

### 8.4 Excess Lifetime Carcinogenic Risk

\[
ELCR=ADD_{ing,cancer}\times CSF.
\]

The conventional comparison range is `10⁻⁶` to `10⁻⁴`; this project uses `10⁻⁴` as the exceedance threshold.

ELCR is a modeled incremental lifetime probability under the assumptions, not an observed cancer rate. For example, `ELCR = 10⁻³` is conventionally read as about one modeled excess case per 1,000 similarly exposed people, but it should not be presented as an epidemiological prediction.

### 8.5 A useful deterministic sanity check

At `C = 0.01 mg/L` (10 µg/L, the WHO guideline), using Bangladesh mean body weights and paper deterministic exposures:

| Population | HI | Primary lifetime ELCR |
|---|---:|---:|
| Adult | 1.199 | `5.37×10⁻⁴` |
| Child | 5.532 | `2.13×10⁻⁴` |

At `C = 0.05 mg/L` (50 µg/L, the Bangladesh standard used in the benchmark):

| Population | HI | Primary lifetime ELCR |
|---|---:|---:|
| Adult | 5.997 | `2.68×10⁻³` |
| Child | 27.66 | `1.06×10⁻³` |

The child HI is much higher because children have far lower body weight while modeled water intake is comparable. Child primary ELCR is not proportionally as high because only six exposure years are averaged over a 70-year lifetime.

---

## 9. Phase 6: equation implementation and validation

Phase 6 contains no random simulation. It establishes that equations and units are correct before any probabilistic result is trusted.

Validation includes:

- hand-calculated adult and child cases;
- SI-unit recomputation of dermal and ingestion dose;
- tests that doubling each input changes only the correct pathway by the expected factor;
- array/scalar identity checks;
- rejection of negative, zero where impossible, NaN, infinity, EF > 366, and ET > 24;
- reproduction of the paper's deterministic Table 4.

The software reproduces 25 of 28 printed paper cells within 0.5%. Three ELCR cells disagree with direct arithmetic and are likely paper typos:

- Ferozepur child: recomputed `4.86×10⁻³`, printed `4.80×10⁻³`;
- Amritsar adult: recomputed `3.56×10⁻³`, printed `3.50×10⁻³`;
- Amritsar child: recomputed `1.49×10⁻²`, printed `1.40×10⁻²`.

This is evidence that the repository preserved the equations rather than blindly copying printed outputs.

---

## 10. Phase 7: Monte Carlo simulation

### 10.1 One simulation run

For each district and population:

1. derive a reproducible random stream from master seed `20260923`;
2. draw six independent uniform arrays in fixed order: `C, IR, BW, EF, ET, SA`;
3. transform each with its inverse CDF, `X = F⁻¹(U)`;
4. evaluate all risk equations element-wise;
5. save all six sampled inputs and seven risk outputs for 10,000 iterations;
6. summarize percentiles and threshold counts directly from the iterations.

The seed key is:

\[
(district\ index,\ population\ index,\ N\ index,\ replicate).
\]

This gives distinct, reproducible streams and prevents accidental reuse across cells.

### 10.2 Why inverse-CDF sampling?

If `U ~ Uniform(0,1)`, then `F⁻¹(U)` follows CDF `F`. Using this explicitly makes the sampling contract transparent and lets Phase 9 reuse the same uniforms to isolate changes due to fitted parameters or model choices.

### 10.3 Primary non-cancer results

`C̄_RKM` is the censoring-aware mean used only in the deterministic benchmark. P95 is the 95th percentile of the 10,000 simulated risk values. `P(HI>1)` is the direct fraction of iterations exceeding 1.

| District | `C̄_RKM` mg/L | Deterministic HI adult | Deterministic HI child | P95 HI adult | P95 HI child | P(HI>1) adult | P(HI>1) child |
|---|---:|---:|---:|---:|---:|---:|---:|
| Dhaka | 0.0408 | 4.89 | 22.56 | 21.6 | 121.1 | 33.5% | 47.7% |
| Chittagong | 0.0314 | 3.77 | 17.38 | 6.66 | 34.7 | 23.2% | 53.8% |
| Rajshahi ⚠ | 0.00725 | 0.87 | 4.01 | 3.17 | 18.2 | 11.2% | 20.3% |
| Khulna | 0.0349 | 4.19 | 19.31 | 15.7 | 84.3 | 30.0% | 44.9% |
| Barisal | 0.0923 | 11.07 | 51.07 | 32.1 | 175.1 | 41.2% | 55.8% |
| Sylhet | 0.0222 | 2.66 | 12.29 | 6.61 | 35.8 | 32.8% | 60.0% |
| Rangpur | 0.00807 | 0.97 | 4.47 | 1.84 | 9.62 | 10.0% | 29.9% |
| Mymensingh | 0.0158 | 1.90 | 8.74 | 7.19 | 37.9 | 19.1% | 32.0% |
| Comilla | **0.1416** | **16.98** | **78.33** | **51.3** | **259.6** | **72.3%** | **85.1%** |
| Bogra | 0.0176 | 2.11 | 9.72 | 6.57 | 33.3 | 15.5% | 32.4% |

⚠ Rajshahi has 71.8% censored wells.

### 10.4 How to read this table

- **Comilla is clearly the highest-risk district** in the selected set and also has relatively low censoring, so this ranking is not just a high-censoring artifact.
- Barisal, Dhaka, and Khulna also have high upper-tail risks.
- Rangpur and Rajshahi are lowest, but “lowest” does not mean safe: their child P95 HI is still 9.62 and 18.2.
- Child P95 HI is about 5.1–5.7 times adult P95 HI in every district.
- Probabilistic P95 is 1.8–5.4 times the deterministic estimate because arsenic is strongly right-skewed.
- The probabilistic median is often below the deterministic mean-based estimate. A deterministic average can therefore understate the upper tail while overstating the typical case.

### 10.5 Cancer-risk results

Primary P95 ELCR values are:

| District | Adult P95 ELCR | Child P95 ELCR, lifetime AT | Child P95 ELCR, paper AT |
|---|---:|---:|---:|
| Dhaka | `9.68×10⁻³` | `4.67×10⁻³` | `5.45×10⁻²` |
| Chittagong | `2.99×10⁻³` | `1.34×10⁻³` | `1.56×10⁻²` |
| Rajshahi | `1.42×10⁻³` | `6.99×10⁻⁴` | `8.16×10⁻³` |
| Khulna | `7.03×10⁻³` | `3.25×10⁻³` | `3.79×10⁻²` |
| Barisal | `1.44×10⁻²` | `6.75×10⁻³` | `7.87×10⁻²` |
| Sylhet | `2.96×10⁻³` | `1.38×10⁻³` | `1.61×10⁻²` |
| Rangpur | `8.27×10⁻⁴` | `3.71×10⁻⁴` | `4.33×10⁻³` |
| Mymensingh | `3.23×10⁻³` | `1.46×10⁻³` | `1.70×10⁻²` |
| Comilla | `2.30×10⁻²` | `1.00×10⁻²` | `1.17×10⁻¹` |
| Bogra | `2.95×10⁻³` | `1.28×10⁻³` | `1.50×10⁻²` |

All are above `10⁻⁴` at P95. The primary child ELCR can be lower than adult ELCR despite higher child dose per kg because the child's six exposure years are averaged over a 70-year lifetime. Under the paper's exposure-duration averaging convention, child ELCR becomes 5.1–5.7 times adult ELCR instead.

This difference is methodological, not a contradiction in the simulation.

---

## 11. Phase 8: convergence and Monte Carlo error

### 11.1 What convergence means here

Monte Carlo output changes slightly with random seed and iteration count. Phase 8 repeats the full model for:

- 10 districts;
- 2 populations;
- `N = 1,000, 5,000, 10,000, 20,000`;
- 5 independent seeds per N.

That is **400 runs and 3.6 million iterations**.

For statistic `T`, it computes the across-seed mean, between-seed standard deviation, coefficient of variation, and successive relative change:

\[
\Delta_T(N)=\frac{|\overline T_N-\overline T_{N_{prev}}|}{|\overline T_N|}\times100\%.
\]

For an exceedance probability `p̂`, the binomial Monte Carlo standard error is:

\[
MCSE(\hat p)=\sqrt{\frac{\hat p(1-\hat p)}{N}}.
\]

### 11.2 Adequacy rule

- Continuous statistics: successive change ≤ 5% and between-seed CV ≤ 5%.
- Probabilities: change ≤ 0.01 and between-seed SD ≤ 0.01.
- “Smallest adequate N” must pass at that N and every larger tested N.

These tolerances are project decisions, not universal laws.

### 11.3 Conclusions

- `N = 10,000` is adequate for **P(HI>1) in all 20 cells**.
- It is adequate for P95 HI and P95 ELCR in 16 of 20 cells.
- `N = 20,000` is needed for P95 in Bogra adult/child, Chittagong child, and Sylhet adult.
- Bogra's mean HI does not converge even at 20,000 because its heavy Lognormal tail makes rare draws dominate the mean; about 1.4 million iterations would be needed for a 5% mean CV.
- Several medians are unstable because they lie below detection limits in an extrapolated near-zero region. Most are far below HI = 1 and do not affect the risk conclusion.
- The maximum MC standard error for the primary `P(HI>1)` estimates is about 0.5 percentage points.

The key distinction:

> Monte Carlo error is randomness from simulating a fixed fitted model; it shrinks with N. Parameter uncertainty comes from having limited well/survey data; it does not shrink by simply simulating more.

---

## 12. Phase 9: sensitivity and uncertainty

Phase 9 answers three different questions.

### 12.1 Which varying inputs drive risk?

It computes tie-aware Spearman rank correlation between each sampled input and HI/ELCR using the exact saved Phase 7 iterations.

| Input | Spearman ρ with HI across cells | Meaning |
|---|---|---|
| Arsenic `C` | **+0.950 to +0.998** | Rank 1 everywhere; almost all rank variation |
| Ingestion `IR` | +0.054 to +0.231 | Rank 2 everywhere |
| Body weight `BW` | −0.022 to −0.157 | Higher weight lowers dose per kg |
| Exposure frequency `EF` | +0.014 to +0.070 | Small; often not significant |
| Exposure time `ET` and skin area `SA` | −0.020 to +0.024 | Noise-level for HI because dermal risk is tiny |

Fixed ED, Kp, CF, reference doses, and slope factor are excluded because a constant has no rank correlation.

The Bangladesh arsenic correlation is higher than the paper's roughly 0.65–0.90 because fitted Bangladesh district concentrations span orders of magnitude. This supports the intervention conclusion: reducing arsenic in water matters much more than changing dermal-exposure assumptions.

### 12.2 How uncertain are the fitted C and BW distributions?

The project performs two-dimensional Monte Carlo:

1. create about 1,000 fitted parameter sets;
2. refit arsenic using well-level bootstrap samples;
3. refit body weight using a Rao–Wu PSU-within-stratum survey bootstrap;
4. reuse the primary simulation uniforms for every parameter set;
5. measure how results change across fitted parameter sets.

Reusing the same uniforms is called **common random numbers**. It prevents ordinary simulation noise from obscuring differences caused by parameter estimates.

P95 HI parameter-uncertainty intervals:

| District | Adult primary [95% interval] | Child primary [95% interval] |
|---|---|---|
| Dhaka | 21.6 [10.2, 31.2] | 121 [58, 174] |
| Chittagong | 6.7 [1.6, 24.8] | 35 [9, 123] |
| Rajshahi | 3.2 [1.1, 6.2] | 18 [7, 35] |
| Khulna | 15.7 [7.5, 24.6] | 84 [40, 133] |
| Barisal | 32.1 [9.6, 51.8] | 175 [52, 284] |
| Sylhet | 6.6 [4.1, 9.7] | 36 [22, 52] |
| Rangpur | 1.8 [1.0, 3.0] | 9.6 [5.0, 15.6] |
| Mymensingh | 7.2 [2.7, 11.0] | 38 [14, 59] |
| Comilla | 51.3 [40.6, 65.1] | 260 [212, 330] |
| Bogra | 6.6 [0.8, 13.8] | 33 [4, 70] |

These intervals are much wider than the 2–9% Monte Carlo CV. The analysis is **data-limited, not simulation-limited**.

Body-weight parameter uncertainty is negligible relative to arsenic uncertainty because the body-weight samples contain thousands of people while districts have only dozens of wells.

### 12.3 How much do modeling choices matter?

Six scenarios use common random numbers:

| Scenario | Main effect |
|---|---|
| Barisal Gamma → provisional Lognormal | P95 ×3.14; mean ×17–46; P(HI>1) changes little |
| Bogra Lognormal → Gamma | P95 ×0.57–0.61; mean ×0.19–0.25 |
| Rangpur Gamma → Lognormal | P95 ×1.87; mean ×3.9–5.1 |
| Child truncated Normal → Gamma | Child P95 falls about 1.5–5% |
| Extend child BW to 0–71 months | Child P95 falls 5.5–8%; P(HI>1) falls 0–3% |
| Read IR/SA `±` values as log-space | P95 increases 3.2–3.5× |

The important conclusion is:

> Family choice strongly changes means and upper percentiles because they depend on extrapolated tails, but it changes P(HI>1) much less because the threshold lies in the better-observed middle of the concentration distribution.

That is why `P(HI>1)` is the preferred headline metric.

---

## 13. What the figures show

You do not need to memorize every PNG. Know the purpose of each class:

- `arsenic_ecdf_cdf_<district>.png`: reverse-KM empirical CDF versus candidate fitted CDFs; shows central fit and censoring plateau.
- `arsenic_histogram_pdf_<district>.png`: observed detections plus fitted densities; censored observations must not be read as exact histogram values.
- `arsenic_qq_<district>.png`: fitted versus empirical quantiles; deviations reveal poor distributional fit.
- `arsenic_tails_<district>.png`: close view of lower and upper tails, crucial for risk extrapolation.
- `bodyweight_*`: analogous weighted diagnostics for adults and children.
- `simulation_hi_distribution_<district>.png`: adult and child HI on a log scale, with HI = 1 and P95 marked. Near-zero draws are pooled in the first bin.
- `convergence_*_vs_N.png`: across-seed mean and seed range as N grows.
- `convergence_mc_error_decay.png`: error falls roughly as `N⁻¹ᐟ²`, as expected.
- `sensitivity_tornado_<district>.png`: Spearman coefficients; arsenic bars dominate.
- `sensitivity_parameter_uncertainty_*`: primary values plus bootstrap uncertainty intervals.

No reported number is manually read from a plot. Tables are generated from machine-readable iteration results.

---

## 14. What the results do and do not mean

### 14.1 Defensible statements

- Within the ten purposively selected districts and this model, a substantial modeled share of exposures exceeds HI = 1.
- Children have substantially higher modeled non-cancer risk than adults because dose is normalized by much lower body weight.
- Comilla is the clearest high-risk district in this set.
- Arsenic concentration is by far the dominant source of input variability.
- Upper-tail risk estimates are much less certain than the simulation's numerical error.
- More well data and better handling of the exposure joint distribution would improve confidence more than merely increasing N from 10,000.

### 14.2 Statements to avoid

Do not say:

- “72.3% of Comilla adults will become ill.” It means 72.3% of modeled iterations exceeded a reference-dose index, not observed disease incidence.
- “HI = 10 means ten times as many illnesses.” HI is not a probability or linear disease-response measure.
- “ELCR proves that exactly X people will get cancer.” It is a model-derived incremental risk under assumptions.
- “The ten districts represent all Bangladesh.” Selection is purposive, and wells are not person-weighted.
- “The simulation proves the Gamma/Lognormal tail.” Distribution tails are extrapolated from limited data.
- “10,000 runs remove uncertainty.” They reduce numerical error, not uncertainty in data, fitted distributions, or assumptions.
- “Censored `<6` values are 6 or 3.” Only the upper bound is known.

### 14.3 Why children have higher HI but lower primary ELCR

This is a likely evaluation question.

- HI uses a non-cancer averaging time equal to exposure duration, so the six-year ED cancels. Lower child body weight dominates and raises HI.
- Primary ELCR averages six years of child exposure across a 70-year lifetime. That multiplies the child lifetime-average dose by `6/70`.
- If the paper's ED-based averaging convention is used, child ELCR again becomes about 5.1–5.7 times adult ELCR.

---

## 15. Main scientific and reproducibility limitations

1. **Independent inputs.** Real intake, weight, behavior, and water source may be correlated. Direction of bias is unknown.
2. **Under-5 versus six-year exposure.** Quantified as a modest conservative bias of about 5.5–8% in child P95.
3. **Sub-detection concentrations are unresolved.** The fitted curve extrapolates below the smallest detection; medians/lower percentiles are especially model-driven.
4. **High censoring.** Rajshahi is 71.8% censored; Dhaka and Mymensingh are near 50%.
5. **Tail-family dependence.** Barisal, Bogra, Rangpur, and Chittagong show material tail concerns.
6. **Paper exposure assumptions are not Bangladesh-specific.** IR, EF, ET, and SA are inherited or repaired from the Punjab model.
7. **Well sample is not population weighted.** A modeled draw represents a sampled-well concentration distribution, not the fraction of residents actually using each well.
8. **Only arsenic is modeled.** Co-contaminants and arsenic speciation are absent.
9. **Cancer model is simplified.** ELCR uses a linear slope-factor model and only ingestion.
10. **Reviewer sign-off is incomplete.** Several decisions were made in Phase 6–9 but formal cross-review fields are blank.
11. **Environment drift exists.** Phase 4 ran with numpy 2.5.3/scipy 1.18.1, while requirements pin numpy 1.26.4/scipy 1.17.1. Regenerated bootstrap summaries agree to about `1.1×10⁻⁶` relative, not bit for bit.
12. **Phase 10 is missing.** The project has strong handoff reports but no reconciled final scientific report.

---

## 16. Documentation inconsistencies you should know

The machine-readable results and current code should outrank stale planning text.

1. `README.md` says `scripts/` currently only has the audit and `results/` is planned. That is obsolete; phases 1–9 are present.
2. The master plan/design document still says paper Lognormal `±` values are direct log-space parameters. Phase 6 reversed this for strong physical and reproduction reasons. Current config/code use arithmetic mean/SD conversion.
3. Some design text says MICS 2019. The actual source is **MICS7 2025**.
4. `config/run_config.json` still labels itself draft and points to `docs/parameter_decision_register.md` and `docs/phase6_risk_equation_design.md`, which do not exist under those names.
5. There is no dedicated Phase 4 report even though the code, outputs, diagnostics, bootstrap, manifests, and tests exist.
6. `reports/bodyweight_report.md` says the child range is 1.0–36.0 kg. The saved clean data and weighted-summary CSV say **1.6–31.5 kg**. This mismatch predates phases 6–9.
7. Phase 5's selected child family is Gamma, but the final Phase 7 simulation uses a newly justified truncated Normal. This is deliberate and documented, not accidental.
8. Phase 4's Barisal selection is Lognormal, but final simulation uses Gamma. Again, this is a documented tail-safety override.
9. `config_version` remains 1.0.0 even after major risk-model additions; newer manifests compensate by recording the full config SHA-256.

If asked which value/model is “final for the simulation,” use `results/simulation/simulation_input_distributions.csv`, not the earlier provisional selection files.

---

## 17. How to reproduce and inspect the project safely

Use the repository virtual environment; the system `python3` currently lacks pytest.

```bash
.venv/bin/python scripts/build_provenance.py --verify
.venv/bin/python -m pytest -q
```

To regenerate in dependency order:

```bash
.venv/bin/python scripts/preprocess_arsenic.py
.venv/bin/python scripts/run_arsenic_fitting.py
.venv/bin/python scripts/run_bodyweight.py
.venv/bin/python scripts/run_risk_phase6.py
.venv/bin/python scripts/run_simulation.py
.venv/bin/python scripts/run_convergence.py
.venv/bin/python scripts/run_sensitivity.py --reuse-bootstrap
```

The full sensitivity run regenerates expensive bootstrap tables:

```bash
.venv/bin/python scripts/run_sensitivity.py
```

Important: regenerating outputs may expose the known cross-environment numerical drift. Do not overwrite raw data. Manifests and SHA-256 hashes provide the traceability chain.

### Best files to open during revision

1. `understanding.md` — this guide.
2. `reports/simulation_report.md` — main Phase 7 decisions and results.
3. `reports/risk_equations_report.md` — equations and deviations from the paper.
4. `reports/sensitivity_report.md` — what drives risk and how uncertain it is.
5. `reports/convergence_report.md` — why 10,000 was retained.
6. `reports/limitations_register.md` — every caveat.
7. `results/tables/risk_exceedance.csv` — headline probabilities.
8. `results/sensitivity/parameter_uncertainty_summary.csv` — confidence in fitted-model results.
9. `results/simulation/simulation_input_distributions.csv` — actual final distributions.

---

## 18. Likely evaluation questions and strong short answers

### What is the novelty of this repository relative to the Punjab paper?

It applies the same risk-equation framework to Bangladesh-specific well and body-weight data, handles non-detects explicitly, uses survey weights, validates convergence, and separates variability, numerical error, parameter uncertainty, and model-choice uncertainty.

### Why not replace `<6` with 3?

Because `<6` only gives an upper bound. LOD/2 substitution creates invented exact observations and biases the distribution, especially where multiple detection limits exist. Reverse Kaplan–Meier and censored likelihood use the inequality directly.

### Why are Gamma shapes below 1 acceptable?

A Gamma with `k < 1` has high density near zero but remains a valid non-negative distribution. That shape is plausible here because many wells are non-detects. The scale must be positive; shape need not exceed 1.

### Why was Normal rejected for arsenic?

It assigns substantial probability to impossible negative concentrations in every district.

### Why use survey weights for body weight?

STEPS and MICS use complex survey sampling; records represent different numbers of people. Ignoring weights estimates the sample distribution rather than the target population distribution.

### Why was child Gamma replaced?

The untruncated Normal fit the observed weighted CDF about 11 times better but had tiny negative mass. Truncating at the observed minimum removes impossible values, fits tails much better, and reproduces `E[1/BW]`, the quantity directly relevant to dose.

### Why was Barisal Lognormal replaced?

Its modestly better central CDF fit created an unsupported extreme tail—P99 19.6 mg/L and mean 2.71 mg/L against a maximum observation of 0.862 and reverse-KM mean of 0.092. Gamma sacrifices little central fit and avoids tail-driven risk inflation.

### Why is P(HI>1) preferred over mean or P95?

The threshold lies in a better-observed part of the concentration distribution. Means and P95 are much more sensitive to extrapolated tail family. Scenario changes moved means by up to 46× but P(HI>1) by at most about 19%.

### Why keep N = 10,000 if four P95 cells need 20,000?

Ten thousand is the predeclared, paper-comparable primary size and is adequate for the headline P(HI>1) in all cells. The 20,000 results are reported as stability checks; raising N would not solve the much larger fitted-parameter uncertainty.

### What is the biggest uncertainty?

District arsenic distributions. Body-weight data are large, and Monte Carlo error is small. Limited wells, censoring, and tail-family choice dominate.

### What intervention follows from the sensitivity analysis?

Reduce arsenic concentration in drinking water—through testing, safe-source substitution, treatment, and monitoring. Arsenic has Spearman `ρ = 0.95–0.998`; dermal parameters barely affect total HI.

### Is the dermal pathway important?

It is included for fidelity to the model, but contributes only about 0.1–0.3% of mean HI. Ingestion dominates.

### What would you improve next?

Complete Phase 10, obtain formal sign-off, reconcile stale docs/config version/environment pins, collect more representative wells, obtain Bangladesh-specific intake and skin-area inputs, model dependence where joint data permit, and externally validate distribution/tail choices.

---

## 19. A five-minute oral summary

This project adapts a 2024 Punjab arsenic health-risk paper to Bangladesh. It uses 810 groundwater measurements from ten selected BGS districts, 8,013 weighted adult observations from STEPS 2018, and 22,576 weighted under-5 observations from MICS7 2025. Non-detect arsenic values are treated as left-censored, not substituted, and district distributions are fitted with reverse-KM CDF least squares plus censored-MLE and bootstrap checks. Adult weight is Lognormal; final child weight is a Normal truncated at 1.6 kg.

Each district and population gets 10,000 independent inverse-CDF samples of concentration, intake, body weight, exposure frequency/time, and skin area. The paper's ADD, HQ, HI, and ELCR equations are then evaluated. The code reproduces almost all of the paper's deterministic table, with three likely paper typos. It corrects two paper ambiguities: lognormal mean/SD are interpreted as arithmetic moments, and primary cancer risk uses a 70-year averaging time.

Comilla has the highest modeled risk: 72.3% of adult and 85.1% of child iterations exceed HI = 1. Across districts, adult exceedance ranges 10.0–72.3% and child exceedance 20.3–85.1%. Children have P95 HI about five times adults because of lower body weight. Arsenic concentration dominates sensitivity with Spearman correlations of 0.95–0.998. Ten thousand iterations are adequate for the headline exceedance probability, but uncertainty in fitted arsenic distributions is far larger than simulation error. Therefore the safest conclusion is that risk concern is substantial in the selected districts, particularly Comilla, but exact tail magnitudes are uncertain and should not be interpreted as observed disease rates or nationally representative prevalence.

---

## 20. Final mental model

Remember four layers:

1. **Data layer:** wells and survey body weights, with censoring and survey design handled explicitly.
2. **Distribution layer:** smooth models let the project generate plausible concentrations and weights; tails introduce uncertainty.
3. **Risk layer:** dose per kg → quotient relative to reference dose → HI; lifetime-average ingestion dose × slope factor → ELCR.
4. **Uncertainty layer:** Monte Carlo error is small; arsenic parameter and family uncertainty are large.

If you can explain those four layers, the child/adult averaging-time distinction, why censored values are not substituted, and why `P(HI>1)` is the most defensible headline, you understand the core of the project.
