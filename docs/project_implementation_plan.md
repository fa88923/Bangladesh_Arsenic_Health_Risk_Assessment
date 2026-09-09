# Bangladesh Arsenic Health-Risk Model: Master Project Plan

## 1. Purpose

This document organizes the existing dataset audit, arsenic methodology, body-weight methodology, and health-risk model design into one implementation sequence. The approved project stance is to apply the Punjab model of Yadav and Kalkal (2024) to Bangladesh, retaining its model structure while using the approved Bangladesh data and fitted inputs. It defines what must be decided, produced, checked, and approved at each stage.

This is a planning document only. It does not perform preprocessing, distribution fitting, simulation, or result generation.

## 2. Planning basis

The plan is based on these repository artifacts:

- `reports/dataset_audit_report.md`
- `reports/audit_outputs/district_counts.csv`
- `reports/audit_outputs/completeness_tables.csv`
- `docs/arsenic_preprocessing_and_distribution_plan.md`
- `docs/bodyweight_distribution_implementation_plan.md`
- `docs/Bangladesh_Arsenic_Risk_Model_Parameters_Full_Design.md`

The raw inputs remain read-only:

- `data/raw/NationalSurveyData.csv`
- `data/raw/bgd2018.csv`
- the Bangladesh MICS7 `ch.sav` child anthropometry file under `data/raw/`

## 3. Project objective

Build a reproducible district-wise Monte Carlo model for arsenic-related health risk in Bangladesh that:

1. preserves left-censored arsenic observations correctly;
2. fits district arsenic distributions using the required least-squares approach;
3. fits Bangladesh-specific adult and child body-weight distributions;
4. propagates fitted and specified exposure inputs through the full ADD, HQ, HI, and ELCR equations;
5. reports convergence, uncertainty, threshold exceedance, and sensitivity results;
6. keeps every transformation traceable back to an unchanged raw value.

## 4. Fixed methodological principles

The following principles should remain fixed unless the project specification changes:

- Raw data files are never overwritten.
- Original arsenic strings, including values such as `<6`, are preserved.
- Censored arsenic values are not treated as exact measured concentrations.
- Arsenic fitting is performed in mg/L, after an explicit conversion from micrograms per litre.
- Survey weights are used for STEPS and MICS body-weight estimation.
- Special survey codes are not interpreted as physical body weights.
- Legitimate extreme observations are not removed by automatic percentile or IQR trimming.
- Adult and child body-weight distributions are selected independently.
- District arsenic distributions are selected independently.
- Monte Carlo inputs $C$, $IR$, $BW$, $EF$, $ET$, and $SA$ are sampled independently because the source datasets do not provide their joint dependence structure. This approved assumption must be reported as a limitation.
- The primary numerical fit is least-squares based, with a scientifically appropriate robustness fit.
- Monte Carlo sampling uses explicit inverse-CDF sampling.
- The complete risk equations are retained without undocumented algebraic or unit changes.
- Random seeds, software versions, configuration values, and row counts are recorded.

## 5. Approved decisions and remaining hold points

### 5.1 District set

The following maximum set of ten districts is approved. Use these exact `DISTRICT` strings from the BGS file:

1. Dhaka
2. Chittagong
3. Rajshahi
4. Khulna
5. Barisal
6. Sylhet
7. Rangpur
8. Mymensingh
9. Comilla
10. Bogra

This decision is frozen. The decision record must state that the selection is purposive, balancing geographic coverage, importance, sample count, and censoring burden. Rajshahi must carry an explicit high-censoring warning.

### 5.2 Population definitions

- Adult body weight: Bangladesh STEPS 2018 participants aged 18-69 years.
- Child body weight: Bangladesh MICS7 children aged 0-59 months.
- The MICS distribution must be described as under-5 and must not be presented as direct measurement for ages 60-71 months.
- Retain the planned six-year child exposure duration. The mismatch between that duration and the available 0-59-month measured body-weight data remains an explicit model factor and limitation in the methods, sensitivity discussion, and final interpretation.

### 5.3 Censoring policy

Freeze the use of a left-censoring-aware empirical distribution for arsenic, followed by:

- primary CDF least-squares fitting;
- censored maximum-likelihood fitting as a robustness comparison.

No substitution rule such as LOD, LOD/2, or zero should be introduced into the primary analysis.

### 5.4 Candidate distributions

Arsenic candidates:

- Normal
- Lognormal
- Gamma

Adult and child body-weight candidates:

- Normal
- Lognormal
- Gamma
- Triangular

### 5.5 Approved model-input interpretations

- Treat the reported adult and child ingestion-rate Lognormal values `1.26` and `0.66` directly as the Lognormal log-space parameters $\mu_{\log}$ and $\sigma_{\log}$, respectively. Do not convert them from arithmetic mean and standard deviation.
- Treat the reported adult skin-area Lognormal values `1.42` and `0.31` directly as $\mu_{\log}$ and $\sigma_{\log}$, respectively. Do not convert them from arithmetic mean and standard deviation.
- Use a Triangular distribution for child exposed skin surface area. The family is frozen, but its required minimum-mode-maximum triplet must be documented before the final probabilistic child run.
- Maintain independent Monte Carlo sampling of $C$, $IR$, $BW$, $EF$, $ET$, and $SA$. Do not introduce correlations without a supported joint-data source or a later approved dependence model.
- Retain the convergence rule from `docs/Bangladesh_Arsenic_Risk_Model_Parameters_Full_Design.md`: $N\in\{1{,}000,5{,}000,10{,}000,20{,}000\}$, a primary $N=10{,}000$ result, successive relative changes, independent replicate seeds, and Monte Carlo standard error for threshold probabilities.

### 5.6 Remaining model-input hold points

Before simulation, create a signed-off parameter register containing every deterministic and stochastic model input, its unit, source, population, distribution interpretation, and status.

Resolve the remaining ambiguities identified in the full design document, especially:

- the minimum, mode, and maximum for the approved child skin-area Triangular distribution;
- exposure-frequency and exposure-time parameter definitions;
- reference dose and slope-factor units;
- adult and child averaging times;
- how the retained 0-59-month body-weight/six-year exposure mismatch is represented in sensitivity analysis and final limitations.

## 6. Phase 1: Reproducibility and data inventory

### Work

- Freeze the raw-file inventory and record file hashes.
- Record the exact MICS `ch.sav` path and survey release name.
- Define one project-root-aware path convention.
- Define output folders without creating alternate copies of raw data.
- Record Python and dependency versions planned for the final run.
- Create a machine-readable run configuration design covering districts, candidate distributions, seeds, simulation sizes, and policy choices.

### Deliverables

- data inventory and provenance table;
- parameter decision register;
- run-configuration specification;
- output naming convention;
- reproducibility checklist.

### Completion gate

Every source file, model parameter, population definition, and unresolved decision has an owner and status. No analytical transformation begins while a unit-critical model input remains ambiguous.

## 7. Phase 2: Arsenic preprocessing design

### Work

- Read the BGS file using its verified preamble/header structure.
- Preserve `As` as `As_raw`.
- classify each arsenic observation as detected, left-censored, missing, or malformed;
- extract censoring limits without treating them as detections;
- convert detected values and censoring limits from micrograms per litre to mg/L;
- retain identifiers, geography, coordinates, well metadata, and sample date;
- inspect duplicate sample IDs, duplicate full rows, and repeated coordinates;
- filter only after the selected district list is frozen;
- summarize detected, censored, missing, and malformed counts by district.

### Planned outputs

- `results/arsenic/arsenic_selected_districts_preprocessed.csv`
- `results/arsenic/arsenic_district_preprocessing_summary.csv`
- `results/arsenic/arsenic_duplicate_review.csv`
- `results/arsenic/arsenic_malformed_value_review.csv`

### Validation gate

- Raw row counts reconcile with classified row counts.
- Every selected row retains its original arsenic string.
- Unit conversion is verified with fixed examples.
- No censored value appears in the detected-value field.
- Duplicate-coordinate records are reviewed rather than automatically deleted.
- District names exactly match the raw source.

## 8. Phase 3: Body-weight preprocessing design

### 8.1 Adult STEPS dataset

Retain measured weight `m12`, survey weight `wstep2`, age, sex, measured height, participant ID, PSU, stratum, division, and urban/rural status.

Apply the following planned rules in a countable sequence:

- restrict age to 18-69;
- remove missing measured weight;
- remove known non-measurement codes `666` and `888`;
- require positive physical weight;
- require a positive survey weight;
- investigate duplicate participant IDs;
- retain height-derived BMI only as a diagnostic;
- preserve legitimate weight extremes unless source-specific evidence shows a coding error.

### 8.2 Child MICS dataset

Use `ch.sav`, not the birth-history file. Retain `AN8`, `chweight`, `CAGE`, `HL4`, `HH7A`, `AN11`, `WAZFLAG`, `HAZFLAG`, `WHZFLAG`, PSU, and stratum.

Apply the following planned rules in a countable sequence:

- restrict age to 0-59 months;
- remove missing measured weight;
- remove `AN8` codes `99.3`, `99.4`, `99.5`, and `99.6` separately;
- flag unexpected `AN8 >= 90` values for codebook review;
- require positive physical weight and positive `chweight`;
- retain `WAZFLAG == 0` for the primary fit;
- retain height-related flags for diagnostics without using them as automatic weight exclusions;
- preserve valid tails after survey-defined flagging.

### Planned outputs

- `results/bodyweight/adult_bw_clean.csv`
- `results/bodyweight/adult_bw_cleaning_audit.csv`
- `results/bodyweight/child_bw_clean.csv`
- `results/bodyweight/child_bw_cleaning_audit.csv`
- `results/bodyweight/bodyweight_weighted_summary.csv`

### Validation gate

- Starting, excluded, and final row counts reconcile exactly.
- Adult `BW_kg` comes from `m12`, never `wstep2`.
- Child `BW_kg` comes from `AN8`, never `chweight`.
- No documented special code survives as body weight.
- Survey weights are present, finite, and positive in fitting datasets.
- Clean outputs retain design variables needed for later uncertainty work.

## 9. Phase 4: Arsenic distribution fitting

### Work

For each selected district:

- construct a left-censoring-aware empirical CDF;
- fit Normal, Lognormal, and Gamma distributions using the declared CDF least-squares objective;
- fit the same families using censored likelihood as a robustness benchmark;
- report parameter estimates in an unambiguous parameterization;
- calculate comparable fit-error measures;
- examine ECDF/CDF agreement, quantiles, and lower and upper tails;
- assess bootstrap stability when the sample size and censoring pattern permit it;
- document support violations, optimizer failures, and high-censoring uncertainty;
- select one provisional family per district using the predeclared rule;
- require human review of every provisional selection before freezing it.

### Planned outputs

- `results/arsenic/arsenic_distribution_fit_results.csv`
- `results/arsenic/arsenic_bootstrap_stability.csv`
- `results/arsenic/arsenic_selected_distributions.csv`
- district diagnostic figures for empirical/fitted CDFs, densities, quantiles, and tails;
- a short district selection rationale table.

### Validation gate

- Censored observations contribute through a censoring-aware method.
- Parameters reproduce the saved fitted CDFs.
- Lognormal and Gamma fits use zero location unless a separately approved model says otherwise.
- Selection is not based on a single statistic alone.
- High-censoring districts are clearly labelled.
- Sampling from each selected fit produces finite, nonnegative concentrations in mg/L.

## 10. Phase 5: Body-weight distribution fitting

### Work

For adults and children separately:

- normalize survey weights for weighted empirical-CDF calculations;
- aggregate survey mass at tied body weights;
- construct weighted midpoint plotting positions;
- fit Normal, Lognormal, Gamma, and Triangular candidates by weighted CDF least squares;
- fit all candidates with survey-weighted pseudo-likelihood as a robustness benchmark;
- report weighted CDF SSE and maximum weighted CDF discrepancy;
- compare parameters and tail behaviour between fitting routes;
- report Normal negative-weight probability;
- verify positive support for Lognormal and Gamma;
- verify that Triangular support covers the observations carrying survey mass;
- review weighted histograms, fitted CDFs, Q-Q plots, and tails;
- select adult and child families independently;
- require human review before declaring either selection final.

### Planned outputs

- `results/bodyweight/bodyweight_distribution_fit_results.csv`
- `results/bodyweight/bodyweight_selected_distributions.csv`
- weighted histogram and density figures;
- weighted empirical and fitted CDF figures;
- candidate-specific weighted Q-Q and tail figures.

### Validation gate

- Multiplying all survey weights by a constant does not change parameter estimates beyond numerical tolerance.
- Every optimizer reports convergence or an explicit failure state.
- Saved parameters reproduce the reported SSE and discrepancy.
- Selected distributions pass physical-support checks.
- Inverse-CDF samples are finite and positive.
- Adult and child selections are justified independently.

## 11. Phase 6: Risk-equation module design

### Work

- Translate the approved ingestion and dermal ADD equations directly from the design document.
- Implement HQ for each pathway, total HI, and ELCR without dropping factors.
- attach a unit annotation to every equation input and intermediate output;
- define adult and child parameter sets explicitly;
- separate deterministic constants from sampled variables;
- define assertions for dimensional consistency and finite outputs;
- prepare hand-calculated benchmark cases for each equation.

### Planned outputs

- equation specification table;
- parameter register with units and sources;
- deterministic benchmark cases with expected results;
- test specification for ADD, HQ, HI, and ELCR.

### Validation gate

- Independent hand calculations match planned software results within declared tolerance.
- All concentration inputs are mg/L and body weights are kg.
- Pathway-specific factors appear exactly once.
- Adult and child averaging-time assumptions are documented.
- Risk outputs have the intended interpretation and units.

## 12. Phase 7: Monte Carlo simulation design

### Work

- define a master seed and deterministic seed derivation by district, population, and run size;
- generate every stochastic input through Uniform sampling followed by the selected inverse CDF;
- run each selected district for adults and children;
- preserve iteration-level inputs needed for sensitivity analysis;
- compute ADD ingestion, ADD dermal, HQ ingestion, HQ dermal, HI, and ELCR per iteration;
- define a principal simulation size consistent with the project specification;
- define smaller and larger run sizes for convergence analysis;
- decide whether full iterations are stored as Parquet or summarized when storage is constrained.

### Planned outputs

- district/population iteration outputs or approved compact equivalents;
- run manifest containing seeds, sample sizes, distributions, parameters, and software versions;
- simulation summary table;
- deterministic-versus-probabilistic comparison table.

### Validation gate

- Repeating a run with the same seed reproduces results exactly.
- Different district/population streams do not accidentally reuse identical random sequences.
- Sample summaries agree with their target distributions.
- No NaN, infinite, or physically impossible sampled value reaches the risk equations.
- Threshold probabilities are calculated directly from iteration-level results.

## 13. Phase 8: Convergence and numerical error analysis

### Work

- use the planned sequence $N\in\{1{,}000,5{,}000,10{,}000,20{,}000\}$ for every district-population combination;
- track stability of mean, median, upper percentiles, and threshold probabilities;
- repeat selected sizes with independent seeds;
- quantify successive relative change between each simulation size and the immediately preceding size;
- retain $N=10{,}000$ as the primary base-paper-comparable result and $N=20{,}000$ as an additional stability check, not as assumed truth;
- calculate Monte Carlo standard error for estimated threshold probabilities;
- distinguish Monte Carlo error from fitted-parameter uncertainty.

### Planned outputs

- `results/convergence/convergence_summary.csv`
- convergence plots by district, population, metric, and simulation size;
- repeated-seed stability table;
- written conclusion identifying the smallest adequate simulation size.

### Validation gate

- Stabilization is assessed using the successive-relative-change procedure defined in the parameter design.
- Upper-tail results are assessed separately from central estimates.
- Convergence claims are supported by more than one random seed.
- Any metric that fails to stabilize is explicitly reported.

## 14. Phase 9: Sensitivity and uncertainty analysis

### Work

- compute tie-aware Spearman rank correlations between sampled inputs and risk outputs;
- produce tornado plots separately for relevant adult and child outcomes;
- preserve signs as well as magnitudes;
- distinguish input variability from parameter-estimation uncertainty;
- add cluster/stratum-aware resampling for fitted body-weight parameters if required;
- use an appropriate resampling strategy for district arsenic fits;
- avoid claiming independent-row bootstrap as survey-design-correct.

### Planned outputs

- sensitivity coefficient tables;
- tornado plots;
- fitted-parameter uncertainty summaries;
- comparison of variability-only and parameter-uncertainty scenarios, if included.

### Validation gate

- Sensitivity calculations use the exact iterations that produced the reported risks.
- Fixed inputs are not presented as sampled sensitivity variables.
- Correlation signs and labels are verified against the equations.
- Survey-design limitations are stated accurately.

## 15. Phase 10: Final reporting

### Required result classes

- data provenance and preprocessing flow;
- row-count and exclusion audit tables;
- district arsenic descriptive statistics;
- fitted-distribution comparison and selection tables;
- adult and child body-weight fit results;
- deterministic risk benchmark;
- district-wise adult and child risk summaries;
- PDF/CDF-style risk figures;
- probabilities of exceeding health-risk thresholds;
- convergence and repeated-seed analysis;
- sensitivity results;
- limitations and uncertainty discussion.

### Reporting rules

- State whether each parameter is measured, fitted, assumed, or fixed.
- State distribution parameterizations explicitly.
- Keep micrograms per litre and mg/L visibly distinct.
- Describe MICS body weight as under-5 only.
- Describe the BGS data as sampled wells, not a perfectly population-weighted census of wells.
- Do not present weighted pseudo-likelihood information criteria as ordinary iid AIC without a justified survey method.
- Do not hide optimizer failures or unsupported tail extrapolation.
- Keep automated provisional distribution selection separate from final scientific approval.

### Completion gate

Every reported number is traceable to a run manifest, source dataset, preprocessing audit, fitted parameter record, and reproducible seed.

## 16. Testing strategy

Testing should be planned at four levels.

### Unit tests

- arsenic-string classification and limit extraction;
- unit conversions;
- survey-weight normalization;
- weighted ECDF plotting positions;
- each distribution's parameter transformation, CDF, PDF, and PPF;
- each risk equation;
- summary and threshold calculations.

### Data-contract tests

- required columns exist;
- raw files are unchanged;
- row-count reconciliation holds;
- special codes are absent from physical-value fields;
- cleaned weights and concentrations satisfy support rules.

### Numerical tests

- fitted objectives recompute from saved parameters;
- weight-scale invariance holds;
- inverse-CDF samples reproduce target quantiles;
- fixed-seed runs are reproducible;
- optimizer failures are surfaced rather than silently accepted.

### End-to-end tests

- one small district/population run completes from frozen configuration to report table;
- deterministic benchmark values match;
- all expected artifacts are produced once in the documented locations;
- report values reconcile with machine-readable outputs.

## 17. Recommended execution order

1. Approve the district set, population definitions, censoring policy, and unresolved parameter decisions.
2. Freeze data provenance and the run-configuration specification.
3. Complete arsenic preprocessing and its reconciliation review.
4. Complete adult and child body-weight preprocessing and review exclusion counts.
5. Fit and approve district arsenic distributions.
6. Fit and approve adult and child body-weight distributions.
7. Validate the complete risk equations against hand calculations.
8. Run a small end-to-end pilot for one district and both populations.
9. Run the convergence experiment and freeze the principal simulation size.
10. Run all district/population simulations with the frozen configuration.
11. Produce sensitivity and optional parameter-uncertainty analyses.
12. Generate final tables and figures from machine-readable results.
13. Perform an independent reconciliation and methods-language review.

## 18. Hold points requiring approval

No later phase should silently decide any of the following:

- treatment of highly censored districts;
- the minimum-mode-maximum triplet for the approved child skin-area Triangular distribution;
- adult and child averaging-time treatment;
- sensitivity treatment of the retained MICS age-coverage mismatch;
- final arsenic distribution for each district;
- final adult and child body-weight distributions;
- inclusion or exclusion of fitted-parameter uncertainty.

Each hold point should be recorded with the selected option, rationale, reviewer, and date.

## 19. Definition of done

The project is complete only when:

- raw data integrity is demonstrable;
- all preprocessing counts reconcile;
- censoring and survey weights are handled as declared;
- selected distributions have documented numerical and graphical support;
- all equations pass benchmark tests;
- simulations are reproducible from a frozen configuration;
- convergence is demonstrated for the reported metrics;
- sensitivity outputs use the actual simulation inputs and outputs;
- every final table and figure is generated from saved machine-readable results;
- limitations are stated without overstating geographic or age representativeness.
