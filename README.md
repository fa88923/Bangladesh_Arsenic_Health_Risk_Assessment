# Bangladesh Arsenic Health-Risk Assessment

This repository contains a reproducible probabilistic health-risk assessment for arsenic exposure in selected districts of Bangladesh. It adapts the Yadav and Kalkal (2024) Punjab arsenic-risk model to Bangladesh arsenic and body-weight data, preserves the original risk equations, and quantifies exposure risk with deterministic benchmarks, Monte Carlo simulation, convergence checks, sensitivity analysis, and fitted-parameter uncertainty.

The primary outputs are district-level non-carcinogenic hazard index (HI) and excess lifetime cancer risk (ELCR) estimates for adults and children.

The comprehensive report of this project is stored in the B_03.pdf. 


## Project Scope

The analysis covers 10 Bangladesh districts:

- Dhaka
- Chittagong
- Rajshahi
- Khulna
- Barisal
- Sylhet
- Rangpur
- Mymensingh
- Comilla
- Bogra

The model combines:

- District arsenic concentration distributions fitted from detected and left-censored well measurements.
- Adult body-weight distributions from Bangladesh STEPS 2018.
- Child body-weight distributions from Bangladesh MICS7 2025.
- Published exposure-factor distributions and risk equations from the base paper.
- Monte Carlo simulation with reproducible random streams.

Raw data are treated as immutable. Provenance checks verify that the raw inputs have not changed before preprocessing and fitting workflows run.

## Repository Structure

| Path | Purpose |
|---|---|
| `config/` | Run configuration and parameter register |
| `data/raw/` | Immutable raw survey and arsenic source files |
| `docs/` | Methodology plans, design notes, and reference paper |
| `reports/` | Phase handoff reports and interpretation notes |
| `results/` | Generated clean data, fitted distributions, tables, figures, manifests, and simulation iterations |
| `scripts/` | Reproducible command-line workflows |
| `src/arsenic_hra/` | Analysis package: preprocessing, fitting, risk equations, simulation, convergence, and sensitivity modules |
| `tests/` | Unit, data-contract, reproduction, and integration tests |

## Data Sources

The project uses three main data inputs:

1. Bangladesh arsenic well measurements, including detected and left-censored concentration values.
2. Bangladesh STEPS 2018 adult body-weight data.
3. Bangladesh MICS7 2025 child body-weight data.

Arsenic concentrations are modeled in mg/L. Body weights are modeled in kg. Censored arsenic values are not substituted with zero, LOD, or LOD/2; they are represented as censoring bounds and handled with reverse Kaplan-Meier empirical CDFs and censored fitting methods.

## Methods Summary

### Arsenic Preprocessing and Fitting

The arsenic preprocessing workflow classifies all source arsenic strings as detected, left-censored, missing, or malformed; converts units from ug/L to mg/L; preserves the original raw value; and writes a selected-district fitting file.

For each district, a censoring-aware empirical CDF is estimated using reverse Kaplan-Meier. Normal, Lognormal, and Gamma candidate distributions are fit by CDF least squares, with censored-MLE and bootstrap checks used as supporting evidence.

Final arsenic simulation models:

| Family | Districts |
|---|---|
| Gamma | Dhaka, Rajshahi, Khulna, Barisal, Sylhet, Rangpur, Mymensingh, Comilla |
| Lognormal | Chittagong, Bogra |

Barisal was changed from its provisional Lognormal fit to Gamma after tail checks showed unsupported extrapolation. Rajshahi is retained with a high-censoring warning.

### Body-Weight Preprocessing and Fitting

Adult body weight is cleaned from STEPS 2018 and child body weight from MICS7 2025. Survey weights are used in weighted empirical CDFs and fitting.

Final body-weight simulation models:

| Population | Model |
|---|---|
| Adult | Lognormal |
| Child | Truncated Normal, bounded below at 1.6 kg |

The child truncated Normal replaced the provisional Gamma model because it better reproduced the lower tail, which is dose-relevant because exposure dose scales inversely with body weight.

### Risk Equations

The implementation preserves the base paper's equations for:

- Average daily dose by ingestion.
- Average daily dose by dermal absorption.
- Ingestion hazard quotient.
- Dermal hazard quotient.
- Total hazard index.
- Excess lifetime cancer risk.

The project reports lifetime-averaged ELCR as the primary cancer-risk quantity. The paper's exposure-duration averaging convention is also stored for comparison.

### Simulation, Convergence, and Sensitivity

The primary Monte Carlo run uses 10,000 iterations for every district and population. Each run samples arsenic concentration, ingestion rate, body weight, exposure frequency, exposure time, and exposed skin area using explicit inverse-CDF sampling and fixed seed streams.

Convergence is assessed with 400 additional runs over 20 district-population cells, 4 sample sizes, and 5 independent seeds per sample size.

Sensitivity analysis includes:

- Spearman rank correlations between sampled inputs and risk outputs.
- Fitted-parameter uncertainty using bootstrap refits.
- Scenario tests for key modeling choices.

## Main Findings

- Comilla has the highest estimated risk, followed by Barisal, Dhaka, and Khulna.
- Children have higher HI than adults in every district because of lower body weight relative to intake.
- Arsenic concentration dominates variation in HI across all districts and populations.
- Dermal exposure contributes only about 0.1-0.3% of mean HI; ingestion dominates.
- Probabilistic P95 HI is higher than deterministic point estimates because arsenic concentrations are strongly right-skewed.
- The headline exceedance metric `P(HI > 1)` is more stable than mean HI or extreme-tail summaries.
- Fitted-parameter uncertainty is much larger than Monte Carlo numerical error; the limiting factor is the available arsenic data, not the simulation size.

## Reproducing the Analysis

Create a Python environment and install the pinned dependencies:

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

Run the workflows in order:

```bash
python scripts/build_provenance.py --verify
python scripts/preprocess_arsenic.py
python scripts/run_arsenic_fitting.py
python scripts/run_bodyweight.py
python scripts/run_risk_phase6.py
python scripts/run_simulation.py
python scripts/run_convergence.py
python scripts/run_sensitivity.py --reuse-bootstrap
```

For a full sensitivity rerun that regenerates bootstrap tables:

```bash
python scripts/run_sensitivity.py
```

That full run is slower than the reuse path because it regenerates the arsenic and survey-bootstrap parameter sets.

## Key Outputs

| Output | Description |
|---|---|
| `results/arsenic/arsenic_selected_districts_preprocessed.csv` | Arsenic fitting input for the 10 selected districts |
| `results/arsenic/arsenic_selected_distributions.csv` | Final and candidate arsenic distribution parameters |
| `results/bodyweight/bodyweight_selected_distributions.csv` | Final body-weight models |
| `results/tables/risk_equation_specification.csv` | Implemented risk equations and units |
| `results/tables/risk_parameter_table.csv` | Fixed and stochastic risk-model inputs |
| `results/tables/risk_deterministic_district.csv` | Deterministic district benchmark table |
| `results/tables/risk_probabilistic_p95.csv` | Probabilistic P95 HI and ELCR table |
| `results/tables/risk_exceedance.csv` | `P(HI > 1)`, `P(HI > 2)`, and `P(ELCR > 1e-4)` |
| `results/simulation/simulation_summary.csv` | Monte Carlo summary statistics |
| `results/simulation/iterations/` | Saved primary iteration-level Parquet outputs |
| `results/convergence/convergence_summary.csv` | Monte Carlo convergence diagnostics |
| `results/sensitivity/parameter_uncertainty_summary.csv` | Bootstrap uncertainty intervals |
| `results/sensitivity/scenario_results.csv` | Model-choice scenario results |
| `results/figures/` | Diagnostic, simulation, convergence, and sensitivity figures |

## Reports

The main report files are:

- `reports/arsenic_preprocessing_report.md`
- `reports/arsenic_fit_report.md`
- `reports/bodyweight_report.md`
- `reports/bodyweight_fit_report.md`
- `reports/risk_equations_report.md`
- `reports/simulation_report.md`
- `reports/convergence_report.md`
- `reports/sensitivity_report.md`
- `reports/limitations_register.md`

These reports record the methodological decisions, validation gates, outputs, and limitations for each phase.

## Testing

Run the complete test suite with:

```bash
python -m pytest -q
```

Focused test files are available for each major component:

```bash
python -m pytest -q tests/test_arsenic_preprocessing.py
python -m pytest -q tests/test_arsenic_fitting.py
python -m pytest -q tests/test_bodyweight_preprocess.py
python -m pytest -q tests/test_bodyweight_fitting.py
python -m pytest -q tests/test_risk_equations.py
python -m pytest -q tests/test_simulation.py
python -m pytest -q tests/test_convergence.py
python -m pytest -q tests/test_sensitivity.py
```

The tests cover raw-data provenance, row-count reconciliation, unit conversion, distribution parameter reproduction, risk-equation benchmarks, seed reproducibility, sampling validation, convergence summaries, and sensitivity outputs.

## Important Limitations

- The arsenic data are treated as sampled wells, not a population-weighted census of wells.
- Arsenic below detection limits is extrapolated by the fitted distributions; lower-percentile and median estimates should be interpreted with this caveat.
- Rajshahi has high censoring, with 71.8% non-detects.
- Bogra's Lognormal upper tail is unstable; mean HI is not numerically well converged in the tested grid.
- The child body-weight source covers ages 0-59 months, while the child exposure duration is 6 years. Sensitivity analysis shows this slightly overstates child risk.
- Inputs are sampled independently because no joint exposure dataset is available.
- Fitted-parameter uncertainty dominates Monte Carlo error, especially in districts with fewer arsenic observations or heavier tails.

## Citation and Reference

The implemented risk equations and Punjab comparison structure are based on:

Yadav and Kalkal (2024), arsenic health-risk assessment in Punjab, available in `docs/references/Yadav_Kalkal_2024_JWH_arsenic_Punjab.pdf`.

The Bangladesh analysis, fitted distributions, simulation results, uncertainty analysis, and interpretation are generated by this repository.
