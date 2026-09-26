# Phase 6: Risk-Equation Module Handoff

Owner: Member 4. Cross-reviewer: Member 3. Member 5 does QA. Phase 7 (Monte Carlo simulation, Member 4) uses these outputs.

## How to reproduce

```bash
python scripts/run_risk_phase6.py   # writes results/tables/risk_*; exits 1 if paper Table 4 is not reproduced
python -m pytest -q tests/test_risk_equations.py   # 38 equation, unit, parameter and benchmark tests
```

No random numbers are drawn in Phase 6.

## Code

| Module | Contents |
|---|---|
| `src/arsenic_hra/risk_equations.py` | Equations (1)–(6) of the base paper in full, unsimplified form: `add_ingestion`, `add_dermal`, `hazard_quotient`, `hazard_index`, `elcr`. `evaluate_risk(inputs, params)` applies them element-wise to arrays and rejects non-finite or physically impossible inputs and outputs. |
| `src/arsenic_hra/risk_parameters.py` | Adult and child parameter sets built from `risk_model` in `config/run_config.json`. Fixed constants (`FIXED_INPUTS`) are kept separate from sampled inputs (`STOCHASTIC_INPUTS`). Every name carries its unit (for example `IR_L_per_day`), and `UNITS` holds the full unit strings. Frozen SciPy distributions for IR, EF, ET and SA are ready for Phase 7 inverse-CDF sampling (`InputDistribution.ppf`). |

## Outputs (`results/tables/`)

| File | Contents |
|---|---|
| `risk_equation_specification.csv` | Each equation with its paper number, inputs and units, and output unit |
| `risk_parameter_table.csv` | Every input for both populations: role (fixed or sampled), unit, family, explicit parameters, deterministic value, and the implied mean, median, P5 and P95 |
| `risk_paper_table4_reproduction.csv` | The base paper's Table 4, recomputed from its own Tables 1 and 2, with the deviation for each cell |
| `risk_deterministic_benchmark.csv` | Bangladesh point-input benchmark at C = 0.01 and 0.05 mg/L |
| `risk_phase6_manifest.json` | Config SHA-256, the BW summary SHA-256, software versions, the settings used, and reproduction counts |

## Equations

Units: C in mg/L, IR in L/day, EF in days/year, ED in years, BW in kg, AT in days, Kp in cm/h, ET in h/day, SA in m², CF in L/(cm·m²).

| Output | Equation | Unit |
|---|---|---|
| ADD_ing (1) | C·IR·EF·ED / (BW·AT_nc) | mg kg⁻¹ day⁻¹ |
| ADD_dermal (2) | C·Kp·EF·ED·ET·SA·CF / (BW·AT_nc) | mg kg⁻¹ day⁻¹ |
| HQ_ing (3) | ADD_ing / RfD_ing, with RfD_ing = 0.0003 | – |
| HQ_dermal (4) | ADD_dermal / RfD_dermal, with RfD_dermal = 0.000285 | – |
| HI (5) | HQ_ing + HQ_dermal | – |
| ELCR (6) | ADD_ing,cancer · CSF, with CSF = 1.5 (mg kg⁻¹ day⁻¹)⁻¹ | – |

AT_nc = ED × 365. ADD_ing,cancer is equation (1) evaluated with AT_cancer (see decision D2).

**CF check (resolved).** 1 cm × 1 m² = 0.01 m³ = 10 L, so CF = 10 L/(cm·m²). That is the same unit the paper prints as L·m/(m³·cm). The product C·Kp·ET·SA·CF is therefore in mg/day. `test_dermal_dose_matches_si_first_principles` recomputes the dermal dose entirely in SI units and matches to 1e-12.

## Decision log (2026-09-25, Member 4; team sign-off pending)

The base paper was obtained as `docs/references/Yadav_Kalkal_2024_JWH_arsenic_Punjab.pdf` (CC BY 4.0). All values below are checked against its Table 2.

| # | Decision | Options considered | Chosen, and why |
|---|---|---|---|
| D1 | How to read "mean ± SD" for the Lognormal IR (1.26 ± 0.66 L/day) and adult SA (1.42 ± 0.31 m²) | (a) arithmetic mean and SD, converted to log-space; (b) log-space μ and σ used directly, as previously approved | **(a).** Reading (b) gives an adult exposed-skin median of 4.14 m² and a P95 of 6.9 m². That is above adult *total* body surface area; the US EPA Exposure Factors Handbook 2011, Table 7-1, gives about 2.5 m² as the male P95. It also gives children an IR P95 of 10.4 L/day. The paper's Table 2 writes BW as "77.1 + 31.5 kg" in the same notation, which can only be arithmetic. Simulating the paper's Moga adult under (b) gives a P95 HI of at least 2.05 even with C held fixed, but the paper reports 0.78 (Table 5), so (b) cannot reproduce the paper's own results. Reading (a) gives IR a median of 1.12 and a P95 of 2.51 L/day, and SA a median of 1.39 and a P95 of 1.98 m². **This reverses the approved master-plan 5.5 reading, so the design docs need an update after team sign-off.** Switch: `risk_model.plus_minus_interpretation`. |
| D2 | Averaging time for ELCR | (a) 70-year lifetime, 70 × 365 = 25,550 days; (b) ED × 365, as the paper does (Table 4 reproduces only with 2,190 days for children); (c) both | **(a) as the primary result** (US EPA RAGS Part A: ELCR is a *lifetime* risk). Adults are unaffected, since ED = 70. Child ELCR under (a) is 6/70 of the paper-convention value. Switch: `risk_model.cancer_averaging_time`. |
| D3 | Child exposed-skin Triangular triplet | (a) EPA-sourced (0.29, 0.60, 0.95) m²; (b) fixed at 0.6 m²; (c) a team-supplied triplet | **(a).** The paper's "Triangular 0.6800 + 0.600" cannot define a triplet. Min 0.29 = EFH Table 7-1 mean for birth to <1 month. Mode 0.60 = the paper's deterministic child SA, which is also close to the EFH mean for 2 to <3 years (0.61). Max 0.95 = EFH P95 for 3 to <6 years. Together they cover ages 0–6, matching ED = 6. US children are larger than Bangladeshi children, so the maximum is conservative. |
| D4 | EF and ET definitions | – | Verified against the paper: EF "345 (180–365)" is Triangular with min 180, mode 345, max 365; ET "0.20 (0.13–0.33)" is Triangular with min 0.13, mode 0.20, max 0.33. The paper's deterministic ET of 0.58 h/day lies *outside* its own probabilistic range. Both values are kept as published and noted here. |
| D5 | Deterministic benchmark point values | – | The deterministic values in the paper's Table 2 are used for the base-paper inputs: IR 2.0 / 1.8 L/day, EF 365, ET 0.58 h/day, SA 1.8 / 0.6 m². BW is the **Bangladesh** survey-weighted mean (adult 55.89 kg, child 10.87 kg), as design doc section 15 requires. |

## Validation against the base paper (Table 4)

The paper's inputs are C from its Table 1, BW 70 / 15 kg, and AT = ED × 365. With them, **25 of the 28 published cells reproduce within 0.5%**, which is the rounding of values printed to 2–3 significant figures. Three ELCR cells disagree with the paper's own inputs. They look like typos in the paper, and they are checked against the hand arithmetic instead:

| Cell | Recomputed from the paper's inputs | Printed |
|---|---:|---:|
| Ferozepur child ELCR (0.027 × 1.8 / 15 × 1.5) | 4.86×10⁻³ | 4.80×10⁻³ |
| Amritsar adult ELCR (0.083 × 2 / 70 × 1.5) | 3.56×10⁻³ | 3.50×10⁻³ |
| Amritsar child ELCR (0.083 × 1.8 / 15 × 1.5) | 1.49×10⁻² | 1.40×10⁻² |

## Deterministic Bangladesh benchmark (reference concentrations)

| C (mg/L) | Population | HQ_ing | HQ_dermal | HI | ELCR (lifetime AT) |
|---|---|---:|---:|---:|---:|
| 0.01 (WHO guideline) | adult | 1.193 | 0.0066 | 1.199 | 5.37×10⁻⁴ |
| 0.01 | child | 5.521 | 0.0112 | 5.532 | 2.13×10⁻⁴ |
| 0.05 (Bangladesh standard) | adult | 5.965 | 0.0328 | 5.997 | 2.68×10⁻³ |
| 0.05 | child | 27.60 | 0.0562 | 27.66 | 1.06×10⁻³ |

Ingestion dominates: the dermal pathway contributes 0.2–0.6% of HI. Under the paper's AT convention, child ELCR would be 70/6 ≈ 11.7 times higher: 2.48×10⁻³ at 0.01 mg/L, and 1.24×10⁻² at 0.05 mg/L.

## Test specification (`tests/test_risk_equations.py`, 38 tests)

- **Hand-calculated benchmarks.** Round-number adult and child cases; a non-cancelling EF case (EF = 300, where EF·ED ≠ AT); and child ELCR under both AT conventions. The arithmetic is written in comments, and the tolerance is 1e-9 or tighter.
- **Paper reproduction.** Table 4 for all 7 districts × 2 populations × HI and ELCR.
- **Dimensional consistency.** Ingestion and dermal doses recomputed from SI first principles.
- **Each factor appears exactly once.** Doubling C, IR, BW, EF, ET, SA, Kp, CF or AT scales each dose by exactly 2^(±1), or leaves it unchanged where that input does not belong to the pathway.
- **Identities.** HI = HQ_ing + HQ_dermal and ELCR = ADD_ing,cancer × CSF, element-wise on arrays of 1,000 values. Scalar and vector results agree.
- **Physical-support guards.** Negative C; BW, IR or SA ≤ 0; NaN or inf; EF > 366; ET > 24; and missing inputs are all rejected.
- **Parameter sets.** Constants, ED and AT values; the fixed and sampled inputs are disjoint; the Lognormal moments reproduce the reported mean and SD exactly; the Triangular SciPy arguments are correct; inverse-CDF draws are finite and positive; adult SA P95 < 2.5 m²; the log-space switch restores the old parameterization.

## Validation gate (master plan section 11)

| Gate | Status |
|---|---|
| Independent hand calculations match within declared tolerance | pass (1e-9 hand cases; 0.5% against the published Table 4) |
| All concentration inputs are mg/L and body weights are kg | pass (unit-suffixed names; the Phase 4 and 5 samplers return mg/L and kg) |
| Pathway-specific factors appear exactly once | pass (`test_each_pathway_factor_appears_exactly_once`, `test_fixed_constants_enter_once`) |
| Adult and child averaging-time assumptions documented | pass (D2) |
| Risk outputs have the intended interpretation and units | pass (`risk_equation_specification.csv`) |

## Open items for the team

1. **Sign off D1 to D3.** D1 reverses an approved decision. After sign-off, update master plan section 5.5 and design doc sections 2, 5 and 10.3/10.7, all of which still describe the log-space reading.
2. **Barisal arsenic fit (Member 2, Phase 4 hold point).** The selected Lognormal has σ_log = 3.53, which gives a mean of 2.7 mg/L, a P99 of 19.6 mg/L and a P99.9 of 291 mg/L. The largest measured value in Barisal is 0.862 mg/L, and 7.5% of draws exceed it. This will dominate Barisal's upper-tail HI in Phase 7. Please review it before the fits are frozen. Bogra (σ_log = 2.86, P99.9 = 6.2 mg/L) is milder but similar.
3. **Child BW family (Member 3, Phase 5 hold point).** The choice is Gamma versus a truncated Normal. It must be frozen before Phase 7, because the lightest children carry the highest dose per kg.
4. **Phase 4 handoff report.** There is no `reports/` document for Phase 4 yet.
5. **`config_version`.** This phase adds a `risk_model` section to `config/run_config.json` but leaves `config_version` at 1.0.0. Bumping it would make the Phase 2 and 4 manifest tests fail until those pipelines are re-run. The Phase 6 manifest records the config SHA-256 instead. Member 5 should bump the version and re-run all pipelines once at integration.
6. **Documentation typos (base paper).** Table 2 lists child BW as "Triangular 26.1 (6.5–15)": the mode lies outside its own range, and the value is not used here. Its deterministic ET (0.58) also lies outside its probabilistic range (0.13–0.33).

Reviewer sign-off (Member 3): ________  Date: ________
