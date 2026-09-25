# Phase 9: Sensitivity and Uncertainty Analysis Handoff

Owner: Member 4. Cross-reviewer: Member 3. Member 5 does QA and assembles Phase 10.

## How to reproduce

```bash
python scripts/run_simulation.py                       # Phase 7 (the primary iterations are required)
python scripts/run_sensitivity.py                      # full run, about 20 min (regenerates the Phase 4 bootstrap)
python scripts/run_sensitivity.py --reuse-bootstrap    # about 1 min, reusing the saved bootstrap tables
python -m pytest -q tests/test_sensitivity.py
```

## Three questions, three analyses

| Question | Method | Uncertainty type |
|---|---|---|
| Which inputs drive the variation in risk? | Spearman rank correlation of HI (and ELCR) with each sampled input, on the exact Phase 7 iterations | **input variability** (between people and wells) |
| How sure are we of the fitted C and BW distributions? | Two-dimensional MC: 1,000 parameter sets from bootstrap refits, each run through the primary uniforms | **fitted-parameter uncertainty** (finite sample of wells and children) |
| How much do modelling choices matter? | Six scenarios, each rerun with the primary uniforms | **model and selection uncertainty** |

The Monte Carlo (numerical) error was quantified in Phase 8 and is the smallest of the four.

## Outputs (`results/sensitivity/`)

| File | Contents |
|---|---|
| `sensitivity_spearman.csv` | For each district, population, output (HI, ELCR) and input: ρ_s, p-value, rank, expected sign, and sign check |
| `sensitivity_spearman_HI_table.csv` | **Same layout as design doc table 19.3**: ρ_s(C, IR, BW, EF, ET, SA) with HI |
| `arsenic_bootstrap_replicates.csv` | The Phase 4 bootstrap, 10 districts × 1,000 replicates, regenerated with Member 2's code and seeds. **The per-replicate table was never saved before**; it is now |
| `arsenic_bootstrap_reproduction_check.csv` | The regenerated table's summary diffed against Member 2's saved `arsenic_bootstrap_stability.csv` |
| `bw_survey_bootstrap_replicates.csv` | 1,000 Rao–Wu PSU bootstrap refits for each population |
| `parameter_uncertainty_replicates.csv` | 20 cells × about 1,000 outer replicates: parameters used and the resulting risk statistics |
| `parameter_uncertainty_summary.csv` | For each cell and statistic: primary value, and the median and 2.5–97.5% interval across outer replicates |
| `scenario_results.csv` | Six scenarios: risk statistics and their ratio to the primary result |
| `sensitivity_manifest.json` | Methods, seeds, checks, and hashes |

Figures: `results/figures/sensitivity_tornado_<district>.png` (10) and `sensitivity_parameter_uncertainty_{P95_HI,P_HI_gt_1,P95_ELCR}.png`.

## Decision log (2026-09-26, Member 4)

### D20. Sensitivity variables and outputs

Only the six **sampled** inputs are used: C, IR, BW, EF, ET and SA (design doc section 16.2). The fixed ED, Kp, CF, RfDs and CSF have zero variance, so they have no rank correlation and are excluded. That also avoids the base paper's inconsistency of calling ED "least sensitive" while listing it as a point input. HI is the primary output, matching the base paper; ELCR is included as a secondary output. The coefficients come from `scipy.stats.spearmanr`, which is tie-aware, applied to the **exact iterations that produced the reported risks** (the Phase 7 Parquet files).

### D21. Sign check: the sign must match the equations only when ρ is distinguishable from zero (p < 0.001)

The first version required every sign to match. It failed on ET and SA, whose correlations are ρ ≈ ±0.02 with random signs. The dermal pathway is 0.1–0.3% of HI, so the true positive effect of ET and SA is smaller than the sampling SE of ρ, about 1/√N ≈ 0.01. A noise-level ρ with the "wrong" sign is not a defect; a *significant* wrong sign would be. After the fix, every ρ that is distinguishable from zero has the sign the equations imply. C and IR are distinguishable from zero in every cell.

### D22. Fitted-parameter uncertainty is included, as a separate analysis

Master plan section 18 lists this as a hold point, and the user delegated the choice. It is **included**, but kept **separate** from the primary result. The primary tables stay the variability-only, base-paper-comparable estimate, and this analysis puts an uncertainty interval around them. The alternative, reporting a pooled two-dimensional distribution as the result, would have mixed variability with uncertainty and broken comparability with the paper.

### D23. Arsenic parameter sets: the Phase 4 bootstrap, regenerated and verified

Member 2 saved only summaries, not the parameters of each replicate. The bootstrap is deterministic (`SeedSequence(20260923, spawn_key=(district, replicate))`), so `arsenic_bootstrap.run_bootstrap` was re-run with the frozen seeds, replicate count and tie tolerance. **Verification:** every winner frequency matches Member 2's saved summary exactly, the pattern of missing values is identical, and every continuous summary matches within **1.1×10⁻⁶ relative**. The remaining difference is the numerical environment: Member 2's manifest records numpy 2.5.3 / scipy 1.18.1, while this environment uses the pinned numpy 1.26.4 / scipy 1.17.1, so the optimizers stop at slightly different points. The tolerance is therefore 10⁻⁵ for continuous values and exact for winner counts. For each replicate, the **final** family's CDF least-squares refit is used (Gamma for Barisal). Replicates that failed or were inadmissible are dropped, leaving 994 for Chittagong and 1,000 elsewhere.

### D24. BW parameter sets: Rao–Wu survey bootstrap, B = 1,000

Within each stratum (adults: 16 strata × about 31 PSUs; children: 131 strata × 6–47 PSUs), n_h − 1 PSUs are drawn with replacement. Their weights are multiplied by n_h/(n_h − 1) times the number of draws, and the final family is refitted by weighted CDF least squares (truncated Normal for children, with the bound kept at 1.6 kg). Resampling whole PSUs keeps the cluster correlation. The master plan warns that an independent-row bootstrap would not be survey-design-correct, which is why this design is used. All 2,000 refits converged. The stream is `SeedSequence(20260923, spawn_key=(9001, population, replicate))`; its three-element key cannot collide with the simulation keys (4 elements) or the arsenic keys (2 elements).

### D25. Common random numbers

Every outer replicate reuses the primary run's uniforms. The spread across replicates is therefore due only to the parameters, with no fresh Monte Carlo noise mixed in. Verified by test: using the point-estimate parameters reproduces the primary statistics to 10⁻¹². Outer replicate b pairs arsenic replicate b with BW replicate b. The two bootstraps come from independent data sources, so the pairing is arbitrary and harmless.

### D26. IR, EF, ET and SA are held at their specified distributions in the two-dimensional analysis

These are base-paper assumptions with no underlying data, so no sampling uncertainty can be estimated for their parameters. Their *choice* is tested in scenario S6 instead.

### D27. Scenarios

Every scenario uses the primary uniforms, so the ratio to the primary result reflects the modelling change alone.

| ID | What changes | Why |
|---|---|---|
| S1 | Barisal C: Gamma → Lognormal (the Phase 4 provisional) | Quantifies the D7 override |
| S2 | Bogra C: Lognormal → Gamma | Bogra's Lognormal tail runs high (L6) |
| S3 | Rangpur C: Gamma → Lognormal | Rangpur's Gamma extreme tail is thin (L7) |
| S4 | Child BW: truncated Normal → Gamma (the Phase 5 provisional) | Quantifies the D8 override |
| S5 | Child BW covers **0–71 months** instead of 0–59 | The retained ED = 6 years vs under-5 data mismatch (L2, master plan 5.2) |
| S6 | IR and SA "mean ± SD" read as **log-space** (the pre-2026-09-25 reading) | Quantifies the Phase 6 D1 reversal |

**S5 construction.** A 6-year exposure spans ages 0–72 months, and MICS covers 60 of those 72. The 60–71-month group is built from the MICS 48–59-month children, shifted by their observed one-year gain: a weighted mean of 14.50 kg against 12.86 kg at 36–47 months, a gain of 1.64 kg. That gives an extrapolated mean of 15.95 kg, which is refitted as a truncated Normal. Child BW is then the mixture (60/72) × MICS fit + (12/72) × extrapolated fit, with mean 11.72 kg against 10.87 kg. The mixture is sampled by the explicit inverse CDF, interpolated on a dense grid of 200,001 points.

## Results

### 1. Sensitivity: arsenic concentration dominates everywhere

| Input | ρ_s with HI, range over 20 cells | Interpretation |
|---|---|---|
| **Arsenic C** | **+0.950 to +0.998** | Rank 1 in every district and population |
| Ingestion rate IR | +0.054 to +0.231 | Rank 2 everywhere; larger where arsenic is less dispersed (Chittagong, Comilla) |
| Body weight BW | −0.022 to −0.157 | Negative, as expected; stronger in children in 9 of 10 districts (child/adult ratio 0.9–3.0) |
| Exposure frequency EF | +0.014 to +0.070 | Small; not distinguishable from zero in 9 of 20 cells |
| Exposure time ET, skin area SA | −0.020 to +0.024 | **Not distinguishable from zero** (dermal is < 0.3% of HI) |

Compared with the base paper: its arsenic coefficients are 0.65–0.90, and ours are higher. The Bangladesh district arsenic distributions are much more dispersed (Gamma shape k = 0.07–0.44, meaning the concentration spans orders of magnitude), so arsenic explains almost all of the rank variation in HI. Mitigation should target arsenic in water, which is the same conclusion as the paper.

### 2. Fitted-parameter uncertainty is far larger than the Monte Carlo error

95% intervals across 1,000 parameter sets:

| District | P95 HI adult: primary [95% interval] | P95 HI child: primary [95% interval] | P(HI>1) adult | P(HI>1) child |
|---|---|---|---|---|
| Dhaka | 21.6 [10.2, 31.2] | 121 [58, 174] | 33% [21, 47] | 48% [32, 65] |
| Chittagong | 6.7 [1.6, 24.8] | 35 [9, 123] | 23% [9, 37] | 54% [40, 65] |
| Rajshahi ⚠ | 3.2 [1.1, 6.2] | 18 [7, 35] | 11% [6, 17] | 20% [13, 29] |
| Khulna | 15.7 [7.5, 24.6] | 84 [40, 133] | 30% [22, 39] | 45% [35, 57] |
| Barisal | 32.1 [9.6, 51.8] | 175 [52, 284] | 41% [31, 50] | 56% [49, 66] |
| Sylhet | 6.6 [4.1, 9.7] | 36 [22, 52] | 33% [24, 42] | 60% [50, 71] |
| Rangpur | 1.8 [1.0, 3.0] | 9.6 [5.0, 15.6] | 10% [5, 15] | 30% [22, 38] |
| Mymensingh | 7.2 [2.7, 11.0] | 38 [14, 59] | 19% [12, 25] | 32% [26, 41] |
| Comilla | 51.3 [40.6, 65.1] | 260 [212, 330] | 72% [61, 85] | 85% [75, 95] |
| Bogra | 6.6 [0.8, 13.8] | 33 [4, 70] | 16% [4, 22] | 32% [21, 41] |

- The **width of the P95 HI interval is 0.45–3.5× the primary value**, while the Monte Carlo CV at N = 10,000 was 2–9% (Phase 8). **The reported risks are limited by the data (roughly 44–110 wells per district), not by the simulation.**
- **Comilla** (n = 110, 8% censored) is the best-determined district. **Chittagong** (n = 44) and **Bogra** have the widest intervals.
- **The child and adult P95 HI intervals do not overlap** in 8 of 10 districts; they overlap only in Chittagong and Bogra, the two least certain. The relative width of the child and adult intervals is almost the same (child/adult ratio 0.95–1.02). That shows BW parameter uncertainty is negligible (n = 22,576 children and 8,013 adults), and almost all of the width comes from arsenic.
- **For Bogra, the primary P95 (6.6) sits above the bootstrap median (4.1).** The heavy Lognormal tail is not a stable feature of Bogra's data, and many resamples give a lighter tail. This agrees with L6 and S2.
- **Every district's child P(HI>1) interval lies above 12.8%, and Comilla's adult interval above 61%.** The conclusion that a substantial share of the population exceeds HI = 1 is robust to parameter uncertainty in every district.

### 3. Model and selection scenarios

Values are ratios to the primary result, using the same uniforms.

| Scenario | P95 HI | P(HI>1) | Mean HI | P95 ELCR |
|---|---|---|---|---|
| S1 Barisal Lognormal | **×3.14** | ×0.90–0.99 | **×17–46** | ×3.14 |
| S2 Bogra Gamma | ×0.57–0.61 | ×1.06–1.08 | ×0.19–0.25 | ×0.57–0.61 |
| S3 Rangpur Lognormal | ×1.87 | ×0.95–1.19 | ×3.9–5.1 | ×1.87 |
| S4 Child BW Gamma | ×0.95–0.985 | ×0.99–1.00 | ×0.95–0.97 | ×0.95–0.985 |
| **S5 Child BW 0–71 months** | **×0.92–0.945** | ×0.97–1.00 | ×0.94 | ×0.92–0.945 |
| S6 IR/SA log-space | **×3.2–3.5** | ×1.07–2.35 | ×3.4–3.6 | ×3.2–3.5 |

- **P(HI>1) is robust to the family choice; P95 and the mean are not.** Swapping Gamma and Lognormal changes P95 by up to 3.1× and the mean by up to 46×, but P(HI>1) by at most 19%. Exceedance near HI = 1 depends on the well-measured middle of the arsenic distribution, while P95 and the mean depend on the extrapolated tail. **P(HI>1) is therefore the most defensible headline metric.** This supports the D7 override: Barisal-Lognormal's 46× mean comes from a tail that no measured well supports.
- **The age gap (L2) is now quantified.** Covering the full 0–6-year exposure window instead of under-5s only **lowers child P95 HI and P95 ELCR by 5.5–8%** and P(HI>1) by 0–3%. The under-5 body weights therefore **overstate** child risk slightly, a conservative bias far smaller than the parameter uncertainty.
- **The child BW family matters little.** The Phase 5 Gamma would have given a 1.5–5% lower child P95, because it under-represents light children, which confirms the direction argued in D8.
- **The ± reading (Phase 6 D1) is the largest single modelling choice.** The old log-space reading would inflate every risk by 3.2–3.5× at P95. It is also physically impossible for SA and cannot reproduce the paper, so the reversal stands, but it needs the team's sign-off.

## Validation gate (master plan section 14)

| Gate | Status |
|---|---|
| Sensitivity uses the exact iterations that produced the reported risks | pass (`test_saved_sensitivity_uses_exact_primary_iterations...` recomputes from the Parquet files) |
| Fixed inputs are not presented as sampled sensitivity variables | pass (only the 6 sampled inputs; `set(input) == STOCHASTIC_INPUTS`) |
| Correlation signs and labels verified against the equations | pass (D21: every significant ρ has the equation-implied sign; BW negative, all others positive) |
| Survey-design limitations stated accurately | pass (D24 Rao–Wu PSU bootstrap; approximations listed in the limitations register) |

Reviewer sign-off (Member 3): ________  Date: ________
