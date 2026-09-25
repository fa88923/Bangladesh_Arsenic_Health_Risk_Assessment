# Phase 8: Convergence and Monte Carlo Error Handoff

Owner: Member 4. Cross-reviewer: Member 3. Member 5 does QA.

## How to reproduce

```bash
python scripts/run_simulation.py    # Phase 7 (the primary runs must exist first)
python scripts/run_convergence.py   # 400 runs, about 10 s; exits 1 if N=10k r0 != the Phase 7 primary run
python -m pytest -q tests/test_convergence.py
```

## Design (design doc sections 12.1 and 17; master plan section 13)

- **Every** district × population (20 cells) at N ∈ {1,000; 5,000; 10,000; 20,000}, with **5 independent seeds** per N: 400 runs and 3.6 M iterations in total.
- Each run is the complete Phase 7 model (`monte_carlo.simulate`) with its own stream, `SeedSequence(20260923, spawn_key=(district, population, N index, replicate))`. The N = 10,000 replicate-0 run **is** the Phase 7 primary run, and it reproduces bit for bit.
- Tracked statistics: mean HI, P50 HI, P95 HI, P(HI>1), P(HI>2) and P95 ELCR (lifetime AT). The first four are the design-doc minimum plus P50; P(HI>2) is added because it is a reported result.
- For each statistic T over the R = 5 seeds: across-seed mean T̄_N; between-seed SD s_N (the Monte Carlo error of one run); CV = s_N / |T̄_N|; and the successive relative change Δ_T(N) = |T̄_N − T̄_{N_prev}| / |T̄_N| × 100%. No N is treated as the true answer.
- For probabilities, the binomial MC SE √(p(1−p)/N) is reported next to s_N as a cross-check.

## Outputs

| File | Contents |
|---|---|
| `results/convergence/convergence_runs.csv` | One row per run: district, population, N, replicate, spawn key, and the 6 statistics (design doc table 19.4) |
| `results/convergence/convergence_summary.csv` | For each district/population/metric/N: across-seed mean, SD, CV, min, max, successive change, binomial SE, and adequacy flag |
| `results/convergence/repeated_seed_stability.csv` | The 5 seeds at N = 10,000 side by side, with CV |
| `results/convergence/convergence_conclusions.csv` | Smallest adequate N for each cell and metric, and whether N = 10,000 is adequate |
| `results/convergence/convergence_manifest.json` | Seeds, grid, criterion, checks, and the Phase 7 manifest hash |
| `results/figures/convergence_<metric>_vs_N.png` | 6 figures, each with 10 district panels: across-seed mean (line) and seed range (band), adult and child |
| `results/figures/convergence_mc_error_decay.png` | Between-seed error vs N on log-log axes, with the N^−½ reference |

## Decision log (2026-09-25, Member 4)

### D16. Numerical convergence criterion

The plans require successive relative change, replicate seeds and MC SE, but give no numeric threshold. A statistic is **adequate at N** when:

- **for continuous statistics** (mean, P50, P95 HI; P95 ELCR): Δ_T(N) ≤ 5% **and** between-seed CV ≤ 5%;
- **for probabilities** (P(HI>1), P(HI>2)): |T̄_N − T̄_{N_prev}| ≤ 0.01 **and** s_N ≤ 0.01, that is, 1 percentage point.

The **smallest adequate N** is the smallest N from which every larger N in the grid is also adequate.

*Why these values.* Both conditions are needed. Δ alone can look stable by chance, and the plan requires more than one seed. CV alone says nothing about drift. A 5% tolerance is small compared with the fitted-parameter uncertainty that Phase 9 will show, so the reported numbers are limited by the data, not by the simulation. Probabilities use an absolute tolerance because a relative one explodes as p → 0, and 1 point is finer than the precision at which exceedance percentages are reported. N = 1,000 can never be adequate by construction, because it has no preceding size.

### D17. Iterations of the 400 convergence runs are not stored

Only the per-run statistics are kept. Every run is reproducible from its spawn key, and the primary N = 10,000 r0 iterations that Phase 9 needs are already saved by Phase 7. Storing everything would add about 450 MB for no analytical gain.

### D18. N = 10,000 stays the primary run size

The plan fixes N = 10,000 as the base-paper-comparable primary result and N = 20,000 as a stability check (master plan section 13). The evidence supports keeping it:

| Metric | Adequate at N = 10,000 | Stabilises within grid | Worst between-seed error at 10k |
|---|---:|---:|---|
| **P(HI > 1)** (headline) | **20 / 20** | 20 / 20 | SD 0.67 points (Barisal child) |
| P(HI > 2) | 19 / 20 | 20 / 20 | Barisal adult: change of 1.02 points, just over the 1.00 tolerance |
| **P95 HI** (base-paper Table 5) | **16 / 20** | 20 / 20 | CV 8.6% (Bogra adult) |
| P95 ELCR | 16 / 20 | 20 / 20 | CV 8.7% (Bogra adult) |
| Mean HI | 16 / 20 | 17 / 20 | **Bogra: CV 26–27%, not converged at 20k** |
| P50 HI | 10 / 20 | 12 / 20 | Rajshahi adult: CV 17% (median HI 0.0002) |

Recommendation for the final tables: report the N = 10,000 values with their between-seed SD as the Monte Carlo error. For the four P95 cells that pass only at 20k (Bogra adult and child, Chittagong child, Sylhet adult), also state the 20k value; the gap is at most 9.1%, for Bogra child. Raising the primary N was considered and rejected: it would break comparability with the paper and gain only a factor of √2 in precision.

### D19. Statistics that must not be reported as converged

- **Bogra mean HI** (adult and child). The selected Lognormal has σ_log = 2.86. A single draw of C then has a CV of about √(e^{σ²} − 1) ≈ 60, so a 5% CV on the mean would need about 1.4 M iterations. The sample mean is dominated by rare extreme draws (child maximum HI 15,865). **Report Bogra's median and P95, and give its mean only with this caveat.** This is limitation L6 appearing numerically.
- **P50 HI** in the high-censoring districts. It sits in the extrapolated region below the detection limit, where HI changes steeply with C (limitation L3). Most of those medians are ≪ 1 (Rajshahi 0.0002, Mymensingh 0.01), so none of them affects a risk conclusion. The exceptions are Dhaka child (0.71, CV 5.6%) and Barisal child (1.98, CV 9.1% at 10k and 3.8% at 20k); report those with their seed spread.

## Results

### Monte Carlo error behaves as theory predicts

The median between-seed error of every metric falls as N^−½ (`convergence_mc_error_decay.png`). For example, the median P95 HI CV falls from 11.1% at 1k to 3.3% at 10k and 2.8% at 20k. For probabilities, the observed between-seed SD averages **0.83×** the binomial MC SE. That is consistent with 1.0, given that an SD estimated from 5 seeds is itself noisy. So the binomial formula in design doc section 17.3 is a valid precision statement for every reported P(HI>1).

### Monte Carlo precision of the primary (N = 10,000) results

| District | P95 HI CV (adult / child) | SD of P(HI>1), points (adult / child) | P95 HI, r0 at 10k vs 20k (adult / child) |
|---|---|---|---|
| Dhaka | 3.1% / 2.6% | 0.26 / 0.39 | −2.7% / −2.1% |
| Chittagong | 3.1% / 5.1% | 0.28 / 0.43 | −1.1% / −3.4% |
| Rajshahi | 4.1% / 4.0% | 0.23 / 0.14 | −1.7% / −2.2% |
| Khulna | 3.5% / 2.9% | 0.31 / 0.53 | −2.5% / −1.9% |
| Barisal | 4.4% / 3.8% | 0.33 / 0.67 | −1.4% / −3.5% |
| Sylhet | 2.2% / 2.0% | 0.40 / 0.40 | +0.6% / +2.4% |
| Rangpur | 3.2% / 3.3% | 0.21 / 0.30 | +5.9% / −2.1% |
| Mymensingh | 4.6% / 2.9% | 0.53 / 0.13 | −5.9% / −5.7% |
| Comilla | 3.0% / 2.8% | 0.64 / 0.47 | +6.1% / −3.3% |
| Bogra | **8.6% / 5.4%** | 0.47 / 0.36 | +4.4% / **+9.1%** |

### Conclusion: smallest adequate simulation size

**N = 10,000 is adequate** for the headline P(HI > 1) in every district and population, and for P95 HI and P95 ELCR in 16 of 20 cells. **N = 20,000 is needed** for P95 in Bogra (both populations), Chittagong child and Sylhet adult. **No N in the grid is enough** for Bogra mean HI; it is an estimator problem caused by the heavy tail, not a coding problem. The N = 20,000 results are an additional check, not the reference truth.

### Monte Carlo error vs fitted-parameter uncertainty

Everything above is **Monte Carlo (numerical) error**: the spread from re-running the *same* fitted model with different random numbers. It shrinks as N grows. **Fitted-parameter uncertainty** is separate: the spread from refitting the distributions to different samples of wells and children. It does **not** shrink with N. It is assessed in Phase 9, using the Phase 4 bootstrap parameter ranges; for example, Dhaka's Gamma shape k has a 95% interval of [0.089, 0.275]. At N = 10,000 the Monte Carlo error (≤ 9% CV) is expected to be far smaller than the parameter uncertainty.

## Validation gate (master plan section 13)

| Gate | Status |
|---|---|
| Stabilisation assessed by successive relative change | pass (Δ_T(N) for every cell, metric and N) |
| Upper-tail results assessed separately from central estimates | pass (P95 HI and P95 ELCR separate from mean and P50; tail CV is higher, as expected) |
| Convergence claims supported by more than one seed | pass (5 seeds per N; adequacy requires between-seed CV or SD within tolerance) |
| Any metric that fails to stabilise is explicitly reported | pass (D19: Bogra mean HI; P50 in 8 cells) |

Reviewer sign-off (Member 3): ________  Date: ________
