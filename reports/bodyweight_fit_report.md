# Body-Weight Fit Report

## Weighted ECDF and Fit Method

For each population, observations were sorted by measured body weight and tied body-weight values were aggregated. Survey weights were normalized as `q = w / sum(w)`. The empirical CDF used midpoint plotting positions:

$$
p_j
=
\sum_{r<j} q_r
+
\frac{1}{2}q_j
$$

where `q_j` is the normalized survey mass at the distinct body-weight value `x_j`.

Each Phase 5 candidate distribution was fitted by survey-weighted CDF least squares:

$$
\operatorname*{arg\,min}_{\theta}
\sum_j q_j\left[F(x_j;\theta)-p_j\right]^2
$$

This was the primary fitting criterion. Survey-weighted pseudo-MLE fits were also computed as robustness checks. In the Phase 7 final review, a child Normal truncated below at the lightest observed weight (1.6 kg) was fitted and assessed as a physically admissible alternative to the provisional Gamma model.

## Fit Results

| Population | Distribution | LS parameters | SSE_w | D_w | Tail error | Final simulation selection |
|---|---|---|---:|---:|---:|:---:|
| Adult | Normal | &mu;=55.08, &sigma;=10.91 | 2.28e-4 | 0.030 | 12.7% |  |
| Adult | Lognormal | &mu;<sub>log</sub>=4.0024, &sigma;<sub>log</sub>=0.2002 | 1.96e-5 | 0.019 | 1.0% | &#10003; |
| Adult | Gamma | k=25.27, &theta;=2.199 | 4.33e-5 | 0.020 | 3.6% |  |
| Adult | Triangular | a=32.3, m=51.4, b=82.8 | 6.94e-5 | 0.027 | 10.0% |  |
| Child | Normal | &mu;=10.87, &sigma;=3.22 | 2.85e-5 | 0.009 | 5.2% |  |
| Child | Lognormal | &mu;<sub>log</sub>=2.3706, &sigma;<sub>log</sub>=0.2993 | 5.96e-4 | 0.043 | 49.7% |  |
| Child | Gamma | k=11.392, &theta;=0.9725 | 3.15e-4 | 0.032 | 37.2% |  |
| Child | Triangular | a=3.17, m=11.23, b=18.09 | 3.08e-5 | 0.015 | 19.6% |  |
| Child | Truncated Normal | parent &mu;=10.8570, parent &sigma;=3.2281, L=1.6 kg | 2.94e-5 | 0.0096 | 1.7% | &#10003; |

Tail error is the largest relative quantile error over P1, P5, P95, and P99. `D_w` is the maximum weighted CDF discrepancy.

## Diagnostic Figures

### Adult

![Adult body-weight fit diagnostics](../results/figures/bodyweight_fit_diagnostics_adult.png)

### Child

![Child body-weight fit diagnostics](../results/figures/bodyweight_fit_diagnostics_child.png)

## Final Body-Weight Models Used in Simulation

| Population | Final model | Parameters |
|---|---|---|
| Adult | Lognormal | &mu;<sub>log</sub>=4.002396, &sigma;<sub>log</sub>=0.200174 |
| Child | Truncated Normal | parent &mu;=10.857019, parent &sigma;=3.228093, L=1.6 kg |

The child lower bound is the lightest cleaned MICS child weight, so the final model does not extrapolate below the observed data. The truncated Normal replaces the provisional Gamma selection: its weighted CDF SSE is about 11 times smaller, its worst P1-P99 quantile error is 1.7% rather than 36%, and its $E[1/BW]$ is within 0.1% of the weighted empirical value. These dose-relevant checks drove D8 in the simulation report.
