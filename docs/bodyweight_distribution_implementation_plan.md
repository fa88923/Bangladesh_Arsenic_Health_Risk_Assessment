# Body-Weight Distribution Fitting Implementation Plan
## Bangladesh Arsenic Health-Risk Monte Carlo Model

## 1. Purpose and methodological position

This plan prepares Bangladesh-specific body-weight (`BW`) distributions for the adult and child Monte Carlo populations in the CSE 402 arsenic risk project.

The project specification requires **least-squares distribution fitting**. The base paper, Yadav & Kalkal (2024), treats body weight as a stochastic Monte Carlo input and uses a **Lognormal adult BW model** and a **Triangular child BW model**. The Bangladesh extension replaces those assumed BW distributions with distributions fitted to Bangladesh survey microdata.

For this implementation, the candidate BW distributions are treated as genuine alternatives for **both adults and children**:

1. Normal;
2. Lognormal;
3. Gamma;
4. Triangular.

Triangular is not ignored or treated as invalid. It is a legitimate candidate distribution. It has bounded support and therefore must be judged from the data just like the other candidates. Its inclusion is additionally defensible because the base paper itself uses a triangular child-BW model, and Joo & Casella (2001, *Environmetrics*, DOI: 10.1002/env.489) studied triangular-distribution estimation specifically for Monte Carlo risk analysis and reported that least-squares/quantile-based estimation was preferable to MLE for that distribution.

Because STEPS and MICS are **complex survey datasets**, survey weights must be respected. The implementation therefore uses two complementary fitting routes:

- **Primary CSE 402 numerical fit:** survey-weighted CDF least-squares minimum-distance fitting;
- **Scientific robustness fit:** survey-weighted pseudo-maximum likelihood.

The weighted pseudo-likelihood is an established complex-survey method. The exact weighted CDF least-squares objective below is the project implementation of the CSE 402 least-squares requirement. It should not be described in the report as the unique standard complex-survey estimator.

---

# 2. Adult body-weight dataset: Bangladesh STEPS 2018

## 2.1 Exact file

Use:

```text
Datasets/bgd2018.csv
```

The dataset audit reports **8,185 rows × 425 columns**. WHO's Bangladesh STEPS 2018 metadata independently confirms the same public-use dataset structure.

## 2.2 Exact columns to extract

| Column | Meaning | Role in BW pipeline |
|---|---|---|
| `m12` | Measured weight in kilograms | **Body-weight observation to fit** |
| `wstep2` | Final analysis weight for Step 2 physical measures | **Survey weight** |
| `age` | Age for analysis | Define/check adult population |
| `sex` | Sex | Diagnostic/subgroup checking |
| `m11` | Measured height in cm | Anthropometric QC only |
| `pid` | Participant ID | Duplicate/record checking |
| `psu` | Primary sampling unit | Retain for design-aware uncertainty analysis |
| `stratum` | Sampling stratum | Retain for design-aware uncertainty analysis |
| `division` | Division | Optional representativeness diagnostics |
| `urbanrural` | Urban/rural indicator | Optional representativeness diagnostics |

### Critical interpretation

`m12` is the measured body weight. **`wstep2` is not body weight.** WHO defines `wstep2` as the final analysis weight for Step 2 physical measurements. WHO also states that the STEPS analysis weights account for unequal selection probabilities and are adjusted for the age-sex composition of the target population.

The audit reports:

- `m12`: 167 missing values (2.04%);
- `wstep2`: no missing values;
- valid complete combinations of `m12 + age + sex + wstep2`: 8,018 / 8,185;
- `m12 = 888` occurs and is not a physical body weight.

WHO STEPS questionnaire documentation uses:

- `888` = participant refused weight measurement;
- `666` = participant too large for the scale.

The audit shows `888` in this dataset; `666` should also be screened for in code even if absent from the audited frequencies.

## 2.3 Adult population definition

Use the STEPS survey's observed adult age range:

$$
18 \leq age \leq 69.
$$

Therefore the fitted adult distribution must be described as the **Bangladesh STEPS 2018 body-weight distribution for adults aged 18–69**, not as a distribution for every possible adult age.

## 2.4 Adult preprocessing algorithm

### Step A1: Load only required columns

```python
import pandas as pd

cols = [
    "pid", "psu", "stratum", "division", "urbanrural",
    "wstep2", "age", "sex", "m11", "m12"
]
df = pd.read_csv("Datasets/bgd2018.csv", usecols=cols)
```

### Step A2: Restrict to the intended adult population

```python
df = df[df["age"].between(18, 69)]
```

### Step A3: Remove records without usable BW or survey weight

Remove:

- missing `m12`;
- missing/nonpositive `wstep2`;
- known STEPS special weight codes `888` and `666`;
- nonpositive physical weight values.

```python
df = df.dropna(subset=["m12", "wstep2", "age", "sex"])
df = df[~df["m12"].isin([666, 888])]
df = df[(df["m12"] > 0) & (df["wstep2"] > 0)]
```

### Step A4: Check duplicate IDs

```python
duplicate_pid = df["pid"].duplicated(keep=False)
```

Investigate duplicates rather than automatically deleting them.

### Step A5: Anthropometric QC without arbitrary outlier trimming

Do **not** automatically delete the lowest/highest 1%, IQR outliers, or unusually heavy/light adults. Legitimate extremes belong to the population distribution.

Use `m11` only as a diagnostic when available:

$$
BMI_i = \frac{BW_i}{(H_i/100)^2}.
$$

```python
valid_h = df["m11"].notna() & (~df["m11"].isin([666, 888])) & (df["m11"] > 0)
df.loc[valid_h, "bmi_check"] = (
    df.loc[valid_h, "m12"] / (df.loc[valid_h, "m11"] / 100.0) ** 2
)
```

Use this to investigate implausible records, not to silently impose a BMI-based exclusion rule that was not pre-specified.

### Step A6: Preserve survey-design variables

Keep `psu` and `stratum` in the cleaned file. They are not inputs to the Monte Carlo risk equation, but they are useful if design-aware bootstrap or variance estimation of fitted BW parameters is added.

### Step A7: Output adult cleaned file

Recommended columns:

```text
pid, age, sex, BW_kg, survey_weight, psu, stratum, division, urbanrural, height_cm
```

Rename:

```python
df = df.rename(columns={
    "m12": "BW_kg",
    "wstep2": "survey_weight",
    "m11": "height_cm"
})
```

---

# 3. Child body-weight dataset: Bangladesh MICS7

## 3.1 Exact file

The dataset audit shows that the anthropometric child file is:

```text
Datasets/BGD_2025_MICS7_v01_M/BGD_2025_MICS7_v01_M/
BGD_2025_MICS7_Datasets/Bangladesh MICS7 Datasets/
Bangladesh MICS7 Datasets/Bangladesh MICS7 SPSS Datasets/ch.sav
```

The important filename is:

```text
ch.sav
```

Do **not** use `bh.sav` for the child BW fit. `bh.sav` is the birth-history file. The audited `ch.sav` contains **24,680 rows × 445 columns** and contains the child anthropometric measurements.

Note on naming: the audited package is labelled `BGD_2025_MICS7...`. The project should not relabel it as a 2024 dataset in the methods section unless the official survey documentation explicitly supports that wording.

## 3.2 Exact columns to extract

| Column | Meaning | Role in BW pipeline |
|---|---|---|
| `AN8` | Child's measured weight (kg) | **Body-weight observation to fit** |
| `chweight` | Children under 5's sample weight | **Survey weight** |
| `CAGE` | Age in months | Define child population |
| `HL4` | Sex | Diagnostic/subgroup checking |
| `HH7A` | District | Geographic diagnostics; not needed for national BW fit |
| `AN11` | Child length/height (cm) | Anthropometric QC |
| `WAZFLAG` | WHO weight-for-age error flag | **Primary anthropometric plausibility flag for BW** |
| `HAZFLAG` | WHO height-for-age error flag | Supporting QC |
| `WHZFLAG` | WHO weight-for-height error flag | Supporting QC |
| `HH1` / `PSU` | Cluster / primary sampling unit | Retain for design-aware uncertainty analysis |
| `stratum` | Sample stratum | Retain for design-aware uncertainty analysis |

The audit directly identifies `AN8` as **Child's weight (kilograms)** and `chweight` as **Children under 5's sample weight**.

## 3.3 Exact child weight special codes found in the audit

The audit shows the following non-physical `AN8` values:

| `AN8` value | Label in audit | Action |
|---:|---|---|
| `99.3` | CHILD NOT PRESENT AFTER REVISITS | Exclude |
| `99.4` | CHILD REFUSED | Exclude |
| `99.5` | RESPONDENT REFUSED | Exclude |
| `99.6` | OTHER | Exclude |
| missing | no measured value | Exclude |

Do not interpret any of these values as kilograms.

## 3.4 Child population definition

The audited `CAGE` values run from 0 to 59 months. Therefore the BW distribution supported by this file is:

$$
0 \leq CAGE \leq 59\ \text{months}.
$$

The report should call this the **under-5 child BW distribution**.

Important model limitation: the base paper uses a child exposure duration of 6 years, but this MICS file does not provide measured BW for ages 60–71 months. Do not extrapolate the fitted under-5 BW distribution and then claim that it was directly measured for ages 0–6 years. If the project requires a 0–6-year BW population, an additional dataset covering 60–71 months is required.

## 3.5 Child preprocessing algorithm

### Step C1: Read the SPSS file while preserving numeric codes

```python
import pyreadstat

child, meta = pyreadstat.read_sav(
    ".../Bangladesh MICS7 SPSS Datasets/ch.sav",
    apply_value_formats=False
)
```

Keep only the required variables:

```python
cols = [
    "HH1", "PSU", "stratum", "HH7A", "HL4", "CAGE",
    "AN8", "AN11", "WAZFLAG", "HAZFLAG", "WHZFLAG", "chweight"
]
child = child[[c for c in cols if c in child.columns]].copy()
```

### Step C2: Restrict to measured under-5 children

```python
child = child[child["CAGE"].between(0, 59)]
```

### Step C3: Remove missing and coded non-measurements

```python
special_an8 = {99.3, 99.4, 99.5, 99.6}
child = child.dropna(subset=["AN8", "CAGE", "HL4", "chweight"])
child = child[~child["AN8"].isin(special_an8)]
child = child[(child["AN8"] > 0) & (child["chweight"] > 0)]
```

Any unexpected `AN8 >= 90` not already documented above should be **flagged for manual codebook review**, not silently interpreted as weight.

### Step C4: Apply WHO weight-for-age error flag

For MICS, `WAZFLAG` is coded:

- `0` = No error;
- `1` = Error flag.

The primary BW dataset should therefore retain:

```python
child = child[child["WAZFLAG"] == 0]
```

This is preferable to arbitrary percentile trimming because it uses the survey's own WHO anthropometric plausibility processing.

`HAZFLAG` and `WHZFLAG` should be retained for diagnostics, but a child should not automatically be removed from a **weight-only** fit merely because a height-related flag fails, unless the corresponding record investigation shows that the weight measurement itself is problematic.

### Step C5: Do not trim valid tails arbitrarily

Do not automatically remove valid low/high measured weights using z-score cutoffs invented for this project, IQR rules, or 1st/99th percentile trimming after `WAZFLAG` filtering. Genuine low or high child weights affect the population BW distribution and therefore the risk simulation.

### Step C6: Output child cleaned file

Recommended columns:

```text
CAGE_months, sex, BW_kg, survey_weight, district, PSU, stratum,
height_cm, WAZFLAG, HAZFLAG, WHZFLAG
```

Rename:

```python
child = child.rename(columns={
    "CAGE": "CAGE_months",
    "HL4": "sex",
    "AN8": "BW_kg",
    "chweight": "survey_weight",
    "HH7A": "district",
    "AN11": "height_cm"
})
```

---

# 4. Survey-weight preparation

For either population let the cleaned observations be

$$
(x_i,w_i),\qquad i=1,\ldots,n,
$$

where:

- $x_i$ = measured body weight in kg;
- $w_i$ = STEPS `wstep2` or MICS `chweight`.

For CDF calculations, define normalized survey mass:

$$
q_i=\frac{w_i}{\sum_{r=1}^{n}w_r},
\qquad
\sum_i q_i=1.
$$

Multiplying every survey weight by the same constant does not alter the weighted ECDF or the pseudo-MLE parameter estimate. For numerical stability in pseudo-likelihood optimization, weights may additionally be normalized to sum to the sample size:

$$
\tilde w_i
=
\frac{n w_i}{\sum_{r=1}^{n}w_r}.
$$

Python:

```python
q = w / w.sum()
w_pml = len(w) * w / w.sum()
```

---

# 5. Construct the survey-weighted empirical CDF

Sort the observed BW values and aggregate survey weight for tied values. Let the distinct sorted weights be

$$
x_{(1)}<x_{(2)}<\cdots<x_{(J)},
$$

with total survey mass $q_j$ at each distinct value.

Use weighted midpoint plotting positions:

$$
p_j
=
\left(\sum_{r<j}q_r\right)+\frac{q_j}{2}.
$$

This avoids comparing a smooth theoretical CDF only to the upper edge of each empirical CDF jump.

Python implementation tools:

```python
import numpy as np
import pandas as pd

agg = (
    pd.DataFrame({"x": x, "w": w})
      .groupby("x", as_index=False)["w"].sum()
      .sort_values("x")
)

x_u = agg["x"].to_numpy()
q_u = agg["w"].to_numpy(dtype=float)
q_u = q_u / q_u.sum()
p_u = np.cumsum(q_u) - 0.5 * q_u
```

The ordinary survey-weighted ECDF itself is

$$
\hat F_w(x)
=
\frac{\sum_i w_i\mathbf 1(x_i\leq x)}{\sum_iw_i}.
$$

---

# 6. Primary fitting method: survey-weighted CDF least squares

For each candidate distribution with CDF $F(x;\theta)$, estimate parameters by minimizing

$$
\boxed{
\hat\theta_{LS}
=
\arg\min_{\theta}
\sum_{j=1}^{J}
q_j
\left[F(x_{(j)};\theta)-p_j\right]^2
}
$$

Equivalently define residuals

$$
r_j(\theta)
=
\sqrt{q_j}\left[F(x_{(j)};\theta)-p_j\right]
$$

and call:

```python
scipy.optimize.least_squares(residual_function, initial_parameters)
```

The objective value to report is

$$
SSE_w
=
\sum_jq_j\left[F(x_{(j)};\hat\theta)-p_j\right]^2.
$$

Lower $SSE_w$ indicates closer agreement with the survey-weighted empirical distribution under this project criterion.

### Important methodological wording

This is a defensible **survey-weighted minimum-distance implementation of the project's required least-squares fitting**. The standard complex-survey benchmark is the weighted pseudo-likelihood fit in Section 8. Do not claim that this particular weighted LS formula is the only standard estimator used in survey statistics.

---

# 7. Candidate distributions and exact parameterizations

## 7.1 Normal

$$
BW\sim N(\mu,\sigma^2),\qquad \sigma>0.
$$

CDF:

$$
F_N(x)=\Phi\left(\frac{x-\mu}{\sigma}\right).
$$

Python:

```python
from scipy.stats import norm
norm.cdf(x, loc=mu, scale=sigma)
```

Recommended initial values:

$$
\mu_0=\sum_iq_ix_i,
$$

$$
\sigma_0=\sqrt{\sum_iq_i(x_i-\mu_0)^2}.
$$

Physical-support check after fitting:

$$
P(BW<0)=F_N(0;\hat\mu,\hat\sigma).
$$

Report this probability. A Normal model that places meaningful probability on negative BW should not be used as the physical Monte Carlo sampling distribution unless a truncated-normal model is explicitly introduced and justified.

---

## 7.2 Lognormal

$$
\ln(BW)\sim N(\mu_{\log},\sigma_{\log}^2),
\qquad BW>0.
$$

CDF:

$$
F_{LN}(x)
=
\Phi\left(\frac{\ln x-\mu_{\log}}{\sigma_{\log}}\right),
\qquad x>0.
$$

SciPy parameterization:

```python
from scipy.stats import lognorm
lognorm.cdf(
    x,
    s=sigma_log,
    loc=0,
    scale=np.exp(mu_log)
)
```

Initial values are the survey-weighted mean and SD of $\ln(BW)$.

Keep `loc=0`; do not introduce an arbitrary three-parameter shifted Lognormal unless there is a separate scientific justification.

---

## 7.3 Gamma

Use shape-scale parameterization:

$$
BW\sim \operatorname{Gamma}(k,\theta),
\qquad k>0,\;\theta>0.
$$

PDF:

$$
f(x;k,\theta)
=
\frac{x^{k-1}e^{-x/\theta}}
{\Gamma(k)\theta^k},
\qquad x>0.
$$

CDF:

$$
F_{\Gamma}(x;k,\theta)
=
\frac{\gamma(k,x/\theta)}{\Gamma(k)}.
$$

Python:

```python
from scipy.stats import gamma
gamma.cdf(x, a=k, loc=0, scale=theta)
```

Moment-based starting values:

$$
k_0=\frac{\bar x_w^2}{s_w^2},
\qquad
\theta_0=\frac{s_w^2}{\bar x_w}.
$$

Keep `loc=0` because body weight is a positive physical quantity and a free shift is not required by the model.

---

## 7.4 Triangular

The triangular distribution is a **valid fourth candidate**:

$$
BW\sim \operatorname{Triangular}(a,m,b),
$$

with

$$
0<a<m<b,
$$

where:

- $a$ = lower support bound;
- $m$ = mode;
- $b$ = upper support bound.

Its CDF is

$$
F_T(x)=
\begin{cases}
0, & x<a,\\[4pt]
\dfrac{(x-a)^2}{(b-a)(m-a)}, & a\leq x\leq m,\\[10pt]
1-\dfrac{(b-x)^2}{(b-a)(b-m)}, & m<x\leq b,\\[10pt]
1, & x>b.
\end{cases}
$$

SciPy uses

$$
c=\frac{m-a}{b-a},
\qquad
loc=a,
\qquad
scale=b-a.
$$

Python:

```python
from scipy.stats import triang
c = (m - a) / (b - a)
triang.cdf(x, c=c, loc=a, scale=b-a)
```

### Triangular least-squares fit

Use the **same primary CDF least-squares criterion** as the other candidates:

$$
(\hat a,\hat m,\hat b)
=
\arg\min_{0<a<m<b}
\sum_jq_j
\left[F_T(x_{(j)};a,m,b)-p_j\right]^2.
$$

To enforce the ordering during numerical optimization, parameterize, for example,

$$
a=e^{\alpha},
$$

$$
m=a+e^{\beta},
$$

$$
b=m+e^{\gamma}.
$$

Then optimize $(\alpha,\beta,\gamma)$ using `scipy.optimize.least_squares` and transform back to $(a,m,b)`.

The fitted support must also be checked against the observed BW range. A candidate that excludes non-negligible observed survey mass is not acceptable.

### Why triangular is methodologically legitimate here

Triangular models are established in quantitative risk analysis. Joo & Casella (2001) specifically studied Normal, Lognormal, and Triangular predictive distributions in Monte Carlo risk analysis and found **quantile least-squares estimation preferable to MLE for the triangular case**. Thus, a triangular BW candidate fitted through a least-squares distributional criterion is not an invented assumption.

The triangular model should nevertheless win or lose on fit. Its bounded support can be advantageous if the data are clearly bounded and approximately triangular, but it should not be selected merely because the base paper used it for children.

---

# 8. Robustness fit: survey-weighted pseudo-maximum likelihood

Complex-survey literature supports estimation by weighted pseudo-likelihood. For candidate density $f(x;\theta)$, estimate

$$
\boxed{
\hat\theta_{PML}
=
\arg\max_{\theta}
\sum_{i=1}^{n}
\tilde w_i\log f(x_i;\theta)
}
$$

or minimize the negative form:

$$
-
\sum_i\tilde w_i\log f(x_i;\theta).
$$

Python:

```python
from scipy.optimize import minimize
```

Use:

```python
norm.logpdf(...)
lognorm.logpdf(...)
gamma.logpdf(...)
triang.logpdf(...)
```

Example structure:

```python
def negative_weighted_loglik(params, x, w, dist_name):
    logf = ...  # corresponding scipy.stats.<dist>.logpdf
    if np.any(~np.isfinite(logf)):
        return np.inf
    return -np.sum(w * logf)

result = minimize(objective, x0=initial_params, method="L-BFGS-B")
```

### Distribution-specific constraints

- Normal: $\sigma>0$;
- Lognormal: $\sigma_{\log}>0$, `loc=0`;
- Gamma: $k>0$, $\theta>0$, `loc=0`;
- Triangular: $0<a<m<b$ and fitted support must contain the observed values contributing positive survey weight.

For Triangular, treat pseudo-MLE as a **robustness comparison**, not automatically as the preferred estimator. The triangular likelihood has special boundary/support behaviour, and the risk-analysis literature gives direct support for least-squares/quantile-based fitting.

Do not use plain calls such as

```python
scipy.stats.gamma.fit(data)
```

as the main scientific fit because they do not incorporate STEPS/MICS survey weights.

They may be used only as debugging/initialization aids.

---

# 9. Fit evaluation and distribution selection

The adult and child distributions are selected **separately**. Do not force the same distribution family on both populations.

## 9.1 Primary numerical ranking

For every candidate calculate:

$$
SSE_w
=
\sum_jq_j
\left[F(x_{(j)};\hat\theta_{LS})-p_j\right]^2.
$$

Lower is better.

## 9.2 Weighted KS-style discrepancy

Calculate

$$
D_w
=
\max_j
\left|
F(x_{(j)};\hat\theta)-p_j
\right|.
$$

Python:

```python
D_w = np.max(np.abs(fitted_cdf - p_u))
```

Lower is better.

Because ordinary `scipy.stats.kstest()` assumes an ordinary iid empirical sample, do not present its default p-value as a design-correct survey goodness-of-fit test. A directly calculated weighted CDF discrepancy is safer for this project.

## 9.3 Graphical diagnostics

For each candidate produce:

1. weighted histogram + fitted PDF;
2. weighted ECDF + fitted CDF;
3. weighted empirical quantiles vs fitted theoretical quantiles (Q-Q plot);
4. lower-tail and upper-tail zooms.

Python tools:

```python
matplotlib.pyplot.hist(..., weights=survey_weight, density=True)
matplotlib.pyplot.plot(...)
scipy.stats.<distribution>.pdf(...)
scipy.stats.<distribution>.cdf(...)
scipy.stats.<distribution>.ppf(...)
```

## 9.4 Robustness comparison with pseudo-MLE

For each family, compare the LS-fitted parameters with the pseudo-MLE parameters. Large disagreements should trigger investigation of:

- tail behaviour;
- bounded-support effects;
- influential survey weights;
- remaining miscoded observations;
- optimizer convergence.

## 9.5 Physical admissibility

A fitted distribution is not acceptable merely because its numeric SSE is small.

Check:

- sampled BW must be positive;
- Normal negative-tail probability must be reported;
- Triangular fitted support must cover the relevant observed data;
- Lognormal/Gamma must retain `loc=0` unless a shifted model is scientifically justified;
- no fitted model may treat survey special codes as physical BW.

## 9.6 Final selection rule

Use the following predeclared rule:

1. Rank all four candidates by weighted CDF $SSE_w$.
2. Check weighted $D_w$ and ECDF/CDF overlays.
3. Check Q-Q/tail behaviour.
4. Reject physically inadmissible models.
5. Compare LS parameters with weighted pseudo-MLE robustness parameters.
6. Select the distribution with the strongest overall support, documenting any case where the smallest SSE model is rejected for physical or tail-fit reasons.

**Do not use ordinary iid AIC/AICc computed from raw survey-weighted pseudo-likelihood as if it had its usual iid interpretation.** If an information criterion is required later, use a survey-design-appropriate criterion and state its methodology explicitly.

---

# 10. Monte Carlo sampling after the BW distribution is selected

Once the best-supported adult and child distributions are selected, their samples are **body weights in kilograms**.

For each Monte Carlo iteration:

$$
U_i\sim\operatorname{Uniform}(0,1),
$$

$$
BW_i
=
F_{BW}^{-1}(U_i;\hat\theta).
$$

Python:

```python
rng = np.random.default_rng(seed)
u = rng.uniform(0.0, 1.0, size=N)
```

Then:

```python
# Normal
BW = norm.ppf(u, loc=mu, scale=sigma)

# Lognormal
BW = lognorm.ppf(u, s=sigma_log, loc=0, scale=np.exp(mu_log))

# Gamma
BW = gamma.ppf(u, a=k, loc=0, scale=theta)

# Triangular
c = (m - a) / (b - a)
BW = triang.ppf(u, c=c, loc=a, scale=b-a)
```

These sampled values are directly inserted as $BW_i$ in the existing risk equations:

$$
ADD_{\text{ing},i}
=
\frac{C_i\times IR_i\times EF_i\times ED}
{BW_i\times AT},
$$

$$
ADD_{\text{dermal},i}
=
\frac{C_i\times K_p\times EF_i\times ED\times ET_i\times SA_i\times CF}
{BW_i\times AT}.
$$

Thus, fitting a distribution to body weight does **not** transform BW into another physical quantity. Sampling the selected fitted BW distribution produces simulated body weights in **kg**.

---

# 11. Recommended Python function structure

```python
def load_steps_bw(path):
    """Load the exact STEPS columns needed for adult BW."""


def clean_steps_bw(df):
    """Apply age, missing-value, survey-weight, and STEPS-code rules."""


def load_mics_bw(path):
    """Load ch.sav and retain child anthropometry/design columns."""


def clean_mics_bw(df):
    """Apply CAGE, AN8 special-code, chweight, and WAZFLAG rules."""


def weighted_ecdf_midpoints(x, w):
    """Return unique x, normalized survey mass q, and midpoint CDF p."""


def weighted_mean_var(x, w):
    """Provide starting values and descriptive statistics."""


def cdf_ls_residuals(params, x, q, p, dist_name):
    """sqrt(q) * [F_theta(x) - p]."""


def fit_cdf_least_squares(x, w, dist_name):
    """Fit Normal/Lognormal/Gamma/Triangular using scipy.optimize.least_squares."""


def negative_weighted_loglik(params, x, w, dist_name):
    """Negative normalized survey-weighted pseudo-log-likelihood."""


def fit_weighted_pseudo_mle(x, w, dist_name):
    """Robustness fit using scipy.optimize.minimize."""


def evaluate_distribution_fit(x, w, ls_fit, pml_fit, dist_name):
    """SSE_w, D_w, support checks, diagnostics."""


def select_bw_distribution(fit_table):
    """Apply the predeclared selection rule."""


def sample_bw(selected_fit, N, seed):
    """Uniform + scipy.stats.<dist>.ppf inverse-transform sampling."""
```

---

# 12. Required outputs

## 12.1 Adult preprocessing outputs

```text
results/bodyweight/adult_bw_clean.csv
results/bodyweight/adult_bw_cleaning_audit.csv
```

Cleaning audit should report:

- starting $n$;
- missing `m12` removed;
- `888` removed;
- `666` removed if present;
- missing/nonpositive survey weights removed;
- final $n$;
- weighted mean, median, SD, P5, P50, P95.

## 12.2 Child preprocessing outputs

```text
results/bodyweight/child_bw_clean.csv
results/bodyweight/child_bw_cleaning_audit.csv
```

Cleaning audit should separately report counts removed for:

- missing `AN8`;
- `99.3`;
- `99.4`;
- `99.5`;
- `99.6`;
- age outside 0–59 months;
- `chweight <= 0`;
- `WAZFLAG = 1`;
- final $n$.

## 12.3 Fit-result table

For each population:

| Population | Distribution | LS parameters | Weighted CDF SSE | Weighted CDF max discrepancy | Pseudo-MLE parameters | Support valid? | Selected? |
|---|---|---|---:|---:|---|---|---|
| Adult | Normal | | | | | | |
| Adult | Lognormal | | | | | | |
| Adult | Gamma | | | | | | |
| Adult | Triangular | | | | | | |
| Child | Normal | | | | | | |
| Child | Lognormal | | | | | | |
| Child | Gamma | | | | | | |
| Child | Triangular | | | | | | |

## 12.4 Diagnostic figures

For adult and child separately:

```text
weighted_histogram_pdf.png
weighted_ecdf_cdf.png
weighted_qq_normal.png
weighted_qq_lognormal.png
weighted_qq_gamma.png
weighted_qq_triangular.png
```

---

# 13. Optional design-aware parameter uncertainty

The point estimates above use survey weights to represent the target population. If time permits, estimate uncertainty in the fitted parameters using a **cluster/stratum-aware resampling procedure**, using STEPS `psu`/`stratum` and the corresponding MICS cluster/stratum variables.

Do not bootstrap individual rows independently and call it design-correct if the original survey used clustered sampling.

This uncertainty analysis is optional for the CSE 402 BW-fitting stage unless explicitly required, but retaining the design variables now prevents having to reconstruct them later.

---

# 14. Final end-to-end execution sequence

1. Load `bgd2018.csv`.
2. Extract `m12`, `wstep2`, `age`, `sex`, `m11`, `pid`, `psu`, `stratum` plus optional geography.
3. Clean adult BW using STEPS special-code rules and preserve legitimate extremes.
4. Load MICS7 `ch.sav`, **not `bh.sav`**.
5. Extract `AN8`, `chweight`, `CAGE`, `HL4`, `HH7A`, `AN11`, `WAZFLAG`, supporting anthropometry flags, and survey-design variables.
6. Restrict child BW to 0–59 months.
7. Remove MICS AN8 special codes `99.3`, `99.4`, `99.5`, `99.6` and missing values.
8. Require positive `chweight` and retain `WAZFLAG == 0` for the primary child BW fit.
9. Construct the survey-weighted ECDF separately for adult and child BW.
10. Fit **Normal, Lognormal, Gamma, and Triangular** with the same survey-weighted CDF least-squares objective.
11. Fit all four again using survey-weighted pseudo-MLE as a robustness benchmark.
12. Calculate weighted CDF SSE, weighted maximum CDF discrepancy, physical-support checks, ECDF overlays, Q-Q plots, and tail diagnostics.
13. Select the best-supported adult and child distributions separately.
14. Save the fitted parameters and all diagnostics.
15. During Monte Carlo simulation, generate $U\sim U(0,1)$ and use the selected distribution's `.ppf()` to generate `BW` values in kg.
16. Insert the generated `BW_i` values directly into the existing full ADD/HQ/HI/ELCR equations.

---

# 15. Source basis

## Project sources

1. **CSE 402 project slide v3**: requires least-squares distribution fitting within the numerical-methods workflow and Monte Carlo propagation of fitted input distributions.
2. **Yadav, S. & Kalkal, S. (2024)**, *Health risk assessment using Monte-Carlo simulations due to arsenic contamination in groundwater in Punjab*, *Journal of Water and Health*, 22(12), 2304–2319. DOI: 10.2166/wh.2024.188. The paper uses stochastic BW, with Lognormal adult BW and Triangular child BW.
3. **Bangladesh_Arsenic_Risk_Model_Parameters_Full_Design_PATCHED_v3_EF_FIXED.md**: establishes Bangladesh STEPS/MICS BW as the project BW source and inverse-transform Monte Carlo sampling from the selected fitted distribution.
4. **dataset_audit_report.md**: provides the exact audited STEPS/MICS filenames, shapes, columns, missingness, MICS AN8 special-code labels, child age range, survey weights, and anthropometric flags used in this plan.

## External methodological validation

5. **WHO NCD Microdata Repository, Bangladesh STEPS 2018**: confirms `m12 = weight (kg)`, `wstep2 = final analysis weight for Step 2 physical measures`, survey weighting purpose, and the 8,185-case public-use dataset.
6. **WHO STEPS questionnaire documentation**: identifies `888` as refused and `666` as too large for the weight scale.
7. **UNICEF MICS variable conventions / MICS microdata documentation**: supports `AN8` as measured child weight, `chweight` as the under-5 sample weight, and `WAZFLAG` coding (`0 = no error`, `1 = error flag`).
8. **Binder (1983) / survey pseudo-likelihood methodology**, as summarized in modern complex-survey statistical literature: survey-weighted pseudo-log-likelihood has the form $\sum_i w_i\log f(y_i\mid\theta)$.
9. **Joo, Y. & Casella, G. (2001)**, *Predictive distributions in risk analysis and estimation for the triangular distribution*, *Environmetrics*, 12(7), 647–658. DOI: **10.1002/env.489**. Establishes triangular-distribution estimation in Monte Carlo risk analysis and reports quantile least squares preferable to MLE for the triangular case.
10. **US Army Corps of Engineers distribution-fitting documentation**: provides the standard triangular PDF/CDF/quantile parameterization used here.

---

# 16. Methodological statements to use in the final report

A concise defensible description is:

> Adult body weight was obtained from measured weight (`m12`) in Bangladesh STEPS 2018 and weighted using the Step-2 physical-measurement analysis weight (`wstep2`). Under-5 child body weight was obtained from measured weight (`AN8`) in the Bangladesh MICS7 child anthropometry file (`ch.sav`) and weighted using `chweight`. Non-measurement codes and survey-defined anthropometric error records were removed according to the source coding. Normal, Lognormal, Gamma, and Triangular distributions were fitted using a survey-weighted CDF least-squares minimum-distance criterion to satisfy the project's least-squares numerical-method requirement. Survey-weighted pseudo-maximum-likelihood fits were used as robustness checks. The final adult and child distributions were selected separately using weighted CDF error, graphical diagnostics, tail behaviour, physical support, and robustness-fit agreement. Monte Carlo body weights were then generated in kilograms by inverse-CDF sampling from the selected fitted distributions.

