# Arsenic Dataset Preprocessing and District-Selection Implementation Plan

## Project context

This plan is for preprocessing the **DPHE/BGS National Hydrochemical Survey of Bangladesh** arsenic dataset before district-wise distribution fitting and Monte Carlo health-risk simulation.

The final modeled variable is:

\[
C = \text{arsenic concentration in groundwater}
\]

with unit:

\[
\text{mg/L}
\]

The fitted Normal, Lognormal, or Gamma distribution for each selected district is therefore a distribution of **arsenic concentration itself**. Samples generated later from the selected fitted distribution are simulated arsenic concentrations in mg/L and are used directly as \(C\) in the ADD equations.

---

## 1. Freeze the selected 10 districts

Use the exact strings stored in the dataset column `DISTRICT`.

Recommended districts:

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

### Selection rationale

The selection is purposive rather than based only on sample count. It gives priority to:

- major nationally important districts;
- broad geographic coverage;
- divisional-center representation;
- adequate well-level observations for district-wise fitting;
- two additional major regional districts with comparatively strong sample counts.

Rajshahi has particularly high censoring and should therefore be retained with an explicit high-uncertainty warning rather than treated as equivalent to a lightly censored district.

---

## 2. Read the raw dataset without altering it

Source file:

```text
NationalSurveyData.csv
```

Important columns to retain:

```text
SAMPLE_ID
SAMPLE_FIELD_ID
SAMPLE_DATE
LAT_DEG
LONG_DEG
YEAR_CONSTRUCTION
WELL_TYPE
WELL_DEPTH
DIVISION
DISTRICT
THANA
UNION
MOUZA
GEOCODE
As
```

Recommended functions:

```python
pandas.read_csv()
DataFrame.copy()
DataFrame.shape
DataFrame.isna()
pandas.to_datetime()
```

### Rules

- Never overwrite the raw source file.
- Preserve the original `As` field exactly.
- Verify and remove only confirmed metadata/unit rows or other non-sample rows.
- Do not remove observations merely because they look extreme.

---

## 3. Preserve and parse the arsenic field

Keep the original arsenic column and create derived variables.

Recommended columns:

```text
As_raw
As_censored
As_bound_ugL
As_detected_ugL
As_model_value_ugL
As_bound_mgL
As_detected_mgL
As_model_value_mgL
```

### Example: censored observation

If:

```text
As_raw = "<6"
```

then:

```text
As_censored = True
As_bound_ugL = 6
As_detected_ugL = NaN
As_model_value_ugL = 6
```

The value `6` is stored only as the **upper censoring boundary**. It is not treated as the measured arsenic concentration.

### Example: detected observation

If:

```text
As_raw = "42"
```

then:

```text
As_censored = False
As_detected_ugL = 42
As_model_value_ugL = 42
```

Recommended functions:

```python
Series.astype("string")
Series.str.strip()
Series.str.startswith("<")
Series.str.extract(r"([0-9]*\.?[0-9]+)")
pandas.to_numeric(errors="coerce")
numpy.where()
```

---

## 4. Convert arsenic units

The BGS arsenic data are reported in \(\mu g/L\), while the health-risk equations require mg/L.

Convert using:

\[
C_{\mathrm{mg/L}} =
\frac{C_{\mu g/L}}{1000}
\]

Apply this conversion to:

- exact detected concentrations;
- censoring limits.

Example:

```python
df["As_detected_mgL"] = df["As_detected_ugL"] / 1000.0
df["As_bound_mgL"] = df["As_bound_ugL"] / 1000.0
df["As_model_value_mgL"] = df["As_model_value_ugL"] / 1000.0
```

All distribution fitting should then be performed in **mg/L** so that subsequent Monte Carlo samples are directly compatible with the risk equations.

---

## 5. Quality-control checks

Perform the following checks before filtering districts.

### Unique sample identifiers

```python
df["SAMPLE_ID"].duplicated()
```

A repeated `SAMPLE_ID` requires investigation.

### Duplicate full rows

```python
df.duplicated()
```

### Repeated coordinates

```python
df[["LAT_DEG", "LONG_DEG"]].duplicated(keep=False)
```

Do not automatically delete repeated coordinates. Instead inspect them with:

```text
SAMPLE_ID
SAMPLE_FIELD_ID
SAMPLE_DATE
WELL_DEPTH
WELL_TYPE
DISTRICT
```

Two different wells or samples may legitimately share coordinates.

### Numerical validity checks

Use:

```python
numpy.isfinite()
```

Check for:

- missing coordinates;
- missing district;
- missing arsenic field;
- negative exact arsenic values;
- malformed censoring strings.

Do not automatically remove high arsenic values. A high concentration may be a genuine hydrogeochemical observation.

---

## 6. Filter the final 10 districts

Use:

```python
SELECTED_DISTRICTS = [
    "Dhaka",
    "Chittagong",
    "Rajshahi",
    "Khulna",
    "Barisal",
    "Sylhet",
    "Rangpur",
    "Mymensingh",
    "Comilla",
    "Bogra"
]

selected = df[df["DISTRICT"].isin(SELECTED_DISTRICTS)].copy()
```

For each district report:

```text
total_n
detected_n
censored_n
censoring_percent
number_of_distinct_censoring_limits
minimum_detected_As
maximum_detected_As
```

Recommended functions:

```python
DataFrame.groupby()
DataFrame.agg()
Series.value_counts()
Series.nunique()
```

---

## 7. Do not substitute censored values

For observations such as:

```text
<6
<1
<0.5
```

do **not** replace them in the primary analysis with:

```text
0
LOD
LOD / 2
```

Such substitutions invent exact measurements and can distort the lower tail, fitted parameters, goodness-of-fit results, and subsequent Monte Carlo samples.

Instead, preserve the information:

\[
C < L
\]

where \(L\) is the reported detection limit.

---

## 8. Construct a left-censor-aware empirical arsenic CDF

The purpose of this step is to estimate the empirical distribution of arsenic while respecting nondetects.

The empirical distribution is still a distribution of:

\[
C = \text{arsenic concentration}
\]

Use a Reverse Kaplan-Meier / left-censored Kaplan-Meier implementation rather than assigning artificial concentrations.

Recommended Python package:

```python
from lifelines import KaplanMeierFitter
```

For each district:

```python
kmf = KaplanMeierFitter()

kmf.fit_left_censoring(
    durations=district_df["As_model_value_mgL"],
    event_observed=~district_df["As_censored"]
)
```

Obtain the cumulative distribution using:

```python
kmf.cumulative_density_
```

This produces points of the form:

\[
(x_j,\hat F_{\mathrm{RKM}}(x_j))
\]

where \(x_j\) is on the **arsenic concentration scale in mg/L**.

Important:

- Reverse Kaplan-Meier does not create fake arsenic values.
- It only estimates the empirical CDF using both detected values and censoring bounds.
- The final parametric distribution is still fitted to arsenic concentration.

---

## 9. Primary project fitting: Least-Squares CDF fitting

For every selected district, fit these candidates:

```text
Normal
Lognormal
Gamma
```

The primary fitting objective is:

\[
\hat\theta_{LS}
=
\arg\min_{\theta}
\sum_j
\left[
F(x_j;\theta)
-
\hat F_{\mathrm{RKM}}(x_j)
\right]^2
\]

where:

- \(x_j\) are arsenic concentration locations in mg/L;
- \(\hat F_{\mathrm{RKM}}\) is the censor-aware empirical arsenic CDF;
- \(F(x;\theta)\) is the candidate arsenic distribution CDF.

Recommended functions:

```python
scipy.optimize.least_squares
scipy.stats.norm.cdf
scipy.stats.lognorm.cdf
scipy.stats.gamma.cdf
```

For positive variables:

```python
floc = 0
```

should be used for Lognormal and Gamma where appropriate so that the physical support remains nonnegative.

For every fitted candidate save:

```text
estimated parameters
CDF SSE
CDF RMSE
RKM empirical CDF vs fitted CDF plot
```

### Interpretation

If the selected model is:

\[
C_{\text{Dhaka}}
\sim
\operatorname{Lognormal}(\hat\mu,\hat\sigma)
\]

then this is a **Lognormal distribution of Dhaka groundwater arsenic concentration**.

---

## 10. Censored MLE robustness fit

Use censored maximum likelihood as an independent robustness benchmark.

Recommended SciPy class:

```python
from scipy.stats import CensoredData
```

Construct censored data:

```python
cdata = CensoredData.left_censored(
    district_df["As_model_value_mgL"],
    district_df["As_censored"]
)
```

Fit candidates:

```python
scipy.stats.norm.fit(cdata)
scipy.stats.lognorm.fit(cdata, floc=0)
scipy.stats.gamma.fit(cdata, floc=0)
```

The censored likelihood uses:

\[
\log f(x_i)
\]

for detected observations and:

\[
\log F(L_i)
\]

for observations reported as:

\[
C < L_i
\]

Thus the model uses censored samples without pretending their exact concentrations are known.

---

## 11. Compare and select the district distribution

Use two parallel lines of evidence.

### Primary project criterion

Least-squares fit against the Reverse-KM arsenic CDF:

```text
CDF SSE
CDF RMSE
CDF overlay
```

### Robustness criterion

Censored-MLE:

```text
log-likelihood
AIC
AICc
```

Compute:

\[
AIC = 2k - 2\ell
\]

and:

\[
AIC_c
=
AIC
+
\frac{2k(k+1)}
{n-k-1}
\]

Useful functions:

```python
dist.logpdf()
dist.logcdf()
numpy.sum()
```

### Decision rule

Strongest case:

```text
LS-CDF winner = censored-MLE/AICc winner
```

If they disagree:

1. inspect the RKM-vs-fitted CDF plot;
2. inspect whether the disagreement is concentrated in the lower or upper tail;
3. examine bootstrap stability;
4. document the uncertainty rather than silently forcing agreement.

Do not apply ordinary uncensored KS or Anderson-Darling tests to a dataset in which nondetects have been replaced by artificial values.

---

## 12. Bootstrap stability analysis

Bootstrap analysis evaluates how sensitive the fitted distribution is to the particular wells present in the observed sample.

For each district:

1. resample the original district rows with replacement;
2. preserve each row's exact/censored status and its censoring limit;
3. reconstruct the Reverse-KM empirical CDF;
4. refit Normal, Lognormal, and Gamma;
5. record fitted parameters, fit error, and winning distribution;
6. repeat approximately 1,000 times.

Recommended function:

```python
district_df.sample(
    n=len(district_df),
    replace=True,
    random_state=seed
)
```

Record:

```text
parameter median
2.5th percentile
97.5th percentile
CDF RMSE distribution
frequency with which each candidate wins
```

This is particularly important for districts with high censoring, especially Rajshahi.

---

## 13. Save the processed arsenic dataset

Save a clean analytical file without changing the raw source.

Recommended filename:

```text
arsenic_selected_districts_preprocessed.csv
```

Recommended retained fields:

```text
SAMPLE_ID
SAMPLE_FIELD_ID
SAMPLE_DATE
LAT_DEG
LONG_DEG
WELL_TYPE
WELL_DEPTH
DIVISION
DISTRICT
THANA
UNION
MOUZA
As_raw
As_censored
As_bound_ugL
As_detected_ugL
As_model_value_ugL
As_bound_mgL
As_detected_mgL
As_model_value_mgL
```

Also save:

```text
arsenic_district_preprocessing_summary.csv
arsenic_distribution_fit_results.csv
arsenic_bootstrap_stability.csv
```

---

## 14. Final flow from raw data to district distributions

The complete preprocessing and fitting flow is:

\[
\boxed{
\text{Raw BGS well data}
\rightarrow
\text{preserve `As`}
\rightarrow
\text{parse exact vs <LOD}
\rightarrow
\text{convert to mg/L}
}
\]

\[
\boxed{
\rightarrow
\text{quality-control checks}
\rightarrow
\text{select 10 districts}
\rightarrow
\text{Reverse-KM arsenic empirical CDF}
}
\]

\[
\boxed{
\rightarrow
\text{LS fit Normal/Lognormal/Gamma}
\rightarrow
\text{censored-MLE robustness fit}
\rightarrow
\text{bootstrap stability}
\rightarrow
\text{select one arsenic distribution per district}
}
\]

---

## 15. Flow from fitted distribution to Monte Carlo arsenic samples

Suppose Dhaka finally selects:

\[
C_{\text{Dhaka}}
\sim
\operatorname{Lognormal}
(\hat\mu,\hat\sigma)
\]

Generate:

\[
U_i
\sim
\operatorname{Uniform}(0,1)
\]

and then:

\[
C_i
=
F_C^{-1}(U_i)
\]

Recommended implementation:

```python
rng = numpy.random.default_rng(seed)

u = rng.uniform(
    low=0.0,
    high=1.0,
    size=N
)

C = scipy.stats.lognorm.ppf(
    u,
    s=sigma_hat,
    scale=numpy.exp(mu_hat)
)
```

The resulting `C` values are:

\[
\boxed{\text{simulated arsenic concentrations in mg/L}}
\]

They are not probabilities, Kaplan-Meier values, or transformed quantities.

Each sampled \(C_i\) is used directly in:

\[
ADD_{\text{ing}}
=
\frac{
C\times IR\times EF\times ED
}{
BW\times AT
}
\]

and:

\[
ADD_{\text{dermal}}
=
\frac{
C\times K_p\times EF\times ED\times ET\times SA\times CF
}{
BW\times AT
}
\]

---

## 16. Required methodological notes in the report

The final project report should explicitly state:

1. censored arsenic observations were not replaced by arbitrary exact concentrations;
2. Reverse Kaplan-Meier / left-censored Kaplan-Meier was used to construct the empirical arsenic CDF;
3. least-squares CDF fitting was the primary numerical fitting method;
4. censored MLE was used as a robustness benchmark;
5. bootstrap resampling was used to examine fitted-parameter stability where needed;
6. the fitted distribution for each district represents **arsenic concentration \(C\)**;
7. Monte Carlo inverse-CDF sampling from the selected distribution produces arsenic concentrations in **mg/L**;
8. districts with high censoring, especially Rajshahi, carry greater fitted-distribution uncertainty;
9. the BGS survey should be interpreted as a distribution of sampled wells within each selected district, not necessarily a perfectly population-weighted distribution of all groundwater wells.

---

## 17. Recommended Python dependencies

```text
pandas
numpy
scipy
lifelines
matplotlib
```

Core functions/classes:

```python
pandas.read_csv
pandas.to_numeric
DataFrame.groupby
DataFrame.agg
DataFrame.sample
Series.str.startswith
Series.str.extract

numpy.where
numpy.isfinite
numpy.random.default_rng

lifelines.KaplanMeierFitter
KaplanMeierFitter.fit_left_censoring

scipy.optimize.least_squares
scipy.stats.CensoredData
scipy.stats.norm
scipy.stats.lognorm
scipy.stats.gamma
```
