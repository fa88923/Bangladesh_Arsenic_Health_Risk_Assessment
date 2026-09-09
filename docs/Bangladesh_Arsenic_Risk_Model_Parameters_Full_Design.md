# Bangladesh Arsenic Health-Risk Model: Agreed Data Sources, Parameters, and Formulas

**Approved project stance:** apply the Punjab model of Yadav and Kalkal (2024) to Bangladesh. Retain the Punjab model structure and equations while replacing the planned location-specific inputs with the approved Bangladesh arsenic and body-weight data and documenting every retained assumption.

## 1. Arsenic concentration, $C$

**Data source:** DPHE/BGS National Hydrochemical Survey of Bangladesh groundwater.

**Dataset to use:** the well-level CSV containing individual borehole/tubewell arsenic measurements with district and geographic information. Do not use arsenic-threshold maps, smoothed arsenic maps, or district-mean tables as the main input.

For each district:

$$
C_d = \{C_1,C_2,\ldots,C_n\}
$$

Fit and compare:
- Lognormal distribution;
- Gamma distribution;
- Normal distribution.

Select the best-supported fitted distribution for each district, estimate its parameters, and generate Monte Carlo arsenic samples from that selected fitted distribution.

**Unit:** $C$ in mg/L.

---

## 2. Drinking-water ingestion rate, $IR$

For consistency with the base paper, retain the probabilistic ingestion-rate assumption exactly as reported by Yadav & Kalkal (2024).

- Adults: **Lognormal**, reported as $1.26\pm0.66\ \text{L/day}$
- Children: **Lognormal**, reported as $1.26\pm0.66\ \text{L/day}$

No Bangladesh-specific fitting is performed for $IR$. During Monte Carlo simulation, $IR$ is sampled parametrically from the reported Lognormal distribution.

For implementation, generate

$$
U_i\sim\operatorname{Uniform}(0,1)
$$

and sample

$$
IR_i=F_{IR}^{-1}(U_i).
$$

**Approved interpretation:** treat `1.26` and `0.66` directly as the Lognormal log-space parameters $\mu_{\log}$ and $\sigma_{\log}$, respectively, for adults and children. They are not arithmetic-scale mean and standard deviation and must not be converted as such. Under the SciPy parameterization, use `s=0.66`, `loc=0`, and `scale=exp(1.26)`.

---

## 3. Exposure duration, $ED$

For the main non-carcinogenic model, keep $ED$ fixed unless a justified Bangladesh-specific exposure-duration model is introduced.

**Unit:** years.

For non-carcinogenic risk:

$$
AT=ED\times365
$$

where $AT$ is averaging time in days. The ADD equations are retained in their full base-paper form without algebraic simplification.

---

## 4. Body weight, $BW$

### Adults

**Preferred dataset:** WHO Bangladesh STEPS 2018 microdata.

Use measured body weight, measured height, age, sex, and the physical-measurement survey weight.

**Unit:** kg.

Primary approach:
1. construct the appropriately weighted BW dataset;
2. inspect histogram, ECDF, skewness, and Q-Q plots;
3. fit Normal, Lognormal, and Gamma as candidate distributions;
4. select the best-supported fitted distribution;
5. use the estimated parameters of that selected distribution for Monte Carlo sampling.

### Children

**Preferred dataset:** Bangladesh MICS 2019.

Useful variables include child weight, height/length, age in months, sex, district, and survey sampling weight.

**Unit:** kg.

The child BW population is approved as ages 0-59 months, matching the available MICS measurements. The planned child exposure duration remains 6 years. Because the BW source does not directly cover ages 60-71 months, this mismatch remains an explicit model factor and must be carried into the sensitivity discussion and final limitations rather than hidden or described as measured 0-6-year coverage.

---

## 5. Exposed skin surface area, $SA$

For consistency with the base paper, retain the paper's **exposed skin surface area** assumptions rather than deriving total body surface area from height and weight.

### Deterministic values

- Adult: $SA=1.8\ \text{m}^2$
- Child: $SA=0.6\ \text{m}^2$

### Probabilistic values reported in the base paper

- Adult: **Lognormal**, reported as $1.42\pm0.31\ \text{m}^2$
- Child: **Triangular**, reported as $0.6800\pm0.600\ \text{m}^2$

**Approved distribution family:** use a Triangular distribution for child exposed skin surface area. The paper does not provide an explicit minimum-mode-maximum triplet, so the family is fixed but its three required parameters remain a hold point. The reported notation must not be converted into a triplet by inventing a missing value. $SA$ refers to **exposed skin surface area**, not total body surface area.

For Monte Carlo simulation, adult $SA$ is sampled from the reported Lognormal distribution by treating `1.42` and `0.31` directly as $\mu_{\log}$ and $\sigma_{\log}$, respectively. Under the SciPy parameterization, use `s=0.31`, `loc=0`, and `scale=exp(1.42)`. Child $SA$ is sampled from the approved Triangular family after its minimum-mode-maximum triplet is established.

---

## 6. Average Daily Dose equations

### Ingestion pathway

$$
ADD_{\text{ing}}
=
\frac{
C\times IR\times EF\times ED
}{
BW\times AT
}
$$

**Output unit:**

$$
ADD_{\text{ing}}:\ \text{mg kg}^{-1}\text{day}^{-1}
$$

### Dermal pathway

$$
ADD_{\text{dermal}}
=
\frac{
C\times K_p\times EF\times ED\times ET\times SA\times CF
}{
BW\times AT
}
$$

**Output unit:**

$$
ADD_{\text{dermal}}:\ \text{mg kg}^{-1}\text{day}^{-1}
$$

---

## 7. Parameter definitions, values, and units

| Parameter | Meaning | Selected value/source | Unit |
|---|---|---|---|
| $C$ | Arsenic concentration | District-specific DPHE/BGS well data, have to find the distribution | mg/L |
| $IR$ | Drinking-water ingestion rate | Base paper: Lognormal, reported as $1.26\pm0.66$ for adults and children | L/day |
| $EF$ | Exposure frequency | Deterministic: 365; Probabilistic: Triangular with minimum 180, mode 345, maximum 365 | days/year |
| $ED$ | Exposure Duration | 70 years for adults, 6 years for children | years
| $BW$ | Body weight | Bangladesh STEPS/MICS, have to find the distribution | kg |
| $AT$ | Averaging time | $ED\times365$ for non-cancer risk | days |
| $ET$ | Dermal exposure time | Triangular, reported as $0.20\ (0.13-0.33)$ | h/day |
| $SA$ | Exposed skin surface area | Adult: Lognormal $1.42\pm0.31$; Child: Triangular $0.6800\pm0.600$, as reported in base paper | m² |
| $K_p$ | Derma7. Parameter definitions, values, and unitsl permeability coefficient | 0.001 | cm/h |
| $CF$ | Unit conversion factor | 10, as used in base paper | L·m/(m³·cm) |
| $RfD_{\text{ing}}$ | Oral reference dose | 0.0003 | mg kg⁻¹ day⁻¹ |
| $RfD_{\text{dermal}}$ | Dermal reference dose | 0.000285 | mg kg⁻¹ day⁻¹ |
| $CSF$ | Cancer slope factor | 1.5 | (mg kg⁻¹ day⁻¹)⁻¹ |

**Notes:**  
- The base paper reports $CF=10$ using the conversion-factor notation $L\cdot m/(m^3\cdot cm)$. This notation is retained here to stay consistent with the reproduced model.  
- The CSF unit is written here as inverse dose so that $ELCR$ is dimensionless; the base paper prints the CSF unit as a dose unit, which is dimensionally inconsistent with its own ELCR equation.

---

## 8. Risk equations

### Ingestion hazard quotient

$$
HQ_{\text{ing}}
=
\frac{ADD_{\text{ing}}}{RfD_{\text{ing}}}
$$

### Dermal hazard quotient

$$
HQ_{\text{dermal}}
=
\frac{ADD_{\text{dermal}}}{RfD_{\text{dermal}}}
$$

### Hazard index

$$
HI=HQ_{\text{ing}}+HQ_{\text{dermal}}
$$

$HQ$ and $HI$ are dimensionless.

### Excess lifetime carcinogenic risk

$$
ELCR=ADD_{\text{ing}}\times CSF
$$

$ELCR$ is dimensionless and interpreted as an excess lifetime risk probability.

---

## 9. Distribution-fitting strategy

For variables with raw Bangladesh data, fit parametric candidate distributions to the observed data and compare their goodness of fit.

For BW and district-level arsenic:
- Normal;
- Lognormal;
- Gamma;
- Triangular only as a comparison/fallback when justified by min-mode-max information.

Use least-squares CDF fitting as the central numerical method required by the project specification, with MLE as a robustness benchmark. Select the best-supported parametric model using AIC/AICc, Anderson-Darling, KS, Q-Q plots, ECDF comparison, and tail-fit behavior. Monte Carlo sampling is then performed from the selected fitted distribution using its estimated parameters.

### Python functions for fitting and sampling

Use the following functions explicitly in the implementation:

| Distribution | Least-squares CDF fitting | MLE robustness fit | Inverse-CDF sampling |
|---|---|---|---|
| Normal | `scipy.optimize.least_squares` with `scipy.stats.norm.cdf` | `scipy.stats.norm.fit(data)` | `scipy.stats.norm.ppf(u, loc=mu, scale=sigma)` |
| Lognormal | `scipy.optimize.least_squares` with `scipy.stats.lognorm.cdf` | `scipy.stats.lognorm.fit(data, floc=0)` | `scipy.stats.lognorm.ppf(u, s=sigma_log, scale=np.exp(mu_log))` |
| Gamma | `scipy.optimize.least_squares` with `scipy.stats.gamma.cdf` | `scipy.stats.gamma.fit(data, floc=0)` | `scipy.stats.gamma.ppf(u, a=shape, scale=scale)` |
| Triangular | `scipy.optimize.least_squares` with `scipy.stats.triang.cdf`, only when a defensible triangular parameterization exists | `scipy.stats.triang.fit(data)` only when fitting a triangular model is justified | `scipy.stats.triang.ppf(u, c=c, loc=minimum, scale=maximum-minimum)` |

For least-squares CDF fitting, SciPy does not provide a single built-in call equivalent to `fit(method="least_squares_cdf")`. Define the residuals between the candidate fitted CDF and the empirical CDF, then optimize the distribution parameters with `scipy.optimize.least_squares`.

For Monte Carlo inverse-transform sampling, generate uniform random probabilities with NumPy:

```python
rng = numpy.random.default_rng(seed)
u = rng.uniform(0.0, 1.0, size=N)
```

and transform them through the selected distribution's inverse CDF using SciPy's percent-point function (`ppf`):

```python
samples = scipy.stats.<distribution>.ppf(u, ...)
```

This implements:

$$
U_i\sim\operatorname{Uniform}(0,1),
$$

$$
X_i=F_X^{-1}(U_i).
$$

Direct `.rvs()` sampling may be used for diagnostic checks, but the primary project implementation should use `Generator.uniform` + `.ppf()` so that inverse-transform sampling is explicit.

---

# Extended Implementation Design: From Parameter Preparation to Monte Carlo Results Analysis

This extension preserves Sections 1–9 above exactly as agreed. It adds implementation notes, the district-wise Monte Carlo design, the result-generation plan, numerical validation, and sensitivity analysis required for the CSE 402 project.

The design follows three sources of methodological guidance:

1. the CSE 402 project specification in **402 project slide v3**;
2. the base paper, **Yadav & Kalkal (2024)**, especially its Monte Carlo, percentile, threshold-exceedance, and Spearman sensitivity outputs;
3. standard probabilistic risk-assessment guidance, especially **US EPA RAGS Volume III, Part A**, which describes Monte Carlo risk assessment as repeated random sampling of exposure-variable distributions followed by evaluation of the risk equation, producing a full risk distribution that is summarized by percentiles, exceedance probabilities, and sensitivity measures.

The primary implementation language should be **Python** so that the numerical methods are visible and reproducible rather than hidden inside a commercial Monte Carlo package.

---

## 10. Recommended software and tool notes for Sections 1–9

### 10.1 Core environment

Recommended environment:

- **Python 3.11+**
- **pandas**: CSV loading, cleaning, district grouping, result tables
- **NumPy**: vectorized arithmetic, random-number generation, quantiles
- **SciPy**: probability distributions, distribution fitting, statistical tests, Spearman correlation, KDE
- **Matplotlib**: histograms, ECDF/CDF plots, Q-Q plots, district-wise risk-distribution figures, tornado plots, convergence plots
- **statsmodels**: optional ECDF and proportion-confidence-interval utilities
- **PyArrow/Parquet**: optional but recommended for storing full Monte Carlo samples efficiently

Use one project-wide reproducible random-number interface:

```python
rng = numpy.random.default_rng(seed)
```

The seed must be recorded for every reported run.

### 10.2 Section 1: arsenic concentration, $C$

**Tools:** `pandas.read_csv`, `DataFrame.groupby`, `numpy`, `scipy.optimize.least_squares`, `scipy.stats.norm`, `scipy.stats.lognorm`, `scipy.stats.gamma`, `matplotlib`.

Implementation notes:

- Clean arsenic values and convert all observations to mg/L before fitting.
- Group observations by district.
- Fit Normal, Lognormal, and Gamma distributions separately for each district.
- For least-squares CDF fitting, use `scipy.optimize.least_squares` together with `scipy.stats.norm.cdf`, `scipy.stats.lognorm.cdf`, or `scipy.stats.gamma.cdf`.
- For the MLE robustness benchmark, use `scipy.stats.norm.fit(data)`, `scipy.stats.lognorm.fit(data, floc=0)`, and `scipy.stats.gamma.fit(data, floc=0)`.
- Plot histogram + fitted PDF, ECDF + fitted CDF, and Q-Q diagnostics.
- Select the best-supported fitted distribution and generate Monte Carlo arsenic values from that fitted model using `rng.uniform(...)` followed by the corresponding `scipy.stats.<distribution>.ppf(...)`.

For inverse-transform sampling:

$$
U_i\sim\operatorname{Uniform}(0,1),
$$

$$
C_{d,i}=F_{C,d}^{-1}(U_i;\hat\theta_d).
$$

### 10.3 Section 2: drinking-water ingestion rate, $IR$

**Tools:** `numpy.random.default_rng`, `Generator.uniform`, and `scipy.stats.lognorm.ppf`.

Use parametric Monte Carlo sampling from the base-paper Lognormal assumption reported as $1.26\pm0.66\ \text{L/day}$ for both adults and children.

Using inverse-transform notation:

$$
U_i\sim\operatorname{Uniform}(0,1),
$$

$$
IR_i=F_{IR}^{-1}(U_i).
$$

Use the approved direct Lognormal parameterization: $\mu_{\log}=1.26$ and $\sigma_{\log}=0.66$, corresponding to `s=0.66`, `loc=0`, and `scale=exp(1.26)` in SciPy. Do not apply arithmetic-moment conversion.

### 10.4 Exposure frequency, $EF$

**Tools:** `numpy.random.default_rng`, `Generator.uniform`, and `scipy.stats.triang.ppf`.

For the probabilistic model, retain the base-paper Triangular exposure-frequency distribution:

$$
EF\sim\operatorname{Triangular}(180,345,365)\ \text{days/year},
$$

where 180 is the minimum, 345 is the mode, and 365 is the maximum. The deterministic baseline uses $EF=365\ \text{days/year}$.

For SciPy, the triangular shape parameter is

$$
c=\frac{345-180}{365-180},
$$

with `loc=180` and `scale=185`.

### 10.5 Section 3: exposure duration, $ED$, and averaging time, $AT$

**Tools:** ordinary NumPy scalar/vector arithmetic.

In the agreed primary model, $ED$ is fixed and

$$
AT=365ED.
$$

Therefore no random-number generator is needed for $ED$ or $AT$ in the primary analysis.

### 10.6 Section 4: body weight, $BW$

**Tools:** `pandas`, `numpy`, `scipy.optimize.least_squares`, `scipy.stats.norm`, `scipy.stats.lognorm`, `scipy.stats.gamma`, `matplotlib`.

Use the survey weights when preparing/fitting the adult and child BW distributions. Fit the candidate distributions, select the best-supported model, and sample BW from that selected fitted distribution using its estimated parameters.

For least-squares CDF fitting, use `scipy.optimize.least_squares` with the relevant candidate CDF. For the MLE robustness benchmark, use `scipy.stats.norm.fit(data)`, `scipy.stats.lognorm.fit(data, floc=0)`, and `scipy.stats.gamma.fit(data, floc=0)`. For Monte Carlo sampling, generate `u = rng.uniform(...)` and transform it with the selected distribution's `.ppf()` function.

For inverse-transform sampling:

$$
U_i\sim\operatorname{Uniform}(0,1),
$$

$$
BW_i=F_{BW}^{-1}(U_i;\hat\theta_{BW}).
$$

**Child-model requirement:** the final child age range must be frozen before the final child BW distribution is fitted.

### 10.7 Section 5: exposed skin area, $SA$

**Tools:** `numpy.random.default_rng`, `Generator.uniform`, `scipy.stats.lognorm.ppf` for adult $SA$, and `scipy.stats.triang.ppf` only when a defensible $(\text{left},\text{mode},\text{right})$ triplet is available.

For adult SA, use the approved direct Lognormal parameterization $\mu_{\log}=1.42$ and $\sigma_{\log}=0.31$, corresponding to `s=0.31`, `loc=0`, and `scale=exp(1.42)` in SciPy. Do not apply arithmetic-moment conversion.

For child SA, the Triangular family is approved. **Do not invent triangular parameters** from the reported $0.6800\pm0.600$ notation. The final probabilistic child run requires a documented minimum, mode, and maximum from an authoritative interpretation or an explicit project decision.

### 10.8 Sections 6–8: ADD, HQ, HI, and ELCR calculations

**Tools:** NumPy vectorized arithmetic. Avoid row-by-row loops for the simulation core.

For arrays of $N$ sampled inputs, compute the equations element-wise and retain each iteration's inputs and outputs for later sensitivity analysis.

### 10.9 Section 9: distribution fitting

**Tools:**

- `scipy.optimize.least_squares` for least-squares CDF fitting required by the project specification;
- `scipy.stats.norm.cdf`, `scipy.stats.lognorm.cdf`, `scipy.stats.gamma.cdf` for candidate CDF evaluation;
- `scipy.stats.norm.fit(data)`, `scipy.stats.lognorm.fit(data, floc=0)`, `scipy.stats.gamma.fit(data, floc=0)` for MLE robustness checks;
- `scipy.stats.norm.ppf`, `scipy.stats.lognorm.ppf`, `scipy.stats.gamma.ppf` for inverse-CDF sampling from the selected fitted model;
- `scipy.stats.kstest` or direct KS-statistic calculation for CDF distance;
- `matplotlib` for Q-Q and ECDF overlays.

For ordered observations $x_{(1)},\ldots,x_{(n)}$, define the empirical plotting positions as, for example,

$$
\hat F_n(x_{(i)})=\frac{i-0.5}{n}.
$$

For candidate CDF $F(x;\theta)$, the least-squares fit is

$$
\hat\theta_{LS}
=
\arg\min_{\theta}
\sum_{i=1}^{n}
\left[
F(x_{(i)};\theta)-\hat F_n(x_{(i)})
\right]^2.
$$

The MLE benchmark is

$$
\hat\theta_{MLE}
=
\arg\max_{\theta}
\sum_{i=1}^{n}\log f(x_i;\theta).
$$

For each fitted model, report at minimum the fitted parameters, least-squares CDF error, log-likelihood, AIC/AICc, KS statistic, and graphical diagnostics.

For positive variables such as arsenic and body weight, do not silently allow a Normal fit to generate negative physical values. If a Normal model has non-negligible negative probability, it should not be selected as the physical sampling model unless a justified truncated-normal model is explicitly used.

---

## 11. District-wise Monte Carlo model definition

Let:

- $d$ denote a Bangladesh district;
- $g\in\{\text{adult},\text{child}\}$ denote population group;
- $i=1,\ldots,N$ denote Monte Carlo iteration.

For each district and group, define the stochastic input vector

$$
\mathbf X_{d,g,i}
=
\left(
C_{d,i},IR_{g,i},BW_{g,i},ET_i,SA_{g,i}
\right).
$$

In the agreed probabilistic primary model, $EF$ is stochastic, while $ED$ and $AT$ remain fixed:

$$
EF_i\sim \operatorname{Triangular}(180,345,365),
$$

$$
ED=\text{fixed},
\qquad
AT=365ED.
$$

The deterministic baseline reproduction uses $EF=365\ \text{days/year}$. The constants $K_p$, $CF$, $RfD_{ing}$, $RfD_{dermal}$, and $CSF$ remain fixed as specified in Section 7.

### 11.1 Risk equations used in simulation

Use the complete base-paper equations from Sections 6–8 directly in every Monte Carlo iteration:

$$
ADD_{\text{ing},i}
=
\frac{
C_i\times IR_i\times EF\times ED
}{
BW_i\times AT
},
$$

$$
ADD_{\text{dermal},i}
=
\frac{
C_i\times K_p\times EF\times ED\times ET_i\times SA_i\times CF
}{
BW_i\times AT
}.
$$

Then calculate $HQ_{\text{ing}}$, $HQ_{\text{dermal}}$, $HI$, and $ELCR$ using the equations in Section 8. Do not algebraically simplify the ADD equations in the implementation plan.

---

## 12. Monte Carlo sampling design

### 12.1 Primary run and mandatory convergence experiment

Use

$$
N=10{,}000
$$

for the **principal reported Monte Carlo run**, so that the main analysis matches the base paper's stated simulation size. However, the implementation must not rely on a single 10,000-iteration run without a convergence check.

For **every district × population group**, run the simulation at

$$
N\in\{1{,}000,\ 5{,}000,\ 10{,}000,\ 20{,}000\}.
$$

At each iteration count, calculate and store at minimum:

- mean HI;
- P95 HI;
- $P(HI>1)$.

Also retain P50 HI and P95 ELCR because they are useful supporting diagnostics and already form part of the planned result package.

The purpose is to check whether the quantities actually reported in the risk analysis stabilize as $N$ increases. The $N=10{,}000$ result remains the **base-paper-comparable primary result**. The $N=20{,}000$ run is an additional higher-iteration stability check, not an assumed reference truth.

For a reported statistic $T$, assess stabilization using the **successive relative change** between the current and immediately preceding iteration count:

$$
\Delta_T(N)
=
\frac{|T_N-T_{N_{prev}}|}{|T_N|}\times 100\%.
$$

For example, for the 95th percentile:

$$
\Delta_{P95}(N)
=
\frac{|P95_N-P95_{N_{prev}}|}{|P95_N|}\times100\%.
$$

This avoids treating the largest simulation as the exact answer. Convergence is judged from the progressive stabilization of the reported statistics as $N$ increases.

For the threshold probability

$$
\hat p_N=\widehat{P(HI>1)},
$$

its Monte Carlo standard error can additionally be estimated by

$$
SE(\hat p_N)\approx\sqrt{\frac{\hat p_N(1-\hat p_N)}{N}}.
$$

**Tools:** NumPy `Generator` for reproducible sampling; pandas for the convergence table; NumPy `mean` and `quantile` for mean/P95; Boolean counting for $P(HI>1)$; Matplotlib for plots of each statistic against $N$. Use explicitly recorded random seeds.

The base paper reports independent 10,000-iteration runs to assess reliability. The CSE 402 specification explicitly requires convergence/error analysis, so the multi-$N$ experiment above is a required part of this Bangladesh implementation rather than an optional post-processing step.

### 12.2 Arsenic sampling from the selected fitted distribution

For district $d$, fit the candidate distributions to the observed arsenic measurements and select the best-supported model:

$$
C_d\sim F_d(\hat\theta_d).
$$

For each Monte Carlo iteration, generate arsenic from that selected fitted distribution:

$$
U_i\sim\operatorname{Uniform}(0,1),
$$

$$
C_{d,i}=F_d^{-1}(U_i;\hat\theta_d).
$$

**Tools:** `scipy.stats` distribution functions/quantiles or equivalent NumPy/SciPy parametric random-variate methods.

### 12.3 Sampling the remaining stochastic variables

For each Monte Carlo iteration, generate:

$$
IR_i\sim F_{IR,g},
$$

$$
BW_i\sim F_{BW,g},
$$

$$
EF_i\sim \operatorname{Triangular}(180,345,365),
$$

$$
ET_i\sim \operatorname{Triangular}(0.13,0.20,0.33),
$$

and the agreed $SA_g$ model.

$IR$ is sampled from the base-paper Lognormal distribution. $BW$ is sampled from the selected fitted BW distribution. $EF$ follows the base-paper probabilistic Triangular distribution with minimum 180, mode 345, and maximum 365 days/year. $ET$ and $SA$ are sampled from their specified parametric distributions, subject to the unresolved child-$SA$ triangular parameterization noted in Section 5.

For the stochastic parametric inputs, use inverse-transform sampling in the form

$$
X_i=F_X^{-1}(U_i),\qquad U_i\sim\operatorname{Uniform}(0,1).
$$


### 12.4 Dependence assumption

If the inputs are independently sampled, the model assumes

$$
f( C,IR,BW,EF,ET,SA )
=
f_C(C)f_{IR}(IR)f_{BW}(BW)f_{EF}(EF)f_{ET}(ET)f_{SA}(SA).
$$

Independent sampling is approved for the primary model and must be written explicitly in the report. US EPA guidance notes that Monte Carlo inputs are often treated as independent when dependence information is unavailable, but known correlations should be modeled rather than ignored.

Because the selected Bangladesh sources do not provide a common individual-level dataset linking arsenic concentration, IR, BW, EF, ET, and SA, the primary implementation will maintain independence, with this limitation stated clearly.

---

## 13. Monte Carlo implementation algorithm

For every valid Bangladesh district $d$:

1. load and clean the district's observed arsenic measurements;
2. fit Normal, Lognormal, and Gamma candidates;
3. select the primary fitted arsenic model according to Section 9;
4. for adults and children separately:
   - generate $C$, $IR$, $BW$, $EF$, $ET$, and $SA$ from their selected/specified parametric distributions;
   - calculate and save ADD, HQ, HI, and ELCR for every iteration;
   - summarize the simulated distributions;
   - compute sensitivity coefficients for HI;
5. **repeat the full simulation at $N=1{,}000$, $5{,}000$, $10{,}000$, and $20{,}000$ for every district × population group; at each $N$ store mean HI, P95 HI, and $P(HI>1)$, then evaluate successive stabilization and Monte Carlo error using independent replicate seeds;**
6. generate district-wise tables and figures.

### 13.1 Vectorized pseudocode

```text
for district d:
    C_obs = cleaned arsenic values for d
    fitted_C_model = fit_and_select(C_obs)

    for group g in {adult, child}:
        C[N]  = sample fitted_C_model
        IR[N] = sample base-paper Lognormal IR model
        BW[N] = sample selected fitted BW model
        EF[N] = sample Triangular(180, 345, 365)
        ET[N] = sample specified Triangular ET model
        SA[N] = sample agreed SA model

        ADD_ing[N]    = full ingestion equation
        ADD_dermal[N] = full dermal equation
        HQ_ing[N]     = ADD_ing / RfD_ing
        HQ_dermal[N]  = ADD_dermal / RfD_dermal
        HI[N]         = HQ_ing + HQ_dermal
        ELCR[N]       = ADD_ing * CSF

        save inputs + outputs
        summarize outputs
        calculate Spearman sensitivity of HI
```

**Tools:** NumPy for the simulation core, pandas for summaries and storage.

---

## 14. Primary Monte Carlo outputs: reproduce the base-paper result structure

The Bangladesh analysis should produce the **same classes of result reported in the base paper**, district-wise for Bangladesh using the selected fitted input distributions.

### 14.1 Base-paper-style probabilistic results table

For each district, report:

| District | P95 HI adult | P95 ELCR adult | P95 HI child | P95 ELCR child |
|---|---:|---:|---:|---:|

This directly mirrors the structure of Table 5 in Yadav & Kalkal (2024).

The empirical $p$-quantile is estimated from the simulated sample:

$$
\hat Q_p
=
\operatorname{Quantile}(Y_1,\ldots,Y_N;p).
$$

For the base-paper-matched result,

$$
p=0.95.
$$

**Tool:** `numpy.quantile(samples, 0.95)`.

Generate one table using the selected fitted arsenic distribution for each district.

### 14.2 District-wise HI probability-distribution figures

The base paper presents one district figure with adult and child HI distributions for each district. Reproduce the same output class for Bangladesh.

For each district, produce a two-panel figure:

- panel (a): adult HI distribution;
- panel (b): child HI distribution.

Each panel should contain:

- normalized histogram of simulated HI;
- optional KDE curve for visualization only;
- vertical line at $HI=1$;
- vertical line at the simulated 95th percentile;
- legend identifying the simulated HI distribution.

Do not allow the KDE visualization to replace the Monte Carlo estimates used in the result tables.

**Tools:** `matplotlib.pyplot.hist`, `scipy.stats.gaussian_kde`, `matplotlib.pyplot.axvline`.

### 14.3 Threshold-exceedance probabilities

For each district and group, report at minimum

$$
P(HI>1).
$$

Because the base paper also discusses the percentage exceeding $HI=2$, also report

$$
P(HI>2)
$$

as a secondary result.

The standard Monte Carlo estimator is direct exceedance counting:

$$
\widehat P(HI>h_0)
=
\frac{1}{N}
\sum_{i=1}^{N}
\mathbf 1(HI_i>h_0).
$$

Therefore,

$$
\widehat P(HI>1)
=
\frac{\#\{HI_i>1\}}{N}.
$$

**Tool:**

```python
p_hi_gt_1 = numpy.mean(HI > 1.0)
p_hi_gt_2 = numpy.mean(HI > 2.0)
```

Report percentages as $100\hat p\%$.

### 14.4 Additional distribution summaries

Although the base paper emphasizes P95, retaining additional summaries makes the analysis auditable. For every district/group/model, save:

- mean;
- standard deviation;
- median/P50;
- P5;
- P25;
- P75;
- P90;
- P95;
- P99;
- minimum and maximum;
- $P(HI>1)$;
- $P(HI>2)$.

These additional summaries should not replace the base-paper-matched P95 table.

**Tools:** `numpy.mean`, `numpy.std`, `numpy.quantile`.

---

## 15. Optional deterministic benchmark matching the base-paper comparison structure

The base paper also reports a deterministic table before the probabilistic results. A Bangladesh deterministic benchmark is recommended so that the effect of using full distributions can be demonstrated directly.

For each district, use the district mean arsenic concentration and clearly documented point values for the other inputs. Then report:

| District | HI adult | HI child | ELCR adult | ELCR child |
|---|---:|---:|---:|---:|

This mirrors the structure of Table 4 in the base paper.

The point values must come from the Bangladesh parameter design, not be copied from the Punjab study. For variables represented by Bangladesh distributions, use a pre-declared central value such as the weighted mean; do not choose the central value after seeing the risk results.

**Tools:** pandas/NumPy arithmetic.

---

## 16. Sensitivity analysis: base-paper-matched output

### 16.1 Method

The base paper performs sensitivity analysis for **non-carcinogenic risk, HI**, using the **Spearman rank-order correlation coefficient**, displayed as tornado plots for adults and children by district.

For each stochastic input $X_j$, calculate

$$
\rho_{s,j}
=
\operatorname{Corr}
\left(
\operatorname{rank}(X_j),
\operatorname{rank}(HI)
\right).
$$

For data without ties, Spearman's coefficient can also be written as

$$
\rho_s
=
1-
\frac{6\sum_i d_i^2}{N(N^2-1)},
$$

but ties can occur, so the implementation should use a standard tie-aware statistical routine rather than the simplified formula.

**Tool:** `scipy.stats.spearmanr`.

### 16.2 Sensitivity variables for the Bangladesh model

In the agreed Bangladesh primary model, the variables that should enter the HI sensitivity analysis are the variables that are genuinely stochastic:

- $C$: arsenic concentration;
- $IR$: ingestion rate;
- $BW$: body weight;
- $EF$: exposure frequency;
- $ET$: dermal exposure time;
- $SA$: exposed skin area, when stochastic.

Do **not** calculate Spearman sensitivity for fixed $ED$, $K_p$, $CF$, or the RfDs. A constant has zero variance and therefore no meaningful rank correlation with HI.

This is an important clarification relative to the base paper. Its text states that ED was the least-sensitive factor, yet its Table 2 labels ED as a point input, and its tornado-plot images do not show variable labels. The exact sensitivity-variable list in the published paper is therefore internally ambiguous. The Bangladesh implementation should not reproduce this inconsistency.

### 16.3 Tornado-plot output

For every district, create a two-panel tornado plot:

- adult HI sensitivity;
- child HI sensitivity.

Sort bars by absolute magnitude $|\rho_s|$. Preserve the sign:

- positive $\rho_s$: the input tends to increase HI as it increases;
- negative $\rho_s$: the input tends to decrease HI as it increases.

**Tool:** `matplotlib.pyplot.barh`.

The base-paper-matched deliverable is the signed Spearman coefficient, not only $\rho_s^2$.

### 16.4 Sensitivity table

In addition to figures, save:

| District | Group | Arsenic $\rho_s$ | IR $\rho_s$ | BW $\rho_s$ | EF $\rho_s$ | ET $\rho_s$ | SA $\rho_s$ |
|---|---|---:|---:|---:|---:|---:|---:|

This makes the plotted results reproducible.

---

## 17. Convergence and Monte Carlo error analysis

### 17.1 Iteration-count experiment

Run the complete simulation at, at minimum,

$$
N\in\{1{,}000,5{,}000,10{,}000,20{,}000\}.
$$

For each $N$, track:

- mean HI;
- P50 HI;
- P95 HI;
- $P(HI>1)$;
- P95 ELCR.

Use the same model definition at every $N$. The minimum convergence indicators are **mean HI, P95 HI, and $P(HI>1)$** because they represent central tendency, upper-tail risk, and threshold exceedance, respectively.

### 17.2 Independent replicate seeds

Because one Monte Carlo realization can appear stable by chance, run several independent simulations with different random seeds at each $N$. A practical project choice is **five replicate seeds per $N$**, but the number five is not a universal HHRA standard. The methodological requirement is to use independent replicate runs so that apparent convergence is not an artifact of a single random sequence. The base paper states that independent 10,000-iteration runs were performed but does not specify the number of replicates.

For statistic $T$, report across-seed mean and spread:

$$
\bar T_N
=
\frac{1}{R}\sum_{r=1}^{R}T_{N,r}.
$$

Assess numerical stabilization by comparing each iteration level with the immediately preceding level rather than treating the largest run as the exact answer:

$$
\Delta_T(N)
=
\frac{|\bar T_N-\bar T_{N_{prev}}|}{|\bar T_N|}\times100\%.
$$

For the upper-tail statistic specifically:

$$
\Delta_{P95}(N)
=
\frac{|P95_N-P95_{N_{prev}}|}{|P95_N|}\times100\%.
$$

Convergence is supported when the tracked statistics and their between-seed spread progressively stabilize as $N$ increases. No iteration count is treated as an exact reference truth.

### 17.3 Monte Carlo standard error for threshold probability

For the Monte Carlo exceedance probability

$$
\hat p=\frac{k}{N},
$$

the approximate Monte Carlo standard error is

$$
SE(\hat p)
\approx
\sqrt{
\frac{\hat p(1-\hat p)}{N}
}.
$$

This gives a direct numerical precision measure for $P(HI>1)$.

An optional 95% confidence interval for the binomial exceedance proportion can be calculated with a Wilson interval.

**Tool:** `statsmodels.stats.proportion.proportion_confint(..., method='wilson')`.

### 17.4 Convergence plots

For each representative district or for all districts in compact form, plot each tracked statistic against $N$. P95 deserves particular attention because upper-tail estimates usually stabilize more slowly than central statistics.

**Tool:** Matplotlib line plots.

---

## 18. Required Bangladesh result package

The final project should produce the following outputs.

### 18.1 Data and fitting deliverables

1. cleaned district-level arsenic dataset;
2. district sample-count table;
3. descriptive arsenic-statistics table;
4. distribution-fit table for Normal/Lognormal/Gamma by district;
5. histogram + fitted PDF plots;
6. ECDF + fitted CDF plots;
7. Q-Q plots;
8. selected fitted distribution and parameters for each district.

### 18.2 Base-paper-matched risk deliverables

1. **Deterministic district table**: HI adult, HI child, ELCR adult, ELCR child;
2. **Probabilistic P95 table**: P95 HI and P95 ELCR for adults and children;
3. **district-wise HI distribution figures**, adult and child panels;
4. **threshold-exceedance table**, including $P(HI>1)$ and $P(HI>2)$;
5. **district-wise Spearman tornado plots** for HI, adults and children;
6. **numeric Spearman sensitivity table**.

### 18.3 Project-extension deliverables

1. convergence/error table across $N$;
2. convergence plots;
3. reproducibility table containing model version, seed, iteration count, date, and selected distribution parameters.

---

## 19. Suggested result-table schema

### 19.1 Main probabilistic risk table

| District | P95 HI Adult | P95 ELCR Adult | P95 HI Child | P95 ELCR Child |
|---|---:|---:|---:|---:|

### 19.2 Exceedance table

| District | Group | P(HI>1) | P(HI>2) | MC SE for P(HI>1) |
|---|---|---:|---:|---:|

### 19.3 Sensitivity table

| District | Group | $\rho_s(C,HI)$ | $\rho_s(IR,HI)$ | $\rho_s(BW,HI)$ | $\rho_s(EF,HI)$ | $\rho_s(ET,HI)$ | $\rho_s(SA,HI)$ |
|---|---|---:|---:|---:|---:|---:|---:|

### 19.4 Convergence table

| District | Group | N | Seed/replicate | Mean HI | P50 HI | P95 HI | P(HI>1) | P95 ELCR |
|---|---|---:|---|---:|---:|---:|---:|---:|

---

## 20. Recommended code organization

```text
project/
├── data/
│   ├── raw/
│   └── processed/
├── src/
│   ├── preprocess.py
│   ├── fit_distributions.py
│   ├── parameters.py
│   ├── risk_equations.py
│   ├── monte_carlo.py
│   ├── sensitivity.py
│   └── plotting.py
├── results/
│   ├── tables/
│   ├── figures/
│   └── simulations/
├── notebooks/
│   └── analysis.ipynb
├── tests/
│   ├── test_risk_equations.py
│   ├── test_sampling.py
│   └── test_reproducibility.py
└── README.md
```

For full iteration-level outputs, prefer Parquet over CSV because district-wise simulations can produce millions of rows.

---

## 21. Minimum implementation tests before reporting results

### 21.1 Unit and algebra tests

Verify that:

1. all concentrations are in mg/L;
2. all BW values are positive;
3. all sampled IR, ET, and SA values are positive;
4. $HI=HQ_{ing}+HQ_{dermal}$ exactly within floating-point tolerance;
5. $ELCR=ADD_{ing}\times CSF$;
6. repeated execution with the same seed gives identical results.

**Tools:** Python `pytest`, `numpy.testing.assert_allclose`.

### 21.2 Distribution-sampling tests

For a large diagnostic sample from each parametric input, compare the simulated mean/SD or quantiles with the theoretical distribution values. This catches parameterization errors, especially confusion between arithmetic and log-space parameters for Lognormal distributions.

### 21.3 Physical-support checks

Reject a final model configuration if it produces impossible values such as negative arsenic concentration, negative BW, or negative ingestion rate because of an unsuitable sampling distribution.

---

## 22. Interpretation rules for the final report

1. **HI is dimensionless.** The principal non-carcinogenic threshold is $HI=1$.
2. Report the **full probability distribution**, not only a single mean.
3. Use **P95 HI and P95 ELCR** for direct comparability with the base paper.
4. Calculate $P(HI>1)$ directly from the Monte Carlo output as the proportion of simulated HI values exceeding 1.
5. Interpret Spearman coefficients as sensitivity/association measures, not causal effects.
6. Do not include fixed parameters in rank-correlation sensitivity analysis.
7. Explicitly document the independence assumption and any unresolved parameter limitations.
8. Do not conceal districts with small well counts; report $n_d$ and discuss tail uncertainty where data are sparse.

---

## 23. Methodological basis and references

### Project and base-study sources

- **CSE 402 project specification:** *Probabilistic Arsenic Health-Risk Assessment for Bangladesh Groundwater Using Numerical Methods*, 402 project slide v3. The numerical stages used in this implementation are district-wise distribution fitting, Monte Carlo simulation from the selected fitted distributions, and convergence/error analysis.
- **Yadav, S. & Kalkal, S. (2024).** *Health risk assessment using Monte-Carlo simulations due to arsenic contamination in groundwater in Punjab.* Journal of Water and Health, 22(12), 2304–2319. DOI: 10.2166/wh.2024.188. The paper uses 10,000-iteration Monte Carlo simulations, reports district-wise P95 HI and cancer risk for adults and children, discusses percentages above HI thresholds, and uses Spearman rank-order sensitivity displayed as tornado plots.

### Standard probabilistic-risk methodology

- **US EPA (2001).** *Risk Assessment Guidance for Superfund, Volume III, Part A: Process for Conducting Probabilistic Risk Assessment.* EPA describes Monte Carlo risk assessment as repeated random sampling from exposure-variable probability distributions, calculation of a risk value for each sampled input set, and summarization of the resulting risk distribution using statistics such as percentiles and graphical PDF/CDF outputs. EPA also recommends evaluation of simulation stability with different iteration counts and discusses Spearman rank correlation for sensitivity analysis.
- **Burmaster, D.E. & Anderson, P.D. (1994).** *Principles of Good Practice for the Use of Monte Carlo Techniques in Human Health and Ecological Risk Assessments.* Risk Analysis, 14(4), 477–481. DOI: 10.1111/j.1539-6924.1994.tb00265.x.
- **Binkowitz, B.S. & Wartenberg, D. (2001).** *Disparity in Quantitative Risk Assessment: A Review of Input Distributions.* Risk Analysis, 21(1), 75–90. DOI: 10.1111/0272-4332.211091. This paper describes Monte Carlo methods as propagating individual exposure-input distributions into an overall risk distribution.

### Computational tools

- **NumPy random Generator documentation:** standard reproducible sampling for Normal, Lognormal, Gamma, and Triangular distributions.
- **SciPy statistics documentation:** `scipy.stats` distributions, `spearmanr`, `gaussian_kde`, KS tests, and distribution fitting.

---

## 24. Final implementation sequence

The recommended order of work is:

1. freeze the adult and child population definitions;
2. resolve the minimum-mode-maximum triplet for the approved child SA Triangular distribution;
3. preprocess district arsenic data;
4. fit and select district arsenic distributions;
5. prepare and fit adult and child BW distributions;
6. implement and unit-test the full base-paper ADD/HQ/HI/ELCR equations;
7. implement Monte Carlo sampling from the selected/specified parametric distributions;
8. reproduce a deterministic Bangladesh benchmark;
9. **run the mandatory convergence experiment for every district × population group at**
   $$
   N=1{,}000,\ 5{,}000,\ 10{,}000,\ 20{,}000,
   $$
   **tracking at minimum mean HI, P95 HI, and $P(HI>1)$ at every $N$;**
10. quantify stabilization using successive relative changes between iteration levels, independent replicate seeds, and, for $P(HI>1)$, Monte Carlo standard error;
11. designate the $N=10{,}000$ result as the primary base-paper-comparable simulation after confirming acceptable stabilization;
12. generate base-paper-matched P95 risk tables and HI-distribution figures;
13. calculate $P(HI>1)$ and $P(HI>2)$ directly from the simulated HI values by exceedance counting;
14. perform district-wise Spearman HI sensitivity analysis and tornado plots;
15. export the convergence tables/plots plus all other tables, figures, configuration values, distribution parameters, and random seeds needed to reproduce the reported results.

This sequence keeps the health-risk methodology aligned with the base paper while retaining the numerical-analysis components used in the CSE 402 implementation.
