# Arsenic Distribution Fit Report

## District List

The arsenic fitting workflow was run for the following 10 selected districts:

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

## Method

For each district, the arsenic concentration distribution was estimated using the detected and left-censored observations from the preprocessed arsenic dataset.

Reverse Kaplan-Meier was applied separately within each district to estimate a censoring-aware empirical CDF. Least-squares CDF fitting was then used to fit three candidate distributions to that district-specific empirical CDF:

- Normal
- Lognormal
- Gamma

The main fitting metric was the CDF RMSE from the least-squares fit. Bootstrap support reports how often the distribution won the admissible LS-RMSE contest. The final simulation distribution is marked with a single tick.

The Phase 7 review accepted nine provisional selections and overrode one: Barisal changed from Lognormal to Gamma after tail checks showed that the Lognormal fitted P99 was 26 times the empirical P99, its fitted mean was 29 times the reverse-KM mean, and 7.5% of draws exceeded the largest measured concentration. Candidate-fit statistics below are retained unchanged; the selection marker reflects the final model actually used in simulation.

## Fit Results

| District | Distribution | Parameters | CDF RMSE | Bootstrap support | Final simulation selection |
|---|---|---|---:|---|:---:|
| Dhaka | Normal | $\mu=-0.0123$; $\sigma=0.1087$ | 0.0228 | LS win 0.0% (n=1000) |  |
| Dhaka | Lognormal | $\mu_{\log}=-5.9154$; $\sigma_{\log}=3.9033$ | 0.0621 | LS win 0.5% (n=1000) |  |
| Dhaka | Gamma | $k=0.1489$; $\theta=0.4425$ | 0.0463 | LS win 99.5% (n=1000) | &#10003; |
| Chittagong | Normal | $\mu=0.0064$; $\sigma=0.0095$ | 0.1073 | LS win 0.0% (n=994) |  |
| Chittagong | Lognormal | $\mu_{\log}=-5.4876$; $\sigma_{\log}=1.9829$ | 0.0462 | LS win 98.8% (n=994) | &#10003; |
| Chittagong | Gamma | $k=0.4564$; $\theta=0.0253$ | 0.0697 | LS win 1.2% (n=994) |  |
| Rajshahi | Normal | $\mu=-0.0296$; $\sigma=0.0420$ | 0.0193 | LS win 0.0% (n=1000) |  |
| Rajshahi | Lognormal | $\mu_{\log}=-9.7730$; $\sigma_{\log}=4.6229$ | 0.0264 | LS win 2.0% (n=1000) |  |
| Rajshahi | Gamma | $k=0.0682$; $\theta=0.1520$ | 0.0164 | LS win 98.0% (n=1000) | &#10003; |
| Khulna | Normal | $\mu=-0.0058$; $\sigma=0.0766$ | 0.0361 | LS win 0.0% (n=1000) |  |
| Khulna | Lognormal | $\mu_{\log}=-6.3854$; $\sigma_{\log}=3.9430$ | 0.0521 | LS win 6.5% (n=1000) |  |
| Khulna | Gamma | $k=0.1447$; $\theta=0.3278$ | 0.0357 | LS win 93.5% (n=1000) | &#10003; |
| Barisal | Normal | $\mu=0.0158$; $\sigma=0.1348$ | 0.0960 | LS win 0.0% (n=1000) |  |
| Barisal | Lognormal | $\mu_{\log}=-5.2479$; $\sigma_{\log}=3.5345$ | 0.0415 | LS win 84.4% (n=1000) |  |
| Barisal | Gamma | $k=0.1724$; $\theta=0.5818$ | 0.0447 | LS win 15.6% (n=1000) | &#10003; |
| Sylhet | Normal | $\mu=0.0084$; $\sigma=0.0320$ | 0.0419 | LS win 0.0% (n=1000) |  |
| Sylhet | Lognormal | $\mu_{\log}=-5.0215$; $\sigma_{\log}=2.0617$ | 0.0355 | LS win 5.8% (n=1000) |  |
| Sylhet | Gamma | $k=0.3498$; $\theta=0.0664$ | 0.0168 | LS win 94.2% (n=1000) | &#10003; |
| Rangpur | Normal | $\mu=-0.0013$; $\sigma=0.0114$ | 0.0344 | LS win 0.0% (n=1000) |  |
| Rangpur | Lognormal | $\mu_{\log}=-7.2105$; $\sigma_{\log}=2.6081$ | 0.0258 | LS win 15.7% (n=1000) |  |
| Rangpur | Gamma | $k=0.2056$; $\theta=0.0267$ | 0.0114 | LS win 84.3% (n=1000) | &#10003; |
| Mymensingh | Normal | $\mu=-0.0179$; $\sigma=0.0539$ | 0.0418 | LS win 0.0% (n=1000) |  |
| Mymensingh | Lognormal | $\mu_{\log}=-7.7219$; $\sigma_{\log}=3.9658$ | 0.0307 | LS win 44.2% (n=1000) |  |
| Mymensingh | Gamma | $k=0.1100$; $\theta=0.1913$ | 0.0230 | LS win 55.8% (n=1000) | &#10003; |
| Comilla | Normal | $\mu=0.1124$; $\sigma=0.1496$ | 0.0431 | LS win 0.0% (n=1000) |  |
| Comilla | Lognormal | $\mu_{\log}=-2.8795$; $\sigma_{\log}=2.3354$ | 0.0973 | LS win 0.0% (n=1000) |  |
| Comilla | Gamma | $k=0.4364$; $\theta=0.4201$ | 0.0664 | LS win 100.0% (n=1000) | &#10003; |
| Bogra | Normal | $\mu=-0.0036$; $\sigma=0.0238$ | 0.0754 | LS win 0.0% (n=1000) |  |
| Bogra | Lognormal | $\mu_{\log}=-7.0228$; $\sigma_{\log}=2.8619$ | 0.0322 | LS win 99.5% (n=1000) | &#10003; |
| Bogra | Gamma | $k=0.1677$; $\theta=0.0646$ | 0.0438 | LS win 0.5% (n=1000) |  |

## Final Arsenic Models Used in Simulation

| Final family | Districts |
|---|---|
| Gamma | Dhaka, Rajshahi, Khulna, Barisal, Sylhet, Rangpur, Mymensingh, Comilla |
| Lognormal | Chittagong, Bogra |

Rajshahi is retained with a high-censoring warning (71.8%). Bogra's heavy Lognormal upper tail and Rangpur's thin Gamma extreme tail remain explicit model-uncertainty items. Full selection evidence and the Barisal override are documented in D6-D10 of the simulation report.

## ECDF-CDF Diagnostic Plots

![District-wise ECDF-CDF subplot panel](../results/figures/arsenic_ecdf_cdf_panel.png)
