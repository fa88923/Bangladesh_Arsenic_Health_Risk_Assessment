# Simulation at N = 10,000 Report

This short report extracts the primary N = 10,000 Monte Carlo results from `reports/simulation_report.md`.

## Probabilistic P95

| District | P95 HI adult | P95 ELCR adult | P95 HI child | P95 ELCR child | P95 ELCR child (paper AT) |
|---|---:|---:|---:|---:|---:|
| Dhaka | 21.6 | 9.68 x 10^-3 | 121.1 | 4.67 x 10^-3 | 5.45 x 10^-2 |
| Chittagong | 6.66 | 2.99 x 10^-3 | 34.7 | 1.34 x 10^-3 | 1.56 x 10^-2 |
| Rajshahi * | 3.17 | 1.42 x 10^-3 | 18.2 | 6.99 x 10^-4 | 8.16 x 10^-3 |
| Khulna | 15.7 | 7.03 x 10^-3 | 84.3 | 3.25 x 10^-3 | 3.79 x 10^-2 |
| Barisal | 32.1 | 1.44 x 10^-2 | 175.1 | 6.75 x 10^-3 | 7.87 x 10^-2 |
| Sylhet | 6.61 | 2.96 x 10^-3 | 35.8 | 1.38 x 10^-3 | 1.61 x 10^-2 |
| Rangpur | 1.84 | 8.27 x 10^-4 | 9.62 | 3.71 x 10^-4 | 4.33 x 10^-3 |
| Mymensingh | 7.19 | 3.23 x 10^-3 | 37.9 | 1.46 x 10^-3 | 1.70 x 10^-2 |
| Comilla | 51.3 | 2.30 x 10^-2 | 259.6 | 1.00 x 10^-2 | 1.17 x 10^-1 |
| Bogra | 6.57 | 2.95 x 10^-3 | 33.3 | 1.28 x 10^-3 | 1.50 x 10^-2 |

*Rajshahi is the high-censoring district, with 71.8% non-detects.

## Threshold Exceedance

MC SE for P(HI > 1) is at most 0.005, or 0.5 percentage points, for every run.

| District | P(HI > 1) adult | P(HI > 1) child | P(HI > 2) adult | P(HI > 2) child | P(ELCR > 10^-4) adult | P(ELCR > 10^-4) child |
|---|---:|---:|---:|---:|---:|---:|
| Dhaka | 33.5% | 47.7% | 26.8% | 42.2% | 46.0% | 39.9% |
| Chittagong | 23.2% | 53.8% | 14.3% | 40.3% | 50.6% | 35.3% |
| Rajshahi * | 11.2% | 20.3% | 7.1% | 16.4% | 19.2% | 15.0% |
| Khulna | 30.0% | 44.9% | 23.1% | 39.2% | 43.4% | 36.7% |
| Barisal | 41.2% | 55.8% | 33.8% | 49.5% | 54.7% | 47.4% |
| Sylhet | 32.8% | 60.0% | 20.6% | 50.4% | 57.5% | 46.1% |
| Rangpur | 10.0% | 29.9% | 4.4% | 20.9% | 27.6% | 18.0% |
| Mymensingh | 19.1% | 32.0% | 13.8% | 26.7% | 30.5% | 25.0% |
| Comilla | 72.3% | 85.1% | 62.5% | 80.1% | 85.6% | 78.1% |
| Bogra | 15.5% | 32.4% | 10.6% | 24.6% | 30.5% | 21.9% |

## Simulation Distribution Panels

![Simulation HI distribution panels at N = 10,000](../results/figures/simulation_hi_distribution_panel_N10000.png)

## Main Findings

1. Comilla has the highest simulated risk: adult P95 HI is 51.3 and child P95 HI is 259.6.
2. Barisal, Dhaka, and Khulna are the next highest-risk districts by P95 HI.
3. Children have higher HI than adults in every district, with child P95 HI about five to six times the adult value.
4. Rangpur and Rajshahi have the lowest P95 HI values, although Rajshahi should be read with the high-censoring caveat.
5. Threshold exceedance is substantial in high-risk districts: Comilla has P(HI > 1) of 72.3% for adults and 85.1% for children.
