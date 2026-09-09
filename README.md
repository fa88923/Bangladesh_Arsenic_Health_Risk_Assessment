# Team Onboarding and Work Plan

## Project Goal

Apply the Yadav and Kalkal Punjab arsenic-risk model using Bangladesh arsenic and body-weight data. Fit distributions, perform Monte Carlo risk calculations, check convergence, and report uncertainty and sensitivity. Raw data must never be overwritten.

## Repository Map

- `data/raw/`: immutable BGS, STEPS, and MICS source files.
- `scripts/`: reproducible analysis code; currently contains the dataset-audit script.
- `reports/`: generated audit reports and supporting tables.
- `docs/`: approved methodology, decisions, and this onboarding guide.
- `results/`: planned location for processed data, fitted parameters, simulations, tables, and figures.

## Reading Order

Read before contributing:

1. `docs/project_implementation_plan.md` - phases, deliverables, hold points, and definition of done.
2. `docs/Bangladesh_Arsenic_Risk_Model_Parameters_Full_Design.md` - Punjab equations, Bangladesh inputs, Monte Carlo design, and convergence rules.
3. `docs/arsenic_preprocessing_and_distribution_plan.md` - required reading for Members 1, 2, 4, and 5.
4. `docs/bodyweight_distribution_implementation_plan.md` - required reading for Members 3, 4, and 5.

Check the Git diff before starting. Existing changes belong to their authors; do not overwrite them.

## Five-Member Allocation

| Owner | Phases in `project_implementation_plan.md` | Primary responsibility | Required deliverables |
|---|---|---|---|
| Member 1: Data and arsenic preprocessing | Phases 1-2 | Freeze provenance; parse detected/censored arsenic; convert units; review duplicates; filter approved districts. | Preprocessed arsenic file, reconciliation audit, district summaries, data-contract tests. |
| Member 2: Arsenic fitting | Phase 4 | Build censor-aware empirical CDFs; fit Normal, Lognormal, and Gamma; run censored-MLE and bootstrap checks; prepare selection evidence. | Fit table, bootstrap table, diagnostics, provisional district choices. |
| Member 3: Body-weight track | Phases 3 and 5 | Clean weighted STEPS and MICS extracts; fit Normal, Lognormal, Gamma, and Triangular candidates separately for adults and children. | Cleaning audits, weighted summaries, fit table, diagnostic plots, provisional BW choices. |
| Member 4: Risk and simulation | Phases 6-9 | Preserve Punjab equations; implement independent inverse-CDF sampling; run simulations, convergence checks, and sensitivity analysis. | Equation tests, run manifests, risk summaries, convergence and sensitivity outputs. |
| Member 5: QA, integration, and reporting | All-phase QA; Phase 10 | Maintain contracts; verify units, counts, seeds, equations, and saved parameters; assemble final outputs. | Integration tests, benchmarks, reproducibility checklist, final tables/figures, traceability review. |

## Efficient Execution

Member 1 starts Phase 1, then Members 1 and 3 run Phases 2 and 3 in parallel. Member 2 begins Phase 4 after Phase 2; Member 3 continues to Phase 5 after Phase 3. Member 4 can prepare Phase 6 immediately, but Phase 7 waits for approved Phase 4 and 5 fits. Phases 8 and 9 follow the simulation. Member 5 reviews every phase and completes Phase 10 after integration.

Use one focused branch per workstream and small pull requests. Every handoff must include machine-readable output, row-count reconciliation, assumptions, tests, and a reviewer. Members 1 and 2 cross-review; Members 3 and 4 cross-review; Member 5 performs the final integration review. Decisions listed as hold points require team approval before downstream work proceeds.
