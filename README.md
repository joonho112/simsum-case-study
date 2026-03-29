# Replicating a Simulation Study

A mini course on replicating a published Monte Carlo simulation study, with a complete case study based on Pustejovsky & Tipton (2022).

**Author:** JoonHo Lee ([jlee296@ua.edu](mailto:jlee296@ua.edu))

**Read online:** <https://joonho112.github.io/simsum-case-study/>

## Overview

This Quarto Book is the second course in the SimSum series, following [Performance Metrics for Simulation Studies](https://joonho112.github.io/simsum-mini-course/). It guides graduate students through a hands-on replication of the simulation study from:

> Pustejovsky, J. E., & Tipton, E. (2022). Meta-analysis with robust variance estimation: Expanding the range of working models. *Prevention Science*, 23, 425–438.

The course covers:

- The ADEMP framework applied to a real published simulation
- Data-generating processes for meta-analytic data with dependent effect sizes
- Five working models (HE, CE, CHE, SCE, CMVE) with full mathematical derivations
- Line-by-line code walkthroughs of `generate_smd()` and `generate_meta()`
- Estimation via `metafor`, `clubSandwich`, and `robumeta` with CR2 robust inference
- Performance metrics using `simhelpers` (primary) and `rsimsum` (secondary)
- Verification of key findings against the published results
- Visualization gallery and reporting best practices

## Getting Started

### Prerequisites

- [Quarto](https://quarto.org/) (>= 1.4)
- R (>= 4.1) with packages: `simhelpers`, `rsimsum`, `metafor`, `clubSandwich`, `robumeta`, `ggplot2`, `dplyr`, `tidyr`, `purrr`, `knitr`, `ggridges`, `mvtnorm`

### Build the Book

```bash
quarto render
```

The rendered book will be in `_book/`.

### Data Files

The `data/` directory contains the TSL15 empirical dataset and pre-computed simulation results from the [OSF repository](https://osf.io/x8yre/) for Pustejovsky & Tipton (2022). Run `source("R/download-data.R")` to verify all required files are present.

## Book Structure

| Part | Chapter | Title |
|------|---------|-------|
| I | 1 | Introduction: Learning by Replication |
| | 2 | The ADEMP Blueprint |
| | 3 | The Empirical Foundation: TSL15 |
| II | 4 | Meta-Analytic Models and Working Models |
| | 5 | Generating Effect Sizes: `generate_smd()` |
| | 6 | Generating Meta-Analytic Datasets: `generate_meta()` |
| III | 7 | Implementing the Working Models |
| | 8 | From Single Iteration to Full Pipeline |
| IV | 9 | The Factorial Design |
| | 10 | Running the Simulation |
| V | 11 | Performance Metrics with simhelpers and rsimsum |
| | 12 | Verifying the Replication |
| | 13 | Visualizing Simulation Results |
| | 14 | Reporting the Replication |
| App | A | The Real-Data Analysis |
| | B | Mathematical Details |
| | C | Software Comparison Reference |

## Key References

- Pustejovsky, J. E., & Tipton, E. (2022). Meta-analysis with robust variance estimation: Expanding the range of working models. *Prevention Science*, 23, 425–438.
- Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation studies to evaluate statistical methods. *Statistics in Medicine*, 38(11), 2074–2102.
- Miratrix, L. W., & Pustejovsky, J. E. (2024). *Designing Monte Carlo Simulations in R*. Online textbook.
- Tanner-Smith, E. E., & Lipsey, M. W. (2015). Brief alcohol interventions for adolescents and young adults. *Journal of Substance Abuse Treatment*, 51, 1–18.

## SimSum Series

1. [Performance Metrics for Simulation Studies](https://joonho112.github.io/simsum-mini-course/) — Theory, MCSEs, ADEMP, visualization
2. **Replicating a Simulation Study** (this course) — Hands-on replication with simhelpers and rsimsum

## License

This work is licensed under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).
