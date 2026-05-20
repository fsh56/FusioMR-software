# FusioMR

> **F**lexible, **U**nified and ver**S**atile Mendel**I**an Rand**O**mization framework — **dev**elopment version.

`FusioMR` is the R implementation of FusioMR, a Bayesian hierarchical framework for single- and multi-outcome Mendelian randomization (MR) using GWAS summary statistics. 
It is designed primarily for molecular trait exposures (e.g., gene expression), where the number of available cis-QTLs as instrumental variables (IVs) is often limited and horizontal pleiotropy is pervasive. FusioMR is also applicable to complex trait exposures with a moderate-to-large number of IVs. For methodological details, please refer to https://doi.org/10.1016/j.ajhg.2026.03.017.

## Installation

`FusioMR` requires **R >= 4.3.0** and a working C++ compiler:

- **macOS**: `xcode-select --install`
- **Windows**: install [Rtools](https://cran.r-project.org/bin/windows/Rtools/)

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("kangbw702/FusioMR")
library(FusioMR)
```

### Dependencies

The following R packages are pulled in automatically as dependencies:

| Package          | Minimum version | Purpose                                  |
|------------------|-----------------|------------------------------------------|
| `Rcpp`           | >= 1.0.10       | R / C++ interface for the Gibbs samplers |
| `RcppArmadillo`  | >= 0.12.0.0     | Linear algebra in the Gibbs samplers     |
| `invgamma`       | >= 1.1          | Inverse-gamma sampling                   |

If `devtools::install_github()` fails to fetch them, install manually:

```r
install.packages(c("Rcpp", "RcppArmadillo", "invgamma"))
```

## Overview

The main entry point is `fusiomr()`, which takes four vectors of GWAS
summary statistics:

| Argument | Meaning                       |
|----------|-------------------------------|
| `b_exp`  | IV–exposure effect estimates  |
| `se_exp` | Standard errors of `b_exp`    |
| `b_out`  | IV–outcome effect estimates   |
| `se_out` | Standard errors of `b_out`    |

All four must have the same length, and the input data should already be
preprocessed (LD-clumped, IV-selected, harmonized). `FusioMR` does not
perform data preprocessing.

`FusioMR` supports four models via the `model` argument. Pick the one
that matches your data:

| Model            | Exposure | Outcome | Use when                                       |
|------------------|----------|---------|------------------------------------------------|
| `seso_uhp_only`  | 1        | 1       | uncorrelated pleiotropy (UHP) is the concern   |
| `seso_with_chp`  | 1        | 1       | correlated pleiotropy (CHP) is the concern     |
| `semo`           | 1        | 2       | one exposure, two outcomes                     |
| `memo`           | 2        | 2       | two exposures, two outcomes                    |

For `semo`, pass `b_out` / `se_out` as a `K x 2` matrix.
For `memo`, pass `b_exp` / `se_exp` **and** `b_out` / `se_out` as
`K x 2` matrices.

The returned object is a list with the MR estimates:

| Return | Meaning                          |
|--------|----------------------------------|
| `est`  | Causal effect estimate           |
| `se`   | Standard error                   |
| `pval` | Two-sided p-value                |
| `ci`   | 95% credible interval            |

## Quick Start

```r
library(FusioMR)

# Load an example dataset (single exposure, single outcome, UHP only)
d <- readRDS("examples/data/seso_uhp_only_example.rds")

fit <- fusiomr(d$b_exp, d$se_exp, d$b_out, d$se_out,
               model = "seso_uhp_only")

fit$est; fit$se; fit$pval; fit$ci
```

For full examples covering all four models, advanced parameter tuning, and
hybrid empirical-Bayes priors, see the
[tutorial vignette](vignettes/FusioMR-tutorial.Rmd).

## Simulation and Method Validation

Simulation scripts used in the FusioMR paper, including data-generating
functions and benchmarking pipelines, are available at
[kangbw702/FusioMR-analysis](https://github.com/kangbw702/FusioMR-analysis):

- [`dgf/`](https://github.com/kangbw702/FusioMR-analysis/tree/main/dgf) —
  data-generating functions for individual-level and summary-level GWAS
  with controllable pleiotropy.
- [`simulation/`](https://github.com/kangbw702/FusioMR-analysis/tree/main/simulation) —
  end-to-end simulation pipelines reproducing the paper results.

## Getting Help

```r
?fusiomr
?parameter_control
```

## Feedback

This is a pre-release version actively under development. Please report
bugs or suggest features via
[GitHub Issues](https://github.com/kangbw702/FusioMR/issues)
or contact the authors at
[kbw@uchicago.edu](mailto:kbw@uchicago.edu),
[sfeng56@uchicago.edu](mailto:sfeng56@uchicago.edu).

## License

MIT