# FusioMRdev

> **F**lexible, **U**nified and ver**S**atile Mendel**I**an Rand**O**mization framework — **dev**elopment version.

An R package for FusioMR, a flexible Bayesian hierarchical framework for single- and multi-outcome Mendelian randomization,
designed for molecular trait exposures and applicable to complex traits, with robust performance under limited instruments.

> **Note**: This is the development version of FusioMR. For the stable
> release, see [kangbw702/FusioMR](https://github.com/kangbw702/FusioMR).

## Installation

Requires **R >= 4.3.0** and a working C++ compiler:

- **macOS**: `xcode-select --install`
- **Windows**: install [Rtools](https://cran.r-project.org/bin/windows/Rtools/)

```r
# install.packages("devtools")
devtools::install_github("fsh56/FusioMR-software")
library(FusioMRdev)
```

### Dependencies

`FusioMRdev` automatically pulls in the following packages during installation:

- `Rcpp`, `RcppArmadillo` — for the Gibbs samplers
- `invgamma` — for inverse-gamma sampling

If `devtools::install_github()` fails to fetch them, install manually:

```r
install.packages(c("Rcpp", "RcppArmadillo", "invgamma"))
```

## Choosing a Model

`FusioMRdev` supports four models via the `model` argument of `fusiomr()`.
Pick the one that matches your data and concerns about pleiotropy:

| Model            | Exposure | Outcome | Use when                                              |
|------------------|----------|---------|-------------------------------------------------------|
| `seso_uhp_only`  | 1        | 1       | Most IVs valid; uncorrelated pleiotropy is the concern |
| `seso_with_chp`  | 1        | 1       | Some IVs may share genetic effects with the outcome    |
| `semo`           | 1        | 2       | One exposure, two related outcomes (joint analysis)    |
| `memo`           | 2        | 2       | Two exposures and two outcomes; full joint model       |

Input format mirrors that of `TwoSampleMR` and `MendelianRandomization` —
you only need summary statistics and a `model` name.

## Quick Start

## Main Function

```r
fusiomr(b_exp, se_exp, b_out, se_out,
        model   = "seso_uhp_only",
        control = parameter_control(),
        verbose = FALSE)
```

| Argument   | Description                                                                                |
|------------|--------------------------------------------------------------------------------------------|
| `b_exp`    | SNP-exposure effects. Vector for one exposure; K x 2 matrix for `memo`.                    |
| `se_exp`   | Standard errors of `b_exp`. Same shape as `b_exp`.                                         |
| `b_out`    | SNP-outcome effects. Vector for one outcome; K x 2 matrix for `semo` and `memo`.           |
| `se_out`   | Standard errors of `b_out`. Same shape as `b_out`.                                         |
| `model`    | One of `"seso_uhp_only"`, `"seso_with_chp"`, `"semo"`, `"memo"`.                            |
| `control`  | Advanced settings from `parameter_control()`. Defaults are tuned for typical MR settings.  |
| `verbose`  | If `TRUE`, print progress messages and a results summary.                                  |

> **Input note:** `FusioMRdev` assumes you have already performed upstream
> IV selection (e.g. LD clumping, p-value filtering). All input SNPs are
> treated as instrumental variables.

## Advanced Usage

For most users, default settings are sufficient. To tune MCMC length,
IV selection, or empirical-Bayes priors, pass a customized
`parameter_control()`:

```r
fit <- fusiomr(b_exp, se_exp, b_out, se_out,
               model = "seso_uhp_only",
               control = parameter_control(
                 niter = 30000,                  # longer MCMC
                 z_thresh = qnorm(1 - 5e-8 / 2),    # winner's-curse correction
                 rho_ov = 0.2                     # sample overlap
               ))
```

See `?parameter_control` for the full list of advanced options.

### Model 1: `seso_uhp_only`

Single exposure, single outcome. The default model when most IVs are
expected to be valid.

```r
library(FusioMRdev)
# Toy summary statistics: 80 SNPs, true causal effect = 0.3
set.seed(42)
K <- 50
b_exp  <- rnorm(K, 0, 0.1);  
se_exp <- rep(0.03, K)
b_out  <- 0.3 * b_exp + rnorm(K, 0, 0.05); 
se_out <- rep(0.03, K)
model1 <- fusiomr(b_exp, se_exp, b_out, se_out,
               model = "seso_uhp_only")
model1$est    # causal effect estimate
model1$se     # standard error estimate
model1$ci     # 95% empirical credible interval
model1$pval   # two-sided p-value
```

### Model 2: `seso_with_chp`

Use this model when you suspect some IVs exhibit **correlated horizontal
pleiotropy** (their outcome effect is partially mediated by the exposure).
The model jointly estimates the causal effect and the proportion of
pleiotropic IVs.

```r
library(FusioMRdev)
set.seed(42)
K <- 50
b_exp  <- rnorm(K, 0, 0.1) + rnorm(K, 0, 0.05)
b_out  <- 0.3 * b_exp + rnorm(K, 0, 0.05) + rnorm(K, 0, 0.03)
se_exp <- rep(0.03, K)
se_out <- rep(0.03, K)
model2 <- fusiomr(b_exp, se_exp, b_out, se_out,
               model = "seso_with_chp")
model2$est    # causal effect estimate
model2$se     # standard error estimate
model2$ci     # 95% empirical credible interval
model2$pval   # two-sided p-value
```


### Model 3: `semo`

Single exposure, two related outcomes. Useful when two outcomes share
the same exposure (e.g. one biomarker affecting two diseases) and you
want a joint estimate that pools information across outcomes.

```r
# Toy data: shared exposure, two outcomes with different causal effects
library(FusioMRdev)

set.seed(42)
K <- 30

b_exp  <- rnorm(K, 0, 0.1)
se_exp <- rep(0.03, K)

b_out <- cbind(
  0.3 * b_exp + rnorm(K, 0, 0.04),
  0.5 * b_exp + rnorm(K, 0, 0.04)
)

se_out <- matrix(0.04, K, 2)

model3 <- fusiomr(
  b_exp, se_exp,
  b_out, se_out,
  model = "semo"
)

model3$est    # length-2 vector of causal effect estiamte
model3$se     # length-2 vector of standard error estimate
model3$ci     # 2 x 2 matrix; row j is CI for outcome j
model3$pval   # length-2 vector of two-sided p-values
```

### Model 4: `memo`

Two exposures, two outcomes — the full joint model with both UHP and CHP.
Each (exposure j, outcome j) pair gets its own causal effect, but
pleiotropy correlations are modeled jointly across the four cells.

```r
```

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
[GitHub Issues](https://github.com/fsh56/FusioMR-software/issues)
or contact the authors at
[kbw@uchicago.edu](mailto:kbw@uchicago.edu),
[sfeng56@uchicago.edu](mailto:sfeng56@uchicago.edu).

## License

MIT
