# FusioMRdev

> **F**lexible, **U**nified and ver**S**atile Mendel**I**an Rand**O**mization framework — **dev**elopment version.

## Installation

Requires **R >= 4.3.0** and a working C++ compiler:
- **macOS**: `xcode-select --install`
- **Windows**: install [Rtools](https://cran.r-project.org/bin/windows/Rtools/)

```r
# install.packages("devtools")
devtools::install_github("fsh56/FusioMR-software")
library(FusioMRdev)
```

## Quick Start

For the fusiomr function input, you only need four summary-statistics vectors and a `model` name:

### Model1: seso_uhp_only
```r
library(FusioMRdev)

# Toy summary statistics: 80 SNPs, true causal effect = 0.3
set.seed(42)
K <- 80
b_exp  <- rnorm(K, 0, 0.1);  se_exp <- rep(0.01, K)
b_out  <- 0.3 * b_exp + rnorm(K, 0, 0.01); se_out <- rep(0.01, K)

model1 <- fusiomr(b_exp, se_exp, b_out, se_out,
               model = "seso_uhp_only")

model1$est    # causal effect estimate
model1$se     # standard error estimate
model1$ci     # 95% empirical credible interval
model1$pval   # two-sided p-value
```

### Model2: seso_with_chp
```r
library(FusioMRdev)

# Same input format as seso_uhp_only
model2 <- fusiomr(b_exp, se_exp, b_out, se_out,
               model = "seso_with_chp")

model2$est    # causal effect estimate
model2$se     # standard error estimate
model2$ci     # 95% empirical credible interval
model2$pval   # two-sided p-value
```

## Main Function

```r
fusiomr(b_exp, se_exp, b_out, se_out,
        model   = "seso_uhp_only",
        control = parameter_control(),
        verbose = FALSE)
```

**Arguments**

| Argument   | Description                                                   |
|------------|---------------------------------------------------------------|
| `b_exp`    | Numeric vector of SNP-exposure effect estimates               |
| `se_exp`   | Standard errors of `b_exp`                                    |
| `b_out`    | Numeric vector of SNP-outcome effect estimates                |
| `se_out`   | Standard errors of `b_out`                                    |
| `model`    | One of `"seso_uhp_only"`, `"seso_with_chp"`, `"semo"`, `"memo"` |
| `control`  | Advanced settings from `parameter_control()`                  |
| `verbose`  | If `TRUE`, print progress messages                            |

## Advanced Usage

For most users, default settings are sufficient. To tune MCMC length, IV
selection, or empirical-Bayes priors, pass a customized `parameter_control()`:

```r
model <- fusiomr(b_exp, se_exp, b_out, se_out,
               model = "seso_uhp_only",
               control = parameter_control(
                 niter             = 20000,             
                 p_value_threshold = 1e-5,              # stricter IV selection
                 z_thresh          = qnorm(1 - 5e-8/2), # winner's-curse correction
                 rho_ov            = 0.2                # account for sample overlap
               ))
```

See `?parameter_control` for all advanced options.

## Getting Help

```r
?fusiomr
?parameter_control
```

## Feedback

This is a pre-release version actively under development. Please report bugs or
suggest features via [GitHub Issues](https://github.com/fsh56/FusioMR-software/issues)
or contact the authors at `kbw@uchicago.edu`.

## License

MIT