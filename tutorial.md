# FusioMRdev Tutorial

This tutorial walks through a basic Mendelian Randomization (MR) analysis with
**FusioMRdev** using GWAS summary statistics. We will use 
single-exposure single-outcome model (`seso_uhp_only`) as the example to explain the software usage
and how to set up advanced parameter settings.

---

## 1. Installation

Requires **R >= 4.3.0** and a working C++ compiler (Xcode CLT on macOS, Rtools
on Windows).

```r
# install.packages("devtools")
devtools::install_github("fsh56/FusioMR-software")
library(FusioMRdev)
```

---

## 2. Input data

`fusiomr()` takes four vectors of GWAS summary statistics:

| Argument | Meaning                              |
|----------|--------------------------------------|
| `b_exp`  | SNP–exposure effect estimates        |
| `se_exp` | Standard errors of `b_exp`           |
| `b_out`  | SNP–outcome effect estimates         |
| `se_out` | Standard errors of `b_out`           |

All four must have the same length, and the input data should already be preprocessed.
FusioMRdev does not perform data preprocessing.

---

## 3. Example1: `seso_uhp_only`

We already generated the example dataset. You can access them under directory ~examples/data/seso_uhp_only_example.rds

```r
library(FusioMRdev)

# load the example data shipped with the repo
dat <- readRDS("examples/data/seso_uhp_only_example.rds")

# run the model:
fit <- fusiomr(
  b_exp  = dat$b_exp,
  se_exp = dat$se_exp,
  b_out  = dat$b_out,
  se_out = dat$se_out,
  model  = "seso_uhp_only"
)
```

The returned object is a list with the estimation summary:

```r
fit$est    # causal effect estimate
fit$se     # standard error
fit$pval   # two-sided p-value
fit$ci     # 95% empirical credible interval
```
---

## 4. Example2: advanced parameter setting (`hybrid = TRUE`)

You are required to offer the estimate the global mean and explicitly assign TRUE value to hybrid
parameter in parameter_control(). 
```r
ctrl <- parameter_control(
  hybrid = TRUE,
  kappa_hybrid  = 5,           
  global_mean_gamma  = <your estiamte>
  global_mean_theta  = <your estiamte>
)

fit <- fusiomr(
  b_exp  = dat$b_exp,
  se_exp = dat$se_exp,
  b_out  = dat$b_out,
  se_out = dat$se_out,
  model   = "seso_uhp_only",
  control = ctrl
)
```
---

## 5. Getting help

```r
?fusiomr
?parameter_control
```
