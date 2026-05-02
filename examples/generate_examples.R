# ============================================================================
#  FusioMRdev: example data generation
# ============================================================================
#  This script generates the example .rds files used in the README/tutorial.
#  Run it once locally to (re)generate the .rds files under examples/data/.
#  The .rds files themselves are committed to the GitHub repo for users.
#
#  Settings strictly follow the FusioMR paper simulation grids
#  (param_sim1.csv, param_sim2.csv, param_sim4.csv):
#    - seso_uhp_only: param_sim1 grid point (b_alpha = 0.10, theta = 0.3,
#                     p_cutoff = 1e-3, q_chp = 0)
#    - seso_with_chp: param_sim2 grid point (q_chp = 0.10, b_phi = 0.10)
#    - semo:          param_sim4 grid point (theta1 = theta2 = 0.3,
#                     rho_alpha = 0.4)
#    - memo:          TODO (data-generating function not yet finalized)
# ============================================================================

# ---- 0. setup --------------------------------------------------------------

# Working directory should be examples/. Adjust if running from elsewhere.
setwd("~/FusioMR-software/examples")

library(Rcpp)

# Compile fastSigLm once. This is independent of the FusioMRdev package.
Rcpp::sourceCpp("dgm/fastlm.cpp")

# Source the data-generating functions.
source("dgm/dgm4.R")   # provides dgm4(...) for seso
source("dgm/dgm5.R")   # provides dgm5(...) for semo

# Output directory.
dir.create("data", showWarnings = FALSE)


# ---- 1. seso_uhp_only ------------------------------------------------------
#  Paper Scenario 1A (param_sim1): limited IVs (n_x = 300), no CHP.
#  Grid point: b_alpha = 0.10, theta = 0.3, p_cutoff = 1e-3.
# ----------------------------------------------------------------------------

set.seed(2026101)

params_seso_uhp <- list(
  m       = 200,    # number of independent SNPs
  nx      = 300,    # exposure GWAS sample size (limited IV scenario)
  ny      = 20000,  # outcome GWAS sample size
  a_gamma = -0.3,   # IV-to-exposure effect: Unif(-b_gamma, b_gamma)
  b_gamma = 0.3,
  a_f     = 0.1,    # MAF range
  b_f     = 0.3,
  a_alpha = -0.10,  # UHP effect: Unif(-b_alpha, b_alpha) with b_alpha = 0.10
  b_alpha = 0.10,
  a_phi   = -0.10,  # CHP effect range (unused when q_chp = 0)
  b_phi   = 0.10,
  theta   = 0.3,    # true causal effect (paper grid: c(0, 0.1, 0.2, 0.3, 0.4))
  q_uhp   = 1,      # all SNPs have UHP effect (paper assumption)
  q_chp   = 0       # no CHP for the seso_uhp_only example
)

sim_seso_uhp <- do.call(dgm4, params_seso_uhp)

# Select IVs by exposure p-value < 1e-3 (matches paper).
z_exp <- abs(sim_seso_uhp$b_exp / sim_seso_uhp$se_exp)
p_exp <- 2 * pnorm(z_exp, lower.tail = FALSE)
sel   <- p_exp < 1e-3

seso_uhp_data <- list(
  b_exp     = sim_seso_uhp$b_exp[sel],
  se_exp    = sim_seso_uhp$se_exp[sel],
  b_out     = sim_seso_uhp$b_out[sel],
  se_out    = sim_seso_uhp$se_out[sel],
  true_beta = params_seso_uhp$theta,
  setting   = "Paper Scenario 1A (param_sim1): FusioMR_s with limited IVs (n_x = 300), no CHP",
  params    = params_seso_uhp,
  n_iv      = sum(sel),
  seed      = 2026101
)

saveRDS(seso_uhp_data, file = "data/seso_uhp_only_example.rds")
cat(sprintf("Saved seso_uhp_only_example.rds : %d IVs selected, true beta = %.2f\n",
            seso_uhp_data$n_iv, seso_uhp_data$true_beta))


# ---- 2. seso_with_chp ------------------------------------------------------
#  Paper Scenario 1A with CHP (param_sim2): same as above but with
#  q_chp = 0.10 and b_phi = 0.10 (10% of IVs have CHP effects).
# ----------------------------------------------------------------------------

set.seed(42)

params_seso_chp <- params_seso_uhp
params_seso_chp$q_chp <- 0.10   # 10% of IVs have CHP, paper grid c(0.05, 0.10)
params_seso_chp$a_phi <- -0.10
params_seso_chp$b_phi <- 0.10

sim_seso_chp <- do.call(dgm4, params_seso_chp)

z_exp2 <- abs(sim_seso_chp$b_exp / sim_seso_chp$se_exp)
p_exp2 <- 2 * pnorm(z_exp2, lower.tail = FALSE)
sel2   <- p_exp2 < 1e-3

seso_chp_data <- list(
  b_exp     = sim_seso_chp$b_exp[sel2],
  se_exp    = sim_seso_chp$se_exp[sel2],
  b_out     = sim_seso_chp$b_out[sel2],
  se_out    = sim_seso_chp$se_out[sel2],
  true_beta = params_seso_chp$theta,
  setting   = "Paper Scenario 1A (param_sim2): FusioMR_s with limited IVs (n_x = 300), 10% CHP",
  params    = params_seso_chp,
  n_iv      = sum(sel2),
  seed      = 42
)

saveRDS(seso_chp_data, file = "data/seso_with_chp_example.rds")
cat(sprintf("Saved seso_with_chp_example.rds : %d IVs selected, true beta = %.2f\n",
            seso_chp_data$n_iv, seso_chp_data$true_beta))


# ---- 3. semo (single exposure, two outcomes) -------------------------------
#  Paper Scenario 2 (param_sim4): shared molecular trait exposure,
#  two correlated outcomes, no CHP. theta1 = theta2 (matched), rho_alpha = 0.4.
# ----------------------------------------------------------------------------

set.seed(20240101)

params_semo <- list(
  m         = 200,
  nx        = 300,
  ny1       = 20000,
  ny2       = 20000,
  a_gamma   = -0.3,
  b_gamma   = 0.3,
  a_f       = 0.1,
  b_f       = 0.3,
  a_alpha1  = -0.10,
  b_alpha1  = 0.10,
  a_alpha2  = -0.10,
  b_alpha2  = 0.10,
  rho_theta = 0.4,    # UHP correlation, paper grid c(0, 0.4, 0.8)
  theta1    = 0.3,    # true causal effect on outcome 1
  theta2    = 0.3,    # matched theta (theta2 = theta1) per paper
  q_uhp1    = 1,
  q_uhp2    = 1
)

sim_semo <- do.call(dgm5, params_semo)

z_exp3 <- abs(sim_semo$b_exp / sim_semo$se_exp)
p_exp3 <- 2 * pnorm(z_exp3, lower.tail = FALSE)
sel3   <- p_exp3 < 1e-3

# semo expects b_out and se_out as K x 2 matrices.
b_out_mat  <- cbind(sim_semo$b_out_1[sel3],  sim_semo$b_out_2[sel3])
se_out_mat <- cbind(sim_semo$se_out_1[sel3], sim_semo$se_out_2[sel3])

semo_data <- list(
  b_exp     = sim_semo$b_exp[sel3],
  se_exp    = sim_semo$se_exp[sel3],
  b_out     = b_out_mat,
  se_out    = se_out_mat,
  true_beta = c(params_semo$theta1, params_semo$theta2),
  setting   = "Paper Scenario 2 (param_sim4): FusioMR_m with shared exposure, UHP rho = 0.4, matched theta = 0.3",
  params    = params_semo,
  n_iv      = sum(sel3),
  seed      = 20240101
)

saveRDS(semo_data, file = "data/semo_example.rds")
cat(sprintf("Saved semo_example.rds : %d IVs selected, true betas = (%.2f, %.2f)\n",
            semo_data$n_iv, semo_data$true_beta[1], semo_data$true_beta[2]))


# ---- 4. memo (TODO) --------------------------------------------------------
#  Paper Scenario 3 (param_sim5): two exposures, two outcomes (full joint
#  model with CHP). Data-generating function (Equation S6) not yet
#  implemented in this folder. To be added in a future commit.
# ----------------------------------------------------------------------------

cat("\nDone. .rds files written to examples/data/\n")