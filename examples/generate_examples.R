# dgm
# ---- 0. setup --------------------------------------------------------------

# Working directory should be examples/. Adjust if running from elsewhere.
setwd("~/FusioMR-software/examples")

library(Rcpp)
# Compile fastSigLm once. This is independent of the FusioMRdev package.
Rcpp::sourceCpp("dgm/fastlm.cpp")

# Source the data-generating functions.
source("dgm/dgm4.R")
source("dgm/dgm5.R")
source("dgm/dgm6.R")

# Output directory.
dir.create("data", showWarnings = FALSE)

# ---- 1. seso_uhp_only ------------------------------------------------------
P_CUTOFF <- 1e-3

set.seed(13)
params_seso_uhp <- list(
  m       = 200,
  nx      = 300,
  ny      = 20000,
  a_gamma = -0.3,
  b_gamma = 0.3,
  a_f     = 0.1,
  b_f     = 0.3,
  a_alpha = -0.06,
  b_alpha = 0.06,
  a_phi   = -0.03,
  b_phi   = 0.03,
  theta   = 0,        
  q_uhp   = 1,
  q_chp   = 0.05
)

sim_seso_uhp <- do.call(dgm4, params_seso_uhp)
z_exp <- abs(sim_seso_uhp$b_exp / sim_seso_uhp$se_exp)
p_exp <- 2 * pnorm(z_exp, lower.tail = FALSE)
sel   <- p_exp < P_CUTOFF
seso_uhp_data <- list(
  b_exp     = sim_seso_uhp$b_exp[sel],
  se_exp    = sim_seso_uhp$se_exp[sel],
  b_out     = sim_seso_uhp$b_out[sel],
  se_out    = sim_seso_uhp$se_out[sel],
  true_beta = params_seso_uhp$theta,
  setting   = "Paper Scenario 1A (param_sim1): FusioMR_s with limited IVs (n_x = 300), no CHP, p_cutoff = 1e-3, true beta = 0.2",
  params    = params_seso_uhp,
  p_cutoff  = P_CUTOFF,
  n_iv      = sum(sel),
  seed      = 13
)
saveRDS(seso_uhp_data, file = "data/seso_uhp_only_example.rds")
cat(sprintf("Saved seso_uhp_only_example.rds : %d IVs selected, true beta = %.2f\n",
            seso_uhp_data$n_iv, seso_uhp_data$true_beta))

# ---- 2. seso_with_chp ------------------------------------------------------
set.seed(4)
P_CUTOFF_CHP <- 1e-5
params_seso_chp <- list(
  m       = 200,
  nx      = 20000,    # complex-trait setting: large sample
  ny      = 20000,
  a_gamma = -0.3,
  b_gamma = 0.3,
  a_f     = 0.1,
  b_f     = 0.3,
  a_alpha = -0.03,
  b_alpha = 0.03,
  a_phi   = -0.05,
  b_phi   = 0.05,
  theta   = 0,
  q_uhp   = 1,
  q_chp   = 0.5
)
sim_seso_chp <- do.call(dgm4, params_seso_chp)
z_exp2 <- abs(sim_seso_chp$b_exp / sim_seso_chp$se_exp)
p_exp2 <- 2 * pnorm(z_exp2, lower.tail = FALSE)
sel2   <- p_exp2 < P_CUTOFF_CHP
seso_chp_data <- list(
  b_exp     = sim_seso_chp$b_exp[sel2],
  se_exp    = sim_seso_chp$se_exp[sel2],
  b_out     = sim_seso_chp$b_out[sel2],
  se_out    = sim_seso_chp$se_out[sel2],
  true_beta = params_seso_chp$theta,
  setting   = "Complex-trait exposure (n_x = 20000)",
  params    = params_seso_chp,
  p_cutoff  = P_CUTOFF_CHP,
  n_iv      = sum(sel2),
  seed      = 4
)
saveRDS(seso_chp_data, file = "data/seso_with_chp_example.rds")
cat(sprintf("Saved seso_with_chp_example.rds : %d IVs selected, true beta = %.3f\n",
            seso_chp_data$n_iv, seso_chp_data$true_beta))



# ---- 3. semo (single exposure, two outcomes) -------------------------------
set.seed(1)
P_CUTOFF_SEMO <- 1e-3

params_semo <- list(
  m         = 200,
  nx        = 300,
  ny1       = 80000,
  ny2       = 20000,
  a_gamma   = -0.3,
  b_gamma   = 0.3,
  a_f       = 0.1,
  b_f       = 0.3,
  a_alpha1  = -0.1,
  b_alpha1  = 0.1,
  a_alpha2  = -0.1,
  b_alpha2  = 0.1,
  rho_theta = 0.8,
  theta1    = 0.2,
  theta2    = 0.2,
  q_uhp1    = 1,
  q_uhp2    = 1
)

sim_semo <- do.call(dgm5, params_semo)
z_exp3 <- abs(sim_semo$b_exp / sim_semo$se_exp)
p_exp3 <- 2 * pnorm(z_exp3, lower.tail = FALSE)
sel3   <- p_exp3 < P_CUTOFF_SEMO
b_out_mat  <- cbind(sim_semo$b_out_1[sel3],  sim_semo$b_out_2[sel3])
se_out_mat <- cbind(sim_semo$se_out_1[sel3], sim_semo$se_out_2[sel3])

semo_data <- list(
  b_exp     = sim_semo$b_exp[sel3],
  se_exp    = sim_semo$se_exp[sel3],
  b_out     = b_out_mat,
  se_out    = se_out_mat,
  true_beta = c(params_semo$theta1, params_semo$theta2),
  setting   = "shared exposure, two outcomes",
  params    = params_semo,
  p_cutoff  = P_CUTOFF_SEMO,
  n_iv      = sum(sel3),
  seed      = 1
)

saveRDS(semo_data, file = "data/semo_example.rds")
cat(sprintf("Saved semo_example.rds : %d IVs selected, true betas = (%.2f, %.2f)\n",
            semo_data$n_iv, semo_data$true_beta[1], semo_data$true_beta[2]))

# ---- 4. memo (two exposures, two outcomes) ---------------------------------
set.seed(4)
P_CUTOFF_MEMO <- 1e-5

params_memo <- list(
  m         = 200,
  nx1       = 50000,  
  nx2       = 50000,   
  ny1       = 50000,   
  ny2       = 50000,   
  a_gamma1  = -0.3,
  b_gamma1  = 0.3,
  a_gamma2  = -0.3,
  b_gamma2  = 0.3,
  rho_gamma = 0.8,      
  a_f       = 0.1,
  b_f       = 0.3,
  a_alpha1  = -0.02,
  b_alpha1  = 0.02,
  a_alpha2  = -0.02,
  b_alpha2  = 0.02,
  rho_alpha = 0.4,     
  a_phi1    = -0.05,
  b_phi1    = 0.05,
  a_phi2    = -0.05,
  b_phi2    = 0.05,
  rho_eta   = 0.4,     
  theta1    = 0.2,
  theta2    = 0.2,
  q_uhp1    = 1,
  q_uhp2    = 1,
  q_chp1    = 0.1,
  q_chp2    = 0.1
)

sim_memo <- do.call(dgm6, params_memo)

z_exp1 <- abs(sim_memo$b_exp_1 / sim_memo$se_exp_1)
z_exp2 <- abs(sim_memo$b_exp_2 / sim_memo$se_exp_2)
p_exp1 <- 2 * pnorm(z_exp1, lower.tail = FALSE)
p_exp2 <- 2 * pnorm(z_exp2, lower.tail = FALSE)
sel4   <- (p_exp1 < P_CUTOFF_MEMO) & (p_exp2 < P_CUTOFF_MEMO)

b_exp_mat_m  <- cbind(sim_memo$b_exp_1[sel4],  sim_memo$b_exp_2[sel4])
se_exp_mat_m <- cbind(sim_memo$se_exp_1[sel4], sim_memo$se_exp_2[sel4])
b_out_mat_m  <- cbind(sim_memo$b_out_1[sel4],  sim_memo$b_out_2[sel4])
se_out_mat_m <- cbind(sim_memo$se_out_1[sel4], sim_memo$se_out_2[sel4])

memo_data <- list(
  b_exp     = b_exp_mat_m,
  se_exp    = se_exp_mat_m,
  b_out     = b_out_mat_m,
  se_out    = se_out_mat_m,
  true_beta = c(params_memo$theta1, params_memo$theta2),
  setting   = "two exposures, two outcomes; UHP + CHP",
  params    = params_memo,
  p_cutoff  = P_CUTOFF_MEMO,
  n_iv      = sum(sel4),
  seed      = 4
)

saveRDS(memo_data, file = "data/memo_example.rds")
cat(sprintf("Saved memo_example.rds : %d IVs selected, true betas = (%.2f, %.2f)\n",
            memo_data$n_iv, memo_data$true_beta[1], memo_data$true_beta[2]))