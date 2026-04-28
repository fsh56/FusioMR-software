# Helper functions: IV selection, posterior summary.

# Select IVs by |Z| p-value threshold on exposure
get_sel_idx <- function(b_exp, se_exp, p_threshold) {
  z = abs(b_exp / se_exp)
  p = 2 * (1 - stats::pnorm(z))
  p < p_threshold
}

# Empirical credible interval from posterior draws
get_empirical_ci <- function(x, alpha = 0.05) {
  stats::quantile(x, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
}

# Normal-approximation CI from mean + se
get_normal_ci <- function(m, se, alpha = 0.05) {
  z = stats::qnorm(1 - alpha / 2)
  c(m - z * se, m + z * se)
}

# Posterior summary from a scalar parameter's draws
get_summary <- function(draws) {
  m = mean(draws)
  s = stats::sd(draws)
  pval = 2 * (1 - stats::pnorm(abs(m / s)))
  list(beta_est = m, beta_se = s, beta_pval = pval,
       ci_emp = get_empirical_ci(draws),
       ci_normal = get_normal_ci(m, s))
}

# Pretty-print a single-outcome summary
print_summary <- function(s) {
  cat(sprintf("Estimated Causal Effect (Beta): %.4f\n", s$beta_est))
  cat(sprintf("Standard Error: %.4f\n", s$beta_se))
  cat(sprintf("P-value: %.4g\n", s$beta_pval))
  cat(sprintf("95%% Empirical CI: [%.4f, %.4f]\n",
              s$ci_emp[1], s$ci_emp[2]))
  cat(sprintf("95%% Normal CI:    [%.4f, %.4f]\n",
              s$ci_normal[1], s$ci_normal[2]))
}

# Label-switching correction for seso_with_chp.
# When q > 0.5, the eta = 0 / eta = 1 labels have flipped during MCMC,
# so the true causal effect lives on (beta + alpha) rather than beta.
label_flip <- function(niter, res) {
  ids <- (floor(niter / 2) + 1):niter
  qq <- mean(res$q_tk[ids], na.rm = TRUE)
  
  if (qq > 0.5) {
    # labels flipped: true beta is beta + alpha
    samples <- res$beta_tk[ids] + res$alpha_tk[ids]
  } else {
    # labels normal: true beta is beta
    samples <- res$beta_tk[ids]
  }
  
  list(qq = qq,
       b_mean = mean(samples, na.rm = TRUE),
       b_sd   = stats::sd(samples, na.rm = TRUE),
       bci    = stats::quantile(samples, c(0.025, 0.975), na.rm = TRUE),
       draws  = samples)
}

# Label-switching correction for memo (joint 2-outcome).
# When q_j > 0.5, the eta_j labels have flipped during MCMC,
# so the true causal effect for outcome j is (beta_j + alpha_j).
label_flip_joint <- function(niter, res) {
  ids <- (floor(niter / 2) + 1):niter
  p00 <- res$pst_tk[ids, 1]
  p01 <- res$pst_tk[ids, 2]
  p10 <- res$pst_tk[ids, 3]
  p11 <- res$pst_tk[ids, 4]
  qq1 <- mean(p10 + p11, na.rm = TRUE)
  qq2 <- mean(p01 + p11, na.rm = TRUE)
  
  draws_1 <- if (qq1 > 0.5) res$beta_1_tk[ids] + res$alpha_1_tk[ids]
  else           res$beta_1_tk[ids]
  draws_2 <- if (qq2 > 0.5) res$beta_2_tk[ids] + res$alpha_2_tk[ids]
  else           res$beta_2_tk[ids]
  
  list(qq1 = qq1, qq2 = qq2,
       b1_mean = mean(draws_1, na.rm = TRUE),
       b1_sd   = stats::sd(draws_1, na.rm = TRUE),
       b2_mean = mean(draws_2, na.rm = TRUE),
       b2_sd   = stats::sd(draws_2, na.rm = TRUE),
       bci1 = stats::quantile(draws_1, c(0.025, 0.975), na.rm = TRUE),
       bci2 = stats::quantile(draws_2, c(0.025, 0.975), na.rm = TRUE))
}
