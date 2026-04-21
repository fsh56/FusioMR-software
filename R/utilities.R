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
