# MCMC trace initialization for the different FusioMR models.

#' Initialize MCMC trace for seso_uhp_only
#'
#' Allocates empty trace matrices/vectors for each parameter and fills
#' the first row (iteration 1) with starting values.
#'
#' @param niter Number of MCMC iterations.
#' @param K Number of selected IVs.
#' @param beta_init Starting value for the causal effect beta.
#' @param sigma_gamma_init Starting SD for gamma (use sqrt of prior mean).
#' @param sigma_theta_init Starting SD for theta (use sqrt of prior mean).
#'
#' @return A list of pre-allocated traces to pass to the C++ sampler.
#' @keywords internal
init_setup_seso_uhp_only <- function(niter, K, beta_init,
                                     sigma_gamma_init, sigma_theta_init) {
  # starting values (latent variables start at zero)
  theta_init = rep(0, K)
  gamma_init = rep(0, K)

  # pre-allocate MCMC traces
  theta_tk = gamma_tk = matrix(NA_real_, nrow = niter, ncol = K)
  beta_tk = sigma2_gamma_tk = sigma2_theta_tk = rep(NA_real_, niter)

  # fill first row with initial values
  theta_tk[1, ] = theta_init
  gamma_tk[1, ] = gamma_init
  beta_tk[1] = beta_init
  sigma2_gamma_tk[1] = sigma_gamma_init^2
  sigma2_theta_tk[1] = sigma_theta_init^2

  list(theta_tk = theta_tk, gamma_tk = gamma_tk, beta_tk = beta_tk,
       sigma2_gamma_tk = sigma2_gamma_tk, sigma2_theta_tk = sigma2_theta_tk)
}
