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

#' Initialize MCMC trace for seso_with_chp
#'
#' Allocates trace storage for the seso_with_chp model, which adds
#' alpha (CHP magnitude), eta (CHP indicators), and q (CHP proportion)
#' on top of the seso_uhp_only parameters.
#'
#' @param niter Number of MCMC iterations.
#' @param K Number of selected IVs.
#' @param alpha_init Starting value for the CHP magnitude alpha.
#' @param beta_init Starting value for the causal effect beta.
#' @param sigma_gamma_init Starting SD for gamma.
#' @param sigma_theta_init Starting SD for theta.
#' @param q_init Starting value for the CHP proportion q in [0, 1].
#'
#' @return A list of pre-allocated traces to pass to the C++ sampler.
#' @keywords internal
init_setup_seso_with_chp <- function(niter, K, alpha_init, beta_init,
                                     sigma_gamma_init, sigma_theta_init,
                                     q_init = 0.1) {
  # latent variables start at zero
  theta_init = rep(0, K)
  gamma_init = rep(0, K)
  eta_init   = rep(0, K)
  
  # pre-allocate MCMC traces
  theta_tk = gamma_tk = eta_tk = matrix(NA_real_, nrow = niter, ncol = K)
  alpha_tk = beta_tk = q_tk = sigma2_gamma_tk = sigma2_theta_tk =
    rep(NA_real_, niter)
  
  # fill first row with initial values
  theta_tk[1, ] = theta_init
  gamma_tk[1, ] = gamma_init
  eta_tk[1, ]   = eta_init
  alpha_tk[1] = alpha_init
  beta_tk[1]  = beta_init
  q_tk[1]     = q_init
  sigma2_gamma_tk[1] = sigma_gamma_init^2
  sigma2_theta_tk[1] = sigma_theta_init^2
  
  list(theta_tk = theta_tk, gamma_tk = gamma_tk, eta_tk = eta_tk,
       alpha_tk = alpha_tk, beta_tk = beta_tk, q_tk = q_tk,
       sigma2_gamma_tk = sigma2_gamma_tk, sigma2_theta_tk = sigma2_theta_tk)
}

#' Initialize MCMC trace for semo (single-exposure, two-outcome) UHP-only model
#'
#' @param niter Number of MCMC iterations.
#' @param K Number of selected IVs.
#' @param beta_1_init,beta_2_init Starting values for the two causal effects.
#' @param sigma_gamma_init Starting SD for gamma.
#'
#' @return A list of pre-allocated traces to pass to the C++ sampler.
#' @keywords internal
init_setup_semo_uhp_only <- function(niter, K, beta_1_init, beta_2_init,
                                     sigma_gamma_init) {
  theta_1_init = theta_2_init = rep(0, K)
  gamma_init = rep(0, K)
  
  gamma_tk = theta_1_tk = theta_2_tk =
    matrix(NA_real_, nrow = niter, ncol = K)
  beta_1_tk = beta_2_tk = sigma2_gamma_tk = rep(NA_real_, niter)
  
  gamma_tk[1, ]   = gamma_init
  theta_1_tk[1, ] = theta_1_init
  theta_2_tk[1, ] = theta_2_init
  beta_1_tk[1] = beta_1_init
  beta_2_tk[1] = beta_2_init
  sigma2_gamma_tk[1] = sigma_gamma_init^2
  
  list(theta_1_tk = theta_1_tk, theta_2_tk = theta_2_tk,
       gamma_tk = gamma_tk,
       beta_1_tk = beta_1_tk, beta_2_tk = beta_2_tk,
       sigma2_gamma_tk = sigma2_gamma_tk)
}