#' Control advanced parameters for FusioMR models
#'
#' Bundles advanced hyper-parameters for the MCMC settings, IV-selection thresholds,
#' and empirical-Bayes variance priors setting ups. Typically, users never need to 
#' touch this: \code{fusiomr()} uses \code{parameter_control()} with sensible defaults
#' tuned for estimation. Pass a customized value of parameters only when needed, e.g., 
#' to shorten MCMC niter, burnin prop, tighten IV selection, or correct for winner's 
#' curse via \code{z_thresh}, account for sample overlap via \code{rho_ov}, or
#' tune prior strength via \code{c_gamma}/\code{c_theta}.
#'
#' @param niter Number of Gibbs iterations. Default 20000 is calibrated for stable 
#' posterior estimates;
#' @param burnin_prop Burn-in proportion in [0, 1).
#' @param c_gamma Prior weight per IV for sigma^2_gamma (a_gamma = 1 + c_gamma*K/2), 
#'  Larger values make the inverse-gamma prior more concentrated around its mean. Default= 0.5;
#' @param c_theta Prior weight per IV for sigma^2_theta (a_theta = 1 + c_theta*K/2),
#'  Larger values make the inverse-gamma prior more concentrated around its mean. Default = 0.8;
#' @param kappa_gamma Tunning parameter for prior mean of sigma2_gamma; 
#' @param kappa_theta Tunning parameter for prior mean of sigma2_theta;
#' @param Kmin,Kmax Lower and upper bound on K used when computing prior shape.
#' @param rho_ov Sampling correlation between exposure and outcome due
#'  to sample overlap, in [-1, 1]. Default = 0;
#' @param z_thresh Optional |Z_gamma| selection threshold used to pick QTLs (winner’s-curse fix).
#'  Example: for p=5e-8 two-sided, use \code{qnorm(1 - 5e-8/2)}.
#' @param trim Tail probability for winsor.
#' @param hybrid Logical; if TRUE, blend local with a global: prior mean = eta*local + (1-eta)*global 
#'  (eta = K/(K+kappa_hybrid));
#' @param kappa_hybrid Pooling control; larger values shrink more
#'   toward the global.
#' @param global_mean_gamma,global_mean_theta Global EB centers for hybrid mode.
#'
#' @return A named list of parameters for more advanced setting.
#' @export
#'
#' @examples
#' # defaults
#' ctrl <- parameter_control()
#' 
parameter_control <- function(
    niter = 20000,
    burnin_prop = 0.5,
    c_gamma = 0.5,
    c_theta = 0.8,
    kappa_gamma = 1,
    kappa_theta = 1,
    Kmin = 5,
    Kmax = 20,
    rho_ov = 0,
    z_thresh = NULL,
    trim = 0.1,
    hybrid = FALSE,
    kappa_hybrid = 5,
    global_mean_gamma = NULL,
    global_mean_theta = NULL
) {
  list(
    niter = niter,
    burnin_prop = burnin_prop,
    c_gamma = c_gamma, 
    c_theta = c_theta,
    kappa_gamma = kappa_gamma, 
    kappa_theta = kappa_theta,
    Kmin = Kmin, 
    Kmax = Kmax,
    rho_ov = rho_ov, 
    z_thresh = z_thresh, 
    trim = trim,
    hybrid = hybrid, 
    kappa_hybrid = kappa_hybrid,
    global_mean_gamma = global_mean_gamma,
    global_mean_theta = global_mean_theta
  )
}