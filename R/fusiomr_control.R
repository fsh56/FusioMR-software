#' Control parameters for FusioMR variance priors
#'
#' Bundles advanced hyper-parameters for the empirical-Bayes variance
#' priors. Typical users never need to touch this: \code{fusiomr()} uses
#' \code{fusiomr_control()} with sensible defaults. Pass a customized
#' control only if you want to, e.g., correct for winner's curse via
#' \code{z_thresh}, account for sample overlap via \code{rho_ov}, or
#' tune prior strength via \code{c_gamma}/\code{c_theta}.
#'
#' @param c_gamma,c_theta Prior weight per IV. Larger values make the
#'   inverse-gamma prior more concentrated around its mean. Defaults
#'   0.5 and 0.8 match the simulation settings in the FusioMR paper.
#' @param kappa_gamma,kappa_theta Conservative inflation multipliers on
#'   the prior means (>= 1). Leave at 1 unless you want a deliberately
#'   wider prior.
#' @param Kmin,Kmax Lower and upper bound on K used when computing prior shape.
#' @param rho_ov Sampling correlation between exposure and outcome due
#'   to sample overlap, in [-1, 1].
#' @param z_thresh Optional |Z_gamma| selection threshold for
#'   winner's-curse (truncated-normal) correction; e.g.
#'   \code{qnorm(1 - 5e-8/2)} for a genome-wide threshold.
#' @param trim Tail probability for winsor.
#' @param hybrid Logical; if TRUE, blend local MoM with a global EB
#'   center (requires \code{global_mean_gamma} and \code{global_mean_theta}).
#' @param kappa_hybrid Pooling strength; larger values shrink more
#'   toward the global center.
#' @param global_mean_gamma,global_mean_theta Optional global EB
#'   centers for hybrid mode.
#'
#' @return A named list of control parameters.
#' @export
#'
#' @examples
#' # defaults
#' ctrl <- fusiomr_control()
#'
#' # with winner's-curse correction at genome-wide threshold
#' ctrl <- fusiomr_control(z_thresh = qnorm(1 - 5e-8 / 2))
fusiomr_control <- function(
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
    c_gamma = c_gamma, c_theta = c_theta,
    kappa_gamma = kappa_gamma, kappa_theta = kappa_theta,
    Kmin = Kmin, Kmax = Kmax,
    rho_ov = rho_ov, z_thresh = z_thresh, trim = trim,
    hybrid = hybrid, kappa_hybrid = kappa_hybrid,
    global_mean_gamma = global_mean_gamma,
    global_mean_theta = global_mean_theta
  )
}
