#' Compute MoM-based variance priors for semo (single-exposure, two-outcome) models
#'
#' Same MoM/EB strategy as \code{\link{set_variance_priors}}, but the outcome
#' is bivariate. Returns an inverse-Wishart prior for the 2x2 pleiotropy
#' covariance Sigma_theta in addition to the gamma prior.
#'
#' @param ghat,gse Exposure GWAS effects and SEs.
#' @param Ghat_mat,Gse_mat K x 2 matrices of outcome GWAS effects and SEs.
#' @param beta0 Optional length-2 init for (beta1, beta2); if NULL, GLS init.
#' @param K Number of IVs after filtering.
#' @param Kmin,Kmax Floor and cap on the effective K used in prior shape.
#' @param rho12 Sampling correlation between the two outcomes (sample overlap).
#' @param rho1g,rho2g Sampling correlations between exposure and each outcome.
#' @param c_gamma,c_theta Prior weights per IV (see set_variance_priors).
#' @param global_mean_gamma,global_mean_theta Optional EB centers (hybrid).
#' @param hybrid Logical; enable hybrid EB if TRUE.
#' @param kappa_hybrid Pooling strength for hybrid mode.
#' @param z_thresh Optional |Z_gamma| selection threshold (winner's-curse fix).
#' @param trim Tail probability for winsorized moments.
#' @param kappa_gamma,kappa_theta Conservative inflation multipliers (>=1).
#'
#' @return A list with sublists \code{gamma} (a, b, prior_mean) and
#'   \code{theta} (nu, Phi, prior_mean — a 2x2 matrix), plus \code{K},
#'   \code{eta}, \code{beta0}, \code{notes}.
#'
#' @keywords internal
#' @export
set_variance_priors_m2 <- function(ghat, gse, Ghat_mat, Gse_mat,
beta0 = NULL, K = NULL,
Kmin = 5, Kmax = 20,
rho12 = 0, rho1g = 0, rho2g = 0,
c_gamma = 0.5, c_theta = 0.8,
global_mean_gamma = NULL, global_mean_theta = NULL,
hybrid = FALSE, kappa_hybrid = 5,
z_thresh = NULL, trim = 0.1,
kappa_gamma = 1, kappa_theta = 1) {

stopifnot(is.matrix(Ghat_mat), is.matrix(Gse_mat),
          ncol(Ghat_mat) == 2, ncol(Gse_mat) == 2,
          nrow(Ghat_mat) == length(ghat),
          nrow(Gse_mat)  == length(gse))
if (is.null(K)) K = length(ghat)
notes = c()

# truncated-normal inflation (winner's-curse)
kappa_TN <- function(c) {
  if (is.null(c) || c <= 0) return(1)
  tail = 1 - stats::pnorm(c)
  1 + c * stats::dnorm(c) / tail
}
infl_g = kappa_TN(z_thresh)

winsor <- function(x, p = trim) {
  if (length(x) < 5) return(x)
  qs = stats::quantile(x, c(p, 1 - p), na.rm = TRUE)
  pmin(pmax(x, qs[[1]]), qs[[2]])
}
wmean <- function(x, p = trim) mean(winsor(x, p), na.rm = TRUE)
wmed  <- function(x, p = trim) stats::median(winsor(x, p), na.rm = TRUE)

# sigma^2_gamma: same as univariate
zg2 = (ghat / gse)^2
mom_gamma = wmed(gse^2, trim) * max(wmean(zg2, trim) - 1, 0)

# beta0: per-outcome GLS init if not provided
if (is.null(beta0)) {
  beta0 = numeric(2)
  for (j in 1:2) {
    wG = 1 / (Gse_mat[, j]^2)
    xtx = sum((ghat^2) * wG); xty = sum((ghat * Ghat_mat[, j]) * wG)
    beta0[j] = if (xtx > 0) xty / xtx else 0
  }
  notes = c(notes, "beta0 estimated by per-outcome GLS")
}

# Per-outcome residual variance (de-noised, overlap-aware)
mom_theta_diag = numeric(2)
for (j in 1:2) {
  rho_jg = if (j == 1) rho1g else rho2g
  r = Ghat_mat[, j] - beta0[j] * ghat
  sv = (Gse_mat[, j]^2) + (beta0[j]^2) * (gse^2) -
       2 * beta0[j] * rho_jg * (Gse_mat[, j] * gse)
  if (!is.null(z_thresh) && z_thresh > 0)
    sv = sv + (beta0[j]^2) * gse^2 * (infl_g - 1)
  zr2 = (r^2) / sv
  mom_theta_diag[j] = wmed(sv, trim) * max(wmean(zr2, trim) - 1, 0)
}

# Off-diagonal: residual cross-product, de-noised by overlap
r1 = Ghat_mat[, 1] - beta0[1] * ghat
r2 = Ghat_mat[, 2] - beta0[2] * ghat
cov_noise = rho12 * (Gse_mat[, 1] * Gse_mat[, 2]) +
            beta0[1] * beta0[2] * (gse^2) -
            beta0[1] * rho2g * (gse * Gse_mat[, 2]) -
            beta0[2] * rho1g * (gse * Gse_mat[, 1])
mom_theta_off = wmean(r1 * r2 - cov_noise, trim)
# bound to a valid correlation
denom_sd = sqrt(mom_theta_diag[1] * mom_theta_diag[2])
if (denom_sd > 1e-12) {
  rho_theta = max(-0.99, min(0.99, mom_theta_off / denom_sd))
  mom_theta_off = rho_theta * denom_sd
} else {
  mom_theta_off = 0
}

# Hybrid EB on diagonals (off-diagonal not pooled)
if (isTRUE(hybrid)) {
  if (is.null(global_mean_gamma) || is.null(global_mean_theta))
    stop("hybrid=TRUE requires global_mean_gamma and global_mean_theta.")
  eta_w = K / (K + kappa_hybrid)
  m_gamma = eta_w * mom_gamma + (1 - eta_w) * global_mean_gamma
  if (length(global_mean_theta) == 1) global_mean_theta = rep(global_mean_theta, 2)
  m_theta_diag = eta_w * mom_theta_diag + (1 - eta_w) * global_mean_theta
} else {
  eta_w = NA
  m_gamma = mom_gamma
  m_theta_diag = mom_theta_diag
}

# Inflation
m_gamma = m_gamma * kappa_gamma
m_theta_diag = m_theta_diag * kappa_theta

# IG shape for gamma
a_gamma = 1 + c_gamma * min(max(K, Kmin), Kmax) / 2
if (a_gamma <= 1) a_gamma = 1 + 1e-3
b_gamma = max((a_gamma - 1) * m_gamma, 1e-6)

# IW prior for theta:
# Sigma_theta ~ IW(nu, Phi) with E[Sigma] = Phi / (nu - p - 1), p = 2.
# Set nu via c_theta on the K scale, then solve for Phi from desired mean.
nu_theta = (2 + 1) + c_theta * min(max(K, Kmin), Kmax)  # = 3 + c_theta*Keff
if (nu_theta <= 3 + 1e-3) nu_theta = 3 + 1e-3
prior_mean_theta = matrix(c(m_theta_diag[1], mom_theta_off,
                            mom_theta_off,   m_theta_diag[2]), 2, 2)
# floor diagonals
prior_mean_theta[1, 1] = max(prior_mean_theta[1, 1], 1e-6)
prior_mean_theta[2, 2] = max(prior_mean_theta[2, 2], 1e-6)
Phi_theta = prior_mean_theta * (nu_theta - 2 - 1)  # nu - p - 1, p=2

list(
  gamma = list(mom_local = mom_gamma, a = a_gamma, b = b_gamma,
               prior_mean = b_gamma / (a_gamma - 1)),
  theta = list(mom_local = prior_mean_theta,
               nu = nu_theta, Phi = Phi_theta,
               prior_mean = prior_mean_theta),
  K = K, eta = eta_w, beta0 = beta0, notes = notes
)
}
