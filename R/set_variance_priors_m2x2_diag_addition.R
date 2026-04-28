#' Compute MoM-based variance priors for memo (2 exposures, 2 outcomes)
#'
#' Bivariate-exposure version. Returns inverse-Wishart priors for both
#' the SNP-effect covariance Sigma_gamma and the pleiotropy covariance
#' Sigma_theta, used internally by the memo Gibbs sampler.
#'
#' @param ghat_mat,gse_mat K x 2 matrices: per-IV exposure effects and SEs
#'   for the two exposures.
#' @param Ghat_mat,Gse_mat K x 2 matrices: outcome effects and SEs for
#'   the two outcomes.
#' @param B0 Optional length-2 vector of starting beta_j; if NULL, GLS init
#'   per (exposure j, outcome j) pair.
#' @param K Number of IVs after filtering.
#' @param Kmin,Kmax Floor and cap on the effective K used in prior shape.
#' @param rho12 Sampling correlation between the two outcomes.
#' @param rho_gg Sampling correlation between the two exposures.
#' @param rho_gj List of length 2; element j is c(rho_g1_outj, rho_g2_outj),
#'   the sampling correlations between exposure j and the two outcomes.
#' @param c_gamma,c_theta Prior weights per IV.
#' @param global_mean_gamma,global_mean_theta Optional EB centers (hybrid).
#' @param hybrid Logical; enable hybrid EB if TRUE.
#' @param kappa_hybrid Pooling strength for hybrid mode.
#' @param z_thresh Optional |Z| selection threshold (winner's-curse fix).
#' @param trim Tail probability for winsorized moments.
#' @param kappa_gamma,kappa_theta Conservative inflation multipliers (>=1).
#'
#' @return A list with sublists \code{gamma} (nu, Phi, prior_mean — a 2x2
#'   matrix) and \code{theta} (same structure), plus \code{K}, \code{eta},
#'   \code{B0}, \code{notes}.
#'
#' @keywords internal
#' @export
set_variance_priors_m2x2_diag <- function(ghat_mat, gse_mat, Ghat_mat, Gse_mat,
B0 = NULL, K = NULL, Kmin = 5, Kmax = 20,
rho12 = 0, rho_gg = 0,
rho_gj = list(c(0, 0), c(0, 0)),
c_gamma = 0.5, c_theta = 0.8,
global_mean_gamma = NULL, global_mean_theta = NULL,
hybrid = FALSE, kappa_hybrid = 5,
z_thresh = NULL, trim = 0.1,
kappa_gamma = 1, kappa_theta = 1) {

stopifnot(is.matrix(ghat_mat), is.matrix(gse_mat),
          is.matrix(Ghat_mat), is.matrix(Gse_mat),
          ncol(ghat_mat) == 2, ncol(gse_mat) == 2,
          ncol(Ghat_mat) == 2, ncol(Gse_mat) == 2,
          nrow(ghat_mat) == nrow(Ghat_mat))
if (is.null(K)) K = nrow(ghat_mat)
notes = c()

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

# --- Sigma_gamma diagonals: per-exposure MoM (de-noised) ----------------
mom_gamma_diag = numeric(2)
for (j in 1:2) {
  zg2 = (ghat_mat[, j] / gse_mat[, j])^2
  mom_gamma_diag[j] = wmed(gse_mat[, j]^2, trim) *
                     max(wmean(zg2, trim) - 1, 0)
}

# --- Sigma_gamma off-diagonal: cross-exposure (overlap-aware) -----------
cov_noise_gg = rho_gg * (gse_mat[, 1] * gse_mat[, 2])
mom_gamma_off = wmean(ghat_mat[, 1] * ghat_mat[, 2] - cov_noise_gg, trim)
denom = sqrt(mom_gamma_diag[1] * mom_gamma_diag[2])
if (denom > 1e-12) {
  rho_g = max(-0.99, min(0.99, mom_gamma_off / denom))
  mom_gamma_off = rho_g * denom
} else {
  mom_gamma_off = 0
}

# --- B0: per-(exposure j, outcome j) GLS init --------------------------
if (is.null(B0)) {
  B0 = numeric(2)
  for (j in 1:2) {
    wG = 1 / (Gse_mat[, j]^2)
    xtx = sum((ghat_mat[, j]^2) * wG)
    xty = sum((ghat_mat[, j] * Ghat_mat[, j]) * wG)
    B0[j] = if (xtx > 0) xty / xtx else 0
  }
  notes = c(notes, "B0 estimated by per-pair GLS")
}

# --- Sigma_theta diagonals: per-outcome residual variance --------------
mom_theta_diag = numeric(2)
for (j in 1:2) {
  rho_jg = rho_gj[[j]][j]    # exposure j to outcome j sample-overlap
  r = Ghat_mat[, j] - B0[j] * ghat_mat[, j]
  sv = (Gse_mat[, j]^2) + (B0[j]^2) * (gse_mat[, j]^2) -
       2 * B0[j] * rho_jg * (Gse_mat[, j] * gse_mat[, j])
  if (!is.null(z_thresh) && z_thresh > 0)
    sv = sv + (B0[j]^2) * gse_mat[, j]^2 * (infl_g - 1)
  zr2 = (r^2) / sv
  mom_theta_diag[j] = wmed(sv, trim) *
                     max(wmean(zr2, trim) - 1, 0)
}

# --- Sigma_theta off-diagonal ------------------------------------------
r1 = Ghat_mat[, 1] - B0[1] * ghat_mat[, 1]
r2 = Ghat_mat[, 2] - B0[2] * ghat_mat[, 2]
cov_noise_t = rho12 * (Gse_mat[, 1] * Gse_mat[, 2]) +
              B0[1] * B0[2] * rho_gg * (gse_mat[, 1] * gse_mat[, 2])
mom_theta_off = wmean(r1 * r2 - cov_noise_t, trim)
denom_t = sqrt(mom_theta_diag[1] * mom_theta_diag[2])
if (denom_t > 1e-12) {
  rho_t = max(-0.99, min(0.99, mom_theta_off / denom_t))
  mom_theta_off = rho_t * denom_t
} else {
  mom_theta_off = 0
}

# --- Hybrid EB on diagonals --------------------------------------------
if (isTRUE(hybrid)) {
  if (is.null(global_mean_gamma) || is.null(global_mean_theta))
    stop("hybrid=TRUE requires global_mean_gamma and global_mean_theta.")
  eta_w = K / (K + kappa_hybrid)
  if (length(global_mean_gamma) == 1) global_mean_gamma = rep(global_mean_gamma, 2)
  if (length(global_mean_theta) == 1) global_mean_theta = rep(global_mean_theta, 2)
  m_gamma_diag = eta_w * mom_gamma_diag + (1 - eta_w) * global_mean_gamma
  m_theta_diag = eta_w * mom_theta_diag + (1 - eta_w) * global_mean_theta
} else {
  eta_w = NA
  m_gamma_diag = mom_gamma_diag
  m_theta_diag = mom_theta_diag
}
m_gamma_diag = m_gamma_diag * kappa_gamma
m_theta_diag = m_theta_diag * kappa_theta

# --- IW priors: nu, Phi (with E[Sigma] = Phi/(nu - p - 1), p = 2) -------
Keff = min(max(K, Kmin), Kmax)
nu_gamma = (2 + 1) + c_gamma * Keff
nu_theta = (2 + 1) + c_theta * Keff
if (nu_gamma <= 3 + 1e-3) nu_gamma = 3 + 1e-3
if (nu_theta <= 3 + 1e-3) nu_theta = 3 + 1e-3

prior_mean_gamma = matrix(c(max(m_gamma_diag[1], 1e-6), mom_gamma_off,
                            mom_gamma_off, max(m_gamma_diag[2], 1e-6)), 2, 2)
prior_mean_theta = matrix(c(max(m_theta_diag[1], 1e-6), mom_theta_off,
                            mom_theta_off, max(m_theta_diag[2], 1e-6)), 2, 2)
Phi_gamma = prior_mean_gamma * (nu_gamma - 2 - 1)
Phi_theta = prior_mean_theta * (nu_theta - 2 - 1)

list(
  gamma = list(mom_local = prior_mean_gamma,
               nu = nu_gamma, Phi = Phi_gamma,
               prior_mean = prior_mean_gamma),
  theta = list(mom_local = prior_mean_theta,
               nu = nu_theta, Phi = Phi_theta,
               prior_mean = prior_mean_theta),
  K = K, eta = eta_w, B0 = B0, notes = notes
)
}
