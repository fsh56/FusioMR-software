# Set variance priors for FusioMR (sigma^2_gamma and sigma^2_theta)
# LD-free; uses clumped, approximately independent SNPs.

#' Compute MoM-based variance priors for seso models
#'
#' Computes method-of-moments (MoM) empirical-Bayes priors for the
#' SNP-effect variance (sigma^2_gamma) and the horizontal-pleiotropy
#' variance (sigma^2_theta), used internally by \code{\link{fusiomr}}.
#'
#' @param ghat,gse QTL effects and SEs for the exposure.
#' @param Ghat,Gse GWAS effects and SEs for the outcome.
#' @param beta0 Optional initial beta; if NULL, GLS init on outcome SEs.
#' @param K Number of IVs after filtering.
#' @param Kmin,Kmax Floor/cap on K when computing prior shape \code{a}.
#' @param rho_ov Sampling correlation due to sample overlap, in [-1,1].
#' @param c_gamma,c_theta Prior weight per IV (see Details).
#' @param global_mean_gamma,global_mean_theta Optional global EB centers (hybrid mode).
#' @param hybrid Logical; if TRUE, prior_mean = eta*local + (1-eta)*global.
#' @param kappa_hybrid Pooling strength for hybrid mode (>0).
#' @param z_thresh Optional |Z_gamma| selection threshold
#'   (winner's-curse truncated-normal correction).
#' @param trim Tail probability for winsorized moments.
#' @param kappa_gamma,kappa_theta Conservative inflation multipliers (>=1).
#'
#' @return A list with sublists \code{gamma} and \code{theta}, each
#'   containing \code{mom_local}, \code{a}, \code{b}, \code{prior_mean};
#'   plus shared \code{K}, \code{eta}, \code{beta0}, \code{notes}.
#'
#' @keywords internal
#' @export
set_variance_priors <- function(ghat, gse, Ghat, Gse, beta0 = NULL, K = NULL,
Kmin = 5, Kmax = 20, rho_ov = 0, c_gamma = 0.5, c_theta = 0.8,
global_mean_gamma = NULL, global_mean_theta = NULL,
hybrid = FALSE, kappa_hybrid = 5, z_thresh = NULL, trim = 0.1,
kappa_gamma = 1, kappa_theta = 1) {

stopifnot(length(ghat) == length(gse), length(Ghat) == length(Gse),
          length(ghat) == length(Ghat))
if (is.null(K)) K = length(ghat)
notes = c()

# truncated-normal inflation E[Z^2 | |Z|>c] for winner's-curse fix
kappa_TN <- function(c) {
  if (is.null(c) || c <= 0) return(1)
  tail = 1 - stats::pnorm(c)
  1 + c * stats::dnorm(c) / tail
}
infl_g = kappa_TN(z_thresh)

# robust trimmed mean/median via winsorization
winsor <- function(x, p = trim) {
  if (length(x) < 5) return(x)
  qs = stats::quantile(x, c(p, 1 - p), na.rm = TRUE)
  pmin(pmax(x, qs[[1]]), qs[[2]])
}
wmean <- function(x, p = trim) mean(winsor(x, p), na.rm = TRUE)
wmed  <- function(x, p = trim) stats::median(winsor(x, p), na.rm = TRUE)

if (!is.null(z_thresh) && z_thresh > 0)
  notes = c(notes, sprintf("TN correction for gamma with |Z|>%.2f (infl=%.3f)",
                           z_thresh, infl_g))

# sigma^2_gamma: local MoM (de-noised, TN-corrected if requested)
zg2 = (ghat / gse)^2
mom_gamma = wmed(gse^2, trim) * max(wmean(zg2, trim) - 1, 0)

# beta0 for theta prior; GLS init on outcome SEs if not provided
if (is.null(beta0)) {
  wG = 1 / (Gse^2)
  xtx = sum((ghat^2) * wG); xty = sum((ghat * Ghat) * wG)
  beta0 = if (xtx > 0) xty / xtx else 0
  notes = c(notes, "beta0 estimated by GLS on outcome SEs")
} else {
  notes = c(notes, "beta0 provided by user")
}

# sigma^2_theta: residual MoM (de-noised, overlap-aware)
r = Ghat - beta0 * ghat
sv = (Gse^2) + (beta0^2) * (gse^2) - 2 * beta0 * rho_ov * (Gse * gse)
if (!is.null(z_thresh) && z_thresh > 0) {
  sv = sv + (beta0^2) * gse^2 * (infl_g - 1)
  notes = c(notes, "TN correction applied to gamma component inside residual variance")
}
zr2 = (r^2) / sv
mom_theta = wmed(sv, trim) * max(wmean(zr2, trim) - 1, 0)

# shapes (pseudo-df) and hybrid prior means
a_gamma = 1 + c_gamma * min(max(K, Kmin), Kmax) / 2
a_theta = 1 + c_theta * min(max(K, Kmin), Kmax) / 2
if (a_gamma <= 1) a_gamma = 1 + 1e-3
if (a_theta <= 1) a_theta = 1 + 1e-3

if (isTRUE(hybrid)) {
  if (is.null(global_mean_gamma) || is.null(global_mean_theta))
    stop("hybrid=TRUE requires global_mean_gamma and global_mean_theta.")
  eta = K / (K + kappa_hybrid)
  m_gamma = eta * mom_gamma + (1 - eta) * global_mean_gamma
  m_theta = eta * mom_theta + (1 - eta) * global_mean_theta
  notes = c(notes, sprintf("hybrid prior means with eta=K/(K+kappa)=%.3f", eta))
} else {
  eta = NA
  m_gamma = mom_gamma
  m_theta = mom_theta
}

# conservative inflation
m_gamma = m_gamma * kappa_gamma
m_theta = m_theta * kappa_theta

# scales b so that E[sigma^2] = prior mean (valid for a > 1)
b_gamma = (a_gamma - 1) * m_gamma
b_theta = (a_theta - 1) * m_theta

# safeguard: floor scale to avoid degenerate IG draws when MoM -> 0
b_gamma = max(b_gamma, 1e-6)
b_theta = max(b_theta, 1e-6)

list(
  gamma = list(mom_local = mom_gamma, a = a_gamma, b = b_gamma,
               prior_mean = b_gamma / (a_gamma - 1)),
  theta = list(mom_local = mom_theta, a = a_theta, b = b_theta,
               prior_mean = b_theta / (a_theta - 1)),
  K = K, eta = eta, beta0 = beta0, notes = unique(notes)
)
}
