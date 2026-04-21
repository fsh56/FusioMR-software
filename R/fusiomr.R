#' FusioMR: Bayesian Mendelian Randomization
#'
#' Main entry point of the FusioMR package. Estimates the causal effect
#' of an exposure on an outcome from GWAS summary statistics, using a
#' Bayesian model fit with Gibbs sampling. The input format mirrors that
#' of \pkg{TwoSampleMR} and \pkg{MendelianRandomization}.
#'
#' @param b_exp Numeric vector of SNP-exposure effect estimates.
#' @param se_exp Numeric vector of standard errors of \code{b_exp}.
#' @param b_out Numeric vector of SNP-outcome effect estimates.
#' @param se_out Numeric vector of standard errors of \code{b_out}.
#' @param model Character string. One of \code{"seso_uhp_only"},
#'   \code{"seso_with_chp"}, \code{"semo"}, \code{"memo"}.
#'   Currently only \code{"seso_uhp_only"} is implemented.
#' @param niter Number of Gibbs iterations (default 20000).
#' @param burnin_prop Burn-in proportion in [0, 1) (default 0.5).
#' @param p_value_threshold Two-sided p-value threshold on the
#'   exposure Z-scores for IV selection (default 1e-3).
#' @param control A list of advanced prior hyper-parameters returned by
#'   \code{\link{fusiomr_control}}. Defaults are suitable for most uses.
#' @param verbose Logical; if TRUE, print progress messages.
#'
#' @return A list with components
#' \describe{
#'   \item{est}{Posterior mean of the causal effect.}
#'   \item{se}{Posterior SD of the causal effect.}
#'   \item{pval}{Two-sided p-value from a normal approximation.}
#'   \item{ci}{95\% empirical credible interval.}
#'   \item{model}{The model used.}
#'   \item{n_iv}{Number of IVs after selection.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' K <- 50
#' b_exp  <- rnorm(K, 0, 0.1);  se_exp <- rep(0.01, K)
#' b_out  <- 0.3 * b_exp + rnorm(K, 0, 0.01); se_out <- rep(0.01, K)
#' fit <- fusiomr(b_exp, se_exp, b_out, se_out,
#'                model = "seso_uhp_only", niter = 2000)
#' fit$est; fit$ci
#' }
fusiomr <- function(b_exp, se_exp, b_out, se_out,
                    model = c("seso_uhp_only", "seso_with_chp", "semo", "memo"),
                    niter = 20000,
                    burnin_prop = 0.5,
                    p_value_threshold = 1e-3,
                    control = fusiomr_control(),
                    verbose = TRUE) {

  model <- match.arg(model)

  # --- input validation ---------------------------------------------------
  if (!is.numeric(b_exp) || !is.numeric(se_exp) ||
      !is.numeric(b_out) || !is.numeric(se_out))
    stop("b_exp, se_exp, b_out, se_out must all be numeric.")
  b_exp <- as.numeric(b_exp); se_exp <- as.numeric(se_exp)
  b_out <- as.numeric(b_out); se_out <- as.numeric(se_out)
  n <- length(b_exp)
  if (length(se_exp) != n || length(b_out) != n || length(se_out) != n)
    stop("b_exp, se_exp, b_out, se_out must have the same length.")
  if (any(is.na(c(b_exp, se_exp, b_out, se_out))))
    stop("Missing values are not allowed in the summary statistics.")
  if (any(se_exp <= 0) || any(se_out <= 0))
    stop("Standard errors must be strictly positive.")
  if (!is.numeric(niter) || niter < 100)
    stop("niter must be an integer >= 100.")
  if (burnin_prop < 0 || burnin_prop >= 1)
    stop("burnin_prop must be in [0, 1).")
  if (p_value_threshold <= 0 || p_value_threshold >= 1)
    stop("p_value_threshold must be in (0, 1).")

  # --- dispatch ----------------------------------------------------------
  if (model == "seso_uhp_only") {
    return(fit_seso_uhp_only(b_exp, se_exp, b_out, se_out,
                             niter = niter, burnin_prop = burnin_prop,
                             p_value_threshold = p_value_threshold,
                             control = control, verbose = verbose))
  }
  stop(sprintf("Model '%s' is not yet implemented.", model))
}


# Internal worker: seso_uhp_only
fit_seso_uhp_only <- function(b_exp, se_exp, b_out, se_out,
                              niter, burnin_prop, p_value_threshold,
                              control, verbose) {

  if (verbose) message("Running model: seso_uhp_only")

  # --- IV selection ------------------------------------------------------
  # Note: prior is computed on all SNPs (pre-selection) to get a more
  # stable MoM estimate, following the convention in seso_fusio_with_chp.R.
  sel <- get_sel_idx(b_exp, se_exp, p_value_threshold)
  K <- sum(sel)
  if (verbose) message(sprintf("Selected %d of %d IVs (p < %g).",
                               K, length(b_exp), p_value_threshold))
  if (K < 3)
    stop(sprintf("Only %d IVs pass the threshold; need at least 3.", K))
  if (K < 5)
    warning("Fewer than 5 IVs selected; results may be unreliable.")

  # post-selection subsets, passed to the Gibbs sampler
  b_exp_s  <- b_exp[sel];  se_exp_s <- se_exp[sel]
  b_out_s  <- b_out[sel];  se_out_s <- se_out[sel]

  # --- prior: compute on full set, shape-scaled by selected K ------------
  vp <- set_variance_priors(
    ghat = b_exp, gse = se_exp, Ghat = b_out, Gse = se_out,
    beta0 = NULL, K = K,
    Kmin = control$Kmin, Kmax = control$Kmax,
    rho_ov = control$rho_ov,
    c_gamma = control$c_gamma, c_theta = control$c_theta,
    global_mean_gamma = control$global_mean_gamma,
    global_mean_theta = control$global_mean_theta,
    hybrid = control$hybrid, kappa_hybrid = control$kappa_hybrid,
    z_thresh = control$z_thresh, trim = control$trim,
    kappa_gamma = control$kappa_gamma, kappa_theta = control$kappa_theta
  )

  # --- initialize MCMC traces -------------------------------------------
  start_val <- init_setup_seso_uhp_only(
    niter = niter, K = K,
    beta_init = vp$beta0,  # GLS estimate as starting point
    sigma_gamma_init = sqrt(vp$gamma$prior_mean),
    sigma_theta_init = sqrt(vp$theta$prior_mean)
  )

  # --- run Gibbs sampler (C++) ------------------------------------------
  if (verbose) message(sprintf("Gibbs sampling: niter=%d, burn-in=%d",
                               niter, floor(niter * burnin_prop)))
  res <- gibbs_seso_uhp_only_cpp(
    niter = niter, K = K,
    beta_tk = start_val$beta_tk,
    theta_tk = start_val$theta_tk,
    gamma_tk = start_val$gamma_tk,
    sigma2_gamma_tk = start_val$sigma2_gamma_tk,
    sigma2_theta_tk = start_val$sigma2_theta_tk,
    Gamma_hat = b_out_s, gamma_hat = b_exp_s,
    s2_hat_Gamma = se_out_s^2, s2_hat_gamma = se_exp_s^2,
    a_gamma = vp$gamma$a, b_gamma = vp$gamma$b,
    a_theta = vp$theta$a, b_theta = vp$theta$b
  )

  # --- post-burnin summary ----------------------------------------------
  burnin <- floor(niter * burnin_prop)
  draws <- res$beta_tk[(burnin + 1):niter]
  s <- get_summary(draws)

  if (verbose) {
    cat("\n--- Results (seso_uhp_only) ---\n")
    print_summary(s)
  }

  list(est = s$beta_est, se = s$beta_se, pval = s$beta_pval,
       ci = s$ci_emp, model = "seso_uhp_only", n_iv = K)
}
