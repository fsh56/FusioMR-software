#' FusioMR: A Flexible and unified MR framework using summary statistics for 
#' single- and multi-outcome, tailored for molecular trait exposures 
#' while also applicable to complex trait exposures. 
#'
#' Main function of the FusioMR package. Offer robust estimates the causal effect
#' of an exposure on an outcome from GWAS summary statistics, using Bayesian 
#' hierarchical model and uses Gibbs sampling.
#'
#' @param b_exp Numeric vector of SNP-exposure effect estimates.
#' @param se_exp Numeric vector of standard errors of \code{b_exp}.
#' @param b_out Numeric vector of SNP-outcome effect estimates.
#' @param se_out Numeric vector of standard errors of \code{b_out}.
#' @param model Character string. One of \code{"seso_uhp_only"},
#'   \code{"seso_with_chp"}, \code{"semo"}, \code{"memo"}.
#' @param control A list of advanced prior hyper-parameters returned by
#'   \code{\link{parameter_control}}. Defaults are suitable for most uses.
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
#'                model = "seso_uhp_only", 
#'                control = parameter_control(niter=2000))
#' fit$est; fit$ci
#' }
fusiomr <- function(b_exp, se_exp, b_out, se_out,
                    model = c("seso_uhp_only", "seso_with_chp", "semo", "memo"),
                    control = parameter_control(),
                    verbose = FALSE) {
  
  model <- match.arg(model)
  
  # --- input validation ---------------------------------------------------
  if (!is.numeric(b_exp) || !is.numeric(se_exp) ||
      !is.numeric(b_out) || !is.numeric(se_out))
    stop("b_exp, se_exp, b_out, se_out must all be numeric.")
  
  # exposure: vector (single exposure) or K x 2 matrix (two exposures)
  if (is.matrix(b_exp) || is.matrix(se_exp)) {
    if (!is.matrix(b_exp) || !is.matrix(se_exp))
      stop("b_exp and se_exp must both be matrices or both be vectors.")
    if (ncol(b_exp) != 2 || ncol(se_exp) != 2)
      stop("Matrix b_exp / se_exp must have exactly 2 columns.")
    if (nrow(b_exp) != nrow(se_exp))
      stop("b_exp / se_exp rows must match.")
    n <- nrow(b_exp)
  } else {
    b_exp <- as.numeric(b_exp); se_exp <- as.numeric(se_exp)
    n <- length(b_exp)
    if (length(se_exp) != n)
      stop("b_exp and se_exp must have the same length.")
  }
  # outcome: vector (single outcome) or K x 2 matrix (two outcomes)
  if (is.matrix(b_out) || is.matrix(se_out)) {
    if (!is.matrix(b_out) || !is.matrix(se_out))
      stop("b_out and se_out must both be matrices or both be vectors.")
    if (ncol(b_out) != 2 || ncol(se_out) != 2)
      stop("Matrix b_out / se_out must have exactly 2 columns.")
    if (nrow(b_out) != n || nrow(se_out) != n)
      stop("b_out / se_out rows must match length/rows of b_exp.")
  } else {
    b_out <- as.numeric(b_out); se_out <- as.numeric(se_out)
    if (length(b_out) != n || length(se_out) != n)
      stop("b_out and se_out must have the same length as b_exp.")
  }
  if (any(is.na(c(b_exp, se_exp, b_out, se_out))))
    stop("Missing values are not allowed in the summary statistics.")
  if (any(se_exp <= 0) || any(se_out <= 0))
    stop("Standard errors must be strictly positive.")
  
  # --- dispatch ----------------------------------------------------------
  # unpack MCMC/selection settings from control
  niter <- control$niter
  burnin_prop <- control$burnin_prop
  if (!is.numeric(niter) || niter < 100)
    stop("control$niter must be an integer >= 100.")
  if (burnin_prop < 0 || burnin_prop >= 1)
    stop("control$burnin_prop must be in [0, 1).")
  
  # model1
  if (model == "seso_uhp_only") {
    return(fit_seso_uhp_only(b_exp, se_exp, b_out, se_out,
                             control = control, verbose = verbose))
  }
  
  # model2
  if (model == "seso_with_chp") {
    return(fit_seso_with_chp(b_exp, se_exp, b_out, se_out,
                             control = control, verbose = verbose))
  }
  
  # model3
  if (model == "semo") {
    return(fit_semo(b_exp, se_exp, b_out, se_out,
                    control = control, verbose = verbose))
  }
  
  # model4
  if (model == "memo") {
    return(fit_memo(b_exp, se_exp, b_out, se_out,
                    control = control, verbose = verbose))
  }
  
  stop(sprintf("Model '%s' is not yet implemented.", model))
}

  

# Internal worker: seso_uhp_only
fit_seso_uhp_only <- function(b_exp, se_exp, b_out, se_out,
                              control, verbose = FALSE) {
  niter <- control$niter
  burnin_prop <- control$burnin_prop
  
  message("Running model: seso_uhp_only")
  K <- length(b_exp)
  if (K < 5)
    warning("Fewer than 5 IVs selected;")
  
  # --- compute MoM priors on the input data------------
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
    Gamma_hat = b_out, gamma_hat = b_exp,
    s2_hat_Gamma = se_out^2, s2_hat_gamma = se_exp^2,
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

# Internal worker: seso_with_chp
fit_seso_with_chp <- function(b_exp, se_exp, b_out, se_out,
                              control, verbose = FALSE) {
  niter <- control$niter
  burnin_prop <- control$burnin_prop
  
  message("Running model: seso_with_chp")
  
  # all input SNPs are treated as IVs (assume upstream filtering)
  K <- length(b_exp)
  if (K < 5)
    warning("Fewer than 5 IVs selected;")
  
  # --- compute MoM priors ----------------------------------------------
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
  
  # --- initialize MCMC traces ------------------------------------------
  start_val <- init_setup_seso_with_chp(
    niter = niter, K = K,
    alpha_init = 1,                       # default in original script
    beta_init = vp$beta0,                 # GLS estimate as starting point
    sigma_gamma_init = sqrt(vp$gamma$prior_mean),
    sigma_theta_init = sqrt(vp$theta$prior_mean),
    q_init = 0.1
  )
  
  # --- run Gibbs sampler (C++) -----------------------------------------
  if (verbose) message(sprintf("Gibbs sampling: niter=%d, burn-in=%d",
                               niter, floor(niter * burnin_prop)))
  res <- gibbs_seso_with_chp_cpp(
    niter = niter, K = K,
    beta_tk = start_val$beta_tk,
    alpha_tk = start_val$alpha_tk,
    eta_tk = start_val$eta_tk,
    theta_tk = start_val$theta_tk,
    gamma_tk = start_val$gamma_tk,
    q_tk = start_val$q_tk,
    Gamma_hat = b_out, gamma_hat = b_exp,
    s2_hat_Gamma = se_out^2, s2_hat_gamma = se_exp^2,
    sigma2_gamma_tk = start_val$sigma2_gamma_tk,
    sigma2_theta_tk = start_val$sigma2_theta_tk,
    a_gamma = vp$gamma$a, b_gamma = vp$gamma$b,
    a_theta = vp$theta$a, b_theta = vp$theta$b,
    a_q = 1.0, b_q = 1.0     # uniform Beta(1,1) prior on q
  )
  
  # --- post-processing with label-switching correction -----------------
  flip <- label_flip(niter, res)
  pval <- 2 * (1 - stats::pnorm(abs(flip$b_mean / flip$b_sd)))
  
  if (verbose) {
    cat("\n--- Results (seso_with_chp) ---\n")
    cat(sprintf("Estimated Causal Effect (Beta): %.4f\n", flip$b_mean))
    cat(sprintf("Standard Error: %.4f\n", flip$b_sd))
    cat(sprintf("P-value: %.4g\n", pval))
    cat(sprintf("95%% Empirical CI: [%.4f, %.4f]\n", flip$bci[1], flip$bci[2]))
    cat(sprintf("Posterior CHP proportion (q): %.3f%s\n",
                flip$qq, if (flip$qq > 0.5) " (labels flipped)" else ""))
  }
  
  list(est = flip$b_mean, se = flip$b_sd, pval = pval,
       ci = flip$bci, q = flip$qq,
       model = "seso_with_chp", n_iv = K)
}

# Internal worker: semo (single exposure, two outcomes; UHP only)
fit_semo <- function(b_exp, se_exp, b_out, se_out,
                     control, verbose = FALSE) {
  niter <- control$niter
  burnin_prop <- control$burnin_prop
  
  message("Running model: semo")
  
  if (!is.matrix(b_out) || ncol(b_out) != 2)
    stop("Model 'semo' requires b_out and se_out to be K x 2 matrices.")
  
  K <- length(b_exp)
  if (K < 5)
    warning("Fewer than 5 IVs selected;")
  
  # --- compute MoM priors (bivariate outcome version) -------------------
  vp <- set_variance_priors_m2(
    ghat = b_exp, gse = se_exp,
    Ghat_mat = b_out, Gse_mat = se_out,
    beta0 = NULL, K = K,
    Kmin = control$Kmin, Kmax = control$Kmax,
    rho12 = 0, rho1g = control$rho_ov, rho2g = control$rho_ov,
    c_gamma = control$c_gamma, c_theta = control$c_theta,
    global_mean_gamma = control$global_mean_gamma,
    global_mean_theta = control$global_mean_theta,
    hybrid = control$hybrid, kappa_hybrid = control$kappa_hybrid,
    z_thresh = control$z_thresh, trim = control$trim,
    kappa_gamma = control$kappa_gamma, kappa_theta = control$kappa_theta
  )
  
  # --- initialize MCMC traces ------------------------------------------
  start_val <- init_setup_semo_uhp_only(
    niter = niter, K = K,
    beta_1_init = vp$beta0[1],
    beta_2_init = vp$beta0[2],
    sigma_gamma_init = sqrt(vp$gamma$prior_mean)
  )
  
  # --- run Gibbs sampler (C++) -----------------------------------------
  if (verbose) message(sprintf("Gibbs sampling: niter=%d, burn-in=%d",
                               niter, floor(niter * burnin_prop)))
  res <- gibbs_semo_uhp_only_cpp(
    niter = niter, K = K,
    beta_1_tk = start_val$beta_1_tk,
    beta_2_tk = start_val$beta_2_tk,
    theta_1_tk = start_val$theta_1_tk,
    theta_2_tk = start_val$theta_2_tk,
    gamma_tk = start_val$gamma_tk,
    sigma2_gamma_tk = start_val$sigma2_gamma_tk,
    Gamma_hat_1 = b_out[, 1], Gamma_hat_2 = b_out[, 2],
    s2_hat_Gamma_1 = se_out[, 1]^2, s2_hat_Gamma_2 = se_out[, 2]^2,
    gamma_hat = b_exp, s2_hat_gamma = se_exp^2,
    a_gamma = vp$gamma$a, b_gamma = vp$gamma$b,
    Sigma_theta_init = vp$theta$prior_mean,
    nu_theta = vp$theta$nu,
    Phi_theta = vp$theta$Phi
  )
  
  # --- post-burnin summary (no label switching to handle) --------------
  burnin <- floor(niter * burnin_prop)
  ids <- (burnin + 1):niter
  s1 <- get_summary(res$beta_1_tk[ids])
  s2 <- get_summary(res$beta_2_tk[ids])
  
  if (verbose) {
    cat("\n--- Results (semo) ---\n")
    cat("Outcome 1:\n");  print_summary(s1)
    cat("\nOutcome 2:\n"); print_summary(s2)
  }
  
  list(est  = c(s1$beta_est,  s2$beta_est),
       se   = c(s1$beta_se,   s2$beta_se),
       pval = c(s1$beta_pval, s2$beta_pval),
       ci   = rbind(s1$ci_emp, s2$ci_emp),
       model = "semo", n_iv = K)
}

# Internal worker: memo (2 exposures, 2 outcomes; UHP + CHP)
fit_memo <- function(b_exp, se_exp, b_out, se_out,
                     control, verbose = FALSE) {
  niter <- control$niter
  burnin_prop <- control$burnin_prop
  
  message("Running model: memo")
  
  if (!is.matrix(b_exp) || ncol(b_exp) != 2)
    stop("Model 'memo' requires b_exp and se_exp to be K x 2 matrices.")
  if (!is.matrix(b_out) || ncol(b_out) != 2)
    stop("Model 'memo' requires b_out and se_out to be K x 2 matrices.")
  
  K <- nrow(b_exp)
  if (K < 5) warning("Fewer than 5 IVs provided.")
  
  # --- compute MoM priors (bivariate exposure x bivariate outcome) -----
  vp <- set_variance_priors_m2x2_diag(
    ghat_mat = b_exp, gse_mat = se_exp,
    Ghat_mat = b_out, Gse_mat = se_out,
    B0 = NULL, K = K,
    Kmin = control$Kmin, Kmax = control$Kmax,
    rho12 = 0, rho_gg = 0,
    rho_gj = list(c(control$rho_ov, control$rho_ov),
                  c(control$rho_ov, control$rho_ov)),
    c_gamma = control$c_gamma, c_theta = control$c_theta,
    global_mean_gamma = control$global_mean_gamma,
    global_mean_theta = control$global_mean_theta,
    hybrid = control$hybrid, kappa_hybrid = control$kappa_hybrid,
    z_thresh = control$z_thresh, trim = control$trim,
    kappa_gamma = control$kappa_gamma, kappa_theta = control$kappa_theta
  )
  
  # --- initialize MCMC traces ------------------------------------------
  pst_init <- c(0.7, 0.1, 0.1, 0.1)   # mostly valid, small CHP probabilities
  start_val <- init_setup_memo(
    niter = niter, K = K,
    alpha_1_init = 0.1, alpha_2_init = 0.1,
    beta_1_init = vp$B0[1], beta_2_init = vp$B0[2],
    eta_1_init = rep(0, K), eta_2_init = rep(0, K),
    pst_init = pst_init
  )
  
  # --- run Gibbs sampler (C++) -----------------------------------------
  if (verbose) message(sprintf("Gibbs sampling: niter=%d, burn-in=%d",
                               niter, floor(niter * burnin_prop)))
  res <- gibbs_memo_joint_cpp(
    niter = niter, K = K,
    beta_1_tk = start_val$beta_1_tk,
    beta_2_tk = start_val$beta_2_tk,
    alpha_1_tk = start_val$alpha_1_tk,
    alpha_2_tk = start_val$alpha_2_tk,
    eta_1_tk = start_val$eta_1_tk,
    eta_2_tk = start_val$eta_2_tk,
    theta_1_tk = start_val$theta_1_tk,
    theta_2_tk = start_val$theta_2_tk,
    gamma_1_tk = start_val$gamma_1_tk,
    gamma_2_tk = start_val$gamma_2_tk,
    pst_tk = start_val$pst_tk,
    Gamma_hat_1 = b_out[, 1], gamma_hat_1 = b_exp[, 1],
    s2_hat_Gamma_1 = se_out[, 1]^2, s2_hat_gamma_1 = se_exp[, 1]^2,
    Gamma_hat_2 = b_out[, 2], gamma_hat_2 = b_exp[, 2],
    s2_hat_Gamma_2 = se_out[, 2]^2, s2_hat_gamma_2 = se_exp[, 2]^2,
    Sigma_gamma_init = vp$gamma$prior_mean,
    Sigma_theta_init = vp$theta$prior_mean,
    m_gamma = vp$gamma$nu, V_gamma = vp$gamma$Phi,
    m_theta = vp$theta$nu, V_theta = vp$theta$Phi,
    cc = c(1, 1, 1, 1)
  )
  
  # --- post-processing with joint label-switching correction -----------
  flip <- label_flip_joint(niter, res)
  pval1 <- 2 * (1 - stats::pnorm(abs(flip$b1_mean / flip$b1_sd)))
  pval2 <- 2 * (1 - stats::pnorm(abs(flip$b2_mean / flip$b2_sd)))
  
  if (verbose) {
    cat("\n--- Results (memo) ---\n")
    cat("Outcome 1:\n")
    cat(sprintf("  est = %.4f, se = %.4f, pval = %.4g\n",
                flip$b1_mean, flip$b1_sd, pval1))
    cat(sprintf("  95%% CI: [%.4f, %.4f]\n", flip$bci1[1], flip$bci1[2]))
    cat(sprintf("  CHP proportion (q1): %.3f%s\n",
                flip$qq1, if (flip$qq1 > 0.5) " (labels flipped)" else ""))
    cat("Outcome 2:\n")
    cat(sprintf("  est = %.4f, se = %.4f, pval = %.4g\n",
                flip$b2_mean, flip$b2_sd, pval2))
    cat(sprintf("  95%% CI: [%.4f, %.4f]\n", flip$bci2[1], flip$bci2[2]))
    cat(sprintf("  CHP proportion (q2): %.3f%s\n",
                flip$qq2, if (flip$qq2 > 0.5) " (labels flipped)" else ""))
  }
  
  list(est  = c(flip$b1_mean, flip$b2_mean),
       se   = c(flip$b1_sd,   flip$b2_sd),
       pval = c(pval1, pval2),
       ci   = rbind(flip$bci1, flip$bci2),
       q    = c(flip$qq1, flip$qq2),
       model = "memo", n_iv = K)
}