#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Single exposure, single outcome, with both UHP (theta) and CHP (eta, alpha).
// Priors and MCMC trace vectors are pre-computed in R and passed in.

// Inverse-gamma sampler shared with seso_uhp_only (defined in that file).
// Forward-declared here so the linker can find it.
NumericVector rinvgamma_cpp(int n, double shape, double rate);

// [[Rcpp::export]]
List gibbs_seso_with_chp_cpp(int niter, int K,
                             NumericVector beta_tk,
                             NumericVector alpha_tk,
                             NumericMatrix eta_tk,
                             NumericMatrix theta_tk,
                             NumericMatrix gamma_tk,
                             NumericVector q_tk,
                             NumericVector Gamma_hat,
                             NumericVector gamma_hat,
                             NumericVector s2_hat_Gamma,
                             NumericVector s2_hat_gamma,
                             NumericVector sigma2_gamma_tk,
                             NumericVector sigma2_theta_tk,
                             double a_gamma, double b_gamma,
                             double a_theta, double b_theta,
                             double a_q, double b_q) {

  // initialize current values from first row/element of trace
  double beta_cur  = beta_tk[0];
  double alpha_cur = alpha_tk[0];
  NumericVector theta_cur = theta_tk(0, _);
  NumericVector gamma_cur = gamma_tk(0, _);
  NumericVector eta_cur   = eta_tk(0, _);
  double sigma2_gamma_cur = sigma2_gamma_tk[0];
  double sigma2_theta_cur = sigma2_theta_tk[0];
  double q_cur = q_tk[0];

  // posterior probability that eta_k = 1, recorded for diagnostics
  NumericMatrix q_post_tk(niter, K);

  for (int iter = 0; iter < (niter - 1); iter++) {

    // update gamma_k | rest
    for (int k = 0; k < K; k++) {
      double bk = beta_cur + alpha_cur * eta_cur[k];
      double Ak = bk * bk / s2_hat_Gamma[k]
                + 1.0 / s2_hat_gamma[k]
                + 1.0 / sigma2_gamma_cur;
      double Bk = bk * (Gamma_hat[k] - theta_cur[k]) / s2_hat_Gamma[k]
                + gamma_hat[k] / s2_hat_gamma[k];
      gamma_cur[k] = rnorm(1, Bk / Ak, sqrt(1.0 / Ak))[0];
    }
    gamma_tk(iter + 1, _) = gamma_cur;

    // update theta_k | rest (UHP)
    for (int k = 0; k < K; k++) {
      double bk = beta_cur + alpha_cur * eta_cur[k];
      double Ak = 1.0 / s2_hat_Gamma[k] + 1.0 / sigma2_theta_cur;
      double Bk = (Gamma_hat[k] - bk * gamma_cur[k]) / s2_hat_Gamma[k];
      theta_cur[k] = rnorm(1, Bk / Ak, sqrt(1.0 / Ak))[0];
    }
    theta_tk(iter + 1, _) = theta_cur;

    // update sigma^2_gamma and sigma^2_theta | rest
    sigma2_gamma_cur = rinvgamma_cpp(1, a_gamma + K / 2.0,
                        b_gamma + 0.5 * sum(pow(gamma_cur, 2)))[0];
    sigma2_theta_cur = rinvgamma_cpp(1, a_theta + K / 2.0,
                        b_theta + 0.5 * sum(pow(theta_cur, 2)))[0];
    sigma2_gamma_tk[iter + 1] = sigma2_gamma_cur;
    sigma2_theta_tk[iter + 1] = sigma2_theta_cur;

    // update eta_k | rest (CHP indicator) with hard 0.1/0.9 thresholding
    for (int k = 0; k < K; k++) {
      double r1 = Gamma_hat[k] - theta_cur[k] - (beta_cur + alpha_cur) * gamma_cur[k];
      double r0 = Gamma_hat[k] - theta_cur[k] -  beta_cur                * gamma_cur[k];
      double f1 = exp(-0.5 * r1 * r1 / s2_hat_Gamma[k]) * q_cur;
      double f0 = exp(-0.5 * r0 * r0 / s2_hat_Gamma[k]) * (1.0 - q_cur);
      double q_post = f1 / (f1 + f0);
      // hard truncation: stabilizes labels and reduces flipping
      if (q_post <= 0.1) q_post = 0.0;
      if (q_post >= 0.9) q_post = 1.0;
      eta_cur[k] = rbinom(1, 1, q_post)[0];
      q_post_tk(iter + 1, k) = q_post;
    }
    eta_tk(iter + 1, _) = eta_cur;

    // update q | rest (Beta-Binomial conjugate)
    int n0 = sum(eta_cur == 0);
    int n1 = sum(eta_cur == 1);
    q_cur = rbeta(1, n1 + a_q, n0 + b_q)[0];
    q_tk[iter + 1] = q_cur;

    // precision matrix for outcome (diagonal 1/se^2)
    mat Omega(K, K, fill::zeros);
    for (int i = 0; i < K; i++) Omega(i, i) = 1.0 / s2_hat_Gamma[i];

    // update alpha | rest (only IVs with eta=1 contribute)
    vec U_alpha = Gamma_hat - theta_cur - beta_cur * gamma_cur;
    vec W_alpha = eta_cur * gamma_cur;
    if (sum(W_alpha) == 0) {
      alpha_cur = 0.0;
    } else {
      double A_alpha = as_scalar(trans(W_alpha) * Omega * W_alpha);
      double mu_alpha = as_scalar(trans(W_alpha) * Omega * U_alpha) / A_alpha;
      alpha_cur = rnorm(1, mu_alpha, sqrt(1.0 / A_alpha))[0];
    }
    alpha_tk[iter + 1] = alpha_cur;

    // update beta | rest
    vec U_beta = Gamma_hat - theta_cur - alpha_cur * eta_cur * gamma_cur;
    vec W_beta = gamma_cur;
    double A_beta  = std::max(as_scalar(trans(W_beta) * Omega * W_beta), 1e-8);
    double mu_beta = as_scalar(trans(W_beta) * Omega * U_beta) / A_beta;
    beta_cur = rnorm(1, mu_beta, sqrt(1.0 / A_beta))[0];
    beta_tk[iter + 1] = beta_cur;
  }

  return List::create(
    Named("K") = K,
    Named("alpha_tk") = alpha_tk,
    Named("beta_tk") = beta_tk,
    Named("eta_tk") = eta_tk,
    Named("q_tk") = q_tk,
    Named("gamma_tk") = gamma_tk,
    Named("theta_tk") = theta_tk,
    Named("sigma2_gamma_tk") = sigma2_gamma_tk,
    Named("sigma2_theta_tk") = sigma2_theta_tk,
    Named("q_post_tk") = q_post_tk
  );
}
