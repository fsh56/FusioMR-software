#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Single exposure, single outcome, UHP only (no CHP)
// Priors (a_gamma, b_gamma, a_theta, b_theta) and MCMC trace vectors
// are pre-computed in R and passed in.

// Inverse-gamma sampler (uses the R package `invgamma`)
NumericVector rinvgamma_cpp(int n, double shape, double rate) {
  Function ff("rinvgamma", Environment::namespace_env("invgamma"));
  return ff(n, Named("shape") = shape, _["rate"] = rate);
}

// [[Rcpp::export]]
List gibbs_seso_uhp_only_cpp(int niter, int K,
                             NumericVector beta_tk,
                             NumericMatrix theta_tk,
                             NumericMatrix gamma_tk,
                             NumericVector sigma2_gamma_tk,
                             NumericVector sigma2_theta_tk,
                             NumericVector Gamma_hat,
                             NumericVector gamma_hat,
                             NumericVector s2_hat_Gamma,
                             NumericVector s2_hat_gamma,
                             double a_gamma, double b_gamma,
                             double a_theta, double b_theta) {

  // initialize current values from first row/element of trace
  double beta_cur = beta_tk[0];
  NumericVector theta_cur = theta_tk(0, _);
  NumericVector gamma_cur = gamma_tk(0, _);
  double sigma2_gamma_cur = sigma2_gamma_tk[0];
  double sigma2_theta_cur = sigma2_theta_tk[0];

  // Gibbs iterations
  for (int iter = 0; iter < (niter - 1); iter++) {

    // update gamma_k | rest
    for (int k = 0; k < K; k++) {
      double Ak = beta_cur * beta_cur / s2_hat_Gamma[k]
                + 1.0 / s2_hat_gamma[k]
                + 1.0 / sigma2_gamma_cur;
      double Bk = beta_cur * (Gamma_hat[k] - theta_cur[k]) / s2_hat_Gamma[k]
                + gamma_hat[k] / s2_hat_gamma[k];
      gamma_cur[k] = rnorm(1, Bk / Ak, sqrt(1.0 / Ak))[0];
    }
    gamma_tk(iter + 1, _) = gamma_cur;

    // update theta_k | rest
    for (int k = 0; k < K; k++) {
      double Ak = 1.0 / s2_hat_Gamma[k] + 1.0 / sigma2_theta_cur;
      double Bk = (Gamma_hat[k] - beta_cur * gamma_cur[k]) / s2_hat_Gamma[k];
      theta_cur[k] = rnorm(1, Bk / Ak, sqrt(1.0 / Ak))[0];
    }
    theta_tk(iter + 1, _) = theta_cur;

    // update sigma^2_gamma and sigma^2_theta | rest
    sigma2_gamma_cur = rinvgamma_cpp(1,
                        a_gamma + K / 2.0,
                        b_gamma + 0.5 * sum(pow(gamma_cur, 2)))[0];
    sigma2_theta_cur = rinvgamma_cpp(1,
                        a_theta + K / 2.0,
                        b_theta + 0.5 * sum(pow(theta_cur, 2)))[0];
    // bound away from 0 and huge values for numerical stability
    if (sigma2_gamma_cur < 1e-10) sigma2_gamma_cur = 1e-10;
    if (sigma2_gamma_cur > 10)    sigma2_gamma_cur = 10;
    if (sigma2_theta_cur < 1e-10) sigma2_theta_cur = 1e-10;
    if (sigma2_theta_cur > 10)    sigma2_theta_cur = 10;
    sigma2_gamma_tk[iter + 1] = sigma2_gamma_cur;
    sigma2_theta_tk[iter + 1] = sigma2_theta_cur;

    // update beta | rest (GLS on outcome residuals)
    vec U = Gamma_hat - theta_cur;
    vec W = gamma_cur;
    mat Omega(K, K, fill::zeros);
    for (int i = 0; i < K; i++) Omega(i, i) = 1.0 / s2_hat_Gamma[i];

    double A_beta = std::max(as_scalar(trans(W) * Omega * W), 1e-8);
    double mu_beta = as_scalar(trans(W) * Omega * U) / A_beta;
    beta_cur = rnorm(1, mu_beta, sqrt(1.0 / A_beta))[0];
    beta_tk[iter + 1] = beta_cur;
  }

  return List::create(
    Named("K") = K,
    Named("beta_tk") = beta_tk,
    Named("gamma_tk") = gamma_tk,
    Named("theta_tk") = theta_tk,
    Named("sigma2_gamma_tk") = sigma2_gamma_tk,
    Named("sigma2_theta_tk") = sigma2_theta_tk
  );
}
