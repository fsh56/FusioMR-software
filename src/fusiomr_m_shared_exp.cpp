#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Single exposure, multiple outcomes (2). UHP only via covariance Sigma_theta.

// rinvgamma_cpp is defined in fusiomr_s_uhp_only.cpp
NumericVector rinvgamma_cpp(int n, double shape, double rate);

// Robust multivariate-normal sampler with Cholesky + eigen fallback.
static arma::mat mvrnorm_robust(int n, const arma::vec& mu, arma::mat sigma) {
  int p = sigma.n_cols;
  arma::mat Y = arma::randn(n, p);
  sigma = 0.5 * (sigma + sigma.t());                 // symmetrize
  arma::mat L;
  bool ok = arma::chol(L, sigma, "lower");
  if (!ok) {
    arma::vec eigval; arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, sigma);
    eigval = arma::clamp(eigval, 1e-8, arma::datum::inf);
    L = eigvec * arma::diagmat(arma::sqrt(eigval));
  }
  return arma::repmat(mu, 1, n).t() + Y * L.t();
}

// Inverse-Wishart sampler (pure C++ via Armadillo)
static arma::mat rinvwishart_arma(double nu, const arma::mat& S) {
  return arma::iwishrnd(S, nu);
}

// [[Rcpp::export]]
List gibbs_semo_uhp_only_cpp(int niter, int K,
                             arma::vec beta_1_tk,
                             arma::vec beta_2_tk,
                             NumericMatrix theta_1_tk,
                             NumericMatrix theta_2_tk,
                             NumericMatrix gamma_tk,
                             NumericVector sigma2_gamma_tk,
                             arma::vec Gamma_hat_1,
                             arma::vec Gamma_hat_2,
                             arma::vec s2_hat_Gamma_1,
                             arma::vec s2_hat_Gamma_2,
                             arma::vec gamma_hat,
                             arma::vec s2_hat_gamma,
                             double a_gamma, double b_gamma,
                             arma::mat Sigma_theta_init,
                             double nu_theta, arma::mat Phi_theta) {

  // initialize current values from first row of trace
  arma::vec beta_cur(2);
  beta_cur[0] = beta_1_tk[0];
  beta_cur[1] = beta_2_tk[0];
  NumericVector gamma_cur = gamma_tk(0, _);
  arma::mat theta_cur(K, 2);
  for (int k = 0; k < K; k++) {
    theta_cur(k, 0) = theta_1_tk(0, k);
    theta_cur(k, 1) = theta_2_tk(0, k);
  }
  double sigma2_gamma_cur = sigma2_gamma_tk[0];
  arma::mat Sigma_theta_cur = Sigma_theta_init;

  // pre-allocate diagnostic traces for outcome covariance diagonals
  NumericVector sigma2_theta1_tk(niter, 0.0);
  NumericVector sigma2_theta2_tk(niter, 0.0);
  sigma2_theta1_tk[0] = Sigma_theta_init(0, 0);
  sigma2_theta2_tk[0] = Sigma_theta_init(1, 1);

  for (int iter = 0; iter < (niter - 1); iter++) {

    // update gamma_k | rest (one shared exposure effect, two outcomes)
    for (int k = 0; k < K; k++) {
      double Ak = beta_cur[0] * beta_cur[0] / s2_hat_Gamma_1[k]
                + beta_cur[1] * beta_cur[1] / s2_hat_Gamma_2[k]
                + 1.0 / s2_hat_gamma[k]
                + 1.0 / sigma2_gamma_cur;
      double Bk = beta_cur[0] * (Gamma_hat_1[k] - theta_cur(k, 0)) / s2_hat_Gamma_1[k]
                + beta_cur[1] * (Gamma_hat_2[k] - theta_cur(k, 1)) / s2_hat_Gamma_2[k]
                + gamma_hat[k] / s2_hat_gamma[k];
      gamma_cur[k] = R::rnorm(Bk / Ak, sqrt(1.0 / Ak));
    }
    gamma_tk(iter + 1, _) = gamma_cur;

    // update theta_k | rest (2D pleiotropy with covariance Sigma_theta)
    arma::mat Sigma_theta_inv;
    bool ok = arma::inv_sympd(Sigma_theta_inv, Sigma_theta_cur);
    if (!ok) {
      Sigma_theta_cur.diag() += 1e-8;
      Sigma_theta_inv = arma::inv_sympd(Sigma_theta_cur);
    }
    for (int k = 0; k < K; k++) {
      arma::mat S_k_inv(2, 2, arma::fill::zeros);
      S_k_inv(0, 0) = 1.0 / std::max(s2_hat_Gamma_1[k], 1e-12);
      S_k_inv(1, 1) = 1.0 / std::max(s2_hat_Gamma_2[k], 1e-12);
      arma::mat prec_post = S_k_inv + Sigma_theta_inv;
      arma::vec res = {Gamma_hat_1[k] - beta_cur[0] * gamma_cur[k],
                       Gamma_hat_2[k] - beta_cur[1] * gamma_cur[k]};
      arma::mat cov_post = arma::inv_sympd(prec_post);
      arma::vec mu_post = cov_post * S_k_inv * res;
      arma::mat sample = mvrnorm_robust(1, mu_post, cov_post);
      theta_cur(k, 0) = sample(0, 0);
      theta_cur(k, 1) = sample(0, 1);
    }

    // update sigma^2_gamma | rest
    double sumg2 = arma::sum(arma::square(arma::vec(gamma_cur)));
    sigma2_gamma_cur = rinvgamma_cpp(1, a_gamma + K / 2.0,
                                     b_gamma + 0.5 * sumg2)[0];
    if (sigma2_gamma_cur < 1e-4)  sigma2_gamma_cur = 1e-4;
    if (sigma2_gamma_cur > 10.0)  sigma2_gamma_cur = 10.0;
    sigma2_gamma_tk[iter + 1] = sigma2_gamma_cur;

    // update Sigma_theta | rest (Inverse-Wishart conjugate)
    arma::mat S_theta(2, 2, arma::fill::zeros);
    for (int k = 0; k < K; k++) {
      arma::rowvec rk = theta_cur.row(k);
      S_theta += rk.t() * rk;
    }
    double nu_post = nu_theta + K;
    arma::mat Phi_post = S_theta + Phi_theta;
    Phi_post = 0.5 * (Phi_post + Phi_post.t());
    // numerical safety: clip eigenvalues
    arma::vec eigval; arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, Phi_post);
    eigval = arma::clamp(eigval, 1e-8, arma::datum::inf);
    Phi_post = eigvec * arma::diagmat(eigval) * eigvec.t();
    arma::mat Sigma_sample = rinvwishart_arma(nu_post, Phi_post);
    // bound diagonals and constrain correlation
    double s11 = std::max(1e-8, std::min(Sigma_sample(0, 0), 10.0));
    double s22 = std::max(1e-8, std::min(Sigma_sample(1, 1), 10.0));
    double rho = Sigma_sample(0, 1) / std::sqrt(s11 * s22);
    rho = std::max(-0.99, std::min(0.99, rho));
    Sigma_theta_cur(0, 0) = s11;
    Sigma_theta_cur(1, 1) = s22;
    Sigma_theta_cur(0, 1) = Sigma_theta_cur(1, 0) = rho * std::sqrt(s11 * s22);
    sigma2_theta1_tk[iter + 1] = s11;
    sigma2_theta2_tk[iter + 1] = s22;

    // update beta_1 and beta_2 | rest (independent given gamma, theta)
    arma::vec gamma_vec(K);
    for (int k = 0; k < K; k++) gamma_vec[k] = gamma_cur[k];

    double A1 = 0.0, B1 = 0.0, A2 = 0.0, B2 = 0.0;
    for (int k = 0; k < K; k++) {
      double w = gamma_vec[k];
      A1 += w * w / s2_hat_Gamma_1[k];
      A2 += w * w / s2_hat_Gamma_2[k];
      B1 += w * (Gamma_hat_1[k] - theta_cur(k, 0)) / s2_hat_Gamma_1[k];
      B2 += w * (Gamma_hat_2[k] - theta_cur(k, 1)) / s2_hat_Gamma_2[k];
    }
    A1 = std::max(A1, 1e-8);
    A2 = std::max(A2, 1e-8);
    beta_cur[0] = R::rnorm(B1 / A1, sqrt(1.0 / A1));
    beta_cur[1] = R::rnorm(B2 / A2, sqrt(1.0 / A2));
    beta_1_tk[iter + 1] = beta_cur[0];
    beta_2_tk[iter + 1] = beta_cur[1];

    // store theta traces
    for (int k = 0; k < K; k++) {
      theta_1_tk(iter + 1, k) = theta_cur(k, 0);
      theta_2_tk(iter + 1, k) = theta_cur(k, 1);
    }
  }

  return List::create(
    Named("K") = K,
    Named("beta_1_tk") = beta_1_tk,
    Named("beta_2_tk") = beta_2_tk,
    Named("gamma_tk") = gamma_tk,
    Named("theta_1_tk") = theta_1_tk,
    Named("theta_2_tk") = theta_2_tk,
    Named("sigma2_gamma_tk") = sigma2_gamma_tk,
    Named("sigma2_theta1_tk") = sigma2_theta1_tk,
    Named("sigma2_theta2_tk") = sigma2_theta2_tk
  );
}
