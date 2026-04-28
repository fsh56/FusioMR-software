#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Multiple exposures, multiple outcomes (2x2). Joint model with UHP (theta)
// and CHP (eta, alpha) on each of the two outcomes; bivariate gamma and
// theta share inverse-Wishart priors.

// rinvgamma_cpp is defined in fusiomr_s_uhp_only.cpp
NumericVector rinvgamma_cpp(int n, double shape, double rate);

// Robust multivariate-normal sampler (Cholesky + eigen fallback).
static arma::mat mvrnorm_robust_memo(int n, const arma::vec& mu, arma::mat sigma) {
  int p = sigma.n_cols;
  arma::mat Y = arma::randn(n, p);
  sigma = 0.5 * (sigma + sigma.t());
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
static arma::mat rinvwishart_arma_memo(double nu, const arma::mat& S) {
  return arma::iwishrnd(S, nu);
}

// Dirichlet sampler via independent Gammas + normalization (n=1 case).
static arma::vec rdirichlet_arma(const arma::vec& alpha) {
  int K = alpha.n_elem;
  arma::vec x(K);
  for (int i = 0; i < K; i++) x[i] = R::rgamma(alpha[i], 1.0);
  return x / arma::sum(x);
}

// [[Rcpp::export]]
List gibbs_memo_joint_cpp(int niter, int K,
                          arma::vec beta_1_tk, arma::vec beta_2_tk,
                          arma::vec alpha_1_tk, arma::vec alpha_2_tk,
                          NumericMatrix eta_1_tk, NumericMatrix eta_2_tk,
                          NumericMatrix theta_1_tk, NumericMatrix theta_2_tk,
                          NumericMatrix gamma_1_tk, NumericMatrix gamma_2_tk,
                          NumericMatrix pst_tk,
                          arma::vec Gamma_hat_1, arma::vec gamma_hat_1,
                          arma::vec s2_hat_Gamma_1, arma::vec s2_hat_gamma_1,
                          arma::vec Gamma_hat_2, arma::vec gamma_hat_2,
                          arma::vec s2_hat_Gamma_2, arma::vec s2_hat_gamma_2,
                          arma::mat Sigma_gamma_init, arma::mat Sigma_theta_init,
                          double m_gamma, arma::mat V_gamma,
                          double m_theta, arma::mat V_theta,
                          arma::vec cc) {

  // --- initialize current values from first row of trace -----------------
  arma::mat theta_cur(K, 2);
  arma::mat gamma_cur(K, 2);
  arma::mat eta_cur(K, 2);
  for (int k = 0; k < K; k++) {
    theta_cur(k, 0) = theta_1_tk(0, k);
    theta_cur(k, 1) = theta_2_tk(0, k);
    gamma_cur(k, 0) = gamma_1_tk(0, k);
    gamma_cur(k, 1) = gamma_2_tk(0, k);
    eta_cur(k, 0)   = eta_1_tk(0, k);
    eta_cur(k, 1)   = eta_2_tk(0, k);
  }
  arma::vec alpha_cur(2);
  alpha_cur[0] = alpha_1_tk[0]; alpha_cur[1] = alpha_2_tk[0];
  arma::vec beta_cur(2);
  beta_cur[0]  = beta_1_tk[0];  beta_cur[1]  = beta_2_tk[0];
  arma::mat Sigma_gamma_cur = Sigma_gamma_init;
  arma::mat Sigma_theta_cur = Sigma_theta_init;
  arma::vec pst_cur(4);
  for (int j = 0; j < 4; j++) pst_cur[j] = pst_tk(0, j);

  // diagnostic traces
  NumericMatrix q_post_1_tk(niter, K);
  NumericMatrix q_post_2_tk(niter, K);
  NumericVector sigma2_gamma1_tk(niter, 0.0);
  NumericVector sigma2_gamma2_tk(niter, 0.0);
  NumericVector sigma2_theta1_tk(niter, 0.0);
  NumericVector sigma2_theta2_tk(niter, 0.0);
  sigma2_gamma1_tk[0] = Sigma_gamma_init(0, 0);
  sigma2_gamma2_tk[0] = Sigma_gamma_init(1, 1);
  sigma2_theta1_tk[0] = Sigma_theta_init(0, 0);
  sigma2_theta2_tk[0] = Sigma_theta_init(1, 1);

  for (int iter = 0; iter < (niter - 1); iter++) {

    // pre-compute inverses of current covariance matrices
    arma::mat Sigma_gamma_inv, Sigma_theta_inv;
    if (!arma::inv_sympd(Sigma_gamma_inv, Sigma_gamma_cur)) {
      Sigma_gamma_cur.diag() += 1e-8;
      Sigma_gamma_inv = arma::inv_sympd(Sigma_gamma_cur);
    }
    if (!arma::inv_sympd(Sigma_theta_inv, Sigma_theta_cur)) {
      Sigma_theta_cur.diag() += 1e-8;
      Sigma_theta_inv = arma::inv_sympd(Sigma_theta_cur);
    }

    // --- update gamma_k and theta_k for each IV --------------------------
    for (int k = 0; k < K; k++) {
      arma::mat A_k(2, 2, arma::fill::zeros);
      A_k(0, 0) = beta_cur[0] + alpha_cur[0] * eta_cur(k, 0);
      A_k(1, 1) = beta_cur[1] + alpha_cur[1] * eta_cur(k, 1);

      arma::mat S_G_inv(2, 2, arma::fill::zeros);
      S_G_inv(0, 0) = 1.0 / std::max(s2_hat_Gamma_1[k], 1e-12);
      S_G_inv(1, 1) = 1.0 / std::max(s2_hat_Gamma_2[k], 1e-12);

      arma::mat S_g_inv(2, 2, arma::fill::zeros);
      S_g_inv(0, 0) = 1.0 / std::max(s2_hat_gamma_1[k], 1e-12);
      S_g_inv(1, 1) = 1.0 / std::max(s2_hat_gamma_2[k], 1e-12);

      // gamma_k | rest
      arma::mat prec_gamma = A_k.t() * S_G_inv * A_k + S_g_inv + Sigma_gamma_inv;
      arma::vec res_G = {Gamma_hat_1[k] - theta_cur(k, 0),
                         Gamma_hat_2[k] - theta_cur(k, 1)};
      arma::vec res_g = {gamma_hat_1[k], gamma_hat_2[k]};
      arma::vec mu_num_gamma = A_k.t() * S_G_inv * res_G + S_g_inv * res_g;
      arma::mat cov_gamma = arma::inv_sympd(prec_gamma);
      arma::vec mu_gamma  = cov_gamma * mu_num_gamma;
      arma::mat g_sample  = mvrnorm_robust_memo(1, mu_gamma, cov_gamma);
      gamma_cur(k, 0) = g_sample(0, 0);
      gamma_cur(k, 1) = g_sample(0, 1);

      // theta_k | rest
      arma::mat prec_theta = S_G_inv + Sigma_theta_inv;
      arma::vec res_theta = {
        Gamma_hat_1[k] - beta_cur[0] * gamma_cur(k, 0) - alpha_cur[0] * gamma_cur(k, 0) * eta_cur(k, 0),
        Gamma_hat_2[k] - beta_cur[1] * gamma_cur(k, 1) - alpha_cur[1] * gamma_cur(k, 1) * eta_cur(k, 1)
      };
      arma::vec mu_num_theta = S_G_inv * res_theta;
      arma::mat cov_theta = arma::inv_sympd(prec_theta);
      arma::vec mu_theta  = cov_theta * mu_num_theta;
      arma::mat t_sample  = mvrnorm_robust_memo(1, mu_theta, cov_theta);
      theta_cur(k, 0) = t_sample(0, 0);
      theta_cur(k, 1) = t_sample(0, 1);
    }

    // --- update Sigma_gamma | rest (Inverse-Wishart) ---------------------
    arma::mat S_gamma(2, 2, arma::fill::zeros);
    for (int k = 0; k < K; k++) {
      arma::rowvec rk = gamma_cur.row(k);
      S_gamma += rk.t() * rk;
    }
    double m_post_gamma   = m_gamma + K;
    arma::mat V_post_gamma = S_gamma + V_gamma;
    V_post_gamma = 0.5 * (V_post_gamma + V_post_gamma.t());
    {
      arma::vec eigval; arma::mat eigvec;
      arma::eig_sym(eigval, eigvec, V_post_gamma);
      eigval = arma::clamp(eigval, 1e-8, arma::datum::inf);
      V_post_gamma = eigvec * arma::diagmat(eigval) * eigvec.t();
    }
    arma::mat Sgs = rinvwishart_arma_memo(m_post_gamma, V_post_gamma);
    double s11g = std::max(1e-4, std::min(Sgs(0, 0), 10.0));
    double s22g = std::max(1e-4, std::min(Sgs(1, 1), 10.0));
    double rhog = Sgs(0, 1) / std::sqrt(s11g * s22g);
    rhog = std::max(-0.99, std::min(0.99, rhog));
    Sigma_gamma_cur(0, 0) = s11g;
    Sigma_gamma_cur(1, 1) = s22g;
    Sigma_gamma_cur(0, 1) = Sigma_gamma_cur(1, 0) = rhog * std::sqrt(s11g * s22g);
    sigma2_gamma1_tk[iter + 1] = s11g;
    sigma2_gamma2_tk[iter + 1] = s22g;

    // --- update Sigma_theta | rest (Inverse-Wishart) ---------------------
    arma::mat S_theta(2, 2, arma::fill::zeros);
    for (int k = 0; k < K; k++) {
      arma::rowvec rk = theta_cur.row(k);
      S_theta += rk.t() * rk;
    }
    double m_post_theta   = m_theta + K;
    arma::mat V_post_theta = S_theta + V_theta;
    V_post_theta = 0.5 * (V_post_theta + V_post_theta.t());
    {
      arma::vec eigval; arma::mat eigvec;
      arma::eig_sym(eigval, eigvec, V_post_theta);
      eigval = arma::clamp(eigval, 1e-8, arma::datum::inf);
      V_post_theta = eigvec * arma::diagmat(eigval) * eigvec.t();
    }
    arma::mat Sts = rinvwishart_arma_memo(m_post_theta, V_post_theta);
    double s11t = std::max(1e-8, std::min(Sts(0, 0), 10.0));
    double s22t = std::max(1e-8, std::min(Sts(1, 1), 10.0));
    double rhot = Sts(0, 1) / std::sqrt(s11t * s22t);
    rhot = std::max(-0.99, std::min(0.99, rhot));
    Sigma_theta_cur(0, 0) = s11t;
    Sigma_theta_cur(1, 1) = s22t;
    Sigma_theta_cur(0, 1) = Sigma_theta_cur(1, 0) = rhot * std::sqrt(s11t * s22t);
    sigma2_theta1_tk[iter + 1] = s11t;
    sigma2_theta2_tk[iter + 1] = s22t;

    // --- update eta_k (CHP indicators) with hard 0.1/0.9 thresholding ---
    for (int k = 0; k < K; k++) {
      // eta_k1
      double r1_1 = Gamma_hat_1[k] - theta_cur(k, 0)
                  - (beta_cur[0] + alpha_cur[0]) * gamma_cur(k, 0);
      double r1_0 = Gamma_hat_1[k] - theta_cur(k, 0)
                  -  beta_cur[0] * gamma_cur(k, 0);
      // joint with eta_k2 via pst (joint distribution over 4 cells)
      double w1_1 = (eta_cur(k, 1) == 0) ? pst_cur[2] : pst_cur[3];
      double w1_0 = (eta_cur(k, 1) == 0) ? pst_cur[0] : pst_cur[1];
      double f1_1 = std::exp(-0.5 * r1_1 * r1_1 / s2_hat_Gamma_1[k]) * w1_1;
      double f1_0 = std::exp(-0.5 * r1_0 * r1_0 / s2_hat_Gamma_1[k]) * w1_0;
      double q1_post = f1_1 / (f1_1 + f1_0);
      if (q1_post <= 0.1) q1_post = 0.0;
      if (q1_post >= 0.9) q1_post = 1.0;
      eta_cur(k, 0) = R::rbinom(1, q1_post);
      q_post_1_tk(iter + 1, k) = q1_post;
      eta_1_tk(iter + 1, k) = eta_cur(k, 0);

      // eta_k2
      double r2_1 = Gamma_hat_2[k] - theta_cur(k, 1)
                  - (beta_cur[1] + alpha_cur[1]) * gamma_cur(k, 1);
      double r2_0 = Gamma_hat_2[k] - theta_cur(k, 1)
                  -  beta_cur[1] * gamma_cur(k, 1);
      double w2_1 = (eta_cur(k, 0) == 0) ? pst_cur[1] : pst_cur[3];
      double w2_0 = (eta_cur(k, 0) == 0) ? pst_cur[0] : pst_cur[2];
      double f2_1 = std::exp(-0.5 * r2_1 * r2_1 / s2_hat_Gamma_2[k]) * w2_1;
      double f2_0 = std::exp(-0.5 * r2_0 * r2_0 / s2_hat_Gamma_2[k]) * w2_0;
      double q2_post = f2_1 / (f2_1 + f2_0);
      if (q2_post <= 0.1) q2_post = 0.0;
      if (q2_post >= 0.9) q2_post = 1.0;
      eta_cur(k, 1) = R::rbinom(1, q2_post);
      q_post_2_tk(iter + 1, k) = q2_post;
      eta_2_tk(iter + 1, k) = eta_cur(k, 1);
    }

    // --- update pst (Dirichlet conjugate) -------------------------------
    int n00 = 0, n01 = 0, n10 = 0, n11 = 0;
    for (int k = 0; k < K; k++) {
      int e1 = (int) eta_cur(k, 0);
      int e2 = (int) eta_cur(k, 1);
      if      (e1 == 0 && e2 == 0) n00++;
      else if (e1 == 0 && e2 == 1) n01++;
      else if (e1 == 1 && e2 == 0) n10++;
      else if (e1 == 1 && e2 == 1) n11++;
    }
    arma::vec dir_par = {n00 + cc[0], n01 + cc[1], n10 + cc[2], n11 + cc[3]};
    pst_cur = rdirichlet_arma(dir_par);
    for (int j = 0; j < 4; j++) pst_tk(iter + 1, j) = pst_cur[j];

    // --- update alpha_j and beta_j --------------------------------------
    // outcome 1
    double A_a1 = 0.0, B_a1 = 0.0;
    double A_b1 = 0.0, B_b1 = 0.0;
    for (int k = 0; k < K; k++) {
      double Wa = eta_cur(k, 0) * gamma_cur(k, 0);
      double Ua = Gamma_hat_1[k] - theta_cur(k, 0) - beta_cur[0] * gamma_cur(k, 0);
      A_a1 += Wa * Wa / s2_hat_Gamma_1[k];
      B_a1 += Wa * Ua / s2_hat_Gamma_1[k];
      double Wb = gamma_cur(k, 0);
      double Ub = Gamma_hat_1[k] - theta_cur(k, 0) - alpha_cur[0] * eta_cur(k, 0) * gamma_cur(k, 0);
      A_b1 += Wb * Wb / s2_hat_Gamma_1[k];
      B_b1 += Wb * Ub / s2_hat_Gamma_1[k];
    }
    if (A_a1 < 1e-12) alpha_cur[0] = 0.0;
    else              alpha_cur[0] = R::rnorm(B_a1 / A_a1, sqrt(1.0 / A_a1));
    A_b1 = std::max(A_b1, 1e-8);
    beta_cur[0] = R::rnorm(B_b1 / A_b1, sqrt(1.0 / A_b1));

    // outcome 2
    double A_a2 = 0.0, B_a2 = 0.0;
    double A_b2 = 0.0, B_b2 = 0.0;
    for (int k = 0; k < K; k++) {
      double Wa = eta_cur(k, 1) * gamma_cur(k, 1);
      double Ua = Gamma_hat_2[k] - theta_cur(k, 1) - beta_cur[1] * gamma_cur(k, 1);
      A_a2 += Wa * Wa / s2_hat_Gamma_2[k];
      B_a2 += Wa * Ua / s2_hat_Gamma_2[k];
      double Wb = gamma_cur(k, 1);
      double Ub = Gamma_hat_2[k] - theta_cur(k, 1) - alpha_cur[1] * eta_cur(k, 1) * gamma_cur(k, 1);
      A_b2 += Wb * Wb / s2_hat_Gamma_2[k];
      B_b2 += Wb * Ub / s2_hat_Gamma_2[k];
    }
    if (A_a2 < 1e-12) alpha_cur[1] = 0.0;
    else              alpha_cur[1] = R::rnorm(B_a2 / A_a2, sqrt(1.0 / A_a2));
    A_b2 = std::max(A_b2, 1e-8);
    beta_cur[1] = R::rnorm(B_b2 / A_b2, sqrt(1.0 / A_b2));

    alpha_1_tk[iter + 1] = alpha_cur[0];
    alpha_2_tk[iter + 1] = alpha_cur[1];
    beta_1_tk[iter + 1]  = beta_cur[0];
    beta_2_tk[iter + 1]  = beta_cur[1];

    // store gamma/theta traces for this iteration
    for (int k = 0; k < K; k++) {
      gamma_1_tk(iter + 1, k) = gamma_cur(k, 0);
      gamma_2_tk(iter + 1, k) = gamma_cur(k, 1);
      theta_1_tk(iter + 1, k) = theta_cur(k, 0);
      theta_2_tk(iter + 1, k) = theta_cur(k, 1);
    }
  }

  return List::create(
    Named("K") = K,
    Named("alpha_1_tk") = alpha_1_tk,
    Named("alpha_2_tk") = alpha_2_tk,
    Named("beta_1_tk") = beta_1_tk,
    Named("beta_2_tk") = beta_2_tk,
    Named("eta_1_tk") = eta_1_tk,
    Named("eta_2_tk") = eta_2_tk,
    Named("theta_1_tk") = theta_1_tk,
    Named("theta_2_tk") = theta_2_tk,
    Named("gamma_1_tk") = gamma_1_tk,
    Named("gamma_2_tk") = gamma_2_tk,
    Named("pst_tk") = pst_tk,
    Named("q_post_1_tk") = q_post_1_tk,
    Named("q_post_2_tk") = q_post_2_tk,
    Named("sigma2_gamma1_tk") = sigma2_gamma1_tk,
    Named("sigma2_gamma2_tk") = sigma2_gamma2_tk,
    Named("sigma2_theta1_tk") = sigma2_theta1_tk,
    Named("sigma2_theta2_tk") = sigma2_theta2_tk
  );
}
