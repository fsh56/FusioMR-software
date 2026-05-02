#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>


using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
void fastLm(const arma::vec & y, const arma::mat & X, double &coefb, double &stdb) {
  
  int n = X.n_rows, k = X.n_cols;
  
  arma::colvec coef = arma::solve(X, y); 
  arma::colvec resid = y - X*coef; 
  
  double sig2 = arma::as_scalar(arma::trans(resid)*resid/(n-k));
  arma::colvec stderrest = 
    arma::sqrt(sig2 * arma::diagvec( arma::inv(arma::trans(X)*X)) );
  
  coefb = coef[1];
  stdb = stderrest[1];
}

// [[Rcpp::export]]
List fastSigLm(const arma::vec & y, const arma::mat & X) {
  
  // int n = X.n_rows, k = X.n_cols;
  int p = X.n_cols;int n = X.n_rows;
  arma::mat xx = zeros(p, 2);
  double coefb = 0;
  double stdb = 0;
  vec coef = zeros(p, 1);
  vec std = zeros(p, 1);
  
  for( int j = 0; j < p; j = j + 1 )
  {
    xx = join_rows(ones(n, 1), X.col(j));
    fastLm(y, xx, coefb, stdb);
    coef[j] = coefb;
    std[j] = stdb;
  }
  
  return List::create(Named("coef") = coef,
                      Named("std") = std);
}

