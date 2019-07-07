cppInit <- function() {
  library(Rcpp)
  library(RcppArmadillo)
  Rcpp::cppFunction(depends = "RcppArmadillo",code = "
    arma::mat mvrnormArma(int n, arma::mat mu, arma::mat sigma) {
      int ncols = sigma.n_cols;
      arma::mat Y = arma::randn(n, ncols);
      return mu + Y * arma::chol(sigma);
    }",env = .GlobalEnv)
  #   Rcpp::sourceCpp("
  # #include <RcppArmadillo.h>
  # // [[Rcpp::depends(RcppArmadillo)]]
  # 
  # using namespace Rcpp;
  # 
  # // [[Rcpp::export]]
  NULL
}