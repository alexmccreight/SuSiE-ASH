#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Function to compute X'X
// [[Rcpp::export]]
arma::mat compute_XtX(const arma::mat& X) {
  return X.t() * X;
}

// Function to compute X'y
// [[Rcpp::export]]
arma::vec compute_Xty(const arma::mat& X, const arma::vec& y) {
  return X.t() * y;
}

// Function to compute eigendecomposition of a symmetric matrix
// [[Rcpp::export]]
Rcpp::List compute_eigen_sym(const arma::mat& X) {
  arma::vec eigval;
  arma::mat eigvec;

  arma::eig_sym(eigval, eigvec, X);

  return Rcpp::List::create(
    Rcpp::Named("values") = eigval,
    Rcpp::Named("vectors") = eigvec
  );
}

// Function to compute V'X'y
// [[Rcpp::export]]
arma::vec compute_VtXty(const arma::mat& V, const arma::vec& Xty) {
  return V.t() * Xty;
}

// Function to compute diagXtOmegaX
// [[Rcpp::export]]
arma::vec compute_diagXtOmegaX(const arma::mat& V, const arma::vec& Dsq, const arma::vec& var) {
  int p = V.n_rows;
  arma::vec diagXtOmegaX(p);
  arma::vec Dsq_div_var = Dsq / var;

  for (int i = 0; i < p; ++i) {
    diagXtOmegaX(i) = arma::dot(arma::square(V.row(i)), Dsq_div_var);
  }

  return diagXtOmegaX;
}

// Function to compute XtOmegay
// [[Rcpp::export]]
arma::vec compute_XtOmegay(const arma::mat& V, const arma::vec& VtXty, const arma::vec& var) {
  return V * (VtXty / var);
}


