#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List precomputations(const arma::mat& X, const arma::vec& y) {
    int n = X.n_rows;
    int p = X.n_cols;
    
    // Compute z
    arma::vec z = X.t() * y / std::sqrt(n);
    
    // Compute meansq
    double meansq = arma::accu(arma::square(y)) / n;
    
    // Compute XtX and LD
    arma::mat XtX = X.t() * X;
    arma::mat LD = XtX / n;
    
    // Compute eigen decomposition
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, LD);
    
    // Reverse order of eigenvectors and values
    arma::mat V = arma::fliplr(eigvec);
    arma::vec Dsq = arma::reverse(arma::max(n * eigval, arma::zeros(p)));
    
    // Compute Xty, VtXty, yty
    arma::vec Xty = std::sqrt(n) * z;
    arma::vec VtXty = V.t() * Xty;
    double yty = n * meansq;
    
    return Rcpp::List::create(
        Rcpp::Named("z") = z,
        Rcpp::Named("meansq") = meansq,
        Rcpp::Named("XtX") = XtX,
        Rcpp::Named("LD") = LD,
        Rcpp::Named("V") = V,
        Rcpp::Named("Dsq") = Dsq,
        Rcpp::Named("Xty") = Xty,
        Rcpp::Named("VtXty") = VtXty,
        Rcpp::Named("yty") = yty
    );
}