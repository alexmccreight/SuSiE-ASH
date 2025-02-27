// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
SEXP parallelXtOmegaEigen(
    const Eigen::Map<Eigen::MatrixXd> &V,    // (p x p)
    const Eigen::Map<Eigen::MatrixXd> &VtXt, // (p x n)
    const Eigen::Map<Eigen::VectorXd> &var,  // length p
    int n_cores = 1
) {

  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd scaled = VtXt;  // (p x n)

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < scaled.rows(); i++) {
    scaled.row(i) /= var(i);
  }

  Eigen::MatrixXd XtOmega = V * scaled;

  return Rcpp::wrap(XtOmega);
}
