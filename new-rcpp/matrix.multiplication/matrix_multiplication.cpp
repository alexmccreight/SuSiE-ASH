// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
SEXP parallelXtOmegaEigen(
    const Eigen::Map<Eigen::MatrixXd> &V,    // (n x p)
    const Eigen::Map<Eigen::MatrixXd> &VtXt, // (p x m)
    const Eigen::Map<Eigen::VectorXd> &var,  // length p
    int n_cores = 1
) {
  // 1) Set the number of threads for Eigen's internal parallelization
  Eigen::setNbThreads(n_cores);

  // 2) Create a temporary copy of VtXt for row-scaling
  Eigen::MatrixXd scaled = VtXt;  // (p x m)

  // 3) Scale each row i by 1/var(i)
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < scaled.rows(); i++) {
    scaled.row(i) /= var(i);
  }

  // 4) Now multiply: (n x p) * (p x m) => (n x m)
  Eigen::MatrixXd XtOmega = V * scaled;

  // 5) Return result to R
  return Rcpp::wrap(XtOmega);
}
