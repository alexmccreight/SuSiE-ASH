// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
int omp_threads() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return -1; // indicates no openmp
#endif
}
