#include "mr_ash.h"

// FUNCTION DEFINITIONS
// --------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List caisa_rcpp       (const arma::mat& X, const arma::vec& y,
                             const arma::vec& w, const arma::vec& sa2,
                             arma::vec& pi, arma::vec& beta,
                             arma::vec& r, double sigma2, const arma::uvec& o,
                             int maxiter, int miniter,
                             double convtol, double epstol, std::string method_q,
                             bool updatepi, bool updatesigma,
                             bool verbose,
                             const arma::mat& XtOmega, //CHANGED arguments
                             double tausq,
                             double sum_Dsq,
                             const arma::mat& V, //ADDED
                             const arma::vec& Dsq,
                             const arma::mat& VtXt) { //ADDED

  // See "mr_ash.h"

  // ---------------------------------------------------------------------
  // DEFINE SIZES
  // ---------------------------------------------------------------------
  int n                   = X.n_rows;
  int p                   = X.n_cols;
  int K                   = sa2.n_elem;

  // ---------------------------------------------------------------------
  // PREDEFINE LOCAL VARIABLES
  // ---------------------------------------------------------------------
  arma::vec varobj(maxiter);
  int iter               = 0;
  int i                  = 0;
  int j;

  double a1;
  double a2;
  double bjj; //CHANGED
  arma::vec piold;
  arma::vec betaold;
  arma::mat XtOmegaT;


  // ---------------------------------------------------------------------
  // PRECALCULATE
  // ---------------------------------------------------------------------
  arma::mat S2inv        = 1 / outerAddition(sa2, w); //CHAGED 1 / outerAddition(1/sa2, w)
  S2inv.row(0).fill(epstol);

  //CHANGED (added XtOmegat before main loop)
  XtOmegaT = XtOmega.t();

  // Initialize a p x K matrix filled with zeros
  arma::mat prod_mat(p, K, arma::fill::zeros);
  // Initialize a p x K matrix and p x 1 vector filled with zeros;
  arma::vec post_var_vec(p, arma::fill::zeros);//CHANGED
  arma::vec b_vec(p,arma::fill::zeros );//CHANGED

  // ---------------------------------------------------------------------
  // START LOOP : CYCLE THROUGH COORDINATE ASCENT UPDATES
  // ---------------------------------------------------------------------
  for (iter = 0; iter < maxiter; iter++) {

    //if(iter = 1){}

    // (C) Print iteration info
    //Rcpp::Rcout << "\n--- iteration = " << iter << " ---\n";

    // reset parameters
    a1                   = 0;
    a2                   = 0;
    bjj                  = 0; //CHANGED
    piold                = pi;
    pi.fill(0);
    betaold              = beta;

    // initialize b

    // (D) Inspect p, i, and loop logic
    //Rcpp::Rcout << "Looping over j=0..(p-1) with p=" << p << "\n";
    //Rcpp::Rcout << "i (outer index) starts at " << i << "\n";

    arma::vec phij;
    arma::vec muj;
    arma::vec post_var; //CHANGED


    // ---------------------------------------------------------------------
    // RUN COORDINATE ASCENT UPDATES : INDEX 1 - INDEX P CHANGED!  XtOmega.t().col(o(i)) added!
    // ---------------------------------------------------------------------
    for (j = 0; j < p; j++){

      unsigned int oi = o(i); // e.g. 0..(p-1)

      updatebetaj(X.col(o(i)),
                  w(o(i)),
                  beta(o(i)),
                  r,
                  piold,
                  pi,
                  sigma2,
                  sa2,
                  S2inv.col(o(i)),
                  a1,
                  a2,
                  o(i),
                  p,
                  epstol,
                  XtOmegaT.col(o(i)),
                  phij,
                  muj,
                  post_var, //CHANGED
                  bjj); // CHANGED (A)

      // Compute the element-wise product and assign to the j-th row of prod_mat
      prod_mat.row(o(i)) = (phij % muj).t();
      post_var_vec(o(i)) = arma::dot(phij, ((muj % muj) + post_var));
      b_vec(o(i)) = bjj;

      i++;

    }


    if (updatesigma) { //CHANGED

      double term1 = arma::dot(r, r);
      double term2 = arma::dot(beta, beta) * (n - 1) ;
      double term3 = a1 * (n - 1);
      double term4 = tausq * sum_Dsq;
      arma::mat M = arma::square(V.t()); //## CHECK!!
      M.each_col() %= (post_var_vec - (beta % beta));
      arma::vec row_sum = arma::sum(M, /*dim=*/1);
      double m1  = dot(Dsq, row_sum) + term1;
      arma::vec Dsqsq = arma::square(Dsq);
      double term6 = dot(Dsqsq, row_sum);
      arma::vec temp = X.t() * r;
      double rtXXtr = arma::sum(arma::square(temp));
      double sum_Dsq2 = arma::sum(arma::square(Dsq));
      double m2 = term6 + rtXXtr;

      // Compute sigma2
      sigma2 = ((sum_Dsq2 * m1) - (sum_Dsq * m2) )/((n * sum_Dsq2) - (sum_Dsq * sum_Dsq));
      arma::vec var = tausq * Dsq + sigma2;

      // Scale each row of VtXt by 1/var[i]
      arma::mat scaled_VtXt = VtXt;  // Copy VtXt
      #pragma omp parallel for
      for (int i = 0; i < scaled_VtXt.n_rows; ++i) {
        scaled_VtXt.row(i) /= var(i);
      }

      //scaled_VtXt.each_row() /= var.t();  // Armadillo vectorized operation (no OpenMP needed)

      // Compute XtOmega = V * scaled_VtXt
      arma::mat XtOmega = V * scaled_VtXt;

      // Transpose to get XtOmegaT for the next iteration
      XtOmegaT = XtOmega.t();

    }

    if (updatepi) {
      // if updatepi == true, we update pi
      piold               = pi;
    }

    // ---------------------------------------------------------------------
    // CALCULATE VARIATIONAL OBJECTIVE 2: Changed Muted !
    // ---------------------------------------------------------------------
    //varobj(iter)          = varobj(iter) / sigma2 / 2.0 +
    //                        log(2.0 * M_PI * sigma2) / 2.0 * n -
    //                        dot(pi, log(piold + epstol)) * p + a2;

    // for (j = 1; j < K; j++){
    //   varobj(iter)       += pi(j) * log(sa2(j)) * p / 2;
    // }

    if (!updatepi) {
      // if updatepi == false, we do not update pi
      pi                  = piold;
    }

    // ---------------------------------------------------------------------
    // CHECK CONVERGENCE
    // ---------------------------------------------------------------------
    if (iter >= miniter - 1) {
      if (arma::norm(betaold - beta) < convtol * p) {
        iter++;
        break;
      }

     // if (iter > 0) {
    //    if (varobj(iter) > varobj(iter - 1)){
    //      break;
    //    }
    //  }
    }
  }

  if (verbose) {
    Rprintf("Mr.ASH terminated at iteration %d.\n", iter);
  }
  // ---------------------------------------------------------------------
  // RETURN VALUES
  // ---------------------------------------------------------------------
  return Rcpp::List::create(Rcpp::Named("beta")    = beta,
                            Rcpp::Named("sigma2")  = sigma2,
                            Rcpp::Named("pi")      = pi,
                            Rcpp::Named("iter")    = iter,
                            Rcpp::Named("b_vec") = b_vec,
                            Rcpp::Named("post_var") = post_var_vec
                            //Rcpp::Named("varobj")  = varobj.subvec(0,iter-1)
                              );
}
