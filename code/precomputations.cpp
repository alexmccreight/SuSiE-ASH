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

// [[Rcpp::export]]
Rcpp::List initialize_susie(const arma::mat& V, const arma::vec& Dsq, const arma::vec& VtXty,
                            double tausq, double sigmasq, int L,
                            Rcpp::Nullable<Rcpp::NumericVector> ssq_ = R_NilValue,
                            Rcpp::Nullable<Rcpp::NumericMatrix> PIP_ = R_NilValue,
                            Rcpp::Nullable<Rcpp::NumericMatrix> mu_ = R_NilValue,
                            Rcpp::Nullable<Rcpp::NumericVector> pi0_ = R_NilValue) {
  int p = V.n_rows;

  arma::vec var = tausq * Dsq + sigmasq;
  arma::vec diagXtOmegaX = arma::sum(arma::square(V) % arma::repmat(Dsq / var, 1, p).t(), 1);
  arma::vec XtOmegay = V * (VtXty / var);

  arma::vec ssq;
  if (ssq_.isNotNull()) {
    ssq = Rcpp::as<arma::vec>(ssq_);
  } else {
    ssq = arma::vec(L, arma::fill::value(0.2));
  }

  arma::mat PIP;
  if (PIP_.isNotNull()) {
    PIP = Rcpp::as<arma::mat>(PIP_);
  } else {
    PIP = arma::mat(p, L, arma::fill::value(1.0 / p));
  }

  arma::mat mu;
  if (mu_.isNotNull()) {
    mu = Rcpp::as<arma::mat>(mu_);
  } else {
    mu = arma::mat(p, L, arma::fill::zeros);
  }

  arma::mat lbf_variable(p, L, arma::fill::zeros);
  arma::vec lbf(L, arma::fill::zeros);
  arma::mat omega = arma::repmat(diagXtOmegaX, 1, L) + 1.0 / arma::repmat(ssq.t(), p, 1);

  arma::vec logpi0;
  if (pi0_.isNull()) {
    logpi0 = arma::vec(p, arma::fill::value(std::log(1.0 / p)));
  } else {
    Rcpp::NumericVector pi0 = Rcpp::as<Rcpp::NumericVector>(pi0_);
    logpi0 = arma::vec(p, arma::fill::value(-arma::datum::inf));
    for (int i = 0; i < p; ++i) {
      if (pi0[i] > 0) {
        logpi0[i] = std::log(pi0[i]);
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("var") = var,
    Rcpp::Named("diagXtOmegaX") = diagXtOmegaX,
    Rcpp::Named("XtOmegay") = XtOmegay,
    Rcpp::Named("ssq") = ssq,
    Rcpp::Named("PIP") = PIP,
    Rcpp::Named("mu") = mu,
    Rcpp::Named("lbf_variable") = lbf_variable,
    Rcpp::Named("lbf") = lbf,
    Rcpp::Named("omega") = omega,
    Rcpp::Named("logpi0") = logpi0
  );
}

double brent_min(double a, double b, const std::function<double(double)>& f, double tol = 1e-8, int max_iter = 100) {
  const double golden_ratio = 0.3819660;
  double x, w, v, fx, fw, fv;
  x = w = v = a + golden_ratio * (b - a);
  fx = fw = fv = f(x);
  double d = 0.0, e = 0.0;

  for (int iter = 0; iter < max_iter; ++iter) {
    double xm = 0.5 * (a + b);
    double tol1 = tol * std::abs(x) + 1e-10;
    double tol2 = 2.0 * tol1;

    if (std::abs(x - xm) <= (tol2 - 0.5 * (b - a))) {
      return x;
    }

    if (std::abs(e) > tol1) {
      double r = (x - w) * (fx - fv);
      double q = (x - v) * (fx - fw);
      double p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > 0.0) p = -p;
      q = std::abs(q);
      double etemp = e;
      e = d;

      if (std::abs(p) >= std::abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)) {
        e = (x >= xm) ? a - x : b - x;
        d = golden_ratio * e;
      } else {
        d = p / q;
        double u = x + d;
        if (u - a < tol2 || b - u < tol2) {
          d = (xm - x >= 0) ? tol1 : -tol1;
        }
      }
    } else {
      e = (x >= xm) ? a - x : b - x;
      d = golden_ratio * e;
    }

    double u = (std::abs(d) >= tol1) ? x + d : x + ((d >= 0) ? tol1 : -tol1);
    double fu = f(u);

    if (fu <= fx) {
      if (u >= x) a = x; else b = x;
      v = w; w = x; x = u;
      fv = fw; fw = fx; fx = fu;
    } else {
      if (u < x) a = u; else b = u;
      if (fu <= fw || w == x) {
        v = w; w = u;
        fv = fw; fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u;
        fv = fu;
      }
    }
  }
  return x;
}

// Helper function for log_sum_exp
double log_sum_exp(const arma::vec& x) {
  double max_x = x.max();
  return max_x + std::log(arma::sum(arma::exp(x - max_x)));
}

// Helper function for optimization
double f(double x, const arma::vec& diagXtOmegaX, const arma::vec& XtOmegar, const arma::vec& logpi0) {
  arma::vec terms = -0.5 * arma::log(1 + x * diagXtOmegaX) +
    x * arma::square(XtOmegar) / (2 * (1 + x * diagXtOmegaX)) +
    logpi0;
  return -log_sum_exp(terms);
}

// [[Rcpp::export]]
Rcpp::List single_effect_regression(const arma::mat& V, const arma::vec& Dsq, const arma::vec& var,
                                    arma::vec& diagXtOmegaX, arma::vec& XtOmegay,
                                    arma::mat& PIP, arma::mat& mu, arma::vec& ssq,
                                    arma::mat& omega, arma::mat& lbf_variable, arma::vec& lbf,
                                    const arma::vec& logpi0, bool est_ssq,
                                    double ssq_lower, double ssq_upper, bool verbose) {
  int p = V.n_rows;
  int L = PIP.n_cols;

  for (int l = 0; l < L; ++l) {
    // Compute X', Omega r_l for residual r_l
    arma::vec b = arma::sum(mu % PIP, 1) - mu.col(l) % PIP.col(l);
    arma::vec XtOmegaXb = V * ((V.t() * b) % (Dsq / var));
    arma::vec XtOmegar = XtOmegay - XtOmegaXb;

    // Update Prior Variance ssq[l]
    if (est_ssq) {
      auto objective = [&](double x) {
        return f(x, diagXtOmegaX, XtOmegar, logpi0);
      };

      double tol = 1e-8;
      int max_iter = 100;

      // Use our implemented Brent's method for optimization
      double ssq_new = brent_min(ssq_lower, ssq_upper, objective, tol, max_iter);

      ssq[l] = ssq_new;

      if (verbose) {
        Rcpp::Rcout << "Update s^2 for effect " << l + 1 << " to " << ssq[l] << std::endl;
      }
    }

    // Update omega, mu, and PIP
    omega.col(l) = diagXtOmegaX + 1 / ssq[l];
    mu.col(l) = XtOmegar / omega.col(l);
    lbf_variable.col(l) = arma::square(XtOmegar) / (2 * omega.col(l)) - 0.5 * arma::log(omega.col(l) * ssq[l]);
    arma::vec logPIP = lbf_variable.col(l) + logpi0;
    lbf[l] = log_sum_exp(logPIP);
    PIP.col(l) = arma::exp(logPIP - lbf[l]);
  }

  return Rcpp::List::create(
    Rcpp::Named("PIP") = PIP,
    Rcpp::Named("mu") = mu,
    Rcpp::Named("ssq") = ssq,
    Rcpp::Named("omega") = omega,
    Rcpp::Named("lbf_variable") = lbf_variable,
    Rcpp::Named("lbf") = lbf
  );
}

// [[Rcpp::export]]
arma::vec calculate_diagXtOmegaX(const arma::mat& V, const arma::vec& Dsq, const arma::vec& var) {
  int p = V.n_rows;
  arma::vec diagXtOmegaX(p);

  arma::vec Dsq_div_var = Dsq / var;

  for (int i = 0; i < p; ++i) {
    diagXtOmegaX(i) = arma::dot(arma::square(V.row(i)), Dsq_div_var);
  }

  return diagXtOmegaX;
}

// [[Rcpp::export]]
arma::vec calculate_XtOmegay(const arma::mat& V, const arma::vec& VtXty, const arma::vec& var) {
  return V * (VtXty / var);
}


// [[Rcpp::export]]
arma::vec calculate_XtOmegaXb(const arma::mat& V, const arma::vec& b, const arma::vec& Dsq, const arma::vec& var) {
  // Compute t(V) %*% b
  arma::vec VtB = V.t() * b;

  // Element-wise multiplication with Dsq and division by var
  VtB %= Dsq / var;

  // Final matrix multiplication
  return V * VtB;
}

// [[Rcpp::export]]
Rcpp::List MoM_rcpp(const arma::mat& PIP, const arma::mat& mu, const arma::mat& omega,
                    double sigmasq, double tausq, int n, const arma::mat& V,
                    const arma::vec& Dsq, const arma::vec& VtXty, const arma::vec& Xty,
                    double yty, bool est_sigmasq, bool est_tausq, bool verbose) {
  int p = mu.n_rows;
  int L = mu.n_cols;

  // Compute A
  arma::mat A(2, 2, arma::fill::zeros);
  A(0, 0) = n;
  A(0, 1) = A(1, 0) = arma::accu(Dsq);
  A(1, 1) = arma::accu(arma::square(Dsq));

  // Compute diag(V'MV)
  arma::vec b = arma::sum(mu % PIP, 1);
  arma::vec Vtb = V.t() * b;
  arma::vec diagVtMV = arma::square(Vtb);
  arma::vec tmpD(p, arma::fill::zeros);

  for (int l = 0; l < L; ++l) {
    arma::vec bl = mu.col(l) % PIP.col(l);
    arma::vec Vtbl = V.t() * bl;
    diagVtMV -= arma::square(Vtbl);
    tmpD += PIP.col(l) % (arma::square(mu.col(l)) + 1.0 / omega.col(l));
  }

  diagVtMV += arma::sum(arma::square(V.t()) % arma::repmat(tmpD.t(), p, 1), 1);

  // Compute x
  arma::vec x(2);
  x(0) = yty - 2 * arma::dot(b, Xty) + arma::dot(Dsq, diagVtMV);
  x(1) = arma::accu(arma::square(Xty)) - 2 * arma::dot(Vtb % VtXty, Dsq) + arma::dot(arma::square(Dsq), diagVtMV);

  // Solve system of equations
  if (est_tausq) {
    arma::vec sol = arma::solve(A, x);
    if (sol(0) > 0 && sol(1) > 0) {
      sigmasq = sol(0);
      tausq = sol(1);
    } else {
      sigmasq = x(0) / n;
      tausq = 0;
    }
    if (verbose) {
      Rcpp::Rcout << "Update (sigma^2,tau^2) to (" << sigmasq << "," << tausq << ")\n";
    }
  } else if (est_sigmasq) {
    sigmasq = (x(0) - A(0, 1) * tausq) / n;
    if (verbose) {
      Rcpp::Rcout << "Update sigma^2 to " << sigmasq << "\n";
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("sigmasq") = sigmasq,
    Rcpp::Named("tausq") = tausq
  );
}
