#' @title SuSiE-Infinitesimal
#'
#' @description
#' Implements the SuSiE-Infinitesimal model (Cui et al 2024), a mixture of sparse single effects
#' and an infinitesimal component, to estimate posterior inclusion probabilities
#' and effect sizes. This function updates posterior distributions of effect sizes,
#' inclusion probabilities, and variance parameters, either via method-of-moments
#' or maximum likelihood estimation (still needs to be implemented).
#'
#' This function supports two modes of data input:
#' 1) **Individual-level data mode**: Provide `X` (n x p, scaled & centered) and `y` (n-vector, scaled & centered).
#'    The function will compute `z`, `meansq`, `XtX`, `LD`, `V`, and `Dsq` as needed.
#' 2) **Summary-level data mode**: Provide `z` (p-vector), `meansq` (scalar), and `n` (scalar),
#'    along with either `LD` or (`V`, `Dsq`). In this mode, `X` and `y` are not required.
#'
#' @param X A numeric n x p matrix of covariates (scaled & centered) for individual-level mode.
#'
#' @param y A numeric response vector of length n (centered) for individual-level mode.
#'
#' @param z A numeric vector of length p giving summary-level marginal z-scores
#' or effects for summary-level mode.
#'
#' @param meansq A numeric value giving sum(y^2)/n for summary-level mode. If not
#' provided in individual-level mode, it will be computed from `y`.
#'
#' @param n An integer giving the sample size. Required if using summary-level mode.
#'
#' @param L An integer specifying the number of modeled causal effects.
#'
#' @param est_ssq Logical; if TRUE, estimate the prior effect size variances \eqn{s^2}
#'  for each effect using MLE.
#'
#' @param ssq A numeric vector of length L giving initial values of \eqn{s^2} for each effect.
#'   Default is 0.2 for every effect if not provided.
#'
#' @param ssq_range A numeric vector of length 2 providing the lower and upper
#' bounds for \eqn{s^2}.
#'
#' @param pi0 A numeric vector of length p giving the prior causal probability for each SNP.
#'   These must sum to 1. Default is equal probability (1/p) for each SNP.
#'
#' @param est_sigmasq Logical; if TRUE, estimate the residual variance \eqn{\sigma^2}.
#'
#' @param est_tausq Logical; if TRUE, estimate both \eqn{\sigma^2} and \eqn{\tau^2}.
#'
#' @param sigmasq A numeric value giving the initial value for \eqn{\sigma^2}.
#'
#' @param tausq A numeric value giving the initial value for \eqn{\tau^2}.
#'
#' @param method A character string indicating the method used to estimate
#'   (\eqn{\sigma^2}, \eqn{\tau^2}). Options are:
#'   \itemize{
#'     \item \code{"moments"}: Method-of-moments estimation
#'     \item \code{"MLE"}: Maximum likelihood estimation. Not currently implemented.
#'   }
#'
#' @param sigmasq_range A numeric vector of length 2 giving the lower and upper bounds for \eqn{\sigma^2}
#'   if estimated using MLE.
#'
#' @param tausq_range A numeric vector of length 2 giving the lower and upper bounds for \eqn{\tau^2}
#'   if estimated using MLE.
#'
#' @param PIP A p x L matrix of initial posterior inclusion probabilities.
#' Default: 1/p for each SNP-effect pair.
#'
#' @param mu A p x L matrix of initial posterior means conditional on causality.
#' Default: 0 for each SNP-effect pair.
#'
#' @param maxiter An integer specifying the maximum number of SuSiE iterations.
#'
#' @param PIP_tol A numeric threshold for convergence.
#'
#' @param coverage A numeric value between 0 and 1 indicating the coverage level
#'  for credible sets.
#'
#' @param verbose Logical; if TRUE, print progress and updates at each iteration.
#'
#' @param XtX Optional p x p precomputed matrix \eqn{X'X}.
#'
#' @param LD Optional p x p LD matrix, equal to \eqn{X'X/n}. If not provided, it can be computed or
#'   derived from \eqn{V} and \eqn{Dsq}.
#'
#' @param V Optional p x p matrix of eigenvectors of \eqn{X'X}.
#'
#' @param Dsq Optional numeric vector of length p giving the eigenvalues of \eqn{X'X}.
#'   If \code{LD} is not provided, both \code{V} and \code{Dsq} must be provided.
#'
#'
#'
#' @return A list with elements:
#'
#' \item{PIP}{A p x L matrix of posterior inclusion probabilities for each SNP-effect pair.}
#'
#' \item{marginal_PIP}{A length-p vector of marginal posterior inclusion probabilities
#'  computed across all effects.}
#'
#' \item{mu}{A p x L matrix of posterior mean effect sizes conditional on inclusion.}
#'
#' \item{omega}{A p x L matrix of posterior precisions conditional on causality.}
#'
#' \item{lbf}{A length-L vector of log-Bayes-factors for each effect.}
#'
#' \item{lbf_variable}{A p x L matrix of per-variable log-Bayes-factors.}
#'
#' \item{ssq}{A length-L vector of final effect size variances \eqn{s^2}.}
#'
#' \item{sigmasq}{Final estimate of \eqn{\sigma^2}.}
#'
#' \item{tausq}{Final estimate of \eqn{\tau^2}.}
#'
#' \item{alpha}{A length-p vector of posterior means of infinitesimal effects.}
#'
#' \item{fitted}{A length-n vector of fitted values of the trait.}
#'
#' \item{sets}{A list of credible sets achieving the specified coverage.}
#'
#' @seealso \code{\link{MoM}}, \code{\link{susie_inf_get_cs}}
#'
#' @examples
#'
#' set.seed(1)
#' n <- 3000; p <- 3000
#' X <- matrix(rnorm(n*p), n, p)
#' X <- scale(X, center=TRUE, scale=TRUE)
#'
# Sparse effects
#' b_true <- rep(0, p)
#' b_true[sample(p, 5)] <- rnorm(5, 0, sqrt(0.01))
#'
# Add weak infinitesimal effects (tau^2 = 0.0001)
#' tausq_weak <- 0.0001
#' b_inf <- rnorm(p)*sqrt(tausq_weak)
#'
# Phenotype with both sparse and infinitesimal effects
#' y <- X %*% (b_true + b_inf) + rnorm(n)
#' y <- scale(y, scale = TRUE, center = TRUE)
#'
# Example 1: Individual-level data mode (no precomputations)
#' res_ind <- susie_inf(X=X, y=y, L=5, method="moments", est_tausq=TRUE)
#' print(sort(res_ind$marginal_PIP, decreasing = T)[1:10])
#' cat("Estimated sigma^2, tau^2: (", res_ind$sigmasq, ",", res_ind$tausq, ")\n")
#'
# Example 2: Summary-level data mode
# Compute z, meansq, and LD from the generated data, then use only summary stats
#' z <- drop(t(X) %*% y)/sqrt(n)
#' meansq <- sum(y^2)/n
#' LD <- crossprod(X)/n
#'
#' res_sum <- susie_inf(z=z, meansq=meansq, n=n, L=5, LD=LD, method="moments", est_tausq=TRUE)
#' print(sort(res_sum$marginal_PIP, decreasing = T)[1:10])
#' cat("Estimated sigma^2, tau^2: (", res_sum$sigmasq, ",", res_sum$tausq, ")\n")
#'
# Example 3: Individual-level data with precomputations
#' XtX <- crossprod(X)
#' eig <- eigen(LD, symmetric=TRUE)
#' idx <- order(eig$values, decreasing=TRUE)
#' V <- eig$vectors[, idx]
#' Dsq <- pmax(eig$values[idx]*n, 0)
#'
#' res_ind_precomp <- susie_inf(X=X, y=y, L=5, XtX=XtX, LD=LD, V=V, Dsq=Dsq, method="moments", est_tausq=TRUE)
#' print(sort(res_ind_precomp$marginal_PIP, decreasing = T)[1:10])
#' cat("Estimated sigma^2, tau^2: (", res_ind_precomp$sigmasq, ",", res_ind_precomp$tausq, ")\n")
#'
#' @export

susie_inf <- function(X = NULL, y = NULL, z = NULL, meansq = NULL, n = NULL,
                      L,
                      est_ssq = TRUE, ssq = NULL, ssq_range = c(0, 1), pi0 = NULL,
                      est_sigmasq = TRUE, est_tausq = TRUE, sigmasq = 1, tausq = 0,
                      method = "moments", sigmasq_range = NULL, tausq_range = NULL,
                      PIP = NULL, mu = NULL, maxiter = 100, PIP_tol = 1e-3, coverage = 0.9, verbose = TRUE,
                      XtX = NULL, LD = NULL, V = NULL, Dsq = NULL) {

  # Check inputs for mode
  have_individual_data <- !is.null(X) && !is.null(y)
  have_summary_data <- !is.null(z) && !is.null(meansq) && !is.null(n)

  if (have_individual_data && have_summary_data) {
    stop("Please provide either individual-level data (X,y) OR summary-level data (z,meansq,n), not both.")
  }

  if (!have_individual_data && !have_summary_data) {
    stop("Must provide either individual-level data (X,y) OR summary-level data (z,meansq,n).")
  }

  if (have_individual_data) {
    # Individual-level mode
    if (verbose) cat("Using individual-level data mode.\n")
    n <- nrow(X)
    p <- ncol(X)
    # Compute meansq if not given
    if (is.null(meansq)) {
      meansq <- sum(y^2)/n
    }
    # Compute z if not given
    if (is.null(z)) {
      z <- (t(X) %*% y)/sqrt(n)
    }
    # Compute XtX if not given
    if (is.null(XtX)) {
      XtX <- t(X) %*% X
    }
    # Compute LD if not given
    if (is.null(LD)) {
      LD <- XtX / n
    }
    # Compute V, Dsq if not given
    if (is.null(V) || is.null(Dsq)) {
      eig <- eigen(LD, symmetric = TRUE)
      # order decreasing
      idx <- order(eig$values, decreasing = TRUE)
      eig$values <- eig$values[idx]
      eig$vectors <- eig$vectors[, idx]
      V <- eig$vectors
      Dsq <- pmax(eig$values * n, 0)
    }
  } else {
    # Summary-level mode
    if (verbose) cat("Using summary-level data mode.\n")
    p <- length(z)
    # If LD not given, need V,Dsq
    if (is.null(LD) && (is.null(V) || is.null(Dsq))) {
      stop("In summary-data mode, if LD is not provided, V and Dsq must be provided.")
    }
    # If LD given but not V,Dsq, can compute them if needed
    if (!is.null(LD) && (is.null(V) || is.null(Dsq))) {
      eig <- eigen(LD, symmetric = TRUE)
      idx <- order(eig$values, decreasing = TRUE)
      eig$values <- eig$values[idx]
      eig$vectors <- eig$vectors[, idx]
      V <- eig$vectors
      Dsq <- pmax(eig$values * n, 0)
    }
  }

  Xty <- sqrt(n)*z
  VtXty <- t(V) %*% Xty
  yty <- n * meansq

  var <- tausq*Dsq + sigmasq
  diagXtOmegaX <- rowSums(sweep(V^2, 2, (Dsq / var), `*`))
  XtOmegay <- V %*% (VtXty / var)

  if (is.null(ssq)) {ssq <- rep(0.2, L)}
  if (is.null(PIP)) {PIP <- matrix(1/p, nrow = p, ncol = L)}
  if (is.null(mu)) {mu <- matrix(0, nrow = p, ncol = L)}

  lbf_variable <- matrix(0, nrow = p, ncol = L)
  lbf <- rep(0, L)
  omega <- matrix(diagXtOmegaX, nrow = p, ncol = L) + 1 / ssq

  if (is.null(pi0)) {
    logpi0 <- rep(log(1 / p), p)
  } else {
    logpi0 <- rep(-Inf, p)
    inds <- which(pi0 > 0)
    logpi0[inds] <- log(pi0[inds])
  }

  ##### Main SuSiE iteration loop #####
  ### Alternate between updating prior effect size variances, ssq,
  ### and posterior distribution, q(B^(l)), sequentially for each effect
  ### and updating the infinitesimal effect size and noise variances (sigmasq, tausq)

  for(it in seq_len(maxiter)) {
    if (verbose) {
      cat(sprintf("Iteration %d\n", it))
    }
    PIP_prev <- PIP

    # Single Effect Regression for each effect l = 1, ... , L
    for (l in seq_len(L)) {
      b <- rowSums(mu * PIP) - mu[, l] * PIP[, l]
      XtOmegaXb <- V %*% ((t(V) %*% b) * Dsq / var)
      XtOmegar <- XtOmegay - XtOmegaXb

      # Update Prior Variance ssq[l]
      if (est_ssq) {
        f <- function(x) {
          -matrixStats::logSumExp(-0.5 * log(1 + x * diagXtOmegaX) +
                                    x * XtOmegar^2 / (2 * (1 + x * diagXtOmegaX)) +
                                    logpi0)
        }
        res <- optim(par = ssq[l],
                     fn = f,
                     method = "Brent",
                     lower = ssq_range[1],
                     upper = ssq_range[2])

        if (!is.null(res$par) && res$convergence == 0) {
          ssq[l] <- res$par
          if (verbose) {
            cat(sprintf("Update s^2 for effect %d to %f\n", l, ssq[l]))
          }
        } else {
          cat(sprintf("WARNING: s^2 update for iteration %d, effect %d failed to converge; keeping previous parameters\n", it, l))
        }
      }
      # Update omega, mu, lbf, and PIP
      omega[, l] <- diagXtOmegaX + 1 / ssq[l]
      mu[, l] <- XtOmegar / omega[, l]
      lbf_variable[, l] <- XtOmegar^2 / (2 * omega[, l]) - 0.5 * log(omega[, l] * ssq[l])
      logPIP <- lbf_variable[, l] + logpi0
      lbf[l] <- matrixStats::logSumExp(logPIP)
      PIP[, l] <- exp(logPIP - lbf[l])
    } # Single Effect Regression Loop Ends

    # Update variance components
    if (est_sigmasq || est_tausq) {
      if (method == "moments") {
        moments_result <- MoM(PIP, mu, omega, sigmasq, tausq, n, V, Dsq, VtXty, Xty, yty,
                              est_sigmasq, est_tausq, verbose)
        sigmasq <- moments_result$sigmasq
        tausq <- moments_result$tausq
      } else if (method == "MLE") {
        stop("MLE method not yet implemented.")
      } else {
        stop("Unsupported variance estimation method")
      }
      # Update X'OmegaX, X'Omegay
      var <- tausq * Dsq + sigmasq
      diagXtOmegaX <- rowSums(sweep(V^2, 2, Dsq / var, `*`))
      XtOmegay <- V %*% (VtXty / var)
    }

    # Convergence check based on absolute PIP differences
    PIP_diff <- max(abs(PIP_prev - PIP))
    if (verbose) {
      cat(sprintf("Maximum change in PIP: %f\n", PIP_diff))
    }
    if (PIP_diff < PIP_tol) {
      if (verbose) {
        cat("CONVERGED\n")
      }
      break
    }
  } # Main SuSiE Loop

  # Compute posterior means of b and alpha
  b <- rowSums(mu * PIP)
  XtOmegaXb <- V %*% ((t(V) %*% b) * Dsq / var)
  XtOmegar <- XtOmegay - XtOmegaXb
  alpha <- tausq * XtOmegar

  # Compute marginal PIPs
  marginal_PIP <- 1 - apply(1 - PIP, 1, prod)

  # Compute Credible Sets
  cred <- susie_inf_get_cs(PIP = PIP, coverage = coverage, LD = LD, V = V, Dsq = Dsq, n = n)

  # Compute fitted values
  fitted <- if (!is.null(X)) X %*% (rowSums(marginal_PIP * mu) + alpha) else rep(0, n)

  return(list(
    PIP = PIP,
    marginal_PIP = marginal_PIP,
    mu = mu,
    omega = omega,
    lbf = lbf,
    lbf_variable = lbf_variable,
    ssq = ssq,
    sigmasq = sigmasq,
    tausq = tausq,
    alpha = alpha,
    fitted = fitted,
    sets = cred,
    niter = it
  ))
}

#' @title Method-of-Moments Variance Estimation
#'
#' @description
#' A subroutine that estimates the residual variance \eqn{\sigma^2} and the infinitesimal variance \eqn{\tau^2}
#' using method-of-moments (MoM).
#'
#' @return A list containing:
#' \item{sigmasq}{Updated estimate of \eqn{\sigma^2}.}
#' \item{tausq}{Updated estimate of \eqn{\tau^2}.}
#'
#'
#' @export
#'
MoM <- function(PIP, mu, omega, sigmasq, tausq, n, V, Dsq, VtXty, Xty, yty,
                est_sigmasq, est_tausq, verbose){
  # Subroutine to estimate sigma^2, tau^2 using MoM
  p <- nrow(mu)
  L <- ncol(mu)

  ### Compute A. corresponds to the matrix in equation (37) of the supplement:
  ### where Tr(X'X) = sum(Dsq) and Tr(X'X)^2 = sum(Dsq^2)
  A <- matrix(0, nrow = 2, ncol = 2)
  A[1, 1] <- n
  A[1, 2] <- sum(Dsq)
  A[2, 1] <- A[1, 2]
  A[2, 2] <- sum(Dsq^2)

  # Compute diag(V'MV)
  b <- rowSums(mu * PIP) # equation 48
  Vtb <- t(V) %*% b
  diagVtMV <- Vtb^2 # portion of equation 51 + 52
  tmpD <- rep(0, p)

  for (l in seq_len(L)) {
    bl <- mu[, l] * PIP[, l]
    Vtbl <- t(V) %*% bl
    diagVtMV <- diagVtMV - Vtbl^2
    tmpD <- tmpD + PIP[, l] * (mu[, l]^2 + 1 / omega[, l])
  }

  diagVtMV <- diagVtMV + rowSums(sweep(t(V)^2, 2, tmpD, `*`))

  # Compute x
  x <- rep(0, 2)
  x[1] <- yty - 2 * sum(b * Xty) + sum(Dsq * diagVtMV) # equation 51
  x[2] <- sum(Xty^2) - 2 * sum(Vtb * VtXty * Dsq) + sum(Dsq^2 * diagVtMV) # equation 52

  # Solves system of equations from equation 37
  if (est_tausq) {
    sol <- solve(A, x)
    if (sol[1] > 0 && sol[2] > 0) {
      sigmasq <- sol[1]
      tausq <- sol[2]
    } else {
      sigmasq <- x[1] / n
      tausq <- 0
    }
    if (verbose) {
      cat(sprintf("Update (sigma^2,tau^2) to (%f,%e)\n", sigmasq, tausq))
    }
  } else if (est_sigmasq) {
    sigmasq <- (x[1] - A[1, 2] * tausq) / n
    if (verbose) {
      cat(sprintf("Update sigma^2 to %f\n", sigmasq))
    }
  }
  return(list(sigmasq = sigmasq, tausq = tausq))
}

#' @title Credible Sets
#'
#' @description
#' Constructs credible sets of variables (e.g., SNPs) from single-effect posterior
#' inclusion probabilities (PIPs) that achieve a specified coverage level. Each credible set
#' includes the smallest number of variables whose cumulative PIP meets the given coverage.
#' Additionally, this function applies a "purity" filter to remove sets that do not meet
#' a minimum absolute correlation threshold.
#'
#' @param PIP A p x L matrix of posterior inclusion probabilities for each SNP and effect.
#'
#' @param coverage A numeric value (between 0 and 1) specifying the coverage level for each credible set.
#'
#' @param purity A numeric value specifying the minimum absolute correlation required for a
#'   credible set to be retained.
#'
#' @param LD A p x p LD matrix, equal to \eqn{X'X/n}. If not provided, it can be computed
#'   from \code{V}, \code{Dsq}, and \code{n}.
#'
#' @param V A p x p matrix of eigenvectors of \eqn{X'X}. Required if \code{LD} is not provided.
#'
#' @param Dsq A numeric vector of length p giving the eigenvalues of \eqn{X'X}.
#' Required if \code{LD} is not provided.
#'
#' @param n The sample size. Required if \code{LD} is not provided and purity filtering is needed.
#'
#' @param dedup Logical; if TRUE, removes duplicate credible sets.
#'
#' @return A list of credible sets, where each set is a list of variable indices.
#'
#' @export

susie_inf_get_cs = function(PIP, coverage = coverage, purity = 0.5, LD = NULL, V = NULL, Dsq = NULL, n = NULL, dedup = TRUE) {
  if (is.null(V) || is.null(Dsq) || is.null(n) && is.null(LD)) {
    stop("Missing inputs for purity filtering")
  }

  # Compute credible sets
  cred <- list()
  p <- nrow(PIP)
  L <- ncol(PIP)

  for (l in 1:L) {
    sortinds <- order(PIP[, l], decreasing = TRUE)
    cumsums <- cumsum(PIP[sortinds, l])
    ind <- which(cumsums >= coverage)[1]
    credset <- sortinds[1:ind]

    # Filter by purity
    if (length(credset) == 1) {
      cred[[length(cred) + 1]] <- credset
      next
    }

    if (length(credset) < 100) {
      rows <- credset
    } else {
      set.seed(123)
      rows <- sample(credset, 100, replace = FALSE)
    }

    if (!is.null(LD)) {
      LDloc <- LD[rows, rows]
    } else {
      LDloc <- (V[rows, ] %*% diag(Dsq)) %*% t(V[rows, ]) / n
      #LDloc <- (V[rows,] * Dsq) %*% t(V[rows,]) / n
    }

    if (min(abs(LDloc)) > purity) {
      cred[[length(cred) + 1]] <- sort(credset)
    }
  }

  if (dedup) {
    cred <- unique(cred)
  }

  return(cred)
}

