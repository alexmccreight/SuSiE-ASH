susie_inf_revised_no_omega <- function(X, y, L,
                                       est_ssq = TRUE, ssq = NULL, ssq_range = c(0, 1), pi0 = NULL,
                                       est_sigmasq = TRUE, est_tausq = TRUE, sigmasq = 1, tausq = 0,
                                       method = "moments", sigmasq_range = NULL, tausq_range = NULL,
                                       PIP = NULL, mu = NULL, maxiter = 100, PIP_tol = 1e-3, coverage = 0.9, verbose = TRUE) {

  n <- nrow(X)
  p <- ncol(X)
  XtX <- t(X) %*% X
  Xty <- t(X) %*% y

  if (is.null(ssq)) { ssq <- rep(0.2, L) } # Prior effect size variance initialization
  if (is.null(PIP)) { PIP <- matrix(1 / p, nrow = p, ncol = L) } # PIP initialization
  if (is.null(mu)) { mu <- matrix(0, nrow = p, ncol = L) } # Posterior mean initialization

  if (is.null(pi0)) {
    pi0 <- rep(1/p, p)  # Uniform prior
  }

  lbf <- rep(0, L)

  # Main loop for SuSiE iterations
  for(it in seq_len(maxiter)) {
    if (verbose) {
      cat(sprintf("Iteration %d\n", it))
    }

    PIP_prev <- PIP

    for(l in seq_len(L)) {
      # Compute the residual for effect l
      r_l <- y - X %*% (rowSums(mu * PIP) - mu[, l] * PIP[, l])

      # Compute X'X_r_l and X'y
      Xtr_l <- t(X) %*% r_l

      # Update Prior Variance ssq[l] (optimization procedure)
      if (est_ssq) {
        f <- function(x) {
          -matrixStats::logSumExp(-0.5 * log(1 + x * diag(XtX)) +
                                    x * Xtr_l^2 / (2 * (1 + x * diag(XtX))) +
                                    log(pi0))
        }
        res <- optim(par = ssq[l], fn = f, method = "Brent", lower = ssq_range[1], upper = ssq_range[2])
        if (!is.null(res$par) && res$convergence == 0) {
          ssq[l] <- res$par
          if (verbose) {
            cat(sprintf("Update s^2 for effect %d to %f\n", l, ssq[l]))
          }
        } else {
          cat(sprintf("WARNING: s^2 update for iteration %d, effect %d failed to converge; keeping previous parameters\n", it, l))
        }
      }

      # Update omega, mu, and PIP directly using X and X'
      omega <- diag(XtX) + 1 / ssq[l]
      mu[, l] <- Xtr_l / omega
      lbf_variable <- Xtr_l^2 / (2 * omega) - 0.5 * log(omega * ssq[l])
      logPIP <- lbf_variable + log(pi0)
      lbf[l] <- matrixStats::logSumExp(logPIP)
      PIP[, l] <- exp(logPIP - lbf[l])
    }

    # Update variance components if required
    if (est_sigmasq || est_tausq) {
      if (method == "moments") {
        moments_result <- MoM_no_omega(PIP, mu, XtX, Xty, y, est_sigmasq, est_tausq, verbose)
        sigmasq <- moments_result$sigmasq
        tausq <- moments_result$tausq
      } else if (method == "MLE") {
        mle_result <- MLE(PIP, mu, XtX, Xty, y, est_sigmasq, est_tausq, sigmasq_range, tausq_range, it, verbose)
        sigmasq <- mle_result$sigmasq
        tausq <- mle_result$tausq
      } else {
        stop("Unsupported variance estimation method")
      }
    }

    # Convergence check based on PIP differences
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
  }

  # Compute posterior means of b and alpha without Omega
  b <- rowSums(mu * PIP)
  alpha <- tausq * (Xty - XtX %*% b) # No precomputation, direct calculation
  PIP2 <- 1 - apply(1 - PIP, 1, prod)

  # Compute fitted values directly
  #fitted <- X %*% (rowSums(PIP2 * mu) + alpha)
  fitted <- X %*% rowSums(PIP2 * mu) + X %*% alpha

  # Credible sets calculation remains the same
  cred <- susie_inf_get_cs_no_omega(PIP = PIP, coverage = coverage, XtX = XtX/n)

  return(list(PIP = PIP, PIP2 = PIP2, mu = mu, omega = omega, lbf = lbf, ssq = ssq,
              sigmasq = sigmasq, tausq = tausq, alpha = alpha, fitted = fitted, sets = cred))
}


##### METHOD OF MOMENTS #####

MoM_no_omega <- function(PIP, mu, XtX, Xty, y, est_sigmasq, est_tausq, verbose) {
  n <- length(y)
  p <- nrow(mu)
  L <- ncol(mu)

  # Compute residuals based on the current estimates of sparse effects (b)
  b <- rowSums(mu * PIP)
  residuals <- y - X %*% b

  # Compute terms for variance estimation
  X_residuals <- t(X) %*% residuals
  sum_XtX_diag <- sum(diag(XtX))  # Sum of the diagonal elements of XtX

  # Compute A matrix for the linear system
  A <- matrix(0, nrow = 2, ncol = 2)
  A[1, 1] <- n  # Trace of I
  A[1, 2] <- sum_XtX_diag  # Trace of XtX
  A[2, 1] <- A[1, 2]
  A[2, 2] <- sum(diag(XtX %*% XtX))  # Trace of (XtX)^2

  # Compute x vector for the linear system
  x <- rep(0, 2)
  x[1] <- sum(residuals^2)  # y'y - 2*y'Xb + b'X'Xb (simplifies to just sum of squared residuals)
  x[2] <- sum(X_residuals^2)  # (X'y)^2 - 2*X'y*Xb + (Xb)'(XtX)(Xb)

  # Solves the linear system A * [sigmasq, tausq] = x
  if (est_tausq) {
    sol <- solve(A, x)
    if (sol[1] > 0 && sol[2] > 0) {
      sigmasq <- sol[1]
      tausq <- sol[2]
    } else {
      sigmasq <- x[1] / n  # If negative, constrain tausq to 0 and recompute sigmasq
      tausq <- 0
    }
    if (verbose) {
      cat(sprintf("Updated (sigma^2, tau^2) to (%f, %e)\n", sigmasq, tausq))
    }
  } else if (est_sigmasq) {
    sigmasq <- (x[1] - A[1, 2] * tausq) / n
    if (verbose) {
      cat(sprintf("Updated sigma^2 to %f\n", sigmasq))
    }
  }

  return(list(sigmasq = sigmasq, tausq = tausq))
}


##### Credible Set #####

susie_inf_get_cs_no_omega <- function(PIP, coverage = 0.9, purity = 0.5, XtX = NULL, dedup = TRUE) {
  if (is.null(XtX)) {
    stop("Missing XtX matrix for purity filtering")
  }

  cred <- list()  # Initialize list to store credible sets
  p <- nrow(PIP)
  L <- ncol(PIP)

  for (l in seq_len(L)) {
    # Step 1: Sort SNPs by PIP in descending order
    sortinds <- order(PIP[, l], decreasing = TRUE)
    cumsums <- cumsum(PIP[sortinds, l])

    # Step 2: Identify the smallest set that reaches the desired coverage level
    ind <- which(cumsums >= coverage)[1]
    credset <- sortinds[1:ind]

    # Step 3: Filter the credible set based on purity
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

    # Calculate the LD matrix for the selected SNPs using XtX
    LDloc <- XtX[rows, rows] / diag(XtX)[rows]

    # Check if the minimum absolute correlation is greater than the purity threshold
    if (min(abs(LDloc)) > purity) {
      cred[[length(cred) + 1]] <- sort(credset)
    }
  }

  if (dedup) {
    cred <- unique(cred)
  }

  return(cred)
}
