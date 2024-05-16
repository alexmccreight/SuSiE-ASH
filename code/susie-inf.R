######### SuSiE-Infinitesimal #########

susie_inf <- function(X, y, L,
                      est_ssq = TRUE, ssq = NULL, ssq_range = c(0, 1), pi0 = NULL,
                      est_sigmasq = TRUE, est_tausq = TRUE, sigmasq = 1, tausq = 0,
                      method = "moments", sigmasq_range = NULL, tausq_range = NULL,
                      PIP = NULL, mu = NULL, maxiter = 100, PIP_tol = 1e-3, verbose = TRUE) {

  # Compute n, p, z, meansq, XtX, LD, V
  n <- nrow(X)
  p <- ncol(X)
  #z <- crossprod(X, y) / sqrt(n) # Vector of z-scores
  z <- (t(X) %*% y) * (1/sqrt(n)) # Vector of z-scores
  meansq <- mean(y^2) # Mean-squared magnitude of y
  #XtX <- crossprod(X)
  XtX <- t(X) %*% X
  LD <- (XtX) * (1/n) # LD Matrix

  eig <- eigen(LD, symmetric = T)
  V <- eig$vectors # pxp matrix of eigenvectors of XtX
  Dsq <- pmax(n * eig$values, 0) # Precomputed length-p vector of eigen values of XtX

  # Precompute V,D^2 in the SVD X=UDV', and V'X'y and y'y

  if ((is.null(V) || is.null(Dsq)) && is.null(LD)) {stop("Missing LD")}
  else if (is.null(V) || is.null(Dsq)) {
    eigvals <- eigen(LD, symmetric = TRUE, only.values = TRUE)$values
    V <- eigen(LD, symmetric = TRUE, only.vectors = TRUE)$vectors
    Dsq <- pmax(n * eigvals, 0)}
  else{Dsq <- pmax(Dsq, 0)}

  Xty <- sqrt(n) * z
  VtXty <- t(V) %*% Xty
  yty <- n * meansq

  # Initialize diagonal variances, diag(X', Omega X, X' Omega y)
  var <- tausq*Dsq+sigmasq
  diagXtOmegaX <- rowSums(V^2 * (Dsq / var))
  XtOmegay <- V %*% (VtXty / var)

  # DEBUG
  print(diagXtOmegaX)
  try(dim(X), silent = T)
  try(length(X), silent = T)


  # Initialize s_l^2, PIP_j, mu_j
  if (is.null(ssq)) {ssq <- rep(0.2, L)}
  if (is.null(PIP)) {PIP <- matrix(1 / p, nrow = p, ncol = L)}
  if (is.null(mu)) {mu <- matrix(0, nrow = p, ncol = L)}

  # Initialize omega_j
  lbf_variable <- matrix(0, nrow = p, ncol = L)
  lbf <- rep(0, L)
  omega <- matrix(diagXtOmegaX, nrow = p, ncol = L) + 1 / ssq


  # Initialize Prior Causal Probabilities
  if (is.null(pi0)) {
    logpi0 <- rep(log(1 / p), p)
  } else {
    logpi0 <- rep(-Inf, p)
    inds <- which(pi0 > 0)
    logpi0[inds] <- log(pi0[inds])
  }

  ##### Main SuSiE iteration loop #####
  for(it in seq_len(maxiter)){
    if (verbose) {
      cat(sprintf("Iteration %d\n", it))
    }
    PIP_prev <- PIP

    # Single Effect Regression for each effect l = 1, ... , L
    for(l in seq_len(L)){
      # Compute X', Omega r_l for residual r_l
      b <- rowSums(mu * PIP) - mu[, l] * PIP[, l]
      XtOmegaXb <- V %*% ((t(V) %*% b) * Dsq / var)
      XtOmegar <- XtOmegay - XtOmegaXb

      # Update Prior Variance ssq[l]
      if (est_ssq) {
        f <- function(x) {
          -matrixStats::logSumExp(-0.5 * log(1 + x * diagXtOmegaX) +
                       x * XtOmegar^2 / (2 * (1 + x * diagXtOmegaX)) + logpi0)
        }
        res <- optimize(f, lower = ssq_range[1], upper = ssq_range[2], maximum = FALSE)
        if (res$objective < Inf) {
          ssq[l] <- res$minimum
          if (verbose) {
            cat(sprintf("Update s^2 for effect %d to %f\n", l, ssq[l]))
          }
        } else {
          cat(sprintf("WARNING: s^2 update for iteration %d, effect %d failed to converge; keeping previous parameters\n", it, l))
        }
      }

      # Update omega, mu, and PIP
      omega[, l] <- diagXtOmegaX + 1 / ssq[l]
      mu[, l] <- XtOmegar / omega[, l]
      lbf_variable[, l] <- XtOmegar^2 / (2 * omega[, l]) - 0.5 * log(omega[, l] * ssq[l])
      logPIP <- lbf_variable[, l] + logpi0
      lbf[l] <- matrixStats::logSumExp(logPIP)
      PIP[, l] <- exp(logPIP - lbf[l])

    } # Single Effect Regression Loop

    # Update variance components
    if (est_sigmasq || est_tausq) {
      if (method == "moments") {
        moments_result <- MoM(PIP, mu, omega, sigmasq, tausq, n, V, Dsq, VtXty, Xty, yty,
                              est_sigmasq, est_tausq, verbose)
        sigmasq <- moments_result$sigmasq
        tausq <- moments_result$tausq
      } else if (method == "MLE") {
        mle_result <- MLE(PIP, mu, omega, sigmasq, tausq, n, V, Dsq, VtXty, yty,
                          est_sigmasq, est_tausq, sigmasq_range, tausq_range, it, verbose)
        sigmasq <- mle_result$sigmasq
        tausq <- mle_result$tausq
      } else {
        stop("Unsupported variance estimation method")
      }

      # Update X' Omega X, X' Omega y
      var <- tausq * Dsq + sigmasq
      diagXtOmegaX <- rowSums(V^2 * (Dsq / var))
      XtOmegay <- V %*% (VtXty / var)
    }

    # Convergence based from PIP differences
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
  PIP <- 1 - apply(1-PIP, 1, prod)

  return(list(
    PIP = PIP,
    mu = mu,
    omega = omega,
    lbf = lbf,
    lbf_variable = lbf_variable,
    ssq = ssq,
    sigmasq = sigmasq,
    tausq = tausq,
    alpha = alpha
  ))

} # Main Function

######### Method of Moments #########

MoM <- function(PIP, mu, omega, sigmasq, tausq, n, V, Dsq, VtXty, Xty, yty,
                est_sigmasq, est_tausq, verbose){
  # Subroutine to estimate sigma^2, tau^2 using MoM
  p <- nrow(mu)
  L <- ncol(mu)

  # Compute A
  A <- matrix(0, nrow = 2, ncol = 2)
  A[1, 1] <- n
  A[1, 2] <- sum(Dsq)
  A[2, 1] <- A[1, 2]
  A[2, 2] <- sum(Dsq^2)

  # Compute diag(V'MV)
  b <- rowSums(mu * PIP)
  Vtb <- t(V) %*% b
  diagVtMV <- Vtb^2
  tmpD <- rep(0, p)

  for (l in seq_len(L)) {
    bl <- mu[, l] * PIP[, l]
    Vtbl <- t(V) %*% bl
    diagVtMV <- diagVtMV - Vtbl^2
    tmpD <- tmpD + PIP[, l] * (mu[, l]^2 + 1 / omega[, l])
  }

  diagVtMV <- diagVtMV + rowSums((t(V)^2) * tmpD)

  # Compute x
  x <- rep(0, 2)
  x[1] <- yty - 2 * sum(b * Xty) + sum(Dsq * diagVtMV)
  x[2] <- sum(Xty^2) - 2 * sum(Vtb * VtXty * Dsq) + sum(Dsq^2 * diagVtMV)

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
} # Main Method of Moments
