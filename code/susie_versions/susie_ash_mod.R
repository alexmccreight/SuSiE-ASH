######### SuSiE-Infinitesimal #########
# Modified from susie-inf.R

susie_ash_mod <- function(X, y, L,
                      est_ssq = TRUE, ssq = NULL, ssq_range = c(0, 1), pi0 = NULL,
                      est_sigmasq = TRUE, est_tausq = TRUE, sigmasq = 1, tausq = 0,
                      method = "moments", sigmasq_range = NULL, tausq_range = NULL,
                      PIP = NULL, mu = NULL, maxiter = 100, PIP_tol = 1e-3, coverage = 0.9, verbose = TRUE,
                      XtX = NULL, LD = NULL, V = NULL, Dsq = NULL, ash_sd = NULL, ash_sd_method = c("default", "quadratic"),
                      K = 10, upper_bound = 3) {

  mean_y <- mean(y)

  # Compute n, z, p
  n <- nrow(X)
  z <- (t(X) %*% y)/sqrt(n)
  p <- length(z)

  # Compute mean squared magnitude of y
  meansq <- sum(y^2)/n

  # Use precomputed XtX if provided
  if (is.null(XtX)) {
    XtX <- t(X) %*% X
  }

  # Use precomputed LD if provided
  if (is.null(LD)) {
    LD <- XtX / n
  }

  # Use precomputed V and Dsq if provided
  if (is.null(V) || is.null(Dsq)) {
    eig <- eigen(LD, symmetric = TRUE)
    V <- eig$vectors[, ncol(eig$vectors):1]
    Dsq <- pmax(n * sort(eig$values), 0)
  }

  # Ensure non-negative eigenvalues
  Dsq <- pmax(Dsq, 0)

  # Pre-compute X'y, V'X'y, y'y
  Xty <- sqrt(n) * z # p x 1 matrix
  VtXty <- t(V) %*% Xty # p x 1 matrix
  yty <- n * meansq # singular value

  # Initialize diagonal variances, diag(X' Omega X), X' Omega y
  var <- tausq*Dsq+sigmasq # vector of length p
  diagXtOmegaX <- rowSums(sweep(V^2, 2, (Dsq / var), `*`)) # vector of length p (initialized each value at n)
  XtOmegay <- V %*% (VtXty / var) # p x 1 matrix

  # Initialize s_l^2, PIP_j, mu_j
  if (is.null(ssq)) {ssq <- rep(0.2, L)} # vector of length l
  if (is.null(PIP)) {PIP <- matrix(1 / p, nrow = p, ncol = L)} # p x L matrix
  if (is.null(mu)) {mu <- matrix(0, nrow = p, ncol = L)} # p x L matrix

  # Initialize omega_j
  lbf_variable <- matrix(0, nrow = p, ncol = L) # p x L matrix
  lbf <- rep(0, L) # l vector
  omega <- matrix(diagXtOmegaX, nrow = p, ncol = L) + 1 / ssq # Equation 42b omega_j(s_j^2)

  # Initialize Prior Causal Probabilities
  if (is.null(pi0)) {
    logpi0 <- rep(log(1 / p), p) # vector of length p
  } else {
    logpi0 <- rep(-Inf, p)
    inds <- which(pi0 > 0)
    logpi0[inds] <- log(pi0[inds])
  }

  ##### Main SuSiE iteration loop #####
  ### Alternate between updating prior effect size variances, ssq,
  ### and posterior distribution, q(B^(l)), sequentially for each effect
  ### and updating the infinitesimal effect size and noise variances (sigmasq, tausq)

  for(it in seq_len(maxiter)){
    if (verbose) {
      cat(sprintf("Iteration %d\n", it))
    }
    PIP_prev <- PIP

    # Single Effect Regression for each effect l = 1, ... , L
    for(l in seq_len(L)){
      # Compute X', Omega r_l for residual r_l
      b <- rowSums(mu * PIP) - mu[, l] * PIP[, l] # sum_{k: k \neq l} E[\Beta^{(k)}]; vector of length p
      XtOmegaXb <- V %*% ((t(V) %*% b) * Dsq / var) # p x 1 matrix
      XtOmegar <- XtOmegay - XtOmegaXb # Equation 42a

      # Update Prior Variance ssq[l]. Equation 43
      if (est_ssq) {
        f <- function(x) {
          -matrixStats::logSumExp(-0.5 * log(1 + x * diagXtOmegaX) +
                                    x * XtOmegar^2 / (2 * (1 + x * diagXtOmegaX)) +
                                    logpi0)
        }
        #res <- optimize(f, lower = ssq_range[1], upper = ssq_range[2], maximum = FALSE)
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

      # Update omega, mu, and PIP
      omega[, l] <- diagXtOmegaX + 1 / ssq[l] # omega_j(s_l^2)
      mu[, l] <- XtOmegar / omega[, l] # Posterior Mean = z_j^l / omega_j(s_l^2)
      lbf_variable[, l] <- XtOmegar^2 / (2 * omega[, l]) - 0.5 * log(omega[, l] * ssq[l]) # logged numerator equation 44
      logPIP <- lbf_variable[, l] + logpi0 # adding logged prior probability to get full numerator for posterior equation 44
      lbf[l] <- matrixStats::logSumExp(logPIP) # denominator for equation 44
      PIP[, l] <- exp(logPIP - lbf[l]) # combines everything and computes equation 44

    } # Single Effect Regression Loop


    ## Update variance estimation:
    bhat <- apply(mu*PIP, 1, sum)
    y_residuals <- y - X %*% bhat

    if (it == 1) {
      cat("\nUpper bound mulitpler:", upper_bound)
      cat("\nNumber of Variance Components:", K,"\n")
      if (ash_sd_method == "default") {
        # Use init_prior_sd() as is
        ash_sd <- init_prior_sd(X, y_residuals, upper_bound = upper_bound, n = K)
      } else if (ash_sd_method == "quadratic") {
        # Use init_prior_sd() to estimate upper bound 'u'
        ash_sd_sequence <- init_prior_sd(X, y_residuals, upper_bound = upper_bound, n = K)
        u <- max(ash_sd_sequence)
        ash_sd <- u * (seq(0, 1, length.out = K))^2
      } else {
        stop("Unsupported ash variance grid. Choose 'default' or 'quadratic'.")
      }
   }

    mrash_output = mr.ash.alpha::mr.ash(X = X, y = y_residuals, sa2 =ash_sd, intercept = F, standardize = F) #CHECK!!! Initial value for theta?
    mrash_post_output <- get.full.posterior(mrash_output)
    mrash.s2 <- mrash_post_output$s2; mrash.mean <- mrash_post_output$m; mrash.phi <-mrash_post_output$phi # p * K
    # mr.ash posterior values check:
    # plot(mrash_output$pi, apply(mrash.phi,2,mean)); abline(0,1,col="red"); cbind(mrash_output$pi, apply(mrash.phi,2,mean))
    # plot(mrash_output$beta, apply(mrash.phi*mrash.mean,1,sum)); abline(0,1,col="red"); head(cbind(mrash_output$beta, apply(mrash.phi*mrash.mean,1,sum)))

    #Update variance componenets from mrash
    sigmasq <- mrash_output$sigma2
    tausq <- sum(ash_sd*mrash_output$pi)


    # # # Update variance components
    # if (est_sigmasq || est_tausq) {
    #   if (method == "moments") {
    #     moments_result <- MoM(PIP, mu, omega, sigmasq, tausq, n, V, Dsq, VtXty, Xty, yty,
    #                           est_sigmasq, est_tausq, verbose)
    #     sigmasq <- moments_result$sigmasq
    #     tausq <- moments_result$tausq
    #   } else if (method == "MLE") {
    #     mle_result <- MLE(PIP, mu, omega, sigmasq, tausq, n, V, Dsq, VtXty, yty,
    #                       est_sigmasq, est_tausq, sigmasq_range, tausq_range, it, verbose)
    #     sigmasq <- mle_result$sigmasq
    #     tausq <- mle_result$tausq
    #   } else {
    #     stop("Unsupported variance estimation method")
    #   }
    #
      # Update X'OmegaX, X'Omegay
      var <- tausq * Dsq + sigmasq
      diagXtOmegaX <- rowSums(sweep(V^2, 2, Dsq / var, `*`)) # equation 19
      XtOmegay <- V %*% (VtXty / var) # equation 21 but multiply both sides by y
    # }

    # Convergence# Convergence# Convergence based from PIP differences
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
  b <- rowSums(mu * PIP) # posterior mean of the sparse effect for SNP j not conditional on inclusion
  XtOmegaXb <- V %*% ((t(V) %*% b) * Dsq / var)
  XtOmegar <- XtOmegay - XtOmegaXb # X'Omega(y - Xb)
  # alpha <- tausq * XtOmegar # Equation 27
  alpha <- mrash_output$beta

  # Compute Combined PIPs
  PIP2 <- 1 - apply(1-PIP, 1, prod) # PIP for each snp

  # Compute Credible Sets
  cred <- susie_inf_get_cs(PIP = PIP, coverage = coverage, LD = LD, V = V, Dsq = Dsq, n = n)

  # Compute fitted values.
  #fitted <- X %*% (rowSums(PIP2 * mu) + alpha)

  fitted <- mean_y + X %*% (rowSums(PIP2 * mu) + alpha) # reintroduce the mean_y for fitted values because input y is already centered.


  return(list(
    PIP = PIP,
    PIP2 = PIP2,
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
    mrash.s2 = mrash_post_output$s2,
    mrash.mean = mrash_post_output$m,
    mrash.phi = mrash_post_output$phi, # p * K
    mrash.pi = mrash_output$pi,
    mrash.sa2 = mrash_output$data$sa2
  ))

}

######### Method of Moments #########

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

  # Solves system of equations from equation 37 (using more-efficient eigenvalue decomposition values)
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

######### Credible Set Generation #########

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

    #ind <- min(which(cumsum(PIP[sortinds,l]) >= coverage))
    #credset <- sortinds[1:(ind+1)]

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

init_prior_sd <- function(X, y, upper_bound = 3, n = 30) {
  res <- univariate_regression(X, y)
  smax <- upper_bound*max(res$betahat)
  seq(0, smax, length.out = n)
}
