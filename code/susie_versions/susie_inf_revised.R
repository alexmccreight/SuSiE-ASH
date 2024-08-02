susie_inf_revised <- function(X, y, L,
                      est_ssq = TRUE, ssq = NULL, ssq_range = c(0, 1), pi0 = NULL,
                      est_sigmasq = TRUE, est_tausq = TRUE, sigmasq = 1, tausq = 0,
                      method = "moments", sigmasq_range = NULL, tausq_range = NULL,
                      PIP = NULL, mu = NULL, maxiter = 100, PIP_tol = 1e-3, coverage = 0.9, verbose = TRUE) {

  # Initialize values

  # Main Loop

  # within Main Loop run SER for each effect: 1, ... , L

  # within Main Loop update prior variance for each effect and update

  # within Main Loop update posterior mean, etc.

  # within Main Loop use MoM to update residual variance and infinitesimal variance

  # within Main loop check for convergence

  # finalize posterior means, calculate BLUP for infinitesimal posterior means

  # calculate final parameter (credible sets, fitted values, etc.)
}
