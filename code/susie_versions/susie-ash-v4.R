## v4. Normal susie until convergence, ash on residuals

susie_ash_v4 = function(X,y,L = min(10,ncol(X)),
                        scaled_prior_variance = 0.2,
                        residual_variance = NULL,
                        prior_weights = NULL,
                        null_weight = 0,
                        standardize = TRUE,
                        intercept = TRUE,
                        estimate_residual_variance = TRUE,
                        estimate_prior_variance = TRUE,
                        estimate_prior_method = c("optim", "EM", "simple"),
                        check_null_threshold = 0,
                        prior_tol = 1e-9,
                        residual_variance_upperbound = Inf,
                        s_init = NULL,
                        coverage = 0.95,
                        min_abs_corr = 0.5,
                        median_abs_corr = NULL,
                        compute_univariate_zscore = FALSE,
                        na.rm = FALSE,
                        max_iter = 100,
                        tol = 1e-3,
                        verbose = FALSE,
                        track_fit = FALSE,
                        residual_variance_lowerbound = var(drop(y))/1e4,
                        refine = FALSE,
                        n_purity = 100){

  # Run SuSiE until convergence
  susie_output <- susie(X = X, y = y, L = L, intercept = intercept, standardize = standardize)

  # Check for high heritability and update y_residuals
  high_heritability_ls <- which(susie_output$V >= 0.001)

  if (length(high_heritability_ls) > 0) {
    if (length(high_heritability_ls) == 1) {
      y_residuals <- y - X %*% (susie_output$alpha[high_heritability_ls,] * susie_output$mu[high_heritability_ls,])
    } else {
      y_residuals <- y - X %*% colSums(susie_output$alpha[high_heritability_ls,] * susie_output$mu[high_heritability_ls,])
    }
  }

  # Run mr.ash until convergence
  mrash_output <- mr.ash(X = X,
                         y = y_residuals,
                         intercept = intercept,
                         standardize = standardize,
                         sa2 = nrow(X) * (2^((0:19)/20) - 1)^2)

  Xtheta <- X %*% mr.ash.alpha::coef.mr.ash(mrash_output)[-1]
  Xr <- susie_output$Xr

  fitted <- Xr + Xtheta



  output <- list(
    susie_output = susie_output,
    mrash_output = mrash_output,
    y_residuals = y_residuals,
    fitted = fitted,
    Xr = Xr,
    Xtheta = Xtheta
  )

  return(output)
}
