# ================================
# run_methods.R
# ================================
# This file contains wrapper functions for running different methods:
# - run_susie(): Standard SuSiE method.
# - run_susie_ash(): SuSiE.ash (Marginal) using susie_ash_RE_Marg().
# - run_susie_inf(): SuSiE-inf method.
# - run_fineboost(): Fineboost using colocboost().


# Wrapper for SuSiE
run_susie <- function(data, L, intercept = TRUE, standardize = TRUE) {
  cat("Starting SuSiE\n")
  out <- susie(
    X = data$X,
    y = data$y,
    L = L,
    intercept = intercept,
    standardize = standardize
  )
  return(out)
}

# Wrapper for SuSiE.ash (Marginal)
run_susie_ash <- function(data, precomp, L, K.length, upper_bound, intercept = TRUE, standardize = TRUE) {
  cat("Starting SuSiE.ash (Marginal)\n")

  # Manually scale X and y if intercept or standardize is TRUE.
  X_in <- if (intercept || standardize) scale(data$X, center = intercept, scale = standardize) else data$X
  y_in <- if (intercept || standardize) scale(data$y, center = intercept, scale = standardize) else data$y

  out <- susie_ash_RE_Marg(
    X                = X_in,
    y                = y_in,
    L                = L,
    verbose          = FALSE,
    coverage         = 0.95,
    XtX              = precomp$XtX,
    LD               = precomp$LD,
    V                = precomp$V,
    Dsq              = precomp$Dsq,
    VtXt             = precomp$VtXt,
    update_ash_sigma = FALSE,
    K.length         = K.length,
    upper_bound      = upper_bound
  )
  return(out)
}

# Wrapper for SuSiE-inf
run_susie_inf <- function(data, precomp, L, intercept = TRUE, standardize = TRUE) {
  cat("Starting SuSiE-inf\n")

  # Manually scale X and y if needed.
  X_in <- if (intercept || standardize) scale(data$X, center = intercept, scale = standardize) else data$X
  y_in <- if (intercept || standardize) scale(data$y, center = intercept, scale = standardize) else data$y

  out <- susie_inf(
    X       = X_in,
    y       = y_in,
    L       = L,
    verbose = FALSE,
    coverage = 0.95,
    XtX     = precomp$XtX,
    LD      = precomp$LD,
    V       = precomp$V,
    Dsq     = precomp$Dsq
  )
  return(out)
}

# Wrapper for Fineboost
run_fineboost <- function(data, null_max, intercept = TRUE, standardize = TRUE) {
  cat("Starting Fineboost\n")
  out <- colocboost(
    X              = data$X,
    Y              = data$y,
    check_null_max = null_max,
    intercept      = intercept,
    standardize    = standardize
  )
  return(out)
}
