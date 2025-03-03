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
run_susie_ash <- function(data, precomp, L, K.length, upper_bound) {
  cat("Starting SuSiE.ash (Marginal)\n")

  out <- susie_ash_RE_Marg(
    X                = scale(data$X),
    y                = scale(data$y),
    L                = L,
    verbose          = FALSE,
    coverage         = 0.95,
    update_ash_sigma = FALSE,
    K.length         = K.length,
    upper_bound      = upper_bound
  )
  return(out)
}

# Wrapper for SuSiE-inf
run_susie_inf <- function(data, precomp, L) {
  cat("Starting SuSiE-inf\n")

  out <- susie_inf(
    X       = scale(data$X),
    y       = scale(data$y),
    L       = L,
    verbose = FALSE,
    coverage = 0.95
  )
  return(out)
}


# Wrapper for Fineboost
run_fineboost <- function(data, null_max, intercept = TRUE, standardize = TRUE) {
  cat("Starting Fineboost\n")

  # Check that X is a matrix
  if (!is.matrix(data$X)) {
    stop("data$X must be a matrix.")
  }

  # If column names are missing in X, assign default names
  if (is.null(colnames(data$X))) {
    cat("No column names found in data$X. Assigning default column names.\n")
    colnames(data$X) <- paste0("var", seq_len(ncol(data$X)))
  }

  # Get the number of rows in X for consistency check
  nX <- nrow(data$X)

  # Check and process y: it must be a vector or a single-column matrix.
  if (is.vector(data$y)) {
    if (length(data$y) != nX) {
      stop("Mismatch: length of data$y does not equal number of rows in data$X.")
    }
    # Convert vector to matrix then wrap in a list
    data$y <- list(as.matrix(data$y))
  } else if (is.matrix(data$y)) {
    if (nrow(data$y) != nX) {
      stop("Mismatch: number of rows in data$y does not equal number of rows in data$X.")
    }
    if (ncol(data$y) != 1) {
      warning("data$y has more than one column. Only the first column will be used.")
      data$y <- list(as.matrix(data$y[, 1, drop = FALSE]))
    } else {
      data$y <- list(data$y)
    }
  } else {
    stop("data$y must be a vector or a matrix.")
  }

  out <- colocboost(
    X              = data$X,
    Y              = data$y,
    check_null_max = null_max,
    intercept      = intercept,
    standardize    = standardize
  )
  return(out)
}
