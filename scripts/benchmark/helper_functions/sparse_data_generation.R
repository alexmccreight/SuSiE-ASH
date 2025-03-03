# helper function to sample causal indices that satisfy an LD threshold
get_valid_causal <- function(G, ncausal, ld_threshold, max_attempts = 10000, ld_matrix = NULL) {
  snp_indices <- seq_len(ncol(G))
  for (attempt in 1:max_attempts) {
    causal_indices <- sort(sample(snp_indices, ncausal))
    if (!is.null(ld_matrix)) {
      # use the provided LD matrix (assumed to already be absolute)
      corr_mat <- ld_matrix[causal_indices, causal_indices]
    } else {
      # compute correlation matrix and take absolute values
      corr_mat <- abs(cor(G[, causal_indices]))
    }
    # set diagonal and lower triangle to zero (only off-diagonals matter)
    corr_mat[lower.tri(corr_mat, diag = TRUE)] <- 0
    if (max(corr_mat) < ld_threshold) {
      return(causal_indices)
    }
  }
  stop("Could not find a set of causal variants with LD (|r|) below ", ld_threshold, " after ", max_attempts, " attempts.")
}

# helper function to parse number of causals
parse_num_causal_snps <- function(value) {
  is_pct <- FALSE
  if (grepl("[^0-9]+(pct|avg)?$", value, ignore.case = TRUE)) {
    num_tmp <- as.numeric(gsub("[^0-9]+", "", value))
    num_mod <- tolower(gsub("[^a-zA-Z]+", "", value))

    if (num_mod %in% c("pct", "avg")) {
      if (num_mod == "pct") {
        if (num_tmp <= 0 || num_tmp > 100) {
          stop("Percentage of causal SNPs must be in (0, 100].")
        }
        num_tmp <- num_tmp / 100
        is_pct <- TRUE
      } else if (num_mod == "avg") {
        if (num_tmp == 0) {
          stop("Average number of causal SNPs must be at least 1")
        }
        num_smpl <- 0
        while (num_smpl == 0) {
          num_smpl <- rpois(1, num_tmp)
        }
        num_tmp <- num_smpl
      }
    } else {
      if (num_tmp < 1) {
        stop("Number of causal SNPs must be at least 1")
      }
    }
  } else {
    num_tmp = as.numeric(gsub("[^0-9]+", "", value))
  }

  return(list(value = num_tmp, is_pct = is_pct))
}

# helper function to generate phenotype
simulate_polygenic_trait <- function(g, h2g) {
  n <- length(g)

  if (h2g > 0) {
    s2g <- var(g)  # Genetic variance
    s2e <- s2g * ((1 / h2g) - 1)  # Environmental variance
    e <- rnorm(n, mean = 0, sd = sqrt(s2e))
    y <- g + e
  } else {
    y <- rnorm(n, mean = 0, sd = 1)
  }

  # Standardize the phenotype
  y <- scale(y, center = TRUE, scale = FALSE)

  return(y)
}

generate_sparse_eqtl_data <- function(X, K = 10, h2 = 0.3, seed = NULL,
                                      ld_mode = "random", ld_matrix = NULL, max_attempts = 10000) {
  if (!is.null(seed)) set.seed(seed)

  n_samples <- nrow(X)
  n_features <- ncol(X)

  # parse K to get the number of causal SNPs
  causal_info <- parse_num_causal_snps(as.character(K))
  if (causal_info$is_pct) {
    n_causal <- max(1, round(causal_info$value * n_features))
  } else {
    n_causal <- min(causal_info$value, n_features)
  }

  # determine causal indices based on ld_mode
  if (ld_mode == "random") {
    causal_indices <- sort(sample(seq_len(n_features), n_causal, replace = FALSE))
  } else if (ld_mode == "minimal_ld") {
    causal_indices <- get_valid_causal(G = X, ncausal = n_causal,
                                       ld_threshold = 0.05, max_attempts = max_attempts,
                                       ld_matrix = ld_matrix)
  } else if (ld_mode == "low_ld") {
    causal_indices <- get_valid_causal(G = X, ncausal = n_causal,
                                       ld_threshold = 0.30, max_attempts = max_attempts,
                                       ld_matrix = ld_matrix)
  } else {
    stop("Invalid ld_mode. Choose from 'random', 'minimal_ld', or 'low_ld'.")
  }

  # create beta vector and assign effect sizes from N(0, 0.6^2) for the causal SNPs
  beta <- rep(0, n_features)
  beta[causal_indices] <- rnorm(length(causal_indices), mean = 0, sd = 0.6)

  # compute the latent genetic effect
  g <- as.vector(X %*% beta)

  # generate phenotype using the GitHub simulate_polygenic_trait function
  y <- simulate_polygenic_trait(g, h2)

  # compute estimated heritability
  var_g <- var(g)
  var_e <- var(y - g)
  h2_estimated <- var_g / (var_g + var_e)

  return(list(
    X = X,
    y = as.vector(y),
    beta = beta,
    causal_indices = causal_indices,
    var_epsilon = var_e,
    h2_input = h2,
    h2_estimated = h2_estimated
  ))
}
