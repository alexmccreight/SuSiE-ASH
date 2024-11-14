# New simulation Setting Oct 19/2024
init_prior_sd <- function(X, y, n = 30) {
  res <- univariate_regression(X, y)
  smax <- 3 * max(res$betahat)
  seq(0, smax, length.out = n)
}

generate_eqtl_data <- function(X,
                               h2_total = 0.3,            # Total heritability
                               prop_h2_sparse = 0.65,     # Proportion of h2_total explained by sparse effects (including sentinel)
                               prop_h2_oligogenic = 0.20, # Proportion of h2_total explained by oligogenic effects
                               prop_h2_infinitesimal = 0.15, # Proportion of h2_total explained by infinitesimal effects
                               prop_h2_sentinel = 0.7,    # Proportion of h2_sparse explained by sentinel SNP
                               n_oligogenic = 20,
                               mixture_props = c(0.6, 0.4), # Adjusted mixture proportions
                               mixture_sds = c(0.0025, 0.005), # Number of oligogenic SNPs
                               seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  ori.X <- X
  X <- scale(X)
  
  n_samples <- nrow(X)
  n_features <- ncol(X)
  
  # Calculate effect sizes for each component
  h2_sparse <- h2_total * prop_h2_sparse
  h2_sentinel <- h2_sparse * prop_h2_sentinel
  h2_other_sparse <- h2_sparse - h2_sentinel
  h2_oligogenic <- h2_total * prop_h2_oligogenic
  h2_infinitesimal <- h2_total * prop_h2_infinitesimal
  
  # Generate effect sizes
  beta <- rep(0, n_features)
  
  # Sentinel SNP effect (largest effect among sparse effects)
  sentinel_index <- sample(1:n_features, 1)
  beta[sentinel_index] <- rnorm(1, 0, sqrt(h2_sentinel))
  
  # Other sparse effects (large and mappable)
  #n_other_sparse <- max(1, rpois(1, lambda = 2))  # Ensure at least one other sparse effect
  n_other_sparse <- 2
  other_sparse_indices <- sample(setdiff(1:n_features, sentinel_index), n_other_sparse)
  if (n_other_sparse > 0) {
    # Distribute h2_other_sparse equally among other sparse SNPs
    beta[other_sparse_indices] <- rnorm(n_other_sparse, 0, sqrt(h2_other_sparse / n_other_sparse))
    # Ensure the sentinel SNP has the largest effect size
    max_other_sparse_effect <- max(abs(beta[other_sparse_indices]))
    if (abs(beta[sentinel_index]) <= max_other_sparse_effect) {
      beta[sentinel_index] <- sign(beta[sentinel_index]) * (max_other_sparse_effect + 0.01)
    }
  }
  
  # After assigning effect sizes to sparse SNPs
  # Combined sparse effects
  sparse_indices <- c(sentinel_index, other_sparse_indices)
  sparse_effects <- X[, sparse_indices] %*% beta[sparse_indices]
  
  # Scale sparse effects to achieve desired heritability
  scaling_factor_sparse <- sqrt(h2_sparse / var(sparse_effects))
  beta[sparse_indices] <- beta[sparse_indices] * as.vector(scaling_factor_sparse)
  
  # Ensure the sentinel SNP has the largest effect
  max_other_sparse_effect <- max(abs(beta[other_sparse_indices]))
  if (abs(beta[sentinel_index]) <= max_other_sparse_effect) {
    beta[sentinel_index] <- sign(beta[sentinel_index]) * (max_other_sparse_effect + 0.01)
  }
  
  # Oligogenic effects (adjust mixture proportions and sds)
  non_sparse_indices <- setdiff(1:n_features, c(sentinel_index, other_sparse_indices))
  n_oligogenic <- min(n_oligogenic, length(non_sparse_indices))
  oligogenic_indices <- sample(non_sparse_indices, n_oligogenic, replace = FALSE)
  
  mixture_assignments <- sample(1:length(mixture_props), length(oligogenic_indices), replace = TRUE, prob = mixture_props)
  beta[oligogenic_indices] <- rnorm(length(oligogenic_indices), 0, mixture_sds[mixture_assignments])
  
  # Scale oligogenic effects to achieve desired heritability
  oligogenic_effects <- X[, oligogenic_indices] %*% beta[oligogenic_indices]
  scaling_factor <- sqrt(h2_oligogenic / var(oligogenic_effects))
  beta[oligogenic_indices] <- beta[oligogenic_indices] * as.vector(scaling_factor)
  
  # Infinitesimal effects (small effects on remaining SNPs)
  infinitesimal_indices <- setdiff(non_sparse_indices, oligogenic_indices)
  infinitesimal_effects <- rep(0, n_features)
  if (length(infinitesimal_indices) > 0) {
    infinitesimal_effects[infinitesimal_indices] <- rnorm(length(infinitesimal_indices), 0, sqrt(h2_infinitesimal / length(infinitesimal_indices)))
  }
  beta <- beta + as.vector(infinitesimal_effects)
  
  # Generate y
  y <- X %*% beta
  
  # Add noise to achieve desired total heritability
  var_y <- var(as.vector(y))
  var_epsilon <- var_y * (1 - h2_total) / h2_total
  epsilon <- rnorm(n_samples, 0, sqrt(var_epsilon))
  y <- y + epsilon
  
  # Calculate actual heritabilities
  var_y_total <- var(as.vector(y))
  h2_sentinel_actual <- var(X[, sentinel_index] * beta[sentinel_index]) / var_y_total
  
  # Sparse effects heritability
  sparse_indices <- c(sentinel_index, other_sparse_indices)
  h2_sparse_actual <- var(X[, sparse_indices] %*% beta[sparse_indices]) / var_y_total
  
  # Oligogenic effects heritability
  h2_oligogenic_actual <- var(X[, oligogenic_indices] %*% beta[oligogenic_indices]) / var_y_total
  
  # Infinitesimal effects heritability
  h2_infinitesimal_actual <- var(X %*% infinitesimal_effects) / var_y_total
  
  h2_total_actual <- var(as.vector(X %*% beta)) / var_y_total
  
  ori.y <- y
  y <- scale(y, center = TRUE, scale = FALSE)
  
  # Create a full-length mixture_assignments vector
  mixture_assignments_full <- rep(NA, n_features)
  mixture_assignments_full[oligogenic_indices] <- mixture_assignments
  
  return(list(
    ori.X = ori.X,
    X = X,
    ori.y = ori.y,
    y = y,
    beta = beta,
    h2_total = h2_total_actual,
    h2_sparse = h2_sparse_actual,
    h2_sentinel = h2_sentinel_actual,
    h2_oligogenic = h2_oligogenic_actual,
    h2_infinitesimal = h2_infinitesimal_actual,
    sentinel_index = sentinel_index,
    other_sparse_indices = other_sparse_indices,
    oligogenic_indices = oligogenic_indices,
    infinitesimal_indices = infinitesimal_indices,
    mixture_assignments = mixture_assignments_full,
    var_epsilon = var_epsilon,
    sparse_indices = sparse_indices
  ))
}

is_causal <- function(eqtl_data, pve_threshold){
  # Get the beta vector and residual variance for this simulation
  beta <- eqtl_data$beta
  var_epsilon <- eqtl_data$var_epsilon
  
  # Compute variance explained by each SNP (since Var(X_j) = 1)
  variance_explained <- beta^2
  
  # Compute total genetic variance (assuming SNPs are uncorrelated)
  var_g <- sum(variance_explained)
  
  # Compute total variance (genetic variance + residual variance)
  total_variance <- var_g + var_epsilon
  
  # Compute PVE for each SNP
  proportion_var_explained <- variance_explained / total_variance
  
  # Define causal SNPs based on the current PVE threshold
  causal_SNPs <- which(proportion_var_explained > pve_threshold)
  return(causal = causal_SNPs)
}


compute_metrics <- function(test_cs, causal_SNPs) {
  if (length(test_cs) > 0) {
    cs_size <- mean(sapply(test_cs, length))
    
    # Calculate coverage (proportion of credible sets that include at least one causal SNP)
    coverage <- mean(sapply(test_cs, function(cs) any(cs %in% causal_SNPs)))
    
    # CS-based FDR
    TP_fdr <- sum(sapply(test_cs, function(cs) any(cs %in% causal_SNPs)))
    FP_fdr <- length(test_cs) - TP_fdr
    cs_fdr <- FP_fdr / length(test_cs)
    
    # CS-based Recall
    causal_SNPs_detected <- unique(unlist(test_cs))[unique(unlist(test_cs)) %in% causal_SNPs]
    cs_recall <- length(causal_SNPs_detected) / length(causal_SNPs)
  } else {
    cs_size <- NA
    coverage <- NA
    cs_fdr <- NA
    cs_recall <- NA
  }
  
  return(list(
    cs_size = cs_size,
    coverage = coverage,
    cs_fdr = cs_fdr,
    cs_recall = cs_recall
  ))
}
