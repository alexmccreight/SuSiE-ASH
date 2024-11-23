# =====================================
# Updated eQTL Simulation Script
# =====================================

# Load Required Libraries
library(susieR)
library(mr.ash.alpha)
library(dplyr)
library(magrittr)

# Source Custom Functions
source("susie-ash-data/susie_ash_mod.R")
source("susie-ash-data/susie_ash_mod_v2.R")
source("susie-ash-data/susie_inf.R")

# Function to Generate eQTL Data
generate_eqtl_data <- function(X,
                               h2_total = 0.3,            # Total heritability
                               prop_h2_sparse = 0.65,     # Proportion of h2_total explained by sparse effects (including sentinel)
                               prop_h2_oligogenic = 0.20, # Proportion of h2_total explained by oligogenic effects
                               prop_h2_infinitesimal = 0.15, # Proportion of h2_total explained by infinitesimal effects
                               prop_h2_sentinel = 0.7,    # Proportion of h2_sparse explained by sentinel SNP
                               n_oligogenic = 20,
                               mixture_props = c(0.6, 0.4), # Adjusted mixture proportions
                               mixture_sds = c(0.0025, 0.005), # Standard deviations for mixture components
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

# Function to Identify Causal SNPs Based on PVE Threshold
is_causal <- function(eqtl_data, pve_threshold) {
  # Get the beta vector and residual variance for this simulation
  beta <- eqtl_data$beta
  var_epsilon <- eqtl_data$var_epsilon

  # Compute variance explained by each SNP (since Var(X_j) = 1)
  variance_explained <- beta^2

  # Compute total genetic variance
  var_g <- sum(variance_explained)

  # Compute total variance (genetic variance + residual variance)
  total_variance <- var_g + var_epsilon

  # Compute PVE for each SNP
  proportion_var_explained <- variance_explained / total_variance

  # Define causal SNPs based on the current PVE threshold
  causal_SNPs <- which(proportion_var_explained > pve_threshold)
  return(causal = causal_SNPs)
}

# Method and Metrics Function
method_and_score <- function(X,
                             y,
                             beta,
                             causal,
                             L = 10,
                             precomputed_matrices,
                             seed) {

  set.seed(seed)
  XtX <- precomputed_matrices$XtX
  LD <- precomputed_matrices$LD
  V <- precomputed_matrices$V
  Dsq <- precomputed_matrices$Dsq

  #### Run Various Methods ####

  # SuSiE
  cat("Starting SuSiE\n")
  susie_output <- susie(X = X, y = y, L = L, intercept = TRUE, standardize = TRUE)

  # SuSiE-ash (default)
  cat("Starting SuSiE-ash RE (default)\n")
  susie_ash_default_output <- susie_ash_mod(
    X = scale(X),
    y = scale(y, center = TRUE, scale = FALSE),
    L = L,
    verbose = FALSE,
    coverage = 0.95,
    XtX = XtX,
    LD = LD,
    V = V,
    Dsq = Dsq,
    ash_sd_method = "default"
  )

  # SuSiE-ash (quadratic)
  cat("Starting SuSiE-ash RE (quadratic)\n")
  susie_ash_quad_output <- susie_ash_mod(
    X = scale(X),
    y = scale(y, center = TRUE, scale = FALSE),
    L = L,
    verbose = FALSE,
    coverage = 0.95,
    XtX = XtX,
    LD = LD,
    V = V,
    Dsq = Dsq,
    ash_sd_method = "quadratic"
  )

  # SuSiE-ash v2 (default)
  cat("Starting SuSiE-ash RE v2 (default)\n")
  susie_ash_v2_default_output <- susie_ash_mod_v2(
    X = scale(X),
    y = scale(y, center = TRUE, scale = FALSE),
    L = L,
    verbose = FALSE,
    coverage = 0.95,
    XtX = XtX,
    LD = LD,
    V = V,
    Dsq = Dsq,
    ash_sd_method = "default"
  )

  # SuSiE-ash v2 (quadratic)
  cat("Starting SuSiE-ash RE v2 (quadratic)\n")
  susie_ash_v2_quad_output <- susie_ash_mod_v2(
    X = scale(X),
    y = scale(y, center = TRUE, scale = FALSE),
    L = L,
    verbose = FALSE,
    coverage = 0.95,
    XtX = XtX,
    LD = LD,
    V = V,
    Dsq = Dsq,
    ash_sd_method = "quadratic"
  )

  # SuSiE-inf
  cat("Starting SuSiE-inf\n")
  susie_inf_output <- susie_inf(
    X = scale(X),
    y = scale(y, center = TRUE, scale = FALSE),
    L = L,
    verbose = FALSE,
    coverage = 0.95,
    XtX = XtX,
    LD = LD,
    V = V,
    Dsq = Dsq
  )

  #### Function to Calculate Metrics ####
  calc_metrics <- function(mod, X = X, y = y, causal = causal) {
    #### Initialize values ####
    test.cs <- susie_get_cs(mod, X = X, coverage = 0.95)$cs
    coverage <- 0
    cs_fdr <- 0
    cs_recall <- 0
    cs_size <- 0

    if (length(test.cs) > 0) {
      # Calculate Average CS Size
      cs_size <- length(unlist(test.cs)) / length(test.cs)

      # Calculate Coverage (proportion of credible sets with a causal effect)
      coverage <- (lapply(1:length(test.cs), function(cs.l) { ifelse(sum(causal %in% test.cs[[cs.l]]) != 0, TRUE, FALSE) }) %>% unlist() %>% sum()) / (length(test.cs))

      # CS Based FDR
      TP_fdr <- lapply(1:length(test.cs), function(cs.l) { ifelse(sum(test.cs[[cs.l]] %in% causal) != 0, TRUE, FALSE) }) %>% unlist() %>% sum()
      FP_fdr <- length(test.cs) - TP_fdr
      cs_fdr <- ifelse((TP_fdr + FP_fdr) > 0, FP_fdr / (TP_fdr + FP_fdr), NA)

      # CS Based Recall
      TP_recall <- sum(causal %in% unlist(test.cs))
      FN_recall <- length(causal) - TP_recall
      cs_recall <- TP_recall / (TP_recall + FN_recall)
    }

    #### Calculate RMSE ####
    RMSE_y <- sqrt(mean((y - mod$fitted)^2))

    #### Store Results ####
    return(list(
      RMSE_y = RMSE_y,
      cs_size = cs_size,
      coverage = coverage,
      cs_fdr = cs_fdr,
      cs_recall = cs_recall
    ))
  }

  calc_metrics_inf <- function(mod, X = X, y = y, causal = causal) {
    #### Initialize values ####
    test.cs <- mod$sets
    coverage <- 0
    cs_fdr <- 0
    cs_recall <- 0
    cs_size <- 0

    if (length(test.cs) > 0) {
      # Calculate Average CS Size
      cs_size <- length(unlist(test.cs)) / length(test.cs)

      # Calculate Coverage (proportion of credible sets with a causal effect)
      coverage <- (lapply(1:length(test.cs), function(cs.l) { ifelse(sum(causal %in% test.cs[[cs.l]]) != 0, TRUE, FALSE) }) %>% unlist() %>% sum()) / (length(test.cs))

      # CS Based FDR
      TP_fdr <- lapply(1:length(test.cs), function(cs.l) { ifelse(sum(test.cs[[cs.l]] %in% causal) != 0, TRUE, FALSE) }) %>% unlist() %>% sum()
      FP_fdr <- length(test.cs) - TP_fdr
      cs_fdr <- ifelse((TP_fdr + FP_fdr) > 0, FP_fdr / (TP_fdr + FP_fdr), NA)

      # CS Based Recall
      TP_recall <- sum(causal %in% unlist(test.cs))
      FN_recall <- length(causal) - TP_recall
      cs_recall <- TP_recall / (TP_recall + FN_recall)
    }

    #### Calculate RMSE ####
    RMSE_y <- sqrt(mean((y - mod$fitted)^2))

    #### Store Results ####
    return(list(
      RMSE_y = RMSE_y,
      cs_size = cs_size,
      coverage = coverage,
      cs_fdr = cs_fdr,
      cs_recall = cs_recall
    ))
  }

  #### Calculate Metrics for Each Method ####
  susie_metrics <- calc_metrics(susie_output, X, y, causal)
  susie_ash_default_metrics <- calc_metrics_inf(susie_ash_default_output, X, y, causal)
  susie_ash_quad_metrics <- calc_metrics_inf(susie_ash_quad_output, X, y, causal)
  susie_ash_v2_default_metrics <- calc_metrics_inf(susie_ash_v2_default_output, X, y, causal)
  susie_ash_v2_quad_metrics <- calc_metrics_inf(susie_ash_v2_quad_output, X, y, causal)
  susie_inf_metrics <- calc_metrics_inf(susie_inf_output, X, y, causal)

  #### Create a Data Frame with the Results ####
  metrics_table <- data.frame(
    Model = c("SuSiE",
              "SuSiE-ash RE (default)",
              "SuSiE-ash RE (quadratic)",
              "SuSiE-ash RE v2 (default)",
              "SuSiE-ash RE v2 (quadratic)",
              "SuSiE-inf"),
    RMSE_y = c(susie_metrics$RMSE_y,
               susie_ash_default_metrics$RMSE_y,
               susie_ash_quad_metrics$RMSE_y,
               susie_ash_v2_default_metrics$RMSE_y,
               susie_ash_v2_quad_metrics$RMSE_y,
               susie_inf_metrics$RMSE_y),
    CS_FDR = c(susie_metrics$cs_fdr,
               susie_ash_default_metrics$cs_fdr,
               susie_ash_quad_metrics$cs_fdr,
               susie_ash_v2_default_metrics$cs_fdr,
               susie_ash_v2_quad_metrics$cs_fdr,
               susie_inf_metrics$cs_fdr),
    CS_Recall = c(susie_metrics$cs_recall,
                  susie_ash_default_metrics$cs_recall,
                  susie_ash_quad_metrics$cs_recall,
                  susie_ash_v2_default_metrics$cs_recall,
                  susie_ash_v2_quad_metrics$cs_recall,
                  susie_inf_metrics$cs_recall),
    CS_Size = c(susie_metrics$cs_size,
                susie_ash_default_metrics$cs_size,
                susie_ash_quad_metrics$cs_size,
                susie_ash_v2_default_metrics$cs_size,
                susie_ash_v2_quad_metrics$cs_size,
                susie_inf_metrics$cs_size),
    Coverage = c(susie_metrics$coverage,
                 susie_ash_default_metrics$coverage,
                 susie_ash_quad_metrics$coverage,
                 susie_ash_v2_default_metrics$coverage,
                 susie_ash_v2_quad_metrics$coverage,
                 susie_inf_metrics$coverage)
  )

  #### Return the Results Table ####
  return(list(
    metrics = metrics_table
  ))
}

# Simulation Function
simulation <- function(num_simulations = NULL,
                       h2_total = NULL,
                       prop_h2_sentinel = NULL,
                       L = NULL,
                       n_oligogenic = NULL,
                       pve_threshold = NULL,
                       mixture_small = NULL,
                       LD_blocks_dir = "LD_blocks_precomputed") {

  # Set Default Values
  if (is.null(num_simulations)) num_simulations <- 200
  if (is.null(h2_total)) h2_total <- 0.3
  if (is.null(prop_h2_sentinel)) prop_h2_sentinel <- 0.7
  if (is.null(L)) L <- 10
  if (is.null(n_oligogenic)) n_oligogenic <- 20
  if (is.null(pve_threshold)) pve_threshold <- 0.005
  if (is.null(mixture_small)) mixture_small <- 0.4

  # Parse Command-Line Arguments (if any)
  for (arg in commandArgs(trailingOnly = TRUE)) {
    eval(parse(text = arg))
  }

  # Generate Seeds: Seed i for replicate i
  all_seeds <- 1:num_simulations

  # Prepare List of LD Block Filenames
  ld_block_files <- list.files(path = LD_blocks_dir, pattern = "\\.rds$", full.names = TRUE)

  # Check if the number of LD blocks matches num_simulations
  if (length(ld_block_files) < num_simulations) {
    stop("Number of LD block files is less than num_simulations.")
  }

  # If more than needed, truncate the list
  if (length(ld_block_files) > num_simulations) {
    ld_block_files <- ld_block_files[1:num_simulations]
  }

  # Initialize Lists to Store Results
  all_metrics <- vector("list", num_simulations)
  all_betas <- vector("list", num_simulations)
  all_causal_indices <- vector("list", num_simulations)
  all_epsilons <- numeric(num_simulations)
  all_h2_estimated <- numeric(num_simulations)
  ld_block_names <- vector("character", num_simulations)

  # Loop Over Each Simulation Replicate
  for (i in 1:num_simulations) {
    cat("Running simulation", i, "out of", num_simulations, "\n")

    # Set Seed for Current Simulation
    seed <- all_seeds[i]

    # Read in Precomputed LD Block
    ld_block_file <- ld_block_files[i]
    cat("Processing LD block file:", ld_block_file, "\n")

    ld_block_names[i] <- basename(ld_block_file)

    # Load the Precomputed Data
    precomputed_data <- readRDS(ld_block_file)

    # Extract the Genotype Matrix and Precomputed Matrices
    X <- precomputed_data$X  # Unscaled, imputed genotype matrix
    XtX <- precomputed_data$XtX  # Computed using scaled X
    LD <- precomputed_data$LD
    V <- precomputed_data$V
    Dsq <- precomputed_data$Dsq

    # Store Precomputed Matrices
    precomputed_matrices <- list(
      XtX = XtX,
      LD = LD,
      V = V,
      Dsq = Dsq
    )

    mixture_props <- c(mixture_small, 1 - mixture_small)

    # Generate Data Using the eQTL Data Generation Function
    data <- generate_eqtl_data(
      X = X,
      h2_total = h2_total,
      prop_h2_sparse = 0.65,
      prop_h2_oligogenic = 0.20,
      prop_h2_infinitesimal = 0.15,
      prop_h2_sentinel = prop_h2_sentinel,
      n_oligogenic = n_oligogenic,
      mixture_props = mixture_props,
      mixture_sds = c(0.0025, 0.005),
      seed = seed
    )

    # Identify Causal SNPs Based on PVE Threshold
    data$causal <- is_causal(data, pve_threshold)

    # Run Methods and Calculate Metrics
    results <- method_and_score(
      X = data$ori.X,
      y = data$ori.y,
      beta = data$beta,
      causal = data$causal,
      L = L,
      precomputed_matrices = precomputed_matrices,
      seed = seed
    )

    # Store Results
    all_metrics[[i]] <- results$metrics
    all_betas[[i]] <- data$beta
    all_causal_indices[[i]] <- data$causal
    all_epsilons[i] <- data$var_epsilon
    all_h2_estimated[i] <- data$h2_total

    # Remove Large Objects to Free Memory
    rm(precomputed_data, X, XtX, LD, V, Dsq, precomputed_matrices, data, results)
    gc()
  }

  #### Calculate Average Metrics ####

  avg_metrics <- data.frame(
    Model = unique(all_metrics[[1]]$Model),
    RMSE_y = Reduce("+", lapply(all_metrics, function(x) x$RMSE_y)) / num_simulations,
    CS_FDR = Reduce("+", lapply(all_metrics, function(x) x$CS_FDR)) / num_simulations,
    CS_Recall = Reduce("+", lapply(all_metrics, function(x) x$CS_Recall)) / num_simulations,
    CS_Size = Reduce("+", lapply(all_metrics, function(x) x$CS_Size)) / num_simulations,
    Coverage = Reduce("+", lapply(all_metrics, function(x) x$Coverage)) / num_simulations
  )

  #### Save Simulation Results as RDS File ####
  output_dir <- "/home/apm2217/output"

  # Compile All Results into a List
  simulation_results <- list(
    avg_metrics = avg_metrics,
    all_metrics = all_metrics,
    all_betas = all_betas,
    all_causal_indices = all_causal_indices,
    all_epsilons = all_epsilons,
    all_h2_estimated = all_h2_estimated,
    ld_block_names = ld_block_names
  )

  # Create a Descriptive Filename
  file_name <- paste0("numIter", num_simulations,
                      "_h2total", h2_total,
                      "_h2sentinel", prop_h2_sentinel,
                      "_L", L,
                      "_numOligogenic", n_oligogenic,
                      "_pvethreshold", pve_threshold,
                      "_mixturesmall", mixture_small)

  # Save the Results
  saveRDS(simulation_results, file.path(output_dir, paste0(file_name, ".rds")))

  # Return All Results
  return(simulation_results)
}

# =====================================
# Run the Simulation
# =====================================

simulation_results <- simulation(
  num_simulations = NULL,  # Defaults to 200
  h2_total = NULL,         # Defaults to 0.3
  prop_h2_sentinel = NULL, # Defaults to 0.7
  L = NULL,                # Defaults to 10
  n_oligogenic = NULL,     # Defaults to 20
  pve_threshold = NULL,    # Defaults to 0.005
  mixture_small = NULL,    # Defaults to 0.4
  LD_blocks_dir = "LD_blocks_precomputed"  # Directory containing precomputed LD block files
)
