# ================================
# Simulation Script for Sparse eQTL Data
# ================================

# Load Required Libraries
library(susieR)
library(mr.ash.alpha)
library(dplyr)
library(magrittr)

# Source Custom Functions
source("susie-ash-data/susie_ash_mod.R")
source("susie-ash-data/susie_ash_mod_v2.R")
source("susie-ash-data/susie_inf.R")

# Function for Mean Imputation
mean_impute <- function(geno) {
  # Compute column-wise means excluding NAs
  col_means <- apply(geno, 2, function(x) mean(x, na.rm = TRUE))

  # Replace NAs with the corresponding column mean
  for (i in seq_along(col_means)) {
    na_indices <- which(is.na(geno[, i]))
    if (length(na_indices) > 0) {
      geno[na_indices, i] <- col_means[i]
    }
  }

  return(geno)
}

# Updated Data Generation Function for Sparse Setting
generate_sparse_eqtl_data <- function(X,
                                      K = 10,      # Number of non-zero effect SNPs
                                      h2 = 0.3,    # Total heritability
                                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_samples <- nrow(X)
  n_features <- ncol(X)

  # Scale the genotype matrix X
  X_scaled <- scale(X)

  # Initialize beta (effect sizes) to zero
  beta <- rep(0, n_features)

  # Randomly select K SNPs to have non-zero effect sizes
  causal_indices <- sample(1:n_features, K, replace = FALSE)

  # Assign effect sizes to the causal SNPs
  beta[causal_indices] <- rnorm(K, mean = 0, sd = 0.6)

  # Generate the genetic component of the phenotype
  genetic_effect <- X_scaled %*% beta

  # Scale genetic effect to achieve desired heritability
  var_g <- var(genetic_effect)
  var_e <- var_g * (1 - h2) / h2  # Residual variance
  noise <- rnorm(n_samples, mean = 0, sd = sqrt(var_e))

  # Generate the phenotype y
  y <- genetic_effect + noise

  # Center y
  ori.y <- y
  y <- scale(y, center = TRUE, scale = FALSE)

  # Recalculate variances after centering
  var_y <- var(y)
  var_g_estimated <- var(genetic_effect)
  h2_estimated <- var_g_estimated / var_y

  # Return a list containing the data and causal indices
  return(list(
    X = X_scaled,
    ori.X = X,
    y = y,
    ori.y = ori.y,
    beta = beta,
    causal_indices = causal_indices,
    var_epsilon = var_e,
    h2_input = h2,
    h2_estimated = h2_estimated
  ))
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
  cat("Starting SuSiE-ash (default)\n")
  susie_ash_default_output <- susie_ash_mod(
    X = scale(X),
    y = scale(y, center = T, scale = F),
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
  cat("Starting SuSiE-ash (quadratic)\n")
  susie_ash_quad_output <- susie_ash_mod(
    X = scale(X),
    y = scale(y, center = T, scale = F),
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
  cat("Starting SuSiE-ash v2 (default)\n")
  susie_ash_v2_default_output <- susie_ash_mod_v2(
    X = scale(X),
    y = scale(y, center = T, scale = F),
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
  cat("Starting SuSiE-ash v2 (quadratic)\n")
  susie_ash_v2_quad_output <- susie_ash_mod_v2(
    X = scale(X),
    y = scale(y, center = T, scale = F),
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
    y = scale(y, center = T, scale = F),
    L = L,
    verbose = FALSE,
    coverage = 0.95,
    XtX = XtX,
    LD = LD,
    V = V,
    Dsq = Dsq
  )

  #### Function to Calculate Metrics ####
  calc_metrics <- function(mod, X = X, y = y, causal = causal){
    #### Initialize values ####
    test.cs <- susie_get_cs(mod, X = X, coverage = 0.95)$cs
    coverage <- 0
    cs_fdr <- 0
    cs_recall <- 0
    cs_size <- 0

    if(length(test.cs) > 0){
      # Calculate Average CS Size
      cs_size <- length(unlist(test.cs)) / length(test.cs)

      # Calculate Coverage (proportion of credible sets with a causal effect)
      coverage <- (lapply(1:length(test.cs), function(cs.l){ ifelse(sum(causal %in% test.cs[[cs.l]]) != 0, T, F)}) %>% unlist(.) %>% sum(.)) / (length(test.cs))

      # CS Based FDR
      TP_fdr = lapply(1:length(test.cs), function(cs.l){ ifelse(sum(test.cs[[cs.l]] %in% causal)!=0,T,F)}) %>% unlist(.) %>% sum(.)
      FP_fdr = length(test.cs) - TP_fdr
      FN_fdr = length(causal) - sum(causal %in% unlist(test.cs))
      cs_fdr = FP_fdr/(TP_fdr+FP_fdr)

      # CS Based Recall
      TP_recall = sum(causal %in% unlist(test.cs))
      FN_recall = length(causal) - TP_recall
      cs_recall = TP_recall/(TP_recall+FN_recall)
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

  calc_metrics_inf <- function(mod, X = X, y = y, causal = causal){

    #### Initialize values ####
    test.cs <- mod$sets
    coverage <- 0
    cs_fdr <- 0
    cs_recall <- 0
    cs_size <- 0

    if(length(test.cs) > 0){
      # Calculate Average CS Size
      cs_size <- length(unlist(test.cs)) / length(test.cs)

      # Calculate Coverage (proportion of credible sets with a causal effect)
      coverage <- (lapply(1:length(test.cs), function(cs.l){ ifelse(sum(causal %in% test.cs[[cs.l]]) != 0, T, F)}) %>% unlist(.) %>% sum(.)) / (length(test.cs))

      # CS Based FDR
      TP_fdr = lapply(1:length(test.cs), function(cs.l){ ifelse(sum(test.cs[[cs.l]] %in% causal)!=0,T,F)}) %>% unlist(.) %>% sum(.)
      FP_fdr = length(test.cs) - TP_fdr
      FN_fdr = length(causal) - sum(causal %in% unlist(test.cs))
      cs_fdr = FP_fdr/(TP_fdr+FP_fdr)

      # CS Based Recall
      TP_recall = sum(causal %in% unlist(test.cs))
      FN_recall = length(causal) - TP_recall
      cs_recall = TP_recall/(TP_recall+FN_recall)
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
  metrics_table  <- data.frame(
    Model = c("SuSiE",
              "SuSiE-ash (default)",
              "SuSiE-ash (quadratic)",
              "SuSiE-ash v2 (default)",
              "SuSiE-ash v2 (quadratic)",
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
                       K = NULL,
                       L = NULL,
                       LD_blocks_dir = "LD_blocks") {

  # Set Default Values
  if (is.null(num_simulations)) num_simulations <- 200
  if (is.null(h2_total)) h2_total <- 0.3
  if (is.null(K)) K <- 10
  if (is.null(L)) L <- 10

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
  ld_block_names <- vector("character", num_simulations)  # New line to store LD block names


  # Loop Over Each Simulation Replicate
  for (i in 1:num_simulations) {
    cat("Running simulation", i, "out of", num_simulations, "\n")

    # Set Seed for Current Simulation
    seed <- all_seeds[i]

    cat("Processing LD block file:", ld_block_file, "\n")
    ld_block_names[i] <- basename(ld_block_file)

    # Load the LD Block
    ld_block <- readRDS(ld_block_files[i])

    # Extract and Impute Genotype Matrix
    X <- mean_impute(ld_block$genotypes)

    # Precompute Matrices for Current Genotype Matrix
    n_samples <- nrow(X)
    n_features <- ncol(X)

    X_scaled <- scale(X)
    XtX <- t(X_scaled) %*% X_scaled
    LD <- XtX / n_samples
    eig <- eigen(LD, symmetric = TRUE)
    V <- eig$vectors
    Dsq <- pmax(n_samples * eig$values, 0)

    # Store Precomputed Matrices
    precomputed_matrices <- list(
      XtX = XtX,
      LD = LD,
      V = V,
      Dsq = Dsq
    )

    # Generate Data Using the Updated Data Generation Function
    data <- generate_sparse_eqtl_data(
      X = X,
      K = K,
      h2 = h2_total,
      seed = seed
    )

    causal_indices <- data$causal_indices

    # Run Methods and Calculate Metrics
    results <- method_and_score(
      X = data$ori.X,
      y = data$ori.y,
      beta = data$beta,
      causal = causal_indices,
      L = L,
      precomputed_matrices = precomputed_matrices,
      seed = seed
    )

    # Store Results
    all_metrics[[i]] <- results$metrics
    all_betas[[i]] <- data$beta
    all_causal_indices[[i]] <- causal_indices
    all_epsilons[i] <- data$var_epsilon
    all_h2_estimated[i] <- data$h2_estimated

    # Remove Large Objects to Free Memory
    rm(ld_block, X, X_scaled, XtX, LD, eig, V, Dsq, precomputed_matrices, data, results)
    gc()
  }

  #### Calculate Average Metrics ####
  # Initialize a Data Frame with Model Names
  model_names <- all_metrics[[1]]$Model
  avg_metrics <- data.frame(
    Model = model_names,
    RMSE_y = numeric(length(model_names)),
    CS_FDR = numeric(length(model_names)),
    CS_Recall = numeric(length(model_names)),
    CS_Size = numeric(length(model_names)),
    Coverage = numeric(length(model_names))
  )

  # Aggregate Metrics Across All Simulations
  for (metric in c("RMSE_y", "CS_FDR", "CS_Recall", "CS_Size", "Coverage")) {
    avg_metrics[[metric]] <- rowMeans(
      do.call(rbind, lapply(all_metrics, function(x) x[[metric]]))
    )
  }

  #### Save Simulation Results as RDS File ####
  output_dir <- "/home/apm2217/output"  # Adjust if necessary

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
                      "_K", K,
                      "_L", L)

  # Save the Results
  saveRDS(simulation_results, file.path(output_dir, paste0(file_name, ".rds")))

  # Return All Results
  return(simulation_results)
}

# ================================
# Run the Simulation
# ================================

simulation_results <- simulation(
  num_simulations = NULL,  # Defaults to 200
  h2_total = NULL,         # Defaults to 0.3
  K = NULL,                # Defaults to 10
  L = NULL,                # Defaults to 10
  LD_blocks_dir = "LD_blocks"  # Directory containing LD block files
)
