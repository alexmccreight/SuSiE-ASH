# Set Path to Mounted S3 Bucket
data_path <- "/home/apm2217/data/"

# Libraries
library(susieR)
library(mr.ash.alpha)
library(dplyr)
library(magrittr)
source("susie_ash_joint_ELBO_v3.R")
source("susie_inf.R")

# Annotation Matrix (from S3 Bucket)
X_full <- readRDS("X20")
all_seeds <- sample(1:1e9, 100, replace = FALSE)

# New Data Generation Function for Sparse Setting
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

  # Center y (as in your original code)
  ori.y <- y
  y <- scale(y, center = TRUE, scale = FALSE)

  # Return a list containing the data and causal indices
  return(list(
    X = X_scaled,
    ori.X = X,
    y = y,
    ori.y = ori.y,
    beta = beta,
    causal_indices = causal_indices,
    var_epsilon = var_e,
    h2 = h2
  ))
}

# Method and Metrics Function
method_and_score <- function(X, y, beta, causal, L = 10, v_threshold, precomputed_matrices, seed) {
  set.seed(seed)

  XtX <- precomputed_matrices$XtX
  LD <- precomputed_matrices$LD
  V <- precomputed_matrices$V
  Dsq <- precomputed_matrices$Dsq

  #### Run various methods ####
  cat("Starting SuSiE\n")
  susie_output <- susie(X = X, y = y, L = L, intercept = TRUE, standardize = FALSE, track_fit = TRUE)

  cat("Starting mr.ash\n")
  mrash_output <- mr.ash(X = X, y = y, sa2 = nrow(X) * (2^((0:19)/20) - 1)^2, intercept = TRUE, standardize = FALSE)

  cat("Starting SuSiE-ash (MLE)\n")
  susie_ash_mle_output <- susie_ash(X = X, y = y, L = L, tol = 0.001, intercept = TRUE, standardize = FALSE, est_var = "cal_v", true_var_res = NULL, v_threshold = v_threshold, track_fit = TRUE)

  cat("Starting SuSiE-ash (MoM)\n")
  susie_ash_mom_output <- susie_ash(X = X, y = y, L = L, tol = 0.001, intercept = TRUE, standardize = FALSE, est_var = "mom", true_var_res = NULL, v_threshold = v_threshold, track_fit = TRUE)

  cat("Starting SuSiE-inf\n")
  susie_inf_output <- susie_inf(X = X,
                                y = y,
                                L = L,
                                verbose = FALSE,
                                coverage = 0.95,
                                XtX = XtX,
                                LD = LD,
                                V = V,
                                Dsq = Dsq)

  calc_metrics <- function(mod, X, y, causal) {
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
      coverage <- (lapply(1:length(test.cs), function(cs.l) { ifelse(sum(causal %in% test.cs[[cs.l]]) != 0, TRUE, FALSE) }) %>% unlist() %>% sum()) / length(test.cs)

      # CS Based FDR
      TP_fdr <- lapply(1:length(test.cs), function(cs.l) { ifelse(sum(test.cs[[cs.l]] %in% causal) != 0, TRUE, FALSE) }) %>% unlist() %>% sum()
      FP_fdr <- length(test.cs) - TP_fdr
      denominator <- TP_fdr + FP_fdr
      if (denominator == 0) {
        cs_fdr <- NA
      } else {
        cs_fdr <- FP_fdr / denominator
      }

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

  calc_metrics_ash <- function(mod, X, y, causal) {
    #### Calculate RMSE ####
    RMSE_y <- sqrt(mean((y - predict(mod, X))^2))

    #### Store Results ####
    return(list(
      RMSE_y = RMSE_y,
      cs_size = NA,
      coverage = NA,
      cs_fdr = NA,
      cs_recall = NA
    ))
  }

  calc_metrics_inf <- function(mod, X, y, causal) {
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
      coverage <- (lapply(1:length(test.cs), function(cs.l) { ifelse(sum(causal %in% test.cs[[cs.l]]) != 0, TRUE, FALSE) }) %>% unlist() %>% sum()) / length(test.cs)

      # CS Based FDR
      TP_fdr <- lapply(1:length(test.cs), function(cs.l) { ifelse(sum(test.cs[[cs.l]] %in% causal) != 0, TRUE, FALSE) }) %>% unlist() %>% sum()
      FP_fdr <- length(test.cs) - TP_fdr
      denominator <- TP_fdr + FP_fdr
      if (denominator == 0) {
        cs_fdr <- NA
      } else {
        cs_fdr <- FP_fdr / denominator
      }

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

  #############
  # Calculate Metrics for each method
  susie_metrics <- calc_metrics(susie_output, X, y, causal)
  mrash_metrics <- calc_metrics_ash(mrash_output, X, y, causal)
  susie_ash_mle_metrics <- calc_metrics(susie_ash_mle_output, X, y, causal)
  susie_ash_mom_metrics <- calc_metrics(susie_ash_mom_output, X, y, causal)
  susie_inf_metrics <- calc_metrics_inf(susie_inf_output, X, y, causal)

  # Create a data frame with the results
  metrics_table <- data.frame(
    Model = c("SuSiE",
              "mr.ash",
              "SuSiE-ash (MLE)",
              "SuSiE-ash (MoM)",
              "SuSiE-inf"),
    RMSE_y = c(susie_metrics$RMSE_y,
               mrash_metrics$RMSE_y,
               susie_ash_mle_metrics$RMSE_y,
               susie_ash_mom_metrics$RMSE_y,
               susie_inf_metrics$RMSE_y),
    CS_FDR = c(susie_metrics$cs_fdr,
               mrash_metrics$cs_fdr,
               susie_ash_mle_metrics$cs_fdr,
               susie_ash_mom_metrics$cs_fdr,
               susie_inf_metrics$cs_fdr),
    CS_Recall = c(susie_metrics$cs_recall,
                  mrash_metrics$cs_recall,
                  susie_ash_mle_metrics$cs_recall,
                  susie_ash_mom_metrics$cs_recall,
                  susie_inf_metrics$cs_recall),
    CS_Size = c(susie_metrics$cs_size,
                mrash_metrics$cs_size,
                susie_ash_mle_metrics$cs_size,
                susie_ash_mom_metrics$cs_size,
                susie_inf_metrics$cs_size),
    Coverage = c(susie_metrics$coverage,
                 mrash_metrics$coverage,
                 susie_ash_mle_metrics$coverage,
                 susie_ash_mom_metrics$coverage,
                 susie_inf_metrics$coverage)
  )

  # Return the results table and outputs
  return(list(
    metrics = metrics_table,
    susie_output = susie_output,
    mrash_output = mrash_output,
    susie_ash_mle_output = susie_ash_mle_output,
    susie_ash_mom_output = susie_ash_mom_output,
    susie_inf_output = susie_inf_output,
    causal = causal,
    betas = beta
  ))
}

# Simulation Function
simulation <- function(num_simulations = NULL,
                       h2_total = NULL,
                       K = NULL,
                       L = NULL,
                       v_threshold = NULL,
                       sample_size = NULL) {

  # Set default values
  num_simulations <- 2
  h2_total <- 0.3
  K <- 10
  L <- 10
  v_threshold <- 0.005
  sample_size <- 5000

  # Parse command-line arguments (if any)
  for (arg in commandArgs(trailingOnly = TRUE)) {
    eval(parse(text = arg))
  }

  # If sample_size is not provided, default to the number of rows in X_full
  if (is.null(sample_size)) {
    sample_size <- nrow(X_full)
  }

  # Sample n rows from X_full without replacement
  if (sample_size > nrow(X_full)) {
    stop("sample_size cannot be greater than the number of rows in X_full")
  }
  sample_indices <- sample(1:nrow(X_full), sample_size, replace = FALSE)
  X_sampled <- X_full[sample_indices, ]

  # Scale the genotype matrix
  X_scaled <- scale(X_sampled)
  n_samples <- nrow(X_scaled)
  n_features <- ncol(X_scaled)

  # Precompute values for susie-inf
  XtX <- t(X_scaled) %*% X_scaled
  LD <- XtX / n_samples
  eig <- eigen(LD, symmetric = TRUE)
  V <- eig$vectors
  Dsq <- pmax(n_samples * eig$values, 0)

  # Store precomputed matrices
  precomputed_matrices <- list(
    XtX = XtX,
    LD = LD,
    V = V,
    Dsq = Dsq
  )

  # Initialize lists to store results
  all_metrics <- list()
  all_betas <- list()
  all_causal_indices <- list()
  all_susie_outputs <- list()
  all_mrash_outputs <- list()
  all_susie_ash_mle_outputs <- list()
  all_susie_ash_mom_outputs <- list()
  all_susie_inf_outputs <- list()
  all_epsilons <- numeric(num_simulations)

  for (i in 1:num_simulations) {
    cat("Running simulation", i, "out of", num_simulations, "\n")

    # Set random seed for each simulation
    seed <- all_seeds[i]

    # Generate data using the new function
    data <- generate_sparse_eqtl_data(
      X = X_sampled,
      K = K,
      h2 = h2_total,
      seed = seed
    )

    causal_indices <- data$causal_indices

    # Run methods and calculate metrics
    results <- method_and_score(
      X = data$X,
      y = data$y,
      beta = data$beta,
      causal = causal_indices,
      L = L,
      v_threshold = v_threshold,
      precomputed_matrices = precomputed_matrices,
      seed = seed
    )

    # Store results
    all_metrics[[i]] <- results$metrics
    all_betas[[i]] <- results$betas
    all_causal_indices[[i]] <- causal_indices
    all_susie_outputs[[i]] <- results$susie_output
    all_mrash_outputs[[i]] <- results$mrash_output
    all_susie_ash_mle_outputs[[i]] <- results$susie_ash_mle_output
    all_susie_ash_mom_outputs[[i]] <- results$susie_ash_mom_output
    all_susie_inf_outputs[[i]] <- results$susie_inf_output
    all_epsilons[i] <- data$var_epsilon
  }

  # Calculate average metrics
  avg_metrics <- data.frame(
    Model = unique(all_metrics[[1]]$Model),
    RMSE_y = Reduce("+", lapply(all_metrics, function(x) x$RMSE_y)) / num_simulations,
    CS_FDR = Reduce("+", lapply(all_metrics, function(x) x$CS_FDR)) / num_simulations,
    CS_Recall = Reduce("+", lapply(all_metrics, function(x) x$CS_Recall)) / num_simulations,
    CS_Size = Reduce("+", lapply(all_metrics, function(x) x$CS_Size)) / num_simulations,
    Coverage = Reduce("+", lapply(all_metrics, function(x) x$Coverage)) / num_simulations
  )

  # Save simulation results as Rds file
  output_dir <- "/home/apm2217/output"
  #output_dir <- "analysis"
  simulation_results <- list(
    avg_metrics = avg_metrics,
    all_metrics = all_metrics,
    all_betas = all_betas,
    all_causal_indices = all_causal_indices,
    all_susie_outputs = all_susie_outputs,
    all_mrash_outputs = all_mrash_outputs,
    all_susie_ash_mle_outputs = all_susie_ash_mle_outputs,
    all_susie_ash_mom_outputs = all_susie_ash_mom_outputs,
    all_susie_inf_outputs = all_susie_inf_outputs,
    all_seeds = all_seeds,
    all_epsilons = all_epsilons
  )

  file_name <- paste0("numIter", num_simulations,
                      "_h2total", h2_total,
                      "_K", K,
                      "_L", L,
                      "_vthreshold", v_threshold,
                      "_samplesize", sample_size)

  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Save the results
  saveRDS(simulation_results, file.path(output_dir, paste0(file_name, ".rds")))

  # Return all results
  return(simulation_results)
}

# Run the simulation
simulation_results <- simulation(
  num_simulations = NULL,
  h2_total = NULL,
  K = NULL,
  L = NULL,
  v_threshold = NULL,
  sample_size = NULL
)
