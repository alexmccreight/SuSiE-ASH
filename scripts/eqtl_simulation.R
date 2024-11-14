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

# Generate Data Function
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

# Method and Metrics
method_and_score <- function(X = data$ori.X, y = data$ori.y, beta = data$beta, causal = data$causal, L = 10, v_threshold = v_threshold, precomputed_matrices = precomputed_matrices, seed) {

  set.seed(seed)

  XtX <- precomputed_matrices$XtX
  LD <- precomputed_matrices$LD
  V <- precomputed_matrices$V
  Dsq <- precomputed_matrices$Dsq

  #### Run various methods ####
  cat("Starting SuSiE\n")
  susie_output <- susie(X = X, y = y, L = L, intercept = T, standardize = T, track_fit = T)

  cat("Starting mr.ash\n")
  mrash_output <- mr.ash(X = X, y = y, sa2 = nrow(X) * (2^((0:19)/20) - 1)^2, intercept = T, standardize = T)

  cat("Starting SuSiE-ash (MLE)\n")
  susie_ash_mle_output <- susie_ash(X = X, y = y, L = L, tol = 0.001, intercept = T, standardize = T, est_var = "cal_v", true_var_res = NULL, v_threshold = v_threshold, track_fit = T)

  cat("Starting SuSiE-ash (MoM)\n")
  susie_ash_mom_output <- susie_ash(X = X, y = y, L = L, tol = 0.001, intercept = T, standardize = T, est_var = "mom", true_var_res = NULL, v_threshold = v_threshold, track_fit = T)

  cat("Starting SuSiE-inf\n")
  #susie_inf_output <- susie_inf(X = scale(X), y = scale(y, center = T, scale = F), L = L, verbose = F, coverage = 0.95)
  susie_inf_output <- susie_inf(X = scale(X),
                                y = scale(y, center = T, scale = F),
                                L = L,
                                verbose = F,
                                coverage = 0.95,
                                XtX = XtX,
                                LD = LD,
                                V = V,
                                Dsq = Dsq
                                )

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

  calc_metrics_ash <- function(mod, X = X, y = y, causal = causal){
    #### Set up truth ####
    causal <- causal

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

  #############
  # Calculate Metrics for each method
  susie_metrics <- calc_metrics(susie_output, X, y, causal)
  mrash_metrics <- calc_metrics_ash(mrash_output, X, y, causal)
  susie_ash_mle_metrics <- calc_metrics(susie_ash_mle_output, X, y, causal)
  susie_ash_mom_metrics <- calc_metrics(susie_ash_mom_output, X, y, causal)
  susie_inf_metrics <- calc_metrics_inf(susie_inf_output, X, y, causal)



  #Create a data frame with the results
  metrics_table  <- data.frame(
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
  # Return the results table
  return(list(
    metrics = metrics_table,
    susie_output = susie_output,
    mrash_output = mrash_output,
    susie_ash_mle_output = susie_ash_mle_output,
    susie_ash_mom_output = susie_ash_mom_output,
    susie_inf_output = susie_inf_output,
    causal = causal,
    betas = beta)
  )
}


# Simulation Function
simulation <- function(num_simulations = NULL,
                       h2_total = NULL,
                       prop_h2_sentinel = NULL,
                       L = NULL,
                       n_oligogenic = NULL,
                       v_threshold = NULL,
                      # sample_size = NULL,
                       pve_threshold = NULL) {

  # Parse command-line arguments
  num_simulations = 2
  h2_total = 0.3
  prop_h2_sentinel = 0.7
  L = 10
  n_oligogenic = 20
  v_threshold = 0.005
  #sample_size = 1000
  pve_threshold = 0.005

  for (arg in commandArgs(trailingOnly = TRUE)) {
    eval(parse(text=arg))
  }

  # If sample_size is not provided, default to the number of rows in X_full
  # if (is.null(sample_size)) {
  #   sample_size <- nrow(X_full)
  # }
  #
  # # Sample n rows from X_full without replacement
  # if (sample_size > nrow(X_full)) {
  #  stop("n cannot be greater than the number of rows in X_full")
  # }
  # sample_indices <- sample(1:nrow(X_full), sample_size, replace = FALSE)
  # X_sampled <- X_full[sample_indices, ]

  # Precompute values for susie-inf
  scaled_X_sampled <- scale(X_sampled)
  n_samples <- nrow(scaled_X_sampled)
  XtX <- t(scaled_X_sampled) %*% scaled_X_sampled
  LD <- XtX / n_samples
  eig <- eigen(LD, symmetric = TRUE)
  V <- (eig$vectors[, ncol(eig$vectors):1])  # pxp matrix of eigenvectors of XtX
  Dsq <- pmax(n_samples * sort(eig$values), 0)

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
  #all_seeds <- numeric(num_simulations)
  all_epsilons <- numeric(num_simulations)

  for (i in 1:num_simulations) {
    cat("Running simulation", i, "out of", num_simulations, "\n")

    # Set random seed for each simulation
    seed <- all_seeds[i]

    # Generate data
    data <- generate_eqtl_data(X_sampled,
                               h2_total = h2_total,       # Total heritability
                               prop_h2_sparse = 0.65,     # Proportion of h2_total explained by sparse effects (including sentinel)
                               prop_h2_oligogenic = 0.20, # Proportion of h2_total explained by oligogenic effects
                               prop_h2_infinitesimal = 0.15, # Proportion of h2_total explained by infinitesimal effects
                               prop_h2_sentinel = prop_h2_sentinel,    # Proportion of h2_sparse explained by sentinel SNP
                               n_oligogenic = n_oligogenic,
                               mixture_props = c(0.6, 0.4), # Adjusted mixture proportions
                               mixture_sds = c(0.0025, 0.005), # Standard deviations for mixture components
                               seed = seed)

    data$causal <- is_causal(data, pve_threshold)


    # Run methods and calculate metrics
    results <- method_and_score(X = data$ori.X, y = data$ori.y, beta = data$beta, causal = data$causal, L = L, v_threshold = v_threshold, precomputed_matrices = precomputed_matrices, seed = seed)

    # Store results
    all_metrics[[i]] <- results$metrics
    all_betas[[i]] <- results$beta
    all_causal_indices[[i]] <- results$causal
    all_susie_outputs[[i]] <- results$susie_output
    all_mrash_outputs[[i]] <- results$mrash_output
    all_susie_ash_mle_outputs[[i]] <- results$susie_ash_mle_output
    all_susie_ash_mom_outputs[[i]] <- results$susie_ash_mom_output
    all_susie_inf_outputs[[i]] <- results$susie_inf_output
    all_seeds[i] <- seed
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
  #output_dir <- "/home/apm2217/output"
  output_dir <- "analysis"
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
                      "_h2sentinel", prop_h2_sentinel,
                      "_L", L,
                      "_numOligogenic", n_oligogenic,
                      "_vthreshold", v_threshold,
                    #  "_samplesize", sample_size,
                      "_pvethreshold", pve_threshold)

  saveRDS(simulation_results, file.path(output_dir, file_name))

  # Return all results
  return(simulation_results)
}

# Run the simulation
simulation_results <- simulation(num_simulations = NULL,
                                 h2_total = NULL,
                                 prop_h2_sentinel = NULL,
                                 L = NULL,
                                 n_oligogenic = NULL,
                                 v_threshold = NULL,
                               #  sample_size = NULL,
                                 pve_threshold = NULL)
