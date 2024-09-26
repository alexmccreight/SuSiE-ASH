# Set Path to Mounted S3 Bucket
data_path <- "/home/apm2217/data/"

# Libraries
library(susieR)
library(mr.ash.alpha)
library(dplyr)
library(magrittr)
#source("susie-ash.R")
source("susie_inf.R")

# Annotation Matrix (from S3 Bucket)
X <- readRDS("X20")

# Generate Data Function
generate_eqtl_data <- function(X,
                               h2_total = 0.3,            # Total heritability
                               prop_h2_sparse = 0.65,     # Proportion of h2_total explained by sparse effects (including sentinel)
                               prop_h2_oligogenic = 0.20, # Proportion of h2_total explained by oligogenic effects
                               prop_h2_infinitesimal = 0.15, # Proportion of h2_total explained by infinitesimal effects
                               prop_h2_sentinel = 0.7,    # Proportion of h2_sparse explained by sentinel SNP
                               mixture_props = c(0.7, 0.2, 0.1),  # Mixture proportions for oligogenic effects
                               mixture_sds = c(0.01, 0.03, 0.05), # Standard deviations for mixture components
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

  # Sentinel SNP effect
  sentinel_index <- sample(1:n_features, 1)
  beta[sentinel_index] <- rnorm(1, 0, sqrt(h2_sentinel))

  # Other sparse effects
  n_other_sparse <- rpois(1, lambda = 2)  # Random number of other sparse effects
  other_sparse_indices <- sample((1:n_features)[-sentinel_index], n_other_sparse)
  beta[other_sparse_indices] <- rnorm(n_other_sparse, 0, as.vector(sqrt(h2_other_sparse / n_other_sparse)))

  # Oligogenic effects (mixture of normals)
  oligogenic_indices <- setdiff(1:n_features, c(sentinel_index, other_sparse_indices))
  mixture_assignments <- sample(1:length(mixture_props), length(oligogenic_indices), replace = TRUE, prob = mixture_props)
  beta[oligogenic_indices] <- rnorm(length(oligogenic_indices), 0, mixture_sds[mixture_assignments])

  # Scale oligogenic effects to achieve desired heritability
  oligogenic_effects <- X[, oligogenic_indices] %*% beta[oligogenic_indices]
  scaling_factor <- sqrt(h2_oligogenic / var(oligogenic_effects))
  beta[oligogenic_indices] <- beta[oligogenic_indices] * as.vector(scaling_factor)

  # Infinitesimal effects (small effects on all SNPs)
  infinitesimal_effects <- rnorm(n_features, 0, sqrt(h2_infinitesimal / n_features))
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

  # Check if other_sparse_indices is not empty
  if (length(other_sparse_indices) > 0) {
    sparse_indices <- c(sentinel_index, other_sparse_indices)
    h2_sparse_actual <- var(X[, sparse_indices] %*% beta[sparse_indices]) / var_y_total
  } else {
    h2_sparse_actual <- h2_sentinel_actual  # If no other sparse effects, sparse is just the sentinel
  }

  h2_oligogenic_actual <- var(X[, oligogenic_indices] %*% beta[oligogenic_indices]) / var_y_total
  h2_infinitesimal_actual <- var(X %*% infinitesimal_effects) / var_y_total
  h2_total_actual <- var(as.vector(X %*% beta)) / var_y_total

  ori.y <- y
  y <- scale(y, center = T, scale = F)

  return(list(
    ori.X = ori.X,
    X = X,
    ori.y = y,
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
    mixture_assignments = mixture_assignments,
    var_epsilon = var_epsilon,
    sparse_indices = c(sentinel_index, other_sparse_indices)
  ))
}

# Method and Metrics
method_and_score <- function(X = data$X, y = data$y, beta = data$beta, sparse = data$sparse_indices, oligo = data$oligogenic_indices, L = 10) {
  #### Run various methods ####
  cat("Starting SuSiE\n")
  susie_output <- susie(X = X, y = y, L = L, intercept = F, standardize = F)

  cat("Starting mr.ash\n")
  mrash_output <- mr.ash(X = X, y = y, sa2 = nrow(X) * (2^((0:19)/20) - 1)^2, intercept = F, standardize = F)

  #cat("Starting SuSiE-ash\n")
  #susie_ash_output <- susie_ash(X = X, y = y, L = L, warm_start = 2, tol = 0.001, intercept = F, standardize = F)

  cat("Starting SuSiE-inf\n")
  susie_inf_output <- susie_inf(X = X, y = y, L = L, verbose = F, coverage = 0.95)


  calc_metrics <- function(mod, X = X, y = y, sparse = sparse, oligo = oligo){
    #### Set truth ####
    sparse_causal <- sparse
    sparse_oligo_causal <- c(sparse, oligo)


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
      coverage <- (lapply(1:length(test.cs), function(cs.l){ ifelse(sum(sparse_oligo_causal %in% test.cs[[cs.l]]) != 0, T, F)}) %>% unlist(.) %>% sum(.)) / (length(test.cs))

      # CS Based FDR
      TP_fdr = sum(sparse_oligo_causal %in% unlist(test.cs))
      FN_fdr = length(sparse_oligo_causal) - TP_fdr
      FP_fdr = length(test.cs) - lapply(1:length(test.cs), function(cs.l){ ifelse(sum(test.cs[[cs.l]] %in% sparse_oligo_causal)!=0,T,F)}) %>% unlist(.) %>% sum(.)
      cs_fdr = FP_fdr/(TP_fdr+FP_fdr)

      # CS Based Recall
      TP_recall = sum(sparse_causal %in% unlist(test.cs))
      FN_recall = length(sparse_causal) - TP_recall
      FP_recall = length(test.cs) - lapply(1:length(test.cs), function(cs.l){ ifelse(sum(test.cs[[cs.l]] %in% sparse_causal)!=0,T,F)}) %>% unlist(.) %>% sum(.)
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

  calc_metrics_ash <- function(mod, X = X, y = y, sparse = sparse, oligo = oligo){
    #### Set up truth ####
    sparse_causal <- sparse
    sparse_oligo_causal <- c(sparse, oligo)

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

  calc_metrics_inf <- function(mod, X = X, y = y, sparse = sparse, oligo = oligo){
    #### Set truth ####
    sparse_causal <- sparse
    sparse_oligo_causal <- c(sparse, oligo)

    #### Initialize values ####
    test.cs <- mod$sets
    coverage <- 0
    cs_fdr <- 0
    cs_recall <- 0
    cs_size <- 0

    if(length(test.cs) > 0){
      # Calculate Average CS Size
      cs_size <- length(unlist(mod$sets)) / length(mod$sets)

      # Calculate Coverage (proportion of credible sets with a causal effect)
      coverage <- (lapply(1:length(test.cs), function(cs.l){ ifelse(sum(sparse_oligo_causal %in% test.cs[[cs.l]]) != 0, T, F)}) %>% unlist(.) %>% sum(.)) / length(mod$sets)

      # CS Based FDR
      TP_fdr = sum(sparse_oligo_causal %in% unlist(test.cs))
      FN_fdr = length(sparse_oligo_causal) - TP_fdr
      FP_fdr = length(test.cs) - lapply(1:length(test.cs), function(cs.l){ ifelse(sum(test.cs[[cs.l]] %in% sparse_oligo_causal)!=0,T,F)}) %>% unlist(.) %>% sum(.)
      cs_fdr = FP_fdr/(TP_fdr+FP_fdr)

      # CS Based Recall
      TP_recall = sum(sparse_causal %in% unlist(test.cs))
      FN_recall = length(sparse_causal) - TP_recall
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
  susie_metrics <- calc_metrics(susie_output, X, y, sparse, oligo)
  mrash_metrics <- calc_metrics_ash(mrash_output, X, y, sparse, oligo)
  #susie_ash_metrics <- calc_metrics(susie_ash_output, X, y, sparse, oligo)
  susie_inf_metrics <- calc_metrics_inf(susie_inf_output, X, y, sparse, oligo)

  #Create a data frame with the results
  metrics_table  <- data.frame(
    Model = c("SuSiE",
              "mr.ash",
              #"SuSiE-ash",
              "SuSiE-inf"),
    RMSE_y = c(susie_metrics$RMSE_y,
               mrash_metrics$RMSE_y,
               #susie_ash_metrics$RMSE_y,
               susie_inf_metrics$RMSE_y),
    CS_FDR = c(susie_metrics$cs_fdr,
               mrash_metrics$cs_fdr,
               #susie_ash_metrics$cs_fdr,
               susie_inf_metrics$cs_fdr),
    CS_Recall = c(susie_metrics$cs_recall,
                  mrash_metrics$cs_recall,
                  #susie_ash_metrics$cs_recall,
                  susie_inf_metrics$cs_recall),
    CS_Size = c(susie_metrics$cs_size,
                mrash_metrics$cs_size,
                #susie_ash_metrics$cs_size,
                susie_inf_metrics$cs_size),
    Coverage = c(susie_metrics$coverage,
                 mrash_metrics$coverage,
                 #susie_ash_metrics$coverage,
                 susie_inf_metrics$coverage)
  )
  # Return the results table
  return(list(
    metrics = metrics_table,
    susie_output = susie_output,
    mrash_output = mrash_output,
    #susie_ash_output = susie_ash_output,
    susie_inf_output = susie_inf_output,
    sparse_indices = sparse,
    oligogenic_indices = oligo,
    betas = beta)
  )
}


# Simulation Function
simulation <- function(num_simulations = NULL,
                       h2_total = NULL,
                       prop_h2_sentinel = NULL,
                       L = NULL) {
  # Parse command-line arguments
  num_simulations = 2
  h2_total = 0.3
  prop_h2_sentinel = 0.7
  L = 10

  for (arg in commandArgs(trailingOnly = TRUE)) {
    eval(parse(text=arg))
  }

  # Initialize lists to store results
  all_metrics <- list()
  all_betas <- list()
  all_sparse_indices <- list()
  all_oligogenic_indices <- list()
  all_susie_outputs <- list()
  all_mrash_outputs <- list()
  #all_susie_ash_outputs <- list()
  all_susie_inf_outputs <- list()
  all_seeds <- numeric(num_simulations)
  all_epsilons <- numeric(num_simulations)

  for (i in 1:num_simulations) {
    cat("Running simulation", i, "out of", num_simulations, "\n")

    # Set random seed for each simulation
    seed <- round(runif(1, min = 1, max = 1000000000))
    set.seed(seed)

    # Generate data
    data <- generate_eqtl_data(X,
                               h2_total = h2_total,       # Total heritability
                               prop_h2_sparse = 0.65,     # Proportion of h2_total explained by sparse effects (including sentinel)
                               prop_h2_oligogenic = 0.20, # Proportion of h2_total explained by oligogenic effects
                               prop_h2_infinitesimal = 0.15, # Proportion of h2_total explained by infinitesimal effects
                               prop_h2_sentinel = prop_h2_sentinel,    # Proportion of h2_sparse explained by sentinel SNP
                               mixture_props = c(0.7, 0.2, 0.1),  # Mixture proportions for oligogenic effects
                               mixture_sds = c(0.01, 0.03, 0.05), # Standard deviations for mixture components
                               seed = seed)

    # Run methods and calculate metrics
    results <- method_and_score(X = data$X, y = data$y, beta = data$beta, sparse = data$sparse_indices, oligo = data$oligogenic_indices, L = L)

    # Store results
    all_metrics[[i]] <- results$metrics
    all_betas[[i]] <- results$beta
    all_sparse_indices[[i]] <- results$sparse_indices
    all_oligogenic_indices[[i]] <- results$oligogenic_indices
    all_susie_outputs[[i]] <- results$susie_output
    all_mrash_outputs[[i]] <- results$mrash_output
    #all_susie_ash_outputs[[i]] <- results$susie_ash_output
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
    all_sparse_indices = all_sparse_indices,
    all_oligogenic_indices = all_oligogenic_indices,
    all_susie_outputs = all_susie_outputs,
    all_mrash_outputs = all_mrash_outputs,
    #all_susie_ash_outputs = all_susie_ash_outputs,
    all_susie_inf_outputs = all_susie_inf_outputs,
    all_seeds = all_seeds,
    all_epsilons = all_epsilons
  )

  file_name <- paste0("numIter", num_simulations,
                      "_h2total", h2_total,
                      "_h2sentinel", prop_h2_sentinel,
                      "_L", L)

  saveRDS(simulation_results, file.path(output_dir, file_name))

  # Return all results
  return(simulation_results)
}

# Run the simulation
simulation_results <- simulation(num_simulations = NULL,
                                 h2_total = NULL,
                                 prop_h2_sentinel = NULL,
                                 L = NULL)
