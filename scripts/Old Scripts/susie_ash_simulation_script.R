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
X <- readRDS("X4")

# Generate Data Function
generate_data <- function(X,
                          n_large_effects = 3,
                          n_medium_effects = 3,
                          total_pve = 0.4,
                          seed = NULL) {

  # Set seed for reproducibility
  if (!is.null(seed)) set.seed(seed)

  # Center and scale X
  X <- scale(X)

  n_samples <- nrow(X)
  n_features <- ncol(X)

  # Calculate effect sizes
  large_pve <- 0.6 * total_pve
  medium_pve <- 0.15 * total_pve
  small_pve <- 0.25 * total_pve

  # Generate effect sizes
  beta <- rep(0, n_features)
  phi <- rep(0, n_features)
  theta <- rnorm(n_features, 0, sqrt(small_pve / n_features))

  # Set large effect sizes
  large_indices <- sample(1:n_features, n_large_effects)
  beta[large_indices] <- rnorm(n_large_effects, 0, 0.25)

  # Set medium effect sizes
  medium_indices <- sample(setdiff(1:n_features, large_indices), n_medium_effects)
  phi[medium_indices] <- rnorm(n_medium_effects, 0, 0.15)

  # Generate components of y
  Xbeta <- X %*% beta
  Xphi <- X %*% phi
  Xtheta <- X %*% theta

  # Calculate initial variances
  var_Xbeta <- var(as.vector(Xbeta))
  var_Xphi <- var(as.vector(Xphi))
  var_Xtheta <- var(as.vector(Xtheta))

  # Scale components to match desired PVE
  Xbeta <- Xbeta * sqrt(large_pve / var_Xbeta)
  Xphi <- Xphi * sqrt(medium_pve / var_Xphi)
  Xtheta <- Xtheta * sqrt(small_pve / var_Xtheta)

  # Combine components
  y <- Xbeta + Xphi + Xtheta

  # Calculate actual PVE for each component
  var_y <- var(as.vector(y))
  pve_beta <- var(as.vector(Xbeta)) / var_y
  pve_phi <- var(as.vector(Xphi)) / var_y
  pve_theta <- var(as.vector(Xtheta)) / var_y

  # Add noise to achieve desired total PVE
  var_epsilon <- var_y * (1 - total_pve) / total_pve
  epsilon <- rnorm(n_samples, 0, sqrt(var_epsilon))
  y <- y + epsilon

  # Recalculate total variance and PVEs after adding noise
  var_y_total <- var(as.vector(y))
  pve_beta <- var(as.vector(Xbeta)) / var_y_total
  pve_phi <- var(as.vector(Xphi)) / var_y_total
  pve_theta <- var(as.vector(Xtheta)) / var_y_total
  pve_total <- pve_beta + pve_phi + pve_theta

  # Scale + center y (for no intercept)
  y <- scale(y, center = T, scale = F)

  # Return results
  return(list(
    X = X,
    y = y,
    var_epsilon = var_epsilon,
    beta = beta,
    phi = phi,
    theta = theta,
    pve_beta = pve_beta,
    pve_phi = pve_phi,
    pve_theta = pve_theta,
    pve_total = pve_total
  ))
}


# Method and Metrics
method_and_score <- function(X = data$X, y = data$y, beta = data$beta, phi = data$phi, L = 10) {
  #### Run various methods ####
  cat("Starting SuSiE\n")
  susie_output <- susie(X = X, y = y, L = L, intercept = F, standardize = F)

  cat("Starting mr.ash\n")
  mrash_output <- mr.ash(X = X, y = y, sa2 = nrow(X) * (2^((0:19)/20) - 1)^2, intercept = F, standardize = F)

  #cat("Starting SuSiE-ash\n")
  #susie_ash_output <- susie_ash(X = X, y = y, L = L, warm_start = 2, tol = 0.001, intercept = F, standardize = F)

  cat("Starting SuSiE-inf\n")
  susie_inf_output <- susie_inf(X = X, y = y, L = L, verbose = F, coverage = 0.95)


  calc_metrics <- function(mod, X = X, y = y, beta = beta, phi = phi){
    #### Set truth ####
    all_causal <- c(which(beta != 0), which(phi != 0))

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
      coverage <- (lapply(1:length(test.cs), function(cs.l){ ifelse(sum(all_causal %in% test.cs[[cs.l]]) != 0, T, F)}) %>% unlist(.) %>% sum(.)) / (length(test.cs))

      # CS Based FDR
      TP = sum(all_causal %in% unlist(test.cs))
      FN = length(all_causal) - TP
      FP = length(test.cs) - lapply(1:length(test.cs), function(cs.l){ ifelse(sum(test.cs[[cs.l]] %in% all_causal)!=0,T,F)}) %>% unlist(.) %>% sum(.)
      cs_fdr = FP/(TP+FP)

      # CS Based Recall
      TP = sum(all_causal %in% unlist(test.cs))
      FN = length(all_causal) - TP
      FP = length(test.cs) - lapply(1:length(test.cs), function(cs.l){ ifelse(sum(test.cs[[cs.l]] %in% all_causal)!=0,T,F)}) %>% unlist(.) %>% sum(.)
      cs_recall = TP/(TP+FN)
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

  calc_metrics_ash <- function(mod, X = X, y = y, beta = beta, phi = phi){
    #### Set up truth ####
    Xbeta <- X %*% matrix(beta, ncol = 1)
    Xphi <- X %*% matrix(phi, ncol = 1)

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

  calc_metrics_inf <- function(mod, X = X, y = y, beta = beta, phi = phi){
    #### Set truth ####
    all_causal <- c(which(beta != 0), which(phi != 0))

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
      coverage <- (lapply(1:length(test.cs), function(cs.l){ ifelse(sum(all_causal %in% test.cs[[cs.l]]) != 0, T, F)}) %>% unlist(.) %>% sum(.)) / length(mod$sets)

      # CS Based FDR
      TP = sum(all_causal %in% unlist(test.cs))
      FN = length(all_causal) - TP
      FP = length(test.cs) - lapply(1:length(test.cs), function(cs.l){ ifelse(sum(test.cs[[cs.l]] %in% all_causal)!=0,T,F)}) %>% unlist(.) %>% sum(.)
      cs_fdr = FP/(TP+FP)

      # CS Based Recall
      TP = sum(all_causal %in% unlist(test.cs))
      FN = length(all_causal) - TP
      cs_recall = TP/(TP+FN)
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
  susie_metrics <- calc_metrics(susie_output, X, y, beta, phi)
  mrash_metrics <- calc_metrics_ash(mrash_output, X, y, beta, phi)
  #susie_ash_metrics <- calc_metrics(susie_ash_output, X, y, beta, phi)
  susie_inf_metrics <- calc_metrics_inf(susie_inf_output, X, y, beta, phi)

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
    betas = beta,
    phis = phi)
  )
}


# Simulation Function
simulation <- function(num_simulations = NULL,
                       n_large_effects = NULL,
                       n_medium_effects = NULL,
                       total_pve = NULL,
                       L = NULL) {
  # Parse command-line arguments
  num_simulations = 2
  n_large_effects = 3
  n_medium_effects = 3
  total_pve = 0.4
  L = 10

  for (arg in commandArgs(trailingOnly = TRUE)) {
    eval(parse(text=arg))
  }

  # Initialize lists to store results
  all_metrics <- list()
  all_betas <- list()
  all_phis <- list()
  all_susie_outputs <- list()
  all_mrash_outputs <- list()
  #all_susie_ash_outputs <- list()
  all_susie_inf_outputs <- list()
  all_seeds <- numeric(num_simulations)

  for (i in 1:num_simulations) {
    cat("Running simulation", i, "out of", num_simulations, "\n")

    # Set random seed for each simulation
    seed <- abs(round(rnorm(1, mean = 0, sd = 10000)))
    set.seed(seed)

    # Generate data
    data <- generate_data(X, n_large_effects, n_medium_effects, total_pve, seed = seed)

    # Run methods and calculate metrics
    results <- method_and_score(X = data$X, y = data$y, beta = data$beta, phi = data$phi, L = L)

    # Store results + betas/thetas
    all_metrics[[i]] <- results$metrics
    all_betas[[i]] <- data$beta
    all_phis[[i]] <- data$phi
    all_susie_outputs[[i]] <- results$susie_output
    all_mrash_outputs[[i]] <- results$mrash_output
    #all_susie_ash_outputs[[i]] <- results$susie_ash_output
    all_susie_inf_outputs[[i]] <- results$susie_inf_output
    all_seeds[i] <- seed
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
    all_phis = all_phis,
    all_susie_outputs = all_susie_outputs,
    all_mrash_outputs = all_mrash_outputs,
    #all_susie_ash_outputs = all_susie_ash_outputs,
    all_susie_inf_outputs = all_susie_inf_outputs,
    all_seeds = all_seeds
  )

  file_name <- paste0("numIter", num_simulations,
                      "_largeEffects", n_large_effects,
                      "_mediumEffects", n_medium_effects,
                      "_totalPVE", total_pve,
                      "_L", L)

  saveRDS(simulation_results, file.path(output_dir, file_name))

  # Return all results
  return(simulation_results)
}

# Run the simulation
simulation_results <- simulation(num_simulations = NULL,
                                 n_large_effects = NULL,
                                 n_medium_effects = NULL,
                                 total_pve = NULL,
                                 L = NULL)
