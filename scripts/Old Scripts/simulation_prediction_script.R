# Set Path to Mounted S3 Bucket
data_path <- "/home/apm2217/data/"

# Libraries
library(susieR)
library(mr.ash.alpha)
library(dplyr)
library(magrittr)
source("susie-ash.R")
source("susie_inf.R")
source("susie_inf_ash.R")
source("susie-ash-v4.R")
source("susie-ash-v5.R")
source("susie-ash-v10.R")
source("susie-ash-v11.R")


# Annotation Matrix (from S3 Bucket)
X <- readRDS("X4")

# Generating Data
generate_data <- function(X, total_heritability, sparse_effects, nonsparse_coverage, theta_beta_ratio) {
  n <- nrow(X)
  p <- ncol(X)

  # Generate sparse effects (beta.true)
  beta.true <- rep(0, p)
  beta.true[sample(p, sparse_effects)] <- rnorm(sparse_effects, mean = 0, sd = .5)

  # Generate nonsparse effects (theta.true)
  num_nonsparse <- round(p * nonsparse_coverage)
  theta.true <- rep(0, p)
  theta.true[sample(p, num_nonsparse)] <- rnorm(num_nonsparse, mean = 0 ,sd = 0.0005)

  # Scale Annotation Matrix
  X <- scale(X, center = TRUE, scale = TRUE)

  # Calculate the variance of the sparse and nonsparse effects
  var_beta <- var(X %*% beta.true)
  var_theta <- var(X %*% theta.true)

  # Adjust the effect sizes based on the theta_beta_ratio
  ratio_factor <- as.numeric((theta_beta_ratio * var_beta) / var_theta)
  theta.true <- theta.true * sqrt(ratio_factor)

  # Recalculate the variance of the adjusted nonsparse effects
  var_theta_adjusted <- var(X %*% theta.true)

  # Calculate the residual variance based on the total heritability
  sigmasq_error <- (var_beta + var_theta_adjusted) * (1 - total_heritability) / total_heritability

  # Create Outcomes
  y <- X %*% matrix(beta.true, ncol = 1) + X %*% matrix(theta.true, ncol = 1) + rnorm(n, 0, sqrt(sigmasq_error))
  y <- scale(y, center = TRUE, scale = FALSE)

  # Store Information
  return(list(X = X, y = y, error = sigmasq_error, beta = beta.true, theta = theta.true))
}

# Method and Scoring
method_and_score <- function(X = data$X, y = data$y, beta = data$beta, theta = data$theta, L = 10) {
  #### Run various methods ####
  cat("Starting SuSiE\n")
  susie_output <- susie(X = X, y = y, L = L, intercept = F, standardize = F)

  cat("Starting mr.ash\n")
  mrash_output <- mr.ash(X = X, y = y, sa2 = nrow(X) * (2^((0:19)/20) - 1)^2, intercept = F, standardize = F)

  cat("Starting SuSiE-ash\n")
  susie_ash_output <- susie_ash(X = X, y = y, L = L, warm_start = 2, tol = 0.001, intercept = F, standardize = F)

  cat("Starting SuSiE-ash (v4)\n")
  susie_ash_output_v4 <- susie_ash_v4(X = X, y = y, L = L, tol = 0.001, intercept = F, standardize = F)

  cat("Starting SuSiE-ash (v5)\n")
  susie_ash_output_v5 <- susie_ash_v5(X = X, y = y, L = L, tol = 0.001, intercept = F, standardize = F)

  cat("Starting SuSiE-ash (v10)\n")
  susie_ash_output_v10 <- susie_ash_v10(X = X, y = y, L = L, tol = 0.001, intercept = F, standardize = F, warm_start = 2)

  cat("Starting SuSiE-ash (v11)\n")
  susie_ash_output_v11 <- susie_ash_v11(X = X, y = y, L = L, tol = 0.001, intercept = F, standardize = F, warm_start = 2)

  cat("Starting SuSiE-inf\n")
  susie_inf_output <- susie_inf(X = X, y = y, L = L, verbose = F, coverage = 0.95)

  cat("Starting SuSiE-inf-ash\n")
  susie_inf_ash_output <- susie_inf_ash(X = X, y = y, L = L, verbose = F, coverage = 0.95)



  calc_metrics_predict_all <- function(mod, X = X, y = y, beta = beta, theta = theta){
    #### Set truth ####
    Xbeta <- X %*% matrix(beta, ncol = 1)
    Xtheta <- X %*% matrix(theta, ncol = 1)
    all_causal <- c(which(beta != 0), which(theta != 0))

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
    RMSE_beta <- sqrt(mean((Xbeta - mod$Xr)^2))
    RMSE_theta <- sqrt(mean((Xtheta - mod$Xtheta)^2))

    #### Store Results ####
    return(list(
      RMSE_y = RMSE_y,
      RMSE_beta = RMSE_beta,
      RMSE_theta = RMSE_theta,
      cs_size = cs_size,
      coverage = coverage,
      cs_fdr = cs_fdr,
      cs_recall = cs_recall
    ))
  }

  calc_metrics_predict_sparse <- function(mod, X = X, y = y, beta = beta, theta = theta){
    #### Set up truth ####
    Xbeta <- X %*% matrix(beta, ncol = 1)
    Xtheta <- X %*% matrix(theta, ncol = 1)
    all_causal <-  c(which(beta != 0), which(theta != 0))

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
    RMSE_beta <- sqrt(mean((Xbeta - mod$Xr)^2))

    #### Store Results ####
    return(list(
      RMSE_y = RMSE_y,
      RMSE_beta = RMSE_beta,
      RMSE_theta = NA,
      cs_size = cs_size,
      coverage = coverage,
      cs_fdr = cs_fdr,
      cs_recall = cs_recall
    ))
  }

  calc_metrics_predict_outcome <- function(mod, X = X, y = y, beta = beta, theta = theta){
    #### Set up truth ####
    Xbeta <- X %*% matrix(beta, ncol = 1)
    Xtheta <- X %*% matrix(theta, ncol = 1)

    #### Calculate RMSE ####
    RMSE_y <- sqrt(mean((y - predict(mod, X))^2))

    #### Store Results ####
    return(list(
      RMSE_y = RMSE_y,
      RMSE_beta = NA,
      RMSE_theta = NA,
      cs_size = NA,
      coverage = NA,
      cs_fdr = NA,
      cs_recall = NA
    ))
  }

  calc_metrics_predict_ash_variants <- function(mod, X = X, y = y, beta = beta, theta = theta){
    #### Set truth ####
    Xbeta <- X %*% matrix(beta, ncol = 1)
    Xtheta <- X %*% matrix(theta, ncol = 1)

    #### Calculate RMSE ####
    RMSE_y <- sqrt(mean((y - mod$fitted)^2))
    RMSE_beta <- sqrt(mean((Xbeta - mod$Xr)^2))
    RMSE_theta <- sqrt(mean((Xtheta - mod$Xtheta)^2))

    #### Store Results ####
    return(list(
      RMSE_y = RMSE_y,
      RMSE_beta = RMSE_beta,
      RMSE_theta = RMSE_theta,
      cs_size = NA,
      coverage = NA,
      cs_fdr = NA,
      cs_recall = NA
    ))
  }

  calc_metrics_inf <- function(mod, X = X, y = y, beta = beta, theta = theta){
    #### Set truth ####
    Xbeta <- X %*% matrix(beta, ncol = 1)
    Xtheta <- X %*% matrix(theta, ncol = 1)
    all_causal <- c(which(beta != 0), which(theta != 0))

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
    RMSE_beta <- sqrt(mean((Xbeta - X %*% rowSums(mod$PIP2 * mod$mu))^2))
    RMSE_theta <- sqrt(mean((Xtheta - X %*% susie_inf_output$alpha)^2))

    #### Store Results ####
    return(list(
      RMSE_y = RMSE_y,
      RMSE_beta = RMSE_beta,
      RMSE_theta = RMSE_theta,
      cs_size = cs_size,
      coverage = coverage,
      cs_fdr = cs_fdr,
      cs_recall = cs_recall
    ))
  }

  #############
  # Calculate Metrics for each method
  susie_metrics <- calc_metrics_predict_sparse(susie_output, X, y, beta, theta)
  mrash_metrics <- calc_metrics_predict_outcome(mrash_output, X, y, beta, theta)
  susie_ash_metrics <- calc_metrics_predict_all(susie_ash_output, X, y, beta, theta)
  susie_ash_v4_metrics <- calc_metrics_predict_ash_variants(susie_ash_output_v4, X, y, beta, theta)
  susie_ash_v5_metrics <- calc_metrics_predict_ash_variants(susie_ash_output_v5, X, y, beta, theta)
  susie_ash_v10_metrics <- calc_metrics_predict_all(susie_ash_output_v10, X, y, beta, theta)
  susie_ash_v11_metrics <- calc_metrics_predict_all(susie_ash_output_v11, X, y, beta, theta)
  susie_inf_metrics <- calc_metrics_inf(susie_inf_output, X, y, beta, theta)
  susie_inf_ash_metrics <- calc_metrics_inf(susie_inf_ash_output, X, y, beta, theta)

  #Create a data frame with the results
  metrics_table  <- data.frame(
    Model = c("SuSiE",
              "mr.ash",
              "SuSiE-ash",
              "SuSiE-ash (v4)",
              "SuSiE-ash (v5)",
              "SuSiE-ash (v10)",
              "SuSiE-ash (v11)",
              "SuSiE-inf",
              "SuSiE-inf-ash"),
    RMSE_y = c(susie_metrics$RMSE_y,
               mrash_metrics$RMSE_y,
               susie_ash_metrics$RMSE_y,
               susie_ash_v4_metrics$RMSE_y,
               susie_ash_v5_metrics$RMSE_y,
               susie_ash_v10_metrics$RMSE_y,
               susie_ash_v11_metrics$RMSE_y,
               susie_inf_metrics$RMSE_y,
               susie_inf_ash_metrics$RMSE_y),
    RMSE_beta = c(susie_metrics$RMSE_beta,
                  mrash_metrics$RMSE_beta,
                  susie_ash_metrics$RMSE_beta,
                  susie_ash_v4_metrics$RMSE_beta,
                  susie_ash_v5_metrics$RMSE_beta,
                  susie_ash_v10_metrics$RMSE_beta,
                  susie_ash_v11_metrics$RMSE_beta,
                  susie_inf_metrics$RMSE_beta,
                  susie_inf_ash_metrics$RMSE_beta),
    RMSE_theta = c(susie_metrics$RMSE_theta,
                   mrash_metrics$RMSE_theta,
                   susie_ash_metrics$RMSE_theta,
                   susie_ash_v4_metrics$RMSE_theta,
                   susie_ash_v5_metrics$RMSE_theta,
                   susie_ash_v10_metrics$RMSE_theta,
                   susie_ash_v11_metrics$RMSE_theta,
                   susie_inf_metrics$RMSE_theta,
                   susie_inf_ash_metrics$RMSE_theta),
    CS_FDR = c(susie_metrics$cs_fdr,
               mrash_metrics$cs_fdr,
               susie_ash_metrics$cs_fdr,
               susie_ash_v4_metrics$cs_fdr,
               susie_ash_v5_metrics$cs_fdr,
               susie_ash_v10_metrics$cs_fdr,
               susie_ash_v11_metrics$cs_fdr,
               susie_inf_metrics$cs_fdr,
               susie_inf_ash_metrics$cs_fdr),
    CS_Recall = c(susie_metrics$cs_recall,
                  mrash_metrics$cs_recall,
                  susie_ash_metrics$cs_recall,
                  susie_ash_v4_metrics$cs_recall,
                  susie_ash_v5_metrics$cs_recall,
                  susie_ash_v10_metrics$cs_recall,
                  susie_ash_v11_metrics$cs_recall,
                  susie_inf_metrics$cs_recall,
                  susie_inf_ash_metrics$cs_recall),
    CS_Size = c(susie_metrics$cs_size,
                mrash_metrics$cs_size,
                susie_ash_metrics$cs_size,
                susie_ash_v4_metrics$cs_size,
                susie_ash_v5_metrics$cs_size,
                susie_ash_v10_metrics$cs_size,
                susie_ash_v11_metrics$cs_size,
                susie_inf_metrics$cs_size,
                susie_inf_ash_metrics$cs_size),
    Coverage = c(susie_metrics$coverage,
                 mrash_metrics$coverage,
                 susie_ash_metrics$coverage,
                 susie_ash_v4_metrics$coverage,
                 susie_ash_v5_metrics$coverage,
                 susie_ash_v10_metrics$coverage,
                 susie_ash_v11_metrics$coverage,
                 susie_inf_metrics$coverage,
                 susie_inf_ash_metrics$coverage)
  )
  # Return the results table
  return(list(
    metrics = metrics_table,
    susie_output = susie_output,
    mrash_output = mrash_output,
    susie_ash_output = susie_ash_output,
    susie_ash_output_v4 = susie_ash_output_v4,
    susie_ash_output_v5 = susie_ash_output_v5,
    susie_ash_output_v10 = susie_ash_output_v10,
    susie_ash_output_v11 = susie_ash_output_v11,
    susie_inf_output = susie_inf_output,
    susie_inf_ash_output = susie_inf_ash_output,
    betas = beta,
    thetas = theta)
  )
}

# Simulate
simulation <- function(num_simulations = NULL, total_heritability = NULL, sparse_effects = NULL, nonsparse_coverage = NULL, theta_beta_ratio = NULL, L = NULL) {
  # Parse command-line arguments
  num_simulations = 2
  total_heritability = 0.5
  sparse_effects = 2
  nonsparse_coverage = 0.05
  theta_beta_ratio = 0.75
  L = 10

  for (arg in commandArgs(trailingOnly = TRUE)) {
    eval(parse(text=arg))
  }

  # Initialize lists to store results
  all_metrics <- list()
  all_betas <- list()
  all_thetas <- list()
  all_susie_outputs <- list()
  all_mrash_outputs <- list()
  all_susie_ash_outputs <- list()
  all_susie_ash_outputs_v4 <- list()
  all_susie_ash_outputs_v5 <- list()
  all_susie_ash_outputs_v10 <- list()
  all_susie_ash_outputs_v11 <- list()
  all_susie_inf_outputs <- list()
  all_susie_inf_ash_outputs <- list()
  all_seeds <- numeric(num_simulations)

  for (i in 1:num_simulations) {
    cat("Running simulation", i, "out of", num_simulations, "\n")

    # Set random seed for each simulation
    seed <- abs(round(rnorm(1, mean = 0, sd = 10000)))
    set.seed(seed)

    # Generate data
    data <- generate_data(X = X, total_heritability = total_heritability, sparse_effects = sparse_effects, nonsparse_coverage = nonsparse_coverage, theta_beta_ratio)

    # Run methods and calculate metrics
    results <- method_and_score(X = data$X, y = data$y, beta = data$beta, theta = data$theta, L = L)

    # Store results + betas/thetas
    all_metrics[[i]] <- results$metrics
    all_betas[[i]] <- data$beta
    all_thetas[[i]] <- data$theta
    all_susie_outputs[[i]] <- results$susie_output
    all_mrash_outputs[[i]] <- results$mrash_output
    all_susie_ash_outputs[[i]] <- results$susie_ash_output
    all_susie_ash_outputs_v4[[i]] <- results$susie_ash_output_v4
    all_susie_ash_outputs_v5[[i]] <- results$susie_ash_output_v5
    all_susie_ash_outputs_v10[[i]] <- results$susie_ash_output_v10
    all_susie_ash_outputs_v11[[i]] <- results$susie_ash_output_v11
    all_susie_inf_outputs[[i]] <- results$susie_inf_output
    all_susie_inf_ash_outputs[[i]] <- results$susie_inf_ash_output
    all_seeds[i] <- seed
  }

  # Calculate average metrics
  avg_metrics <- data.frame(
    Model = unique(all_metrics[[1]]$Model),
    RMSE_y = Reduce("+", lapply(all_metrics, function(x) x$RMSE_y)) / num_simulations,
    RMSE_beta = Reduce("+", lapply(all_metrics, function(x) x$RMSE_beta)) / num_simulations,
    RMSE_theta = Reduce("+", lapply(all_metrics, function(x) x$RMSE_theta)) / num_simulations,
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
    all_thetas = all_thetas,
    all_susie_outputs = all_susie_outputs,
    all_mrash_outputs = all_mrash_outputs,
    all_susie_ash_outputs = all_susie_ash_outputs,
    all_susie_ash_outputs_v4 = all_susie_ash_outputs_v4,
    all_susie_ash_outputs_v5 = all_susie_ash_outputs_v5,
    all_susie_ash_outputs_v10 = all_susie_ash_outputs_v10,
    all_susie_ash_outputs_v11 = all_susie_ash_outputs_v11,
    all_susie_inf_outputs = all_susie_inf_outputs,
    all_susie_inf_ash_outputs = all_susie_inf_ash_outputs,
    all_seeds = all_seeds
  )

  file_name <- paste0("numIter", num_simulations,
                      "_totHeritability", total_heritability,
                      "_sparseEffect", sparse_effects,
                      "_nonsparse", nonsparse_coverage,
                      "_ratio", theta_beta_ratio,
                      "_L", L)

  saveRDS(simulation_results, file.path(output_dir, file_name))

  # Return all results
  return(simulation_results)
}

# Run the simulation
simulation_results <- simulation(num_simulations = NULL,
                                 total_heritability = NULL,
                                 sparse_effects = NULL,
                                 nonsparse_coverage = NULL,
                                 theta_beta_ratio = NULL,
                                 L = NULL)
