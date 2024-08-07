# Set Path to Mounted S3 Bucket
data_path <- "/home/apm2217/data/"

# Libraries
library(susieR)
library(mr.ash.alpha)
library(dplyr)
library(magrittr)
source("code/susie_versions/susie-ash.R")
source("code/susie_versions/susie_inf.R")

# Annotation Matrix (from S3 Bucket)
generate_data <- function(Ltrue, ssq, sigmasq, tausq){

  # Generate genotype matrix X
  n = 5000
  p = 500
  MAF = 0.1
  X <- matrix(rbinom(n * p, size = 2, prob = MAF), nrow = n, ncol = p)
  X <- scale(X, center = TRUE, scale = TRUE)

  # Sparse Causal Effects
  beta <- rep(0, p)
  inds <- sample(p, Ltrue, replace = FALSE)
  beta[inds] <- rnorm(Ltrue) * sqrt(ssq)
  order <- order(inds)

  # Generate data with strong infinitesimal effects, tau^2 = 1e-3
  effects <- X %*% beta + X %*% (rnorm(p) * sqrt(tausq))
  y <- effects + rnorm(n) * sqrt(sigmasq)
  total_variance_explained <- var(effects) / var(y)

  # Store Information
  return(list(X = X, y = y, heritability = total_variance_explained, beta = beta))
}

# Run Methods and Metrics
method_and_score <- function(X = data$X, y = data$y, beta = data$beta, Ltrue = Ltrue, threshold = 0.90) {
  cat("Starting susie\n")
  susie_output <- susie(X = X, y = y, L = Ltrue, intercept = F, standardize = F)

  cat("Starting susie-ash\n")
  susie_ash_output <- susie_ash(X = X, y = y, L = Ltrue, warm_start = 2, tol = 0.001, intercept = F, standardize = F, max_iter = 20)

  cat("Starting susie-inf\n")
  susie_inf_output <- susie_inf(X = X, y = y, L = Ltrue, coverage = threshold, verbose = F)

  cat("Starting susie-ash v4\n")
  susie_ash_output_v4 <- susie_ash_v4(X = X, y = y, L = Ltrue, tol = 0.001, intercept = F, standardize = F, max_iter = 20)

  cat("Starting susie-ash v5\n")
  susie_ash_output_v5 <- susie_ash_v5(X = X, y = y, L = Ltrue, tol = 0.001, intercept = F, standardize = F, max_iter = 20)

  #cat("Starting susie-ash v6\n")
  #susie_ash_output <- susie_ash_v6(X = X, y = y, L = Ltrue, warm_start = 2, tol = 0.001, intercept = F, standardize = F, max_iter = 20)

  #cat("Starting susie-ash v7\n")
  #susie_ash_output <- susie_ash_v7(X = X, y = y, L = Ltrue, warm_start = 2, tol = 0.001, intercept = F, standardize = F, max_iter = 20)


  calc_metrics <- function(mod, beta = beta, threshold = threshold) {
    all_causal <-  c(which(beta != 0))
    #all_causal <-  which(abs(beta) >= 0.075)
    test.cs <- susie_get_cs(mod, X = X, coverage = threshold)$cs
    coverage <- 0
    cs_fdr <- 0
    cs_recall <- 0
    cs_size <- 0

    if(length(test.cs) > 0){

    # Calculate Average CS Size
    cs_size <- length(unlist(susie_get_cs(mod, X = X, coverage = threshold)$cs)) / length(susie_get_cs(mod, X = X, coverage = threshold)$cs)

    # Calculate Coverage (proportion of credible sets with a causal effect)
    coverage <- (lapply(1:length(test.cs), function(cs.l){ ifelse(sum(all_causal %in% test.cs[[cs.l]]) != 0, T, F)}) %>% unlist(.) %>% sum(.)) / (length(susie_get_cs(mod, X = X, coverage = threshold)$cs))

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

    # RMSE for outcome
    rmse = sqrt(mean((y - mod$fitted)^2))

}
    return(list(cs_fdr = cs_fdr, cs_recall = cs_recall, cs_size = cs_size, coverage = coverage, rmse = rmse))
  }

  calc_metrics_infinitesimal <- function(mod, beta = beta, threshold = threshold){
    all_causal <- c(which(beta != 0))
    #all_causal <-  which(abs(beta) >= 0.075)
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

    # RMSE for outcome
    rmse = sqrt(mean((y - mod$fitted)^2))
}
    return(list(cs_fdr = cs_fdr, cs_recall = cs_recall, cs_size = cs_size, coverage = coverage, rmse = rmse))
  }

  # Calculate Metrics for each method
  susie_metrics <- calc_metrics(susie_output, beta, threshold)
  susie_ash_metrics <- calc_metrics(susie_ash_output, beta, threshold)
  susie_inf_metrics <- calc_metrics_infinitesimal(susie_inf_output, beta, threshold)

  #Create a data frame with the results
  metrics_table  <- data.frame(
    Model = c("susie",
              "susie-ash",
              "susie-inf"),
    CS_FDR = c(susie_metrics$cs_fdr,
               susie_ash_metrics$cs_fdr,
               susie_inf_metrics$cs_fdr),
    CS_Recall = c(susie_metrics$cs_recall,
                  susie_ash_metrics$cs_recall,
                  susie_inf_metrics$cs_recall),
    CS_Size = c(susie_metrics$cs_size,
                susie_ash_metrics$cs_size,
                susie_inf_metrics$cs_size),
    Coverage = c(susie_metrics$coverage,
                 susie_ash_metrics$coverage,
                 susie_inf_metrics$coverage),
    RMSE = c(susie_metrics$rmse,
             susie_ash_metrics$rmse,
             susie_inf_metrics$rmse)
  )

  # Return the results table
  return(list(
    metrics = metrics_table,
    susie_output = susie_output,
    susie_ash_output = susie_ash_output,
    susie_inf_output = susie_inf_output,
    betas = beta)
  )

}

# Main Simulation Command
simulation <- function(num_simulations = NULL, Ltrue = NULL, ssq = NULL, sigmasq = NULL, tausq = NULL, threshold = NULL) {
  # Parse command-line arguments
  num_simulations = 2
  Ltrue = 10
  ssq = 0.01
  sigmasq = 1
  tausq = .0001
  threshold = 0.90

  for (arg in commandArgs(trailingOnly = TRUE)) {
    eval(parse(text=arg))
  }

  # Initialize lists to store results
  all_metrics <- list()
  all_betas <- list()
  all_heritabilities <- list()
  all_susie_outputs <- list()
  all_susie_ash_outputs <- list()
  all_susie_inf_outputs <- list()
  all_seeds <- numeric(num_simulations)

  for (i in 1:num_simulations) {
    cat("Running simulation", i, "out of", num_simulations, "\n")

    # Set random seed for each simulation
    seed <- as.integer(Sys.time()) + i
    set.seed(seed)

    # Generate data
    data <- generate_data(Ltrue = Ltrue, ssq = ssq, sigmasq = sigmasq, tausq = tausq)

    # Run methods and calculate metrics
    results <- method_and_score(X = data$X, y = data$y, beta = data$beta, L = Ltrue, threshold = threshold)

    # Store results + betas
    all_metrics[[i]] <- results$metrics
    all_betas[[i]] <- data$beta
    all_heritabilities[[i]] <- data$heritability
    all_susie_outputs[[i]] <- results$susie_output
    all_susie_ash_outputs[[i]] <- results$susie_ash_output
    all_susie_inf_outputs[[i]] <- results$susie_inf_output
    all_seeds[i] <- seed
  }

  # Calculate average metrics
  avg_metrics <- data.frame(
    Model = unique(all_metrics[[1]]$Model),
    CS_FDR = Reduce("+", lapply(all_metrics, function(x) x$CS_FDR)) / num_simulations,
    CS_Recall = Reduce("+", lapply(all_metrics, function(x) x$CS_Recall)) / num_simulations,
    CS_Size = Reduce("+", lapply(all_metrics, function(x) x$CS_Size)) / num_simulations,
    Coverage = Reduce("+", lapply(all_metrics, function(x) x$Coverage)) / num_simulations,
    RMSE = Reduce("+", lapply(all_metrics, function(x) x$RMSE)) / num_simulations
  )

  # Save simulation results as Rds file
  #output_dir <- "/home/apm2217/output"
  output_dir <- "analysis"
  simulation_results <- list(
    avg_metrics = avg_metrics,
    all_metrics = all_metrics,
    all_betas = all_betas,
    all_heritabilities = all_heritabilities,
    all_susie_outputs = all_susie_outputs,
    all_susie_ash_outputs = all_susie_ash_outputs,
    all_susie_inf_outputs = all_susie_inf_outputs,
    all_seeds = all_seeds
  )

  file_name <- paste0("numIter", num_simulations,
                      "_Ltrue", Ltrue,
                      "_ssq", ssq,
                      "_sigmasq", sigmasq,
                      "_tausq", tausq)

  saveRDS(simulation_results, file.path(output_dir, file_name))

  # Return all results
  return(simulation_results)
}

# Run the simulation
simulation_results <- simulation(num_simulations = NULL,
                                 Ltrue = NULL,
                                 ssq = NULL,
                                 sigmasq = NULL,
                                 tausq = NULL,
                                 threshold = NULL)
