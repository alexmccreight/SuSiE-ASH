# Set Path to Mounted S3 Bucket
data_path <- "/home/apm2217/data/"

# Libraries
library(susieR)
library(mr.ash.alpha)
library(dplyr)
library(magrittr)
source("susie-ash.R")
source("susie-inf.R")

# Annotation Matrix (from S3 Bucket)
X <- readRDS("X4")

generate_data <- function(X, Ltrue, ssq, sigmasq, tausq){

  # Generate genotype matrix X
  # X <- matrix(rbinom(n * p, size = 2, prob = MAF), nrow = n, ncol = p)

  # Real X matrix input
  X <- scale(X, center = TRUE, scale = TRUE)
  n <- nrow(X)
  p <- ncol(X)

  # Sparse Causal Effects
  beta <- rep(0, p)
  inds <- sample(p, Ltrue, replace = FALSE)
  beta[inds] <- rnorm(Ltrue) * sqrt(ssq)
  order <- order(inds)

  # Generate data with strong infinitesimal effects, tau^2 = 1e-3
  theta <- rnorm(p) * sqrt(tausq)
  effects <- X %*% beta + X %*% theta
  y <- effects + rnorm(n) * sqrt(sigmasq)
  total_variance_explained <- var(effects) / var(y)

  # Store Information
  return(list(X = X, y = y, heritability = total_variance_explained, beta = beta, theta = theta))
}

# Run Methods and Metrics
method_and_score <- function(X = data$X, y = data$y, beta = data$beta, theta = data$theta, L = Ltrue, threshold = 0.90) {
  susie_output <- susie(X = X, y = y, L = L, intercept = F, standardize = F)
  susie_ash_output <- susie_ash(X = X, y = y, L = L, warm_start = 5, tol = 0.001, intercept = F, standardize = F)
  susie_inf_output <- susie_inf(X = X, y = y, L = L, coverage = 0.9, verbose = F)

  calc_metrics <- function(mod, beta = beta, theta = theta, threshold = threshold) {
    all_causal <-  c(which(beta != 0))
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
}
    return(list(cs_fdr = cs_fdr, cs_recall = cs_recall, cs_size = cs_size, coverage = coverage))
  }

  calc_metrics_infinitesimal <- function(mod, beta = beta, theta = theta, threshold = threshold){
    all_causal <- c(which(beta != 0))
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
    return(list(cs_fdr = cs_fdr, cs_recall = cs_recall, cs_size = cs_size, coverage = coverage))
  }

  # Calculate Metrics for each method
  susie_metrics <- calc_metrics(susie_output, beta, theta, threshold)
  susie_ash_metrics <- calc_metrics(susie_ash_output, beta, theta, threshold)
  susie_inf_metrics <- calc_metrics_infinitesimal(susie_inf_output, beta, theta, threshold)

  #Create a data frame with the results
  metrics_table  <- data.frame(
    Model = c("SuSiE","SuSiE-ASH", "SuSiE-Inf"),
    CS_FDR = c(susie_metrics$cs_fdr, susie_ash_metrics$cs_fdr, susie_inf_metrics$cs_fdr),
    CS_Recall = c(susie_metrics$cs_recall, susie_ash_metrics$cs_recall, susie_inf_metrics$cs_recall),
    CS_Size = c(susie_metrics$cs_size, susie_ash_metrics$cs_size, susie_inf_metrics$cs_size),
    Coverage = c(susie_metrics$coverage, susie_ash_metrics$coverage, susie_inf_metrics$coverage)
  )

  # Return the results table
  return(list(
    metrics = metrics_table,
    susie_output = susie_output,
    susie_ash_output = susie_ash_output,
    susie_inf_output = susie_inf_output,
    betas = beta,
    thetas = theta)
  )

}

# Main Simulation Command
simulation <- function(num_simulations = NULL, Ltrue = NULL, ssq = NULL, sigmasq = NULL, tausq = NULL, threshold = NULL) {
  # Parse command-line arguments
  num_simulations = 5
  Ltrue = 5
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
  all_thetas <- list()
  all_heritabilities <- list()
  all_susie_outputs <- list()
  all_susie_ash_outputs <- list()
  all_susie_inf_outputs <- list()
  all_seeds <- numeric(num_simulations)

  for (i in 1:num_simulations) {
    cat("Running simulation", i, "out of", num_simulations, "\n")

    # Set random seed for each simulation
    seed <- abs(round(rnorm(1, mean = 0, sd = 10000)))
    set.seed(seed)

    # Generate data
    data <- generate_data(X = X, Ltrue = Ltrue, ssq = ssq, sigmasq = sigmasq, tausq = tausq)

    # Run methods and calculate metrics
    results <- method_and_score(X = data$X, y = data$y, beta = data$beta, theta = data$theta, L = Ltrue, threshold = threshold)

    # Store results + betas/thetas
    all_metrics[[i]] <- results$metrics
    all_betas[[i]] <- data$beta
    all_thetas[[i]] <- data$theta
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
    Coverage = Reduce("+", lapply(all_metrics, function(x) x$Coverage)) / num_simulations
  )

  # Save simulation results as Rds file
  #output_dir <- "/home/apm2217/output"
  output_dir <- "analysis"
  simulation_results <- list(
    avg_metrics = avg_metrics,
    all_metrics = all_metrics,
    all_betas = all_betas,
    all_thetas = all_thetas,
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
