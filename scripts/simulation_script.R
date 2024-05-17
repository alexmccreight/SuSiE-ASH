# Annotation Matrix (from S3 Bucket)
data_path <- "/home/apm2217/data/X20"
X <- readRDS(data_path)

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
method_and_score <- function(X = data$X, y = data$y, beta = data$beta, theta = data$theta, L = 10, threshold = 0.95) {
  # Run various methods
  susie_output <- susie(X = X, y = y, L = L, intercept = F, standardize = F)
  susie_ash_output <- susie_ash(X = X, y = y, L = L, warm_start = 5, tol = 0.001, intercept = F, standardize = F)

  calc_metrics <- function(mod, beta = beta, theta = theta, threshold = threshold) {
    all_causal <-  c(which(beta != 0), which(theta != 0))

    # Calculate Average CS Size
    cs_size <- length(unlist(susie_get_cs(mod, X = X, coverage = threshold)$cs)) / length(susie_get_cs(mod, X = X, coverage = threshold)$cs)

    test.cs <- susie_get_cs(mod, X = X, coverage = threshold)$cs
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


    return(list(cs_fdr = cs_fdr, cs_recall = cs_recall, cs_size = cs_size, coverage = coverage))
  }

  # Calculate Metrics for each method
  susie_metrics <- calc_metrics(susie_output, beta, theta, threshold)
  susie_ash_metrics <- calc_metrics(susie_ash_output, beta, theta, threshold)

  #Create a data frame with the results
  metrics_table  <- data.frame(
    Model = c("SuSiE","SuSiE-ASH"),
    CS_FDR = c(susie_metrics$cs_fdr, susie_ash_metrics$cs_fdr),
    CS_Recall = c(susie_metrics$cs_recall, susie_ash_metrics$cs_recall),
    CS_Size = c(susie_metrics$cs_size, susie_ash_metrics$cs_size),
    Coverage = c(susie_metrics$coverage, susie_ash_metrics$coverage)
  )
  # Return the results table
  return(list(
    metrics = metrics_table,
    susie_output = susie_output,
    susie_ash_output = susie_ash_output,
    betas = beta,
    thetas = theta)
  )
}

# Simulate
simulation <- function(num_simulations = NULL, total_heritability = NULL, sparse_effects = NULL, nonsparse_coverage = NULL, theta_beta_ratio = NULL, L = NULL, threshold = NULL) {
  # Parse command-line arguments
  args <- commandArgs(trailingOnly = TRUE)
  num_simulations <- as.integer(args[1])
  total_heritability <- as.numeric(args[2])
  sparse_effects <- as.integer(args[3])
  nonsparse_coverage <- as.numeric(args[4])
  theta_beta_ratio <- as.numeric(args[5])
  L <- as.integer(args[6])
  threshold <- as.numeric(args[7])

  # Initialize lists to store results
  all_metrics <- list()
  all_betas <- list()
  all_thetas <- list()
  all_susie_outputs <- list()
  all_susie_ash_outputs <- list()
  all_seeds <- numeric(num_simulations)

  for (i in 1:num_simulations) {
    cat("Running simulation", i, "out of", num_simulations, "\n")

    # Set random seed for each simulation
    seed <- abs(round(rnorm(1, mean = 0, sd = 10000)))
    set.seed(seed)

    # Generate data
    data <- generate_data(X = X, total_heritability = total_heritability, sparse_effects = sparse_effects, nonsparse_coverage = nonsparse_coverage, theta_beta_ratio)

    # Run methods and calculate metrics
    results <- method_and_score(X = X, y = data$y, beta = data$beta, theta = data$theta, L = L, threshold = threshold)

    # Store results + betas/thetas
    all_metrics[[i]] <- results$metrics
    all_betas[[i]] <- data$beta
    all_thetas[[i]] <- data$theta
    all_susie_outputs[[i]] <- results$susie_output
    all_susie_ash_outputs[[i]] <- results$susie_ash_output
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
  output_dir <- "/home/apm2217/output"
  simulation_results <- list(
    avg_metrics = avg_metrics,
    all_metrics = all_metrics,
    all_betas = all_betas,
    all_thetas = all_thetas,
    all_susie_outputs = all_susie_outputs,
    all_susie_ash_outputs = all_susie_ash_outputs,
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
                                 L = NULL,
                                 threshold = NULL)
