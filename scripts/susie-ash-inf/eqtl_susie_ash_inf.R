# Set Path to Mounted S3 Bucket
data_path <- "/home/apm2217/data/"

# Libraries
library(susieR)
library(mr.ash.alpha)
library(dplyr)
library(magrittr)
source("susie_ash_inf.R")


# Annotation Matrix (from S3 Bucket)
X_full <- readRDS("X20")
#all_seeds <- sample(1:1e9, 100, replace = FALSE)

all_seeds <- c(
  643222872, 642370445, 980619597, 592742195, 812940569, 874493715, 283685745, 272970643, 141777831, 966547239,
  336635874,  87934614, 218490824, 828480440, 160984108, 138238651, 393276358, 261576669, 857535310, 786529725,
  53585008, 670486776, 687961954, 214522377, 564120958,  64094188, 610307235, 289130629, 179824304, 115034695,
  706509872, 250609098,  40608125, 834773543, 895149010, 625671285, 804197336, 410647268, 953414311, 188486097,
  138702256, 980597446, 651248859, 284882193, 949728722, 135086103, 282482480, 655007900, 911828052,  34534519,
  734510236, 395523674, 283749441, 747995237, 738276315, 792296938, 864040833, 214341810, 562898550, 546636977,
  389332802, 889281319, 985937258, 756423001, 834950196, 491958827, 586230129, 118288578, 370853254, 474254047,
  318475291, 498923779, 597127078, 175922978, 267345675, 745538518, 981681389, 254303966, 824205461, 655357804,
  388218441, 882713019, 343395511, 325928029, 498971931, 373901722, 630480032, 982936899, 878924985, 203400282,
  409586743,  22441724, 235542027, 284456998, 376685291, 424044515, 832585303, 634070142, 569987153,  56124815
)

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
method_and_score <- function(X = data$ori.X,
                             y = data$ori.y,
                             beta = data$beta,
                             causal = data$causal,
                             L = 10,
                             v_threshold = v_threshold,
                             precomputed_matrices = precomputed_matrices,
                             seed) {

  set.seed(seed)

  XtX <- precomputed_matrices$XtX
  LD <- precomputed_matrices$LD
  V <- precomputed_matrices$V
  Dsq <- precomputed_matrices$Dsq

  #### Run the susie_ash_inf method ####
  cat("Starting SuSiE-ash-inf\n")
  susie_ash_inf_output <- susie_ash_inf(
    X = X,
    y = y,
    L = L,
    intercept = TRUE,
    standardize = TRUE,
    v_threshold = v_threshold,
    XtX = XtX,
    LD = LD,
    V = V,
    Dsq = Dsq
  )

  #### Calculate Metrics ####
  calc_metrics <- function(mod, X = X, y = y, causal = causal) {
    #### Initialize values ####
    test.cs <- mod$sets_inf  # Use the credible sets from susie_ash_inf
    coverage <- 0
    cs_fdr <- 0
    cs_recall <- 0
    cs_size <- 0

    if (length(test.cs) > 0) {
      # Calculate Average CS Size
      cs_size <- length(unlist(test.cs)) / length(test.cs)

      # Calculate Coverage (proportion of credible sets with a causal effect)
      coverage <- sum(sapply(test.cs, function(cs) any(causal %in% cs))) / length(test.cs)

      # CS Based FDR
      TP_fdr <- sum(sapply(test.cs, function(cs) any(cs %in% causal)))
      FP_fdr <- length(test.cs) - TP_fdr
      cs_fdr <- FP_fdr / (TP_fdr + FP_fdr)

      # CS Based Recall
      TP_recall <- sum(causal %in% unlist(test.cs))
      cs_recall <- TP_recall / length(causal)
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

  # Calculate Metrics for susie_ash_inf
  susie_ash_inf_metrics <- calc_metrics(susie_ash_inf_output, X, y, causal)

  # Create a data frame with the results
  metrics_table <- data.frame(
    Model = "SuSiE-ash-inf",
    RMSE_y = susie_ash_inf_metrics$RMSE_y,
    CS_FDR = susie_ash_inf_metrics$cs_fdr,
    CS_Recall = susie_ash_inf_metrics$cs_recall,
    CS_Size = susie_ash_inf_metrics$cs_size,
    Coverage = susie_ash_inf_metrics$coverage
  )

  # Return the results
  return(list(
    metrics = metrics_table,
    susie_ash_inf_output = susie_ash_inf_output,
    causal = causal,
    betas = beta
  ))
}


# Simulation Function
simulation <- function(num_simulations = NULL,
                       h2_total = NULL,
                       prop_h2_sentinel = NULL,
                       L = NULL,
                       n_oligogenic = NULL,
                       v_threshold = NULL,
                       pve_threshold = NULL) {

  # Parse command-line arguments
  num_simulations = 2
  h2_total = 0.3
  prop_h2_sentinel = 0.7
  L = 10
  n_oligogenic = 20
  v_threshold = 0.005
  pve_threshold = 0.005

  for (arg in commandArgs(trailingOnly = TRUE)) {
    eval(parse(text=arg))
  }

  # Precompute values for susie-inf
  cat("Precomputing matrices for susie_ash_inf\n")
  n <- nrow(X_full)
  p <- ncol(X_full)
  out = susieR:::compute_colstats(X_full, center = T, scale = T)
  attr(X_full, "scaled:center") = out$cm
  attr(X_full, "scaled:scale") = out$csd
  attr(X_full, "d") = out$d
  X_sc = t((t(X_full) - out$cm) / out$csd)  # Centered and scaled X
  X_c = t((t(X_full) - out$cm))             # Centered X

  XtX <- t(X_c) %*% X_c
  LD <- XtX / n
  eig <- eigen(LD, symmetric = TRUE)
  V <- eig$vectors[, ncol(eig$vectors):1]
  Dsq <- pmax(n * sort(eig$values), 0)
  Dsq <- pmax(Dsq, 0)


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
  all_susie_ash_inf_outputs <- list()
  all_epsilons <- numeric(num_simulations)

  for (i in 1:num_simulations) {
    cat("Running simulation", i, "out of", num_simulations, "\n")

    # Set random seed for each simulation
    seed <- all_seeds[i]

    # Generate data
    data <- generate_eqtl_data(X_full,
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
    results <- method_and_score(X = data$ori.X,
                                y = data$ori.y,
                                beta = data$beta,
                                causal = data$causal,
                                L = L,
                                v_threshold = v_threshold,
                                precomputed_matrices = precomputed_matrices,
                                seed = seed)

    # Store results
    all_metrics[[i]] <- results$metrics
    all_betas[[i]] <- results$beta
    all_causal_indices[[i]] <- results$causal
    all_susie_ash_inf_outputs[[i]] <- results$susie_ash_inf_output
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
  output_dir <- "/home/apm2217/output"
  #output_dir <- "analysis"
  simulation_results <- list(
    avg_metrics = avg_metrics,
    all_metrics = all_metrics,
    all_betas = all_betas,
    all_causal_indices = all_causal_indices,
    all_susie_ash_inf_outputs = all_susie_ash_inf_outputs,
    all_seeds = all_seeds,
    all_epsilons = all_epsilons
  )

  file_name <- paste0("numIter", num_simulations,
                      "_h2total", h2_total,
                      "_h2sentinel", prop_h2_sentinel,
                      "_L", L,
                      "_numOligogenic", n_oligogenic,
                      "_vthreshold", v_threshold,
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
                                 pve_threshold = NULL)
