# ================================
# Main Simulation Script
# ================================

# Load required libraries
library(susieR)
library(dplyr)
library(magrittr)

# Load Fineboost source files
for (file in list.files("fineboost", pattern = "\\.R$", full.names = TRUE)) {
  source(file)
}

# Source helper functions from the helper_functions folder
source("helper_functions/sparse_data_generation.R")
source("helper_functions/run_methods.R")
source("helper_functions/evaluate_method_performance.R")

##### STILL NEED simxQTL and susieash #####

# Define the simulation function
simulation <- function(num_simulations = NULL,
                       h2_total        = NULL,
                       K               = NULL,
                       L               = NULL,
                       null_max        = NULL,
                       ld_mode         = NULL,
                       K.length        = NULL,
                       upper_bound     = NULL,
                       LD_blocks_dir   = "LD_blocks_precomputed") {

  # Set Default Values if parameters are missing
  if (is.null(num_simulations)) num_simulations <- 2
  if (is.null(h2_total))        h2_total        <- 0.3
  if (is.null(K))               K               <- 5
  if (is.null(L))               L               <- 10
  if (is.null(null_max))        null_max        <- 0.02
  if (is.null(ld_mode))         ld_mode         <- "random"
  if (is.null(K.length))        K.length        <- 20
  if (is.null(upper_bound))     upper_bound     <- 2

  # Parse Command-Line Arguments (if any)
  for (arg in commandArgs(trailingOnly = TRUE)) {
    eval(parse(text = arg))
  }

  # List available LD block files and check count
  ld_block_files <- list.files(path = LD_blocks_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(ld_block_files) < num_simulations) {
    stop("Not enough LD block files.")
  }
  ld_block_files <- ld_block_files[1:num_simulations]

  # Container for simulation metrics
  all_metrics <- vector("list", num_simulations)

  # Loop over simulation replicates
  for (i in seq_len(num_simulations)) {
    cat("Running simulation", i, "of", num_simulations, "\n")

    # Load precomputed LD block data
    ld_block <- readRDS(ld_block_files[i])

    # Store precomputed matrices (these can be large)
    precomp <- list(
      XtX  = ld_block$XtX,
      LD   = ld_block$LD,
      V    = ld_block$V,
      Dsq  = ld_block$Dsq,
      VtXt = ld_block$VtXt
    )

    # Generate simulation data with seed equal to the replicate index
    data <- generate_sparse_eqtl_data(
      X         = ld_block$X,
      K         = K,
      h2        = h2_total,
      ld_mode   = ld_mode,
      ld_matrix = ld_block$LD,
      seed      = i
    )

    # Run methods via the wrapper functions
    susie_out      <- run_susie(data, L, intercept = TRUE, standardize = FALSE)
    susie_ash_out  <- run_susie_ash(data, precomp, L, K.length = K.length, upper_bound = upper_bound,
                                    intercept = TRUE, standardize = FALSE)
    susie_inf_out  <- run_susie_inf(data, precomp, L, intercept = TRUE, standardize = FALSE)
    fineboost_out  <- run_fineboost(data, null_max, intercept = TRUE, standardize = FALSE)

    # Evaluate metrics from all methods using the evaluation function
    metrics <- evaluate_method_performance(
      susie_out     = susie_out,
      susie_ash_out = susie_ash_out,
      susie_inf_out = susie_inf_out,
      fineboost_out = fineboost_out,
      causal        = data$causal_indices,
      data          = data
    )

    # Save the metrics, data, and fine-mapping outputs
    all_metrics[[i]] <- metrics

    # Remove large objects no longer needed to free memory
    rm(ld_block, precomp, susie_out, susie_ash_out, susie_inf_out, fineboost_out, data)
    gc()
  }

  # Compile all results into a list including simulation parameters
  simulation_results <- list(
    all_metrics = all_metrics,
    parameters  = list(
      num_simulations = num_simulations,
      h2_total        = h2_total,
      K               = K,
      L               = L,
      null_max        = null_max,
      ld_mode         = ld_mode,
      K.length        = K.length,
      upper_bound     = upper_bound,
      LD_blocks_dir   = LD_blocks_dir
    )
  )

  # Save Simulation Results as RDS File
  output_dir <- "/home/apm2217/output"
  file_name <- paste0("numIter", num_simulations,
                      "_h2total", h2_total,
                      "_K", K,
                      "_L", L,
                      "_nullMax", null_max,
                      "_ldMode", ld_mode,
                      "_Klength", K.length,
                      "_upperBound", upper_bound)

  saveRDS(simulation_results, file.path(output_dir, paste0(file_name, ".rds")))

  # Return all results
  return(simulation_results)
}

# Run the simulation and save the results
results <- simulation(
  num_simulations = NULL,
  h2_total        = NULL,
  K               = NULL,
  L               = NULL,
  null_max        = NULL,
  ld_mode         = NULL,
  LD_blocks_dir   = "LD_blocks_precomputed"
)
