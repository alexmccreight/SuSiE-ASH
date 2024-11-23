# process_ld_blocks.R

# Load necessary libraries
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)
ld_block_names <- args

# Specify the input directory and output directory within the script
input_dir <- "/home/apm2217/data"   # Replace with your actual input directory
output_dir <- "/home/apm2217/output" # Replace with your actual output directory

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# Function to mean-impute missing values
mean_impute <- function(geno) {
  # Compute column-wise means excluding NAs
  col_means <- apply(geno, 2, function(x) mean(x, na.rm = TRUE))

  # Replace NAs with the corresponding column mean
  for (i in seq_along(col_means)) {
    na_indices <- which(is.na(geno[, i]))
    if (length(na_indices) > 0) {
      geno[na_indices, i] <- col_means[i]
    }
  }

  return(geno)
}

# Function to process a single LD block
process_ld_block <- function(ld_block_name) {
  # Construct the full path to the LD block file
  input_file <- file.path(input_dir, paste0(ld_block_name, "_matrix.rds"))

  # Check if the input file exists
  if (!file.exists(input_file)) {
    cat("Input file does not exist:", input_file, "\n")
    return(NULL)
  }

  cat("Processing LD block file:", input_file, "\n")

  # Load the LD block from the local file system
  ld_block_obj <- readRDS(input_file)

  # Check if 'genotypes' exist in the LD block object
  if (!"genotypes" %in% names(ld_block_obj)) {
    cat("The LD block object does not contain 'genotypes'. Skipping.\n")
    return(NULL)
  }

  # Extract and mean-impute the genotype matrix
  X <- mean_impute(ld_block_obj$genotypes)

  # Remove columns with zero or near-zero variance
  sd_X <- apply(X, 2, sd)
  threshold <- 1e-8  # Adjust as needed
  columns_to_remove <- which(sd_X < threshold)

  if (length(columns_to_remove) > 0) {
    cat("Removing", length(columns_to_remove), "columns with zero or near-zero variance.\n")
    X <- X[, -columns_to_remove, drop = FALSE]
  }

  if (ncol(X) == 0) {
    cat("No columns left after removing zero or near-zero variance columns. Skipping.\n")
    return(NULL)
  }

  # Scale the genotype matrix
  X_scaled <- scale(X)

  # Compute XtX
  cat("Computing XtX\n")
  XtX <- crossprod(X_scaled)

  # Compute LD matrix
  n_samples <- nrow(X)
  LD <- XtX / n_samples

  # Compute eigenvalues and eigenvectors
  cat("Computing eigenvalues of LD matrix\n")
  eig <- eigen(LD, symmetric = TRUE)
  V <- eig$vectors
  Dsq <- pmax(n_samples * eig$values, 0)

  # Prepare the result list (including X)
  result <- list(
    X = X,        # Include the unscaled processed genotype matrix
    XtX = XtX,
    LD = LD,
    V = V,
    Dsq = Dsq
  )

  # Construct the output file path
  output_file <- file.path(output_dir, paste0(ld_block_name, "_processed.rds"))
  cat("Saving processed matrices to:", output_file, "\n")

  # Save the result to the output directory
  saveRDS(result, file = output_file)
}

# Process each LD block
for (ld_block_name in ld_block_names) {
  process_ld_block(ld_block_name)
}

cat("All LD blocks have been processed.\n")
