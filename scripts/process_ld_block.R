# Load necessary libraries
library(bigsnpr)

# Get the LD block name from command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No LD block name provided. Usage: Rscript process_ld_block.R ld_block_name")
}
ld_block_name <- args[1]

# Define your data and output paths
data_path <- "/home/apm2217/data/"
output_path <- "/home/apm2217/output/"

# Function to process a single LD block
process_ld_block <- function(ld_block_name) {
  message("Processing LD block: ", ld_block_name)

  # Define file paths
  bedfile <- file.path(data_path, paste0(ld_block_name, ".bed"))
  bimfile <- file.path(data_path, paste0(ld_block_name, ".bim"))
  famfile <- file.path(data_path, paste0(ld_block_name, ".fam"))

  # Check if files exist
  if (!file.exists(bedfile) || !file.exists(bimfile) || !file.exists(famfile)) {
    stop("One or more PLINK files for ", ld_block_name, " are missing.")
  }

  # Create a temporary directory for processing
  temp_dir <- tempfile(pattern = "temp_ld_block_")
  dir.create(temp_dir)

  # Copy PLINK files to temporary directory
  file.copy(bedfile, temp_dir)
  file.copy(bimfile, temp_dir)
  file.copy(famfile, temp_dir)

  # Define paths to copied files
  temp_bedfile <- file.path(temp_dir, basename(bedfile))
  temp_bimfile <- file.path(temp_dir, basename(bimfile))
  temp_famfile <- file.path(temp_dir, basename(famfile))

  # Read PLINK files using bigsnpr
  backingfile <- file.path(temp_dir, paste0(ld_block_name, "_bk"))

  # Convert PLINK files to bigSNP format (creates .rds and .bk files)
  snp_readBed(temp_bedfile, backingfile = backingfile)

  # Attach the bigSNP object
  obj.bigSNP <- snp_attach(paste0(backingfile, ".rds"))
  G <- obj.bigSNP$genotypes
  sample_ids <- obj.bigSNP$fam$sample.ID
  variant_ids <- obj.bigSNP$map$marker.ID

  # Check the dimensions
  num_individuals <- nrow(G)
  num_variants <- ncol(G)

  message("Number of individuals: ", num_individuals)
  message("Number of variants: ", num_variants)

  # Extract the first 10,000 individuals
  num_individuals_to_sample <- min(10000, num_individuals)
  sampled_individuals <- 1:num_individuals_to_sample

  # Use all variants
  sampled_variants <- 1:num_variants

  # Subset the genotype matrix
  G_sub <- G[sampled_individuals, sampled_variants]

  # Create a list to save
  genotype_data <- list(
    genotypes = G_sub,
    sample_ids = sample_ids[sampled_individuals],
    variant_ids = variant_ids[sampled_variants],
    chromosome = obj.bigSNP$map$chromosome[sampled_variants],
    physical_pos = obj.bigSNP$map$physical.pos[sampled_variants],
    alleles = obj.bigSNP$map[, c("allele1", "allele2")][sampled_variants, ]
  )

  # Save the genotype data to an RDS file in the output directory
  output_file <- file.path(output_path, paste0(ld_block_name, "_matrix.rds"))
  saveRDS(genotype_data, file = output_file)

  message("Saved processed data to ", output_file)

  # Clean up temporary files
  unlink(temp_dir, recursive = TRUE)

  message("Completed processing for ", ld_block_name)
}

# Call the function with the provided LD block name
process_ld_block(ld_block_name)
