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

  # Read PLINK files using bigsnpr
  temp_rdsfile <- tempfile(fileext = ".rds")
  snp_readBed(bedfile, backingfile = temp_rdsfile)
  obj.bigSNP <- snp_attach(temp_rdsfile)
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
  file.remove(temp_rdsfile)

  message("Completed processing for ", ld_block_name)
}

# Call the function with the provided LD block name
process_ld_block(ld_block_name)
