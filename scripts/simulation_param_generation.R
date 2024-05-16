# Define the parameter grid (FULL: 90 settings (w/ 50 iterations))
# parameter_grid <- expand.grid(
#   num_simulations = c(50),
#   total_heritability = c(0.25, 0.50),
#   sparse_effects = c(2),
#   nonsparse_coverage = c(0.01, 0.025, 0.05),
#   theta_beta_ratio = c(1.4, 3, 5),
#   L = c(10, 15, 20, 25, 30),
#   threshold = c(0.90),
#   stringsAsFactors = FALSE
# )

# Define the parameter grid (SUBSET: 8 settings (w/ 5 iterations))
parameter_grid <- expand.grid(
  num_simulations = c(5),
  total_heritability = c(0.25, 0.50),
  sparse_effects = c(2),
  nonsparse_coverage = c(0.01, 0.05),
  theta_beta_ratio = c(1.4, 3),
  L = c(10),
  threshold = c(0.90),
  stringsAsFactors = FALSE
)

# Create a text file to store the parameter combinations
file_name <- "simulation_params.txt"
file_conn <- file(file_name, open = "w")

# Iterate over each row of the parameter grid and write to the file
for (i in 1:nrow(parameter_grid)) {
  params <- parameter_grid[i, ]

  # Extract parameter values
  num_simulations <- params["num_simulations"]
  total_heritability <- params["total_heritability"]
  sparse_effects <- params["sparse_effects"]
  nonsparse_coverage <- params["nonsparse_coverage"]
  theta_beta_ratio <- params["theta_beta_ratio"]
  L <- params["L"]
  threshold <- params["threshold"]

  # Write the parameter values as a single line in the file
  param_line <- paste(num_simulations, total_heritability, sparse_effects, nonsparse_coverage,
                      theta_beta_ratio, L, threshold, sep = " ")
  writeLines(param_line, file_conn)
}

# Close the file connection
close(file_conn)

cat("Parameter file 'simulation_params.txt' created successfully.\n")
