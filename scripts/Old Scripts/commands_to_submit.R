# Define the parameter grid (SUBSET: 8 settings (w/ 5 iterations))
parameter_grid <- expand.grid(
  num_simulations = c(30),
  total_heritability = c(0.25, 0.33, 0.50),
  sparse_effects = c(3),
  nonsparse_coverage = c(0.005, 0.01, 0.025, 0.05),
  theta_beta_ratio = c(0.5, 0.75, 1.4, 2.25, 3, 5),
  L = c(10),
  stringsAsFactors = FALSE
)

# Create the commands_to_submit.txt file
commands_file <- "commands_to_submit.txt"
file_conn <- file(commands_file, open = "w")

# Iterate over each row of the parameter grid and write commands to the file
for (i in 1:nrow(parameter_grid)) {
  params <- parameter_grid[i, ]

  # Extract parameter values
  num_simulations <- params["num_simulations"]
  total_heritability <- params["total_heritability"]
  sparse_effects <- params["sparse_effects"]
  nonsparse_coverage <- params["nonsparse_coverage"]
  theta_beta_ratio <- params["theta_beta_ratio"]
  L <- params["L"]

  # Create the command
  command <- paste0("Rscript /home/apm2217/data/simulation_prediction_script.R",
                    " num_simulations=",num_simulations,
                    " total_heritability=", total_heritability,
                    " sparse_effects=", sparse_effects,
                    " nonsparse_coverage=", nonsparse_coverage,
                    " theta_beta_ratio=", theta_beta_ratio,
                    " L=", L)

  # Write the command to the file
  writeLines(command, file_conn)
}

# Close the file connection
close(file_conn)

cat("Commands file 'commands_to_submit.txt' created successfully.\n")
