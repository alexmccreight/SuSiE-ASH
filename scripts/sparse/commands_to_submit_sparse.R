# Define the parameter grid
# parameter_grid <- expand.grid(
#   num_simulations = c(200),
#   h2_total = c(0.05, 0.1, 0.2, 0.4),
#   K = c(1,2,3,4,5),
#   L = c(10),
#   stringsAsFactors = FALSE
# )

parameter_grid <- expand.grid(
  num_simulations = c(200),
  h2_total = c(0.3),
  K = c(10),
  L = c(20),
  stringsAsFactors = FALSE
)

# parameter_grid <- expand.grid(
#   num_simulations = c(200),
#   h2_total = c(0.4),
#   K = c(1,2,3,4,5),
#   L = c(10),
#   stringsAsFactors = FALSE
# )

# Create the commands_to_submit.txt file
commands_file <- "commands_to_submit.txt"
file_conn <- file(commands_file, open = "w")

# Iterate over each row of the parameter grid and write commands to the file
for (i in 1:nrow(parameter_grid)) {
  params <- parameter_grid[i, ]

  # Extract parameter values
  num_simulations <- params["num_simulations"]
  h2_total <- params["h2_total"]
  K <- params["K"]
  L <- params["L"]

  # Create the command
  command <- paste0("Rscript /home/apm2217/data/susie-ash-data/sparse_eqtl_simulation.R",
                    " num_simulations=", num_simulations,
                    " h2_total=", h2_total,
                    " K=", K,
                    " L=", L)

  # Write the command to the file
  writeLines(command, file_conn)
}

# Close the file connection
close(file_conn)

cat("Commands file 'commands_to_submit.txt' created successfully.\n")
