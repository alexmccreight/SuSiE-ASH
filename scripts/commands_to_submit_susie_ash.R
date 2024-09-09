# Define the parameter grid
parameter_grid <- expand.grid(
  num_simulations = c(30),
  n_large_effects = c(1,2,3,4,5),
  n_medium_effects = c(1,2,3,4,5),
  total_pve = c(0.25,0.3,0.4,0.5),
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
  n_large_effects <- params["n_large_effects"]
  n_medium_effects <- params["n_medium_effects"]
  total_pve <- params["total_pve"]
  L <- params["L"]

  # Create the command
  command <- paste0("Rscript /home/apm2217/data/susie_ash_simulation_script.R",
                    " num_simulations=",num_simulations,
                    " n_large_effects=", n_large_effects,
                    " n_medium_effects=", n_medium_effects,
                    " total_pve=", total_pve,
                    " L=", L)

  # Write the command to the file
  writeLines(command, file_conn)
}

# Close the file connection
close(file_conn)

cat("Commands file 'commands_to_submit.txt' created successfully.\n")

