# Define the parameter grid
parameter_grid <- expand.grid(
  num_simulations = c(100),
  h2_total = c(0.1, 0.2, 0.3, 0.4, 0.5),
  prop_h2_sentinel = c(0.6, 0.7, 0.8),
  L = c(10),
  n_oligogenic = c(50, 100, 150),
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
  h2_total <- params["h2_total"]
  prop_h2_sentinel <- params["prop_h2_sentinel"]
  L <- params["L"]
  n_oligogenic <- params["n_oligogenic"]

  # Create the command
  command <- paste0("Rscript /home/apm2217/data/eqtl_simulation.R",
                    " num_simulations=", num_simulations,
                    " h2_total=", h2_total,
                    " prop_h2_sentinel=", prop_h2_sentinel,
                    " L=", L,
                    " n_oligogenic=", n_oligogenic)

  # Write the command to the file
  writeLines(command, file_conn)
}

# Close the file connection
close(file_conn)

cat("Commands file 'commands_to_submit.txt' created successfully.\n")
