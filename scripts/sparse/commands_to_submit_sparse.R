# Define the parameter grid
parameter_grid <- expand.grid(
  num_simulations = c(100),
  h2_total = c(0.3),
  K = c(10),
  L = c(20),
  v_threshold = c(0.0025, 0.005),
  sample_size = c(5000),
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
  K <- params["K"]
  L <- params["L"]
  v_threshold <- params["v_threshold"]
  sample_size <- params["sample_size"]

  # Create the command
  command <- paste0("Rscript /home/apm2217/data/sparse_eqtl_simulation.R",
                    " num_simulations=", num_simulations,
                    " h2_total=", h2_total,
                    " K=", K,
                    " L=", L,
                    " v_threshold=", v_threshold,
                    " sample_size=", sample_size)

  # Write the command to the file
  writeLines(command, file_conn)
}

# Close the file connection
close(file_conn)

cat("Commands file 'commands_to_submit.txt' created successfully.\n")
