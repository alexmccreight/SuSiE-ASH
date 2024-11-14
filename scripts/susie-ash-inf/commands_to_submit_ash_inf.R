# Define the parameter grid
parameter_grid <- expand.grid(
  num_simulations = c(100),
  h2_total = c(0.3),
  prop_h2_sentinel = c(0.7),
  L = c(10),
  n_oligogenic = c(20),
  v_threshold = c(0.005),
  pve_threshold = c(0.005),
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
  v_threshold <- params["v_threshold"]
  pve_threshold <- params["pve_threshold"]

  # Create the command
  command <- paste0("Rscript /home/apm2217/data/eqtl_susie_ash_inf.R",
                    " num_simulations=", num_simulations,
                    " h2_total=", h2_total,
                    " prop_h2_sentinel=", prop_h2_sentinel,
                    " L=", L,
                    " n_oligogenic=", n_oligogenic,
                    " v_threshold=", v_threshold,
                    " pve_threshold=", pve_threshold)

  # Write the command to the file
  writeLines(command, file_conn)
}

# Close the file connection
close(file_conn)

cat("Commands file 'commands_to_submit.txt' created successfully.\n")
