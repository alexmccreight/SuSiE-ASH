# Define the parameter grid (SUBSET: 4 settings (w/ 25 iterations))
parameter_grid <- expand.grid(
  num_simulations = c(25),
  Ltrue = c(5, 10),
  ssq = c(0.01, 0.05, 0.1, 0.125),
  tausq = c(0.0001, 0.00025, 0.00033, 0.0005),
  threshold = c(0.90),
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
  Ltrue <- params["Ltrue"]
  ssq <- params["ssq"]
  tausq <- params["tausq"]
  threshold <- params["threshold"]

  # Create the command
  command <- paste0("Rscript /home/apm2217/data/simple_simulation_script.R",
                    " num_simulations=",num_simulations,
                    " Ltrue=", Ltrue,
                    " ssq=", ssq,
                    " tausq=", tausq,
                    " threshold=", threshold)

  # Write the command to the file
  writeLines(command, file_conn)
}

# Close the file connection
close(file_conn)

cat("Commands file 'commands_to_submit.txt' created successfully.\n")
