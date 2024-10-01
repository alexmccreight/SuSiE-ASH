# Define the parameter grid
parameter_grid <- expand.grid(
  num_simulations = c(30),
  Ltrue = c(3, 5, 7, 10),
  ssq = c(0.01),
  sigmasq = c(1),
  tausq = c(1e-4, 2.5e-4, 3.3e-4, 5e-4, 7.5e-4, 1e-3),
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
  sigmasq <- params["sigmasq"]
  tausq <- params["tausq"]

  # Create the command
  command <- paste0("Rscript /home/apm2217/data/simulation_prediction_script_2.R",
                    " num_simulations=",num_simulations,
                    " Ltrue=", Ltrue,
                    " ssq=", ssq,
                    " sigmasq=", sigmasq,
                    " tausq=", tausq)

  # Write the command to the file
  writeLines(command, file_conn)
}

# Close the file connection
close(file_conn)

cat("Commands file 'commands_to_submit.txt' created successfully.\n")

