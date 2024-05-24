# Define the parameter grid (SUBSET: 2 settings (w/ 2 iterations))
parameter_grid <- expand.grid(
  num_simulations = c(2),
  n = c(5000),
  p = c(500),
  MAF = c(0.1),
  Ltrue = c(5),
  ssq = c(0.01),
  tausq = c(1e-3,1e-4),
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
  n <- params["n"]
  p <- params["p"]
  MAF <- params["MAF"]
  Ltrue <- params["Ltrue"]
  ssq <- params["ssq"]
  tausq <- params["tausq"]
  threshold <- params["threshold"]

  # Create the command
  command <- paste0("Rscript /home/apm2217/data/simple_simulation_script.R",
                    " num_simulations=",num_simulations,
                    " n=", n,
                    " p=", p,
                    " MAF=", MAF,
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
