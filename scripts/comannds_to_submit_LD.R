# List of LD block file names
ld_block_files <- c(
  "LD_chr1_110", "LD_chr2_139", "LD_chr3_285", "LD_chr4_462",
  "LD_chr5_542", "LD_chr6_658", "LD_chr7_764", "LD_chr8_877",
  "LD_chr9_944", "LD_chr10_1043", "LD_chr11_1137", "LD_chr12_1186",
  "LD_chr13_1311", "LD_chr14_1335", "LD_chr15_1388", "LD_chr16_1483",
  "LD_chr17_1491", "LD_chr18_1536", "LD_chr19_1582", "LD_chr20_1627",
  "LD_chr21_1673", "LD_chr22_1702"
)

# Create the commands_to_submit.txt file
commands_file <- "commands_to_submit.txt"
file_conn <- file(commands_file, open = "w")

# Iterate over each LD block and write commands to the file
for (ld_block_name in ld_block_files) {
  # Create the command
  command <- paste0("Rscript /home/apm2217/data/process_ld_block.R ", ld_block_name)

  # Write the command to the file
  writeLines(command, file_conn)
}

# Close the file connection
close(file_conn)

cat("Commands file 'commands_to_submit.txt' created successfully.\n")
