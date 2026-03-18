# Create the combined results folder if it doesn't exist
if (!dir.exists("SimResultsBootAll")) {
  dir.create("SimResultsBootAll")
}

# Define folder names and file names
folder_names <- paste0("SimResultsBoot", seq(1, 91, by = 10))
file_names <- c("table_par.RDS", "table_disp.RDS", "table_parsd.RDS", 
                "table_Ti_list.RDS", "table_TiCI.RDS", "boot_sd.RDS")

# Loop through each file name
for (file in file_names) {
  temp_list <- list()
  
  # Read in from each folder
  for (folder in folder_names) {
    filepath <- file.path(folder, file)
    if (file.exists(filepath)) {
      temp_data <- readRDS(filepath)

      temp_list[[folder]] <- temp_data
    } else {
      warning(paste("File not found:", filepath))
    }
  }
  
  # Combine and save
  combined <- do.call(rbind, temp_list)
  saveRDS(combined, file = file.path("SimResultsBootAll", file))
}
