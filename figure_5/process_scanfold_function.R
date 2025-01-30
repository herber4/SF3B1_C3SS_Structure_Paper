library(data.table)
library(dplyr)

process_profile_files <- function(directory, source_value, output_df_name) {
  # Get a list of files with the specified pattern
  profile_files <- list.files(directory, pattern = ".win_200.stp_20.rnd_100.shfl_mono.out$", full.names = TRUE)
  
  # Create an empty data frame to store the final result
  combined_df <- data.frame()
  
  # Loop through the profile files
  for (profile_file in profile_files) {
    # Get the name of the current file
    file_name <- basename(profile_file)
    
    # Extract the directory name from the file name
    subdir_name <- strsplit(file_name, ".win_200.stp_20.rnd_100.shfl_mono.out", fixed = TRUE)[[1]][1]
    
    # Read the file into a data frame
    current_df <- fread(profile_file, sep = "\t", header = TRUE)
    
    # Add a "gene" column with the directory name
    current_df$gene <- subdir_name
    
    # Rename the 11th column to "Score"
    colnames(current_df)[11] <- "Score"
    
    # Append the current data frame to the combined data frame
    combined_df <- rbind(combined_df, current_df)
  }
  
  # Rename columns 5 and 6
  colnames(combined_df)[5] <- "Z_score"
  colnames(combined_df)[6] <- "P_score"
  
  # Clean up the gene names
  combined_df$gene <- sub("_.*", "", combined_df$gene)
  
  # Filter the combined data frame
  red <- combined_df %>%
    filter(Z_score < 10)
  
  # Add the source column
  red$source <- source_value
  
  # Assign the filtered data frame to the desired output name
  assign(output_df_name, red, envir = .GlobalEnv)
  
  # Clean up temporary variables
  rm(current_df, subdir_name, file_name, profile_file, profile_files, directory, combined_df, red)
}
