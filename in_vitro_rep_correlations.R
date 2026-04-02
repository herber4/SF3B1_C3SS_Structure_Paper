library(Biostrings)
library(ggseqlogo)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(ggbreak)
library(tools)
library(tidyr)

# rep one
{# Set the path to the directory containing the *_profile.txt files
  directory <- "/Users/herber4/Desktop/R_Projects/SF3B1_Pubs/master_scripts/shape_master/rep_two/all_profs/"
  
  # Get a list of files with the *_profile.txt extension
  profile_files <- list.files(directory, pattern = "_profile.txt$", full.names = TRUE)
  
  # Create an empty data frame to store the final result
  combined_df <- data.frame()
  
  # Loop through the profile files
  for (profile_file in profile_files) {
    # Get the name of the current file
    file_name <- basename(profile_file)
    
    # Extract the directory name from the file name
    subdir_name <- strsplit(file_name, "_profile.txt", fixed = TRUE)[[1]][1]
    
    # Read the file into a data frame
    current_df <- fread(profile_file, sep = "\t", header = TRUE)
    current_df[, 3:27] <- NULL
    # Extract the required columns
    # Add a "sample" column with the directory name
    current_df$sample <- subdir_name
    values_to_remove <- c(1:43, 945:967)
    # Filter the rows where "nucleotide" is not in values_to_remove
    current_df <- current_df[!(current_df$Nucleotide %in% values_to_remove), ]
    
    # Append the current data frame to the combined data frame
    combined_df <- rbind(combined_df, current_df)
  }
  rm(current_df, subdir_name, file_name, profile_file, profile_files, directory)
  rep_one <- combined_df}

# rep two
{# Set the path to the directory containing the *_profile.txt files
  directory <- "/Users/herber4/Desktop/R_Projects/SF3B1_Pubs/master_scripts/shape_master/rep_one/rep_one_new_reps"
  
  # Get a list of files with the *_profile.txt extension
  profile_files <- list.files(directory, pattern = "_profile.txt$", full.names = TRUE)
  
  # Create an empty data frame to store the final result
  combined_df <- data.frame()
  
  # Loop through the profile files
  for (profile_file in profile_files) {
    # Get the name of the current file
    file_name <- basename(profile_file)
    
    # Extract the directory name from the file name
    subdir_name <- strsplit(file_name, "_profile.txt", fixed = TRUE)[[1]][1]
    
    # Read the file into a data frame
    current_df <- fread(profile_file, sep = "\t", header = TRUE)
    current_df[, 3:27] <- NULL
    # Extract the required columns
    # Add a "sample" column with the directory name
    current_df$sample <- subdir_name
    values_to_remove <- c(1:43, 945:967)
    # Filter the rows where "nucleotide" is not in values_to_remove
    current_df <- current_df[!(current_df$Nucleotide %in% values_to_remove), ]
    
    # Append the current data frame to the combined data frame
    combined_df <- rbind(combined_df, current_df)
  }
  rm(current_df, subdir_name, file_name, profile_file, profile_files, directory)
  rep_two <- combined_df}

rep_one$adjusted_profile <- ifelse(rep_one$Norm_profile > 4, 4, ifelse(rep_one$Norm_profile < 0, 0, rep_one$Norm_profile))
rep_two$adjusted_profile <- ifelse(rep_two$Norm_profile > 4, 4, ifelse(rep_two$Norm_profile < 0, 0, rep_two$Norm_profile))

# Replace NA with 0 in the adjusted_profile column for both data frames
rep_one <- rep_one %>%
  mutate(adjusted_profile = ifelse(is.na(adjusted_profile), 0, adjusted_profile))

rep_two <- rep_two %>%
  mutate(adjusted_profile = ifelse(is.na(adjusted_profile), 0, adjusted_profile))

# Merge the two data frames by sample
merged_df <- merge(rep_one, rep_two, by = c("sample", "Nucleotide", "Sequence"), suffixes = c("_one", "_two"))

# Calculate the correlation for each sample
spearman <- merged_df %>%
  group_by(sample) %>%
  summarise(correlation = cor(adjusted_profile_one, adjusted_profile_two, method = "spearman"))

rep_one$rep <- "Rep_One"
rep_two$rep <- "Rep_Two"

reps <- rbind(rep_one, rep_two)

library(dplyr)
library(tidyr)

library(dplyr)
library(tidyr)

df_wide <- reps %>%
  select(Nucleotide, adjusted_profile, sample, rep) %>%
  pivot_wider(names_from = rep, values_from = adjusted_profile)


library(dplyr)
library(tidyr)
library(purrr)

results <- reps %>%
  mutate(rep = as.character(rep)) %>%
  group_by(sample, Nucleotide, rep) %>%
  summarise(adjusted_profile = mean(adjusted_profile, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = rep, values_from = adjusted_profile) %>%
  group_by(sample) %>%
  group_modify(~{
    dat <- .x %>% drop_na()
    
    rep_cols <- setdiff(names(dat), "Nucleotide")
    stopifnot(length(rep_cols) == 2)
    
    x <- dat[[rep_cols[1]]]
    y <- dat[[rep_cols[2]]]
    
    ct <- cor.test(x, y, method = "spearman", exact = FALSE)
    rho <- unname(ct$estimate)
    p0 <- ct$p.value
    
    n <- sum(complete.cases(x, y))
    se_z <- 1 / sqrt(n - 3)
    z <- atanh(rho)
    se_rho <- (1 - rho^2) * se_z
    
    z_null <- atanh(0.999999)
    p1 <- 2 * pnorm(-abs((z - z_null) / se_z))
    
    tibble(
      rho = rho,
      p_value_vs_0 = p0,
      se_rho = se_rho,
      p_value_vs_1 = p1,
      n = n,
      rep1 = rep_cols[1],
      rep2 = rep_cols[2]
    )
  }) %>%
  ungroup()

write.table(results, file = "/Users/herber4/Desktop/Dissertation/Dissertation_Aggregation/github_data/Supplementary_Table_3.3.10.txt",
            sep = "\t", row.names = F, quote = F)
