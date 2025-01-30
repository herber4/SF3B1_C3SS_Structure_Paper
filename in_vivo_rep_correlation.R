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
# Set the path to the directory containing the *_profile.txt files
  directory <- "Figures/in_vivo/rep1/profs"
  
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
    # Filter the rows where "nucleotide" is not in values_to_remove
    
    # Append the current data frame to the combined data frame
    combined_df <- rbind(combined_df, current_df)
  }
  rm(current_df, subdir_name, file_name, profile_file, profile_files, directory)
  rep_one <- combined_df

# rep two
# Set the path to the directory containing the *_profile.txt files
  directory <- "Figures/in_vivo/rep2_pg114/profs"
  
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
    # Filter the rows where "nucleotide" is not in values_to_remove
    
    # Append the current data frame to the combined data frame
    combined_df <- rbind(combined_df, current_df)
  }
  rm(current_df, subdir_name, file_name, profile_file, profile_files, directory)
  rep_two <- combined_df

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

reps

vitro <- read.table(file = "Figures/Figure_3/data/K700E_Adjusted_SHAPE_R1.txt",
                    sep = "\t", header = TRUE)
vitro <- vitro %>%
  filter(sample %in% c("TTC3", "TPP2", "PHGDH",
                       "PPP6R3", "BUB1B"))
bub <- vitro %>%
  filter(Nucleotide >= 402,
         Nucleotide <= 530,
         sample == "BUB1B")
bub$Position <- rownames(bub)
ph <- vitro %>%
  filter(Nucleotide >= 421,
         Nucleotide <= 561,
         sample == "PHGDH")
ph$Position <- rownames(ph)
pp <- vitro %>%
  filter(Nucleotide >= 418,
         Nucleotide <= 557,
         sample == "PPP6R3")
pp$Position <- rownames(pp)
tpp <- vitro %>%
  filter(Nucleotide >= 421,
         Nucleotide <= 566,
         sample == "TPP2")
tpp$Position <- rownames(tpp)
tt <- vitro %>%
  filter(Nucleotide >= 425,
         Nucleotide <= 560,
         sample == "TTC3")
tt$Position <- rownames(tt)

vitro <- rbind(tt, tpp, ph,pp,bub)
vitro$new_index <- NULL
vitro$Nucleotide <- vitro$Position
vitro$Position <- NULL
vitro$Nucleotide <- as.numeric(vitro$Nucleotide)
rep_one$sample <- sub("_.*", "", rep_one$sample)
rep_one$adjusted_profile_in_vivo <- rep_one$adjusted_profile
rep_one$adjusted_profile <- NULL
rep_one <- rep_one[,c(1,2,5,7)]

tmp <- merge(rep_one, vitro, by = c("sample", "Nucleotide", "Sequence"))

tmp[is.na(tmp)] <- 0
spearman <- tmp %>%
  group_by(sample) %>%
  summarise(correlation = cor(adjusted_profile, adjusted_profile_in_vivo, method = "spearman"))

library(ggplot2)

bub <- tmp %>%
  filter(sample == "BUB1B")

# Create the plot
ggplot(bub, aes(x = Nucleotide)) +
  # Plot adjusted_profile on the positive y-axis
  geom_bar(aes(y = adjusted_profile), stat = "identity", fill = "blue", alpha = 0.7) +
  # Plot adjusted_profile_in_vivo on the negative y-axis
  geom_bar(aes(y = -adjusted_profile_in_vivo), stat = "identity", fill = "red", alpha = 0.7) +
  # Add labels and theme
  labs(
    x = "Nucleotide",
    y = "Profile Values",
    title = "Adjusted Profile vs Adjusted Profile In Vivo"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

bub <- tmp %>%
  filter(sample == "BUB1B",
         Nucleotide >= 17, Nucleotide <= 111)

# Create the plot
ggplot(bub, aes(x = Nucleotide)) +
  # Plot adjusted_profile on the positive y-axis
  geom_bar(aes(y = adjusted_profile), stat = "identity", fill = "blue", alpha = 0.7) +
  # Plot adjusted_profile_in_vivo on the negative y-axis
  geom_bar(aes(y = -adjusted_profile_in_vivo), stat = "identity", fill = "red", alpha = 0.7) +
  # Add labels and theme
  labs(
    x = "Nucleotide",
    y = "Profile Values",
    title = "Adjusted Profile vs Adjusted Profile In Vivo"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

tmp <- tmp %>%
  filter(Nucleotide >= 15,
         Nucleotide <= 115)


library(ggplot2)
library(dplyr)
library(ggpubr) # For ggarrange

# Create a function to generate a plot for a given sample
create_plot <- function(sample_name, data) {
  data_filtered <- data %>% filter(sample == sample_name)
  
  ggplot(data_filtered, aes(x = Nucleotide)) +
    geom_bar(aes(y = adjusted_profile), stat = "identity", fill = "blue", alpha = 0.7) +
    geom_bar(aes(y = -adjusted_profile_in_vivo), stat = "identity", fill = "red", alpha = 0.7) +
    labs(
      x = "Nucleotide",
      y = "Profile Values",
      title = paste("Sample:", sample_name)
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# Get a list of unique sample names
unique_samples <- unique(tmp$sample)

# Generate plots for all samples
plots <- lapply(unique_samples, create_plot, data = tmp)

# Combine all plots using ggarrange
combined_plot <- ggarrange(plotlist = plots, ncol = 2, nrow = ceiling(length(plots) / 2))
combined_plot
# Save the combined plot to a file
ggsave("Figures/manuscript_v1.0/combined_invivo_vitro_correlation_plots.png", combined_plot, width = 14, height = 10)

# Display the combined plot
print(combined_plot)
