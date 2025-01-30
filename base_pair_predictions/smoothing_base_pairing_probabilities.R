library(ggplot2)
library(dplyr)
library(purrr)
library(Biostrings)
library(future.apply)
library(data.table)

setwd("/data2/lackey_lab/DownloadedSequenceData/austin/base_pairs/vienna_rnafold/pub/smooth/nalm_can_long/md_fasta/base_pairs")
file_list <- list.files(pattern = "\\.txt$")

# Use purrr's map_dfr to read and combine the files while applying the transformations
combined_probabilities <- map_dfr(file_list, ~{
  # Read the data
  data <- fread(.x, header = FALSE, sep = " ")
  
  # Assign proper column names
  setnames(data, c("V1", "V2", "V3", "V4"), c("i", "j", "sqrt_p", "ubox"))
  
  # Filter for rows where ubox equals "ubox"
  data <- data[ubox == "ubox"]
  
  # Calculate squared probabilities
  data[, p := sqrt_p^2]
  
  # Extract gene name
  gene <- sub("\\.txt$", "", .x)
  data[, gene := gene]
  
  # Calculate the sum of probabilities for `i`
  pair_probabilities_i <- data[, .(p_sum_i = sum(p)), by = .(i, gene)]
  
  # Calculate the sum of probabilities for `j`
  pair_probabilities_j <- data[, .(p_sum_j = sum(p)), by = .(j, gene)]
  
  # Rename columns to avoid conflict and then merge them
  setnames(pair_probabilities_j, c("j", "p_sum_j"), c("i", "p_sum_j"))
  
  # Merge the results from both `i` and `j`
  combined <- merge(pair_probabilities_i, pair_probabilities_j, by = c("i", "gene"), all = TRUE)
  
  # Replace NA values with 0 in the sum columns
  combined[is.na(combined)] <- 0
  
  # Sum the probabilities for both `i` and `j`
  combined[, total_p_sum := p_sum_i + p_sum_j]
  
  return(combined)
})

# Remove "_dp_extracted" from gene names
combined_probabilities[, gene := sub("_dp_extracted", "", gene)]

# Convert to data.table if not already
setDT(combined_probabilities)

# Extract the max_distance using data.table's fast subsetting
combined_probabilities[, max_distance := sub(".*md_(\\d+).*", "\\1", gene)]
# Using data.table to modify the gene column
# Assuming combined_probabilities is your data frame
combined_probabilities$gene <- gsub("_md_\\d+$", "", combined_probabilities$gene)

# Using data.table for efficient group by and summarization
avs <- combined_probabilities[, .(
  av_sum_p = mean(total_p_sum),
  median = median(total_p_sum),
  sd = sd(total_p_sum)
), by = .(gene, i)]
