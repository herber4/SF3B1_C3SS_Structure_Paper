# Load necessary libraries
library(tidyverse)
library(stringr)
set.seed(5443)
# Step 1: Read the file line by line
{file_path <- "data/Nalm6_C3SS_Max_ENT.txt"  # Change this to your file's path
  lines <- readLines(file_path)
  
  # Step 2: Extract relevant information
  genes <- lines[grepl("^>", lines)] %>% str_replace_all(">", "")  # Extract gene names and remove ">"
  sequences <- lines[!grepl("^>", lines)] %>% str_extract("[A-Z]+")  # Extract sequences (uppercase letters)
  
  # Updated regex to include negative values
  maxent <- lines[!grepl("^>", lines)] %>% str_extract("(?<=MAXENT: )-?\\d+\\.\\d+") %>% as.numeric()  # Extract MAXENT values
  mm <- lines[!grepl("^>", lines)] %>% str_extract("(?<=MM: )-?\\d+\\.\\d+") %>% as.numeric()  # Extract MM values
  wmm <- lines[!grepl("^>", lines)] %>% str_extract("(?<=WMM: )-?\\d+\\.\\d+") %>% as.numeric()  # Extract WMM values
  
  # Step 3: Combine into a data frame
  c3ss <- data.frame(
    gene = genes,
    seq = sequences,
    MAXENT = maxent,
    MM = mm,
    WMM = wmm
  )}

# View the resulting data frame
{library(tidyverse)
  
  # Step 1: Read the file line by line
  file_path <- "data/Nalm6_Canonical_Max_ENT.txt"  # Change this to your file's path
  lines <- readLines(file_path)
  
  # Step 2: Extract relevant information
  genes <- lines[grepl("^>", lines)] %>% str_replace_all(">", "")  # Extract gene names and remove ">"
  sequences <- lines[!grepl("^>", lines)] %>% str_extract("[A-Z]+")  # Extract sequences (uppercase letters)
  
  # Updated regex to include negative values
  maxent <- lines[!grepl("^>", lines)] %>% str_extract("(?<=MAXENT: )-?\\d+\\.\\d+") %>% as.numeric()  # Extract MAXENT values
  mm <- lines[!grepl("^>", lines)] %>% str_extract("(?<=MM: )-?\\d+\\.\\d+") %>% as.numeric()  # Extract MM values
  wmm <- lines[!grepl("^>", lines)] %>% str_extract("(?<=WMM: )-?\\d+\\.\\d+") %>% as.numeric()  # Extract WMM values
  
  # Step 3: Combine into a data frame
  can <- data.frame(
    gene = genes,
    seq = sequences,
    MAXENT = maxent,
    MM = mm,
    WMM = wmm
  )}

can$source <- "Canonical"
c3ss$source <- "C3SS"
rm(mm, wmm, sequences, lines, genes, file_path, maxent)

nalm_mas <- rbind(c3ss, can)

ggplot(nalm_mas, aes(x = source, y = WMM, fill = source)) +
  geom_boxplot()


# Load necessary library
library(dplyr)

# Function to perform one bootstrap sampling and calculate the means
bootstrap_sample <- function(df) {
  sample_df <- sample_n(df, 192, replace = TRUE)
  means <- colMeans(sample_df[, c("MAXENT", "MM", "WMM")], na.rm = TRUE)
  return(means)
}

# Perform the bootstrap sampling 1000 times and store the results
c3ss_bs <- replicate(1000, bootstrap_sample(c3ss), simplify = "matrix")

# Convert the results into a data frame
c3ss_bs <- as.data.frame(t(c3ss_bs))
c3ss_bs$source <- "C3SS"

can_bs <- replicate(1000, bootstrap_sample(can), simplify = "matrix")

# Convert the results into a data frame
can_bs <- as.data.frame(t(can_bs))
can_bs$source <- "Canonical"

bs <- rbind(can_bs, c3ss_bs)
# Check the first few rows of the resulting bootstrap means

ggplot(bs, aes(x = source, y = MAXENT, fill = source)) +
  geom_boxplot()

mas <- read.table(file = "/rMATS_MAXENT_C3SS_Can_Combined.txt",
                  sep = "\t", header = TRUE)
mas$source <- ifelse(mas$source == "C3SS", "K700E_C3SS", "K700E_Canonical")
tmp <- mas
tmp$gene <- NULL
tmp$seq <- NULL

tmp <- rbind(bs, tmp)

ggplot(tmp, aes(x = source, y = WMM, fill = source)) +
  geom_boxplot() 


k7_c3ss <- mas %>%
  filter(source == "K700E_C3SS")
k7_wm <- mean(k7_c3ss$MAXENT)
wt <- mas %>%
  filter(source == "K700E_Canonical")
wt_wm <- mean(wt$MAXENT)

mean(c3ss_bs$MAXENT >= k7_wm)
mean(c3ss_bs$MAXENT <= k7_wm)
mean(can_bs$MAXENT <= wt_wm)
mean(can_bs$MAXENT >= wt_wm)

stats <- mas %>%
  group_by(source) %>%
  summarise(MAXENT = mean(MAXENT),
            WMM = mean(WMM),
            MM = mean(MM))

pdf(file = "/C3SS_MAXENT.pdf",
    width = 8, height = 10, paper = "letter")
ggplot(bs, aes(x = source, y = MAXENT)) +
  geom_boxplot() +
  geom_point(data = stats, aes(x = source, y = MAXENT),
             col = "black", size = 4) +
  theme_bw() +
  labs(x = "Source", y = "MEXENT") +
  theme(legend.position = "")
dev.off()
mean(can_bs$MAXENT)
mean(wt$MAXENT)

mas <- rbind(bs, stats)

level <- c("K700E_C3SS", "C3SS", "K700E_Canonical", "Canonical")

mas %>%
  filter(!source == "K700E_C3SS",
         !source == "K700E_Canonical") %>%
  ggplot(aes(x = source, y = MAXENT, fill = source, group = source)) +
  geom_boxplot() +
  geom_point(data = mas %>% filter(source == "K700E_C3SS" |
                                   source == "K700E_Canonical"), aes(x = source, y = MAXENT),
             size = 3) +
  scale_x_discrete(limits = level) +
  theme_bw() +
  labs(x = "Source", y = "MAXENT")


pdf(file = "/C3SS_MAXENT.pdf",
    width = 8, height = 10, paper = "letter")
mas %>%
  filter(!source == "K700E_C3SS",
         !source == "K700E_Canonical") %>%
  ggplot(aes(x = source, y = WMM, group = source)) +
  geom_boxplot() +
  geom_point(data = mas %>% filter(source == "K700E_C3SS" |
                                     source == "K700E_Canonical"), aes(x = source, y = WMM),
             size = 3) +
  scale_x_discrete(limits = level) +
  theme_bw() +
  labs(x = "Source", y = "MAXENT") +
  theme(legend.position = "none")
dev.off()

mean(can_bs$MAXENT >= 5.756508)
mean(can_bs$MAXENT <= 5.756508)
