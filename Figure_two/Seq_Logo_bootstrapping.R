library(Biostrings)
library(GenomicRanges)
library(ggplot2)
library(ggseqlogo)
library(ggpubr)
library(dplyr)
library(tidyr)

set.seed(5443)
# this is for cryptic sites
nalm <- readDNAStringSet(filepath = "/data/Nalm6_Regtools_C3SS_Junctions.fasta",
                          format = "fasta")
  library(Biostrings)
  
  # Define the positions you're interested in
  pos_A <- 449
  pos_G <- 450
  
  # Filter sequences that have 'A' at position 449 and 'G' at position 450
  filtered_nalm <- nalm[substring(as.character(nalm), pos_A, pos_A) == "A" & 
                          substring(as.character(nalm), pos_G, pos_G) == "G"]
  nalm <- filtered_nalm
  rm(filtered_nalm)
  
{k7 <- readDNAStringSet(filepath = "/figure_two/K700E_C3SS_MAXENT.fasta",
                          format = "fasta")
    
    k7$COMT <- NULL
    k7$BRD9 <- NULL
    k7$UBA1 <- NULL
    k7$RPRD1A <- NULL
    k7$EHMT1 <- NULL}
  
rmats_can <- readDNAStringSet(file = "/figure_two/rMATS_K700E_paired3SS_MAXENT_submission.fasta",
                                format = "fasta")
{tmp <- subseq(rmats_can, start = 400, end = 470)
  tmp <- consensusMatrix(tmp, as.prob = TRUE)
  ggplot() +
    geom_logo(data = tmp, method = "probability", col_scheme = "nucleotide") +
    theme_classic() +
    theme(text = element_text(size = 18),
          axis.text.x = element_blank()) +
    labs(x = "V")}
rmats_c3ss <- readDNAStringSet(filepath = "/figure_two/rMATS_K700E_shared_C3SS900.fasta",
                               format = "fasta")
tmp <- subseq(rmats_c3ss, start = 400, end = 470)
tmp <- consensusMatrix(tmp, as.prob = TRUE)
ggplot() +
  geom_logo(data = tmp, method = "probability", col_scheme = "nucleotide") +
  theme_classic() +
  theme(text = element_text(size = 18),
        axis.text.x = element_blank()) +
  labs(x = "V")
# Initialize an empty list to store consensus matrices
consensus_matrices <- list()

# Define the number of bootstrap iterations
num_bootstraps <- 1000

# Perform bootstrapping
for (i in 1:num_bootstraps) {
  # Sample from vdbc3ss with replacement
  tmp <- sample(nalm, size = 192, replace = TRUE)
  
  # Generate consensus matrix
  pfm <- consensusMatrix(tmp, as.prob = TRUE)
  
  # Convert consensus matrix to dataframe
  pfm_df <- as.data.frame(pfm)
  
  # Add column for nucleotide positions
  pfm_df$nucleotide <- rownames(pfm_df)
  
  # Remove row names
  rownames(pfm_df) <- NULL
  
  # Add bootstrap number
  pfm_df$Bootstrap_Number <- i
  
  # Store the dataframe
  consensus_matrices[[i]] <- pfm_df
}

# Combine all consensus matrices into a single dataframe
combined_df <- bind_rows(consensus_matrices)

# Filter rows where Nucleotide is A, C, T, or G
combined_df <- combined_df %>%
  filter(nucleotide %in% c("A", "C", "T", "G"))

tmp <- pivot_longer(combined_df, -c(nucleotide, Bootstrap_Number),
                    values_to = "probability", names_to = "position")
tmp$position <- sub("V", "", tmp$position)
tmp$position <- as.numeric(tmp$position)

tmp %>%
  filter(
    position > 425, 
    position < 475) %>%
  ggplot(aes(x = as.factor(position), y = probability,
             fill = nucleotide)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90))

pfm1 <- consensusMatrix(rmats_c3ss, as.prob = TRUE)
pfm1 <- as.data.frame(pfm1)
pfm1$nucleotide <- rownames(pfm1)
rownames(pfm1) <- NULL
pfm1 <- pfm1 %>%
  filter(nucleotide %in% c("A", "T", "C", "G"))
tmp2 <- pivot_longer(pfm1, -c(nucleotide), values_to = "probability",
                     names_to = "position")
tmp2$position <- sub("V", "", tmp2$position)
#tmp2$new_position <- as.numeric(tmp2$position)-1
k700e_long <- tmp2

tmp %>%
  filter(nucleotide == "T",
         position > 440,
         position < 460) %>%
  ggplot(aes(x = as.factor(position), y = probability, group = position)) +
  geom_boxplot() +
  geom_point(data = tmp2 %>% filter(nucleotide == "T",
                                    position > 440,
                                    position < 460),
             aes(x = as.factor(position), y = probability), col = "red")
tmp %>%
  filter(nucleotide == "T",
         position > 300,
         position < 600) %>%
  ggplot(aes(x = as.factor(position), y = probability, group = position)) +
  geom_boxplot() +
  geom_point(data = tmp2 %>% filter(nucleotide == "T",
                                    position > 300,
                                    position < 600),
             aes(x = as.factor(position), y = probability), col = "red")
tmp2$source <- "K700E"
tmp2$position <- as.numeric(tmp2$position)
#tmp2$position <- tmp2$position-1
tmp$source <- "Nalm6_WT"
tmp$Bootstrap_Number <- NULL
test <- rbind(tmp, tmp2)
test$position <- as.numeric(test$position)

pdf(file = "/figure_two/C3SS_T_content.pdf",
    width = 9, height = 4, paper = "letter")
test %>%
  filter(source == "Nalm6_WT",
         nucleotide == "T",
         position > 429,
         position < 464) %>%
  ggplot(aes(x = as.factor(position), y = probability)) +
  geom_boxplot() +
  geom_point(data = test %>% filter(source == "K700E",
                                    nucleotide == "T",
                                    position > 429,
                                    position < 464),
             aes(x = as.factor(position), y = probability),
             col = "red", size = 3) +
  theme_bw() +
  labs(x = "Position",
       y = "Probability (T)") +
  ylim(0,0.6)
dev.off()

pdf(file = "/figure_two/C3SS_A_content.pdf",
    width = 9, height = 4, paper = "letter")
test %>%
  filter(source == "Nalm6_WT",
         nucleotide == "A",
         position > 429,
         position < 464) %>%
  ggplot(aes(x = as.factor(position), y = probability)) +
  geom_boxplot() +
  geom_point(data = test %>% filter(source == "K700E",
                                    nucleotide == "A",
                                    position > 429,
                                    position < 464),
             aes(x = as.factor(position), y = probability),
             col = "red", size = 3) +
  theme_bw() +
  labs(x = "Position",
       y = "Probability (A)")
dev.off()

pdf(file = "/figure_two/C3SS_G_content.pdf",
    width = 9, height = 4, paper = "letter")
test %>%
  filter(source == "Nalm6_WT",
         nucleotide == "G",
         position > 429,
         position < 464) %>%
  ggplot(aes(x = as.factor(position), y = probability)) +
  geom_boxplot() +
  geom_point(data = test %>% filter(source == "K700E",
                                    nucleotide == "G",
                                    position > 429,
                                    position < 464),
             aes(x = as.factor(position), y = probability),
             col = "red", size = 3) +
  theme_bw() +
  labs(x = "Position",
       y = "Probability (G)")
dev.off()

pdf(file = "/figure_two/C3SS_C_content.pdf",
    width = 9, height = 4, paper = "letter")
test %>%
  filter(source == "Nalm6_WT",
         nucleotide == "C",
         position > 429,
         position < 464) %>%
  ggplot(aes(x = as.factor(position), y = probability)) +
  geom_boxplot() +
  geom_point(data = test %>% filter(source == "K700E",
                                    nucleotide == "C",
                                    position > 429,
                                    position < 464),
             aes(x = as.factor(position), y = probability),
             col = "red", size = 3) +
  theme_bw() +
  labs(x = "Position",
       y = "Probability (C)")
dev.off()


logo1 <- subseq(nalm, start = 431, end = 464)
logo1 <- consensusMatrix(logo1, as.prob = TRUE)
logo2 <- subseq(rmats_c3ss, start = 431, end = 463)
logo2 <- consensusMatrix(logo2, as.prob = TRUE)

ggplot() +
  geom_logo(data = logo1, method = "probability", col_scheme = "nucleotide") +
  theme_bw()
pdf(file = "/figure_two/C3SS_seqlogo.pdf",
    width = 9, height = 4, paper = "letter")
ggplot() +
  geom_logo(data = logo2, method = "probability", col_scheme = "nucleotide") +
  theme_classic()
dev.off()

###################################################
### here is for canonical splice site paired to C3SS

library(Biostrings)
nalm <- readDNAStringSet(filepath = "/data/Nalm6_Regtools_Canonical_Junctions_with_C3SS.regtools.fasta",
                         format = "fasta")


# Define the positions you're interested in
pos_A <- 449
pos_G <- 450

# Filter sequences that have 'A' at position 449 and 'G' at position 450
filtered_nalm <- nalm[substring(as.character(nalm), pos_A, pos_A) == "A" & 
                        substring(as.character(nalm), pos_G, pos_G) == "G"]
nalm <- filtered_nalm
rm(filtered_nalm)


rmats_c3ss <- readDNAStringSet(file = "/figure_two/rMATS_K700E_shared_C3SS_Canonical900.fasta",
                              format = "fasta")


consensus_matrices <- list()

# Define the number of bootstrap iterations
num_bootstraps <- 1000

# Perform bootstrapping
for (i in 1:num_bootstraps) {
  # Sample from vdbc3ss with replacement
  tmp <- sample(nalm, size = 192, replace = TRUE)
  
  # Generate consensus matrix
  pfm <- consensusMatrix(tmp, as.prob = TRUE)
  
  # Convert consensus matrix to dataframe
  pfm_df <- as.data.frame(pfm)
  
  # Add column for nucleotide positions
  pfm_df$nucleotide <- rownames(pfm_df)
  
  # Remove row names
  rownames(pfm_df) <- NULL
  
  # Add bootstrap number
  pfm_df$Bootstrap_Number <- i
  
  # Store the dataframe
  consensus_matrices[[i]] <- pfm_df
}

# Combine all consensus matrices into a single dataframe
combined_df <- bind_rows(consensus_matrices)

# Filter rows where Nucleotide is A, C, T, or G
combined_df <- combined_df %>%
  filter(nucleotide %in% c("A", "C", "T", "G"))

tmp <- pivot_longer(combined_df, -c(nucleotide, Bootstrap_Number),
                    values_to = "probability", names_to = "position")
tmp$position <- sub("V", "", tmp$position)
tmp$position <- as.numeric(tmp$position)

tmp %>%
  filter(
    position > 425, 
    position < 475) %>%
  ggplot(aes(x = as.factor(position), y = probability,
             fill = nucleotide)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90))

pfm1 <- consensusMatrix(rmats_c3ss, as.prob = TRUE)
pfm1 <- as.data.frame(pfm1)
pfm1$nucleotide <- rownames(pfm1)
rownames(pfm1) <- NULL
pfm1 <- pfm1 %>%
  filter(nucleotide %in% c("A", "T", "C", "G"))
tmp2 <- pivot_longer(pfm1, -c(nucleotide), values_to = "probability",
                     names_to = "position")
tmp2$position <- sub("V", "", tmp2$position)
#tmp2$new_position <- as.numeric(tmp2$position)-1
k700e_long <- tmp2

tmp2$source <- "K700E"
tmp2$position <- as.numeric(tmp2$position)
#tmp2$position <- tmp2$position-1
tmp$source <- "Nalm6_WT"
tmp$Bootstrap_Number <- NULL
test <- rbind(tmp, tmp2)
test$position <- as.numeric(test$position)

pdf(file = "/figure_two//Canonical_T_content.pdf",
    width = 9, height = 4, paper = "letter")
test %>%
  filter(source == "Nalm6_WT",
         nucleotide == "T",
         position > 429,
         position < 464) %>%
  ggplot(aes(x = as.factor(position), y = probability)) +
  geom_boxplot() +
  geom_point(data = test %>% filter(source == "K700E",
                                    nucleotide == "T",
                                    position > 429,
                                    position < 464),
             aes(x = as.factor(position), y = probability),
             col = "red", size = 3) +
  theme_bw() +
  labs(x = "Position",
       y = "Probability (T)")
dev.off()

pdf(file = "/figure_two/Canonical_G_content.pdf",
    width = 9, height = 4, paper = "letter")
test %>%
  filter(source == "Nalm6_WT",
         nucleotide == "G",
         position > 429,
         position < 464) %>%
  ggplot(aes(x = as.factor(position), y = probability)) +
  geom_boxplot() +
  geom_point(data = test %>% filter(source == "K700E",
                                    nucleotide == "G",
                                    position > 429,
                                    position < 464),
             aes(x = as.factor(position), y = probability),
             col = "red", size = 3) +
  theme_bw() +
  labs(x = "Position",
       y = "Probability (G)")
dev.off()

pdf(file = "/figure_two/Canonical_A_content.pdf",
    width = 9, height = 4, paper = "letter")
test %>%
  filter(source == "Nalm6_WT",
         nucleotide == "A",
         position > 429,
         position < 464) %>%
  ggplot(aes(x = as.factor(position), y = probability)) +
  geom_boxplot() +
  geom_point(data = test %>% filter(source == "K700E",
                                    nucleotide == "A",
                                    position > 429,
                                    position < 464),
             aes(x = as.factor(position), y = probability),
             col = "red", size = 3) +
  theme_bw() +
  labs(x = "Position",
       y = "Probability (A)")
dev.off()

pdf(file = "/figure_two/Canonical_C_content.pdf",
    width = 9, height = 4, paper = "letter")
test %>%
  filter(source == "Nalm6_WT",
         nucleotide == "C",
         position > 429,
         position < 464) %>%
  ggplot(aes(x = as.factor(position), y = probability)) +
  geom_boxplot() +
  geom_point(data = test %>% filter(source == "K700E",
                                    nucleotide == "C",
                                    position > 429,
                                    position < 464),
             aes(x = as.factor(position), y = probability),
             col = "red", size = 3) +
  theme_bw() +
  labs(x = "Position",
       y = "Probability (C)")
dev.off()

logo1 <- subseq(nalm, start = 431, end = 464)
logo1 <- consensusMatrix(logo1, as.prob = TRUE)
logo2 <- subseq(rmats_c3ss, start = 431, end = 463)
logo2 <- consensusMatrix(logo2, as.prob = TRUE)

ggplot() +
  geom_logo(data = logo1, method = "probability", col_scheme = "nucleotide") +
  theme_bw()
pdf(file = "/figure_two/K700E_Canonical_seqlogo.pdf",
    width = 9, height = 4, paper = "letter")
ggplot() +
  geom_logo(data = logo2, method = "probability", col_scheme = "nucleotide") +
  theme_classic()
dev.off()
