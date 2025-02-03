library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)

setwd("/")

mas <- read.table(file = "/Users/herber4/Desktop/Figures/Figure_3/data/K700E_Adjusted_SHAPE_R1.txt",
                  sep = "\t", header = TRUE)
ses <- mas
colnames(ses)[c(5,7)] <- c("gene", "Position")

result <- ses %>%
  arrange(gene, Position) %>%
  group_by(gene) %>%
  filter(
    (Sequence == "A" & lead(Sequence) == "G" & lead(Position) == Position + 1) |
      (Sequence == "G" & lag(Sequence) == "A" & lag(Position) == Position - 1)
  )

result$plot <- "Random_AG"
result <- result %>%
  filter(Sequence == "G")

tmp <- result[,c(5,7)]
tmp <- tmp %>%
  filter(Position > 160,
         Position < 900)
tmp <- tmp %>%
  filter(Position > 650 |
           Position < 350)
# Expand the dataframe
tmp$unique_id <- paste(tmp$gene, tmp$Position, sep = ":")

tmp <- tmp %>%
  rowwise() %>%
  mutate(new_positions = list((Position - 10):(Position + 10))) %>%
  unnest(new_positions) %>%
  ungroup()

tmp$Position <- tmp$new_positions
tmp$new_positions <- NULL
tmp <- merge(ses, tmp, by = c("gene", "Position"))


tmp <- tmp %>%
  arrange(unique_id, Position)
# Add a new column 'Position' with values 1 through 60 for each gene
test <- tmp %>%
  group_by(gene, unique_id) %>%
  mutate(new_index = row_number())

short <- read.table(file = "/shape_vs_no_shape_bpp_k700e/K700E_SHAPE_md28_34_avs.txt",
                    sep = "\t", header = TRUE)

colnames(short)[2] <- "Nucleotide"
short$gene <- sub("_.*", "", short$gene)
short <- merge(short, test, by = c("gene", "Nucleotide"))
short$plot <- "Short_MD_OutSJs_randomAG"
long <- read.table(file = "/shape_vs_no_shape_bpp_k700e/K700E_SHAPE_md200_225_avs.txt",
                   sep = "\t", header = TRUE)
colnames(long)[2] <- "Nucleotide"
long$gene <- sub("_.*", "", long$gene)
long <- merge(long, test, by = c("gene", "Nucleotide"))
long$plot <- "Long_MD_OutSJs_randomAG"

mas <- rbind(long, short)
write.table(mas, file = "/All_SHAPE_Random_AG_for_linegraphs.txt",
            sep = "\t", row.names = F, quote = F)
