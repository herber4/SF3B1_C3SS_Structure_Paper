library(dplyr)
library(ggplot2)

df <- read.table(file = "Figures/manuscript_v1.0/scanfold/K700E_per_base_z_scores.txt",
                sep = "\t", header = TRUE)

df$gene <- sub("_.*", "", df$gene)


ses <- read.table(file = "Figures/figure_3_new/K700E_SHAPE_C3SS_Canonical_merged.txt",
                  sep = "\t", header = TRUE)

test <- merge(ses, df, by = c("gene", "Nucleotide"))
write.table(test, file = "Figures/manuscript_v1.0/scanfold/K700E_reindexed_perbase_z_scores.txt",
            sep = "\t", row.names = F, quote = F)

test %>%
  filter(source.y %in% c("long_shape"), 
         Position > 120, Position < 180) %>%
  ggplot(aes(x = as.factor(Position), y = avgZ, fill = source.x)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90))


setwd("Figures/Figure_4/data/")

r <- read.table(file = "Figures/manuscript_v1.0/scanfold/Resistant_per_base_z_scores.txt",
                sep = "\t", header = TRUE)

# Assuming your dataframe is named r and the column is named gene
r$gene <- sub("\\..*", "", r$gene)
# Assuming your dataframe is named r and the column is named gene
r$gene <- sub("_.*", "", r$gene)
meta <- read.table(file = "R_Projects/SF3B1_Pubs/master_scripts/shape_master/VastDB/VastDB_gBlock_Meta_Data.txt",
                   sep = "\t", header = TRUE)
meta$gene <- meta$EVENT
meta$EVENT <- NULL
meta$gene <- sub("_.*", "", meta$gene)
meta <- meta[,c(1, 11)]
tmp <- merge(meta, r, by = "gene")
r <- anti_join(r, meta, by = "gene")
tmp$gene <- tmp$sample
tmp$sample <- NULL
r <- rbind(tmp, r)

mas <- read.table(file = "shape_vs_no_shape_resistant/SHAPE_NOSHAPE_combined_both_MDS.txt",
                  sep = "\t", header = TRUE)

mas[,c(4:6, 13)] <- NULL 
mas$gene <- mas$gene_name
mas$gene_name <- NULL
tmp <- merge(r, mas, by = c("gene", "Nucleotide"))
write.table(tmp, file = "Figures/manuscript_v1.0/scanfold/Resistant_reindexed_perbase_z_scores.txt",
            sep = "\t", row.names = F, quote = F)


mas <- read.table(file = "Figures/figure_3_new/Control_SHAPE_BPP_reindexed.txt",
                  sep = "\t", header = TRUE)
con <- read.table(file = "Figures/manuscript_v1.0/scanfold/Control_per_base_z_scores.txt",
                  sep = "\t", header = TRUE)
con$gene <- sub("_.*", "", con$gene)
mas[,10:17] <- NULL

tmp <- merge(con, mas, by = c("gene", "Nucleotide"))
write.table(tmp, file = "Figures/manuscript_v1.0/scanfold/Control_reindexed_perbase_z_scores.txt",
            sep = "\t", row.names = F, quote = F)
