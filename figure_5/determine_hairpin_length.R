library(ggplot2)
library(dplyr)
library(ggpubr)

mas <- read.table(file = "Can_Cryptic_MT_WT_All_Random_SHAPE.27.mfe",
                  sep = "\t", header = TRUE)
mas <- read.table(file = "Can_Cryptic_MT_WT_All_Random_SHAPE.No_MD.mfe",
                  sep = "\t", header = TRUE)

ses <- read.table(file = "Figures/figure_3_new/K700E_SHAPE_C3SS_Canonical_merged.txt",
                  sep = "\t", header = TRUE)
ses <- ses[,c(1, 8)]
ses <- ses %>%
  filter(c3ss_direction == "Upstream")
# Add the count column based on the structure column
mas <- mas %>%
  mutate(count = ifelse(structure %in% c(")", "("), 1, 0))

cryp <- mas %>%
  filter(source == "MT_C3SS")
can <- mas %>%
  filter(source == "MT_Canonical")

cryp <- semi_join(cryp, ses, by = c("gene"))
can <- semi_join(can, ses, by = c("gene"))
mas <- mas %>%
  filter(!source %in% c("MT_C3SS", "MT_Canonical"))
mas <- rbind(can, mas, cryp)

gene_counts <- mas %>%
  group_by(source) %>%
  summarise(total_genes =  n_distinct(gene))
stats <- mas %>%
  group_by(source, gene) %>%
  summarise(sum = sum(count))

colors <- c("MT_Canonical" = "#ae017e",
            "MT_C3SS" = "#fbb4b9",
            "Random_AG_All_SHAPE" = "#696969",
            "WT_Canonical" = "#2171b5",
            "WT_Cryptic" = "#a6cee3")

ggplot(stats, aes(x = source, y = sum, fill = source)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = colors)

hps <- mas %>%
  group_by(source, gene) %>%
  summarise(sum = sum(count))
hps$hp_length <- hps$sum / 2

ggplot(hps, aes(x = source, y = hp_length, fill = source)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = colors) 

hps$source <- factor(hps$source, levels = c(
  "MT_C3SS",
  "MT_Canonical",
  "WT_Cryptic",
  "WT_Canonical",
  "Random_AG_All_SHAPE"
))

pdf(file = "Figures/Figure_5/HairPin_Length_boxie_upstream_only_noMD.pdf",
    width = 6, height = 8, paper = "letter")
ggplot(hps, aes(x = source, y = hp_length, fill = source)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = colors) +
  labs(x = "", 
       y = "Hair Pin Length")
dev.off()

pos <- mas %>%
  group_by(source, Position) %>%
  summarise(pos_sum = sum(count))


# Assuming 'gene_counts' has 'source' and 'total_genes'
pos <- pos %>%
  left_join(gene_counts, by = "source") %>%        # Join to bring in 'total_genes'
  mutate(percent = pos_sum / total_genes * 100)    # Calculate percentage

averaged_pos <- pos %>%
  mutate(Position_Group = ceiling(Position / 4)) %>% # Create groups of 4 positions
  group_by(source, Position_Group) %>% # Group by source and the new position group
  summarise(avg_percent = mean(percent, na.rm = TRUE), .groups = "drop") # Calculate average percent

# View the resulting dataframe
print(averaged_pos)

values <- c("MT_Canonical" = "#ae017e",
            "MT_C3SS" = "#fbb4b9",
            "Random_AG_All_SHAPE" = "#bababa")

pdf(file = "Figures/Figure_5/MT_percent_linegraph_averaged.pdf",
    width = 8, height = 3, paper = "letter")
averaged_pos %>%
  filter(source %in% c("MT_C3SS",
                       "MT_Canonical",
                       "Random_AG_All_SHAPE")) %>%
  ggplot(aes(x = Position_Group, y = avg_percent, col = source, group = source)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = values) +
  theme_bw() +
  geom_vline(xintercept = 30, linetype = "dashed", color = "black")
dev.off()

values <- c("WT_Canonical" = "#2171b5",
            "WT_Cryptic" = "#a6cee3",
            "Random_AG_All_SHAPE" = "#bababa")
pdf(file = "Figures/Figure_5/WT_Percent_linegraph_averaged.pdf",
    width = 8, height = 3, paper = "letter")
averaged_pos %>%
  filter(source %in% c("WT_Cryptic",
                       "WT_Canonical",
                       "Random_AG_All_SHAPE")) %>%
  ggplot(aes(x = Position_Group, y = avg_percent, col = source, group = source)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = values) +
  theme_bw() +
  geom_vline(xintercept = 30, linetype = "dashed", color = "black")
dev.off()


values <- c("WT_Cryptic" = "#a6cee3",
            "MT_C3SS" = "#fbb4b9",
            "Random_AG_All_SHAPE" = "#bababa")

pdf(file = "Figures/Figure_5/WT_vs_MT_Cryptic_percent_linegraph_averaged.pdf",
    width = 8, height = 3, paper = "letter")
averaged_pos %>%
  filter(source %in% c("WT_Cryptic",
                       "MT_C3SS",
                       "Random_AG_All_SHAPE")) %>%
  ggplot(aes(x = Position_Group, y = avg_percent, col = source, group = source)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = values) +
  theme_bw() +
  geom_vline(xintercept = 30, linetype = "dashed", color = "black")
dev.off()


values <- c("WT_Canonical" = "#2171b5",
            "MT_Canonical" = "#ae017e",
            "Random_AG_All_SHAPE" = "#bababa")
pdf(file = "Figures/Figure_5/WT_vs_MT_Canonical_percent_linegraph_averaged.pdf",
    width = 8, height = 3, paper = "letter")
averaged_pos %>%
  filter(source %in% c("WT_Canonical",
                       "MT_Canonical",
                       "Random_AG_All_SHAPE")) %>%
  ggplot(aes(x = Position_Group, y = avg_percent, col = source, group = source)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = values) +
  theme_bw() +
  geom_vline(xintercept = 30, linetype = "dashed", color = "black")
dev.off()



pairwise.wilcox.test(
  hps$hp_length,
  hps$source,
  p.adjust.method = "none" # Adjust p-values for multiple testing
)
