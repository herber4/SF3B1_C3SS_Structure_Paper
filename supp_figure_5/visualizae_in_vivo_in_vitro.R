library(tidyverse)
library(ggplot2)
library(plotly)
library(gapminder)

setwd("/Users/herber4/Desktop/Figures/in_vivo/full_seq_w_in_vivo_shape/")

vitro <- read.table(file = "../figures/in_vitro_scanfold_master.txt",
                    sep = "\t", header = TRUE)
vitro_stats <- vitro %>%
  group_by(plot, Position) %>%
  summarise(Z = mean(avgZ),
            ED = mean(avgED))

Control_long_shape_Canonical
K700E_long_shape_C3SS
K700E_long_shape_Canonical
Resistant_long_shape_C3SS
Resistant_long_shape_Canonical



vitro_means <- vitro_stats %>%
  group_by(plot) %>%
  filter(Position >= 30, Position <= 166) %>%
  summarise(Z = mean(Z),
            ED = mean(ED))

scan <- read.table(file = "rbound_vivo_vitro_raw_scanfold_indexed.txt",
            sep = "\t", header = TRUE)

scan %>%
  filter(source == "long_shape",
         Position >= 100, Position <= 200) %>%
  ggplot(aes(x = as.factor(Position), y = avgED, fill = interaction(experiment, splice_site))) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, size = 4)) +
  scale_fill_manual(values = c("in_vitro.C3SS" = "#ae017e",
                                "in_vitro.Canonical" = "#2171b5",
                                "in_vivo.C3SS" = "#fbb4b9",
                                "in_vivo.Canonical" = "#a6cee3")) +
  facet_grid(rows = vars(splice_site)) +
  geom_vline(xintercept = as.factor(125), linetype = "dashed", color = "black") +
  geom_hline(yintercept = 42.42, linetype = "dashed", color = "grey")

setwd("/Users/herber4/Desktop/Figures/in_vivo/figures/")
pdf(file = "bub1b_text_ED.pdf",
    width = 12, height = 3, paper = "letter")
scan %>%
  filter(source == "long_shape",
         Position >= 30, Position <= 166,
         gene == "BUB1B") %>%
  ggplot(aes(x = as.factor(Position), y = avgED, col = interaction(experiment, splice_site),
             group = interaction(experiment, splice_site))) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 4),
        legend.position = "none") +
  scale_color_manual(values = c("in_vitro.C3SS" = "#ae017e",
                               "in_vitro.Canonical" = "#2171b5",
                               "in_vivo.C3SS" = "#fbb4b9",
                               "in_vivo.Canonical" = "#a6cee3")) +
  facet_grid(rows = vars(splice_site)) +
  geom_vline(xintercept = as.factor(151), linetype = "dashed", color = "black") +
  labs(x = "Position", y = "Ensemble Diversity", fill = "Group") +
  ggtitle("BUB1B")
dev.off()

275.5378
84.451

unique(scan$gene)

library(dplyr)
library(ggplot2)
library(purrr)

# Get all unique genes
unique_genes <- unique(scan$gene)

library(dplyr)
library(ggplot2)
library(purrr)

unique_genes <- unique(scan$gene)

walk(unique_genes, function(g) {
  df <- scan %>%
    filter(source == "long_shape",
           Position >= 30, Position <= 166,
           gene == g)
  
  if (nrow(df) > 0) {
    pdf(file = paste0(g, "_ED.pdf"), width = 12, height = 3, paper = "letter", useDingbats = FALSE)
    
    print(
      ggplot(df, aes(x = as.factor(Position), y = avgED, 
                     col = interaction(experiment, splice_site),
                     group = interaction(experiment, splice_site))) +
        geom_point() +
        geom_line() +
        theme_bw() +
        theme(axis.text.x = element_blank(),
              legend.position = "none") +
        scale_color_manual(values = c("in_vitro.C3SS" = "#ae017e",
                                      "in_vitro.Canonical" = "#2171b5",
                                      "in_vivo.C3SS" = "#fbb4b9",
                                      "in_vivo.Canonical" = "#a6cee3")) +
        facet_grid(rows = vars(splice_site)) +
        geom_vline(xintercept = as.factor(151), linetype = "dashed", color = "black") +
        labs(y = "Ensemble Diversity", x = "", fill = "Group") +
        ggtitle(g)
    )
    
    dev.off()
  }
})


280.8324
86.0738
stats <- scan %>%
  group_by(source, experiment, splice_site, Position) %>%
  summarise(MFE = mean(avgMFE),
            ED = mean(avgED),
            Z = mean(avgZ))


pdf(file = "ED_means.pdf",
    width = 12, height = 3, paper = "letter")
stats %>%
  filter(source == "long_shape",
         Position >= 30, Position <= 166) %>%
  ggplot(aes(x = as.factor(Position), y = ED, col = interaction(experiment, splice_site), group = interaction(experiment, splice_site))) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = c("in_vitro.C3SS" = "#ae017e",
                                "in_vitro.Canonical" = "#2171b5",
                                "in_vivo.C3SS" = "#fbb4b9",
                                "in_vivo.Canonical" = "#a6cee3")) +
  facet_grid(rows = vars(splice_site)) +
  geom_vline(xintercept = as.factor(151), linetype = "dashed", color = "black") +
  geom_hline(yintercept = 42.3, linetype = "dashed", color = "grey")
dev.off()

pdf(file = "Z_means.pdf",
    width = 12, height = 3, paper = "letter")
stats %>%
  filter(source == "long_shape",
         Position >= 30, Position <= 166) %>%
  ggplot(aes(x = as.factor(Position), y = Z, col = interaction(experiment, splice_site), group = interaction(experiment, splice_site))) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = c("in_vitro.C3SS" = "#ae017e",
                                "in_vitro.Canonical" = "#2171b5",
                                "in_vivo.C3SS" = "#fbb4b9",
                                "in_vivo.Canonical" = "#a6cee3")) +
  facet_grid(rows = vars(splice_site)) +
  geom_vline(xintercept = as.factor(151), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -1.008, linetype = "dashed", color = "grey")
dev.off()




pdf(file = "ED_boxie.pdf",
    width = 9, height = 5, paper = "letter")
scan %>%
  filter(source == "long_shape",
         Position >= 30, Position <= 166) %>%
  ggplot(aes(x = interaction(experiment, splice_site), y = avgED, fill = interaction(experiment, splice_site))) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none") +
  scale_fill_manual(values = c("in_vitro.C3SS" = "#ae017e",
                                "in_vitro.Canonical" = "#2171b5",
                                "in_vivo.C3SS" = "#fbb4b9",
                                "in_vivo.Canonical" = "#a6cee3"))
dev.off()

pdf(file = "Z_boxie.pdf",
    width = 9, height = 5, paper = "letter")
scan %>%
  filter(source == "long_shape",
         Position >= 30, Position <= 166) %>%
  ggplot(aes(x = interaction(experiment, splice_site), y = avgZ, fill = interaction(experiment, splice_site))) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none") +
  scale_fill_manual(values = c("in_vitro.C3SS" = "#ae017e",
                               "in_vitro.Canonical" = "#2171b5",
                               "in_vivo.C3SS" = "#fbb4b9",
                               "in_vivo.Canonical" = "#a6cee3"))
dev.off()

269.6048
131.1763