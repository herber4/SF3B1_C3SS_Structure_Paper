library(dplyr)
library(ggplot2)
library(tidyr)

con <- read.table(file = "../Control_reindexed_perbase_z_scores.txt",
                    sep = "\t", header = TRUE)
con[,c(3, 14:15)] <- NULL
con$i <- con$i.x
con$i.x <- NULL
con$type <- "Canonical"
k <- read.table(file = "../K700E_reindexed_perbase_z_scores.txt",
                sep = "\t", header = T)
colnames(k)[10] <- "type"
colnames(k)[18] <- "source"
k[,c(3, 8, 11:12)] <- NULL

r <- read.table(file = "../Resistant_reindexed_perbase_z_scores.txt",
                sep = "\t", header = TRUE)
r$nt <- NULL
con$group <- "Control"
k$group <- "K700E"
r$group <- "Resistant"

mas <- rbind(con, k, r)




colors <- c("K700E_short_shape_Canonical" = "#ae017e",
            "K700E_short_shape_C3SS" = "#fbb4b9",
            "Control_short_shape_Canonical" = "#696969",
            "Resistant_short_shape_Canonical" = "#2171b5",
            "Resistant_short_shape_C3SS" = "#a6cee3")

cols <- c("K700E_long_shape_Canonical" = "#ae017e",
          "K700E_long_shape_C3SS" = "#fbb4b9",
          "Control_long_shape_Canonical" = "#696969",
          "Resistant_long_shape_Canonical" = "#2171b5",
          "Resistant_long_shape_C3SS" = "#a6cee3")
#############

stats <- mas %>%
  group_by(plot, Position) %>%
  summarise(mean_z = mean(avgZ),
            mean_ED = mean(avgED))
long <- stats %>%
  filter(plot %in% c("Resistant_long_shape_C3SS",
                     "K700E_long_shape_C3SS",
                     "Resistant_long_shape_Canonical",
                     "K700E_long_shape_Canonical",
                     "Control_long_shape_Canonical"))
pdf(file = "200md_perbase_z_score_mean_linegraph_fig5.pdf",
    width = 12, height = 3, paper = "letter")
long %>%
  filter(Position >= 30, Position <= 166) %>%
  mutate(plot = factor(plot, levels = c(
    "Resistant_long_shape_C3SS",
    "K700E_long_shape_C3SS",
    "Resistant_long_shape_Canonical",
    "K700E_long_shape_Canonical",
    "Control_long_shape_Canonical"
  ))) %>%
  ggplot(aes(x = as.factor(Position), y = mean_z, col = plot, group = plot)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.position = "") +
  labs(x = "Position", y = "Avg Z Score", fill = "Group") +
  scale_color_manual(values = cols) +
  geom_vline(xintercept = as.factor(151), linetype = "dashed", color = "black")
dev.off()

pdf(file = "200md_perbase_ED_mean_linegraph_fig5.pdf",
    width = 12, height = 3, paper = "letter")
long %>%
  filter(Position >= 30, Position <= 166) %>%
  mutate(plot = factor(plot, levels = c(
    "Resistant_long_shape_C3SS",
    "K700E_long_shape_C3SS",
    "Resistant_long_shape_Canonical",
    "K700E_long_shape_Canonical",
    "Control_long_shape_Canonical"
  ))) %>%
  ggplot(aes(x = as.factor(Position), y = mean_ED, col = plot, group = plot)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.position = "") +
  labs(x = "Position", y = "Avg Ensemble Diversity", fill = "Group") +
  scale_color_manual(values = cols) +
  geom_vline(xintercept = as.factor(151), linetype = "dashed", color = "black")
dev.off()


#############


colors <- c("K700E_short_shape_Canonical" = "#ae017e",
            "K700E_short_shape_C3SS" = "#fbb4b9",
            "Control_short_shape_Canonical" = "#696969",
            "Resistant_short_shape_Canonical" = "#2171b5",
            "Resistant_short_shape_C3SS" = "#a6cee3")


short <- stats %>%
  filter(plot %in% c("Resistant_short_shape_C3SS",
                     "K700E_short_shape_C3SS",
                     "Resistant_short_shape_Canonical",
                     "K700E_short_shape_Canonical",
                     "Control_short_shape_Canonical"))
pdf(file = "27md_perbase_Zscore_mean_linegraph_Suppfig5.pdf",
    width = 12, height = 3, paper = "letter")
short %>%
  filter(Position >= 30, Position <= 166) %>%
  mutate(plot = factor(plot, levels = c(
    "Resistant_short_shape_C3SS",
    "K700E_short_shape_C3SS",
    "Resistant_short_shape_Canonical",
    "K700E_short_shape_Canonical",
    "Control_short_shape_Canonical"
  ))) %>%
  ggplot(aes(x = as.factor(Position), y = mean_z, col = plot, group = plot)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
      legend.position = "") +
  labs(x = "Position", y = "Avg Z Score", fill = "Group") +
  scale_color_manual(values = colors) +
  geom_vline(xintercept = as.factor(151), linetype = "dashed", color = "black")
dev.off()

pdf(file = "27md_perbase_ED_mean_linegraph_Suppfig5.pdf",
    width = 12, height = 3, paper = "letter")
short %>%
  filter(Position >= 30, Position <= 166) %>%
  mutate(plot = factor(plot, levels = c(
    "Resistant_short_shape_C3SS",
    "K700E_short_shape_C3SS",
    "Resistant_short_shape_Canonical",
    "K700E_short_shape_Canonical",
    "Control_short_shape_Canonical"
  ))) %>%
  ggplot(aes(x = as.factor(Position), y = mean_ED, col = plot, group = plot)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.position = "") +
  labs(x = "Position", y = "Avg Ensemble Diversity", fill = "Group") +
  scale_color_manual(values = colors) +
  geom_vline(xintercept = as.factor(151), linetype = "dashed", color = "black")
dev.off()


pdf(file = "27_window_ED_boxplot_Suppfig5.pdf",
    width = 9, height = 5, paper = "letter")
short %>%
  filter(Position >= 30, Position <= 166) %>%
  mutate(plot = factor(plot, levels = c(
    "K700E_short_shape_C3SS",
    "K700E_short_shape_Canonical",
    "Resistant_short_shape_C3SS",
    "Resistant_short_shape_Canonical",
    "Control_short_shape_Canonical"
  ))) %>%
  filter(Position >= 30, Position <= 160) %>%
  ggplot(aes(x = plot, y = mean_ED, fill = plot)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = colors)
dev.off()

pdf(file = "27_window_Zscore_boxplot_Suppfig5.pdf",
    width = 9, height = 5, paper = "letter")
short %>%
  filter(Position >= 30, Position <= 166) %>%
  mutate(plot = factor(plot, levels = c(
    "K700E_short_shape_C3SS",
    "K700E_short_shape_Canonical",
    "Resistant_short_shape_C3SS",
    "Resistant_short_shape_Canonical",
    "Control_short_shape_Canonical"
  ))) %>%
  filter(Position >= 30, Position <= 160) %>%
  ggplot(aes(x = plot, y = mean_z, fill = plot)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = colors)
dev.off()



pdf(file = "200_window_ED_boxplot_fig5.pdf",
    width = 9, height = 5, paper = "letter")
long %>%
  filter(Position >= 30, Position <= 166) %>%
  mutate(plot = factor(plot, levels = c(
    "K700E_long_shape_C3SS",
    "K700E_long_shape_Canonical",
    "Resistant_long_shape_C3SS",
    "Resistant_long_shape_Canonical",
    "Control_long_shape_Canonical"
  ))) %>%
  ggplot(aes(x = plot, y = mean_ED, fill = plot)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = cols) +
  labs(x = "",
       y = "Avg Ensemble Diversity")
dev.off()

pdf(file = "200_window_Z_score_boxplot_fig5.pdf",
    width = 9, height = 5, paper = "letter")
long %>%
  filter(Position >= 30, Position <= 166) %>%
  mutate(plot = factor(plot, levels = c(
    "K700E_long_shape_C3SS",
    "K700E_long_shape_Canonical",
    "Resistant_long_shape_C3SS",
    "Resistant_long_shape_Canonical",
    "Control_long_shape_Canonical"
  ))) %>%
  ggplot(aes(x = plot, y = mean_z, fill = plot)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = cols) +
  labs(x = "",
       y = "Avg Z Score")
dev.off()
