library(ggplot2)
library(dplyr)
library(ggbeeswarm)
library(ggridges)

con <- read.table(file = "Desktop/Figures/figure_3_new/Control_SHAPE_BPP_reindexed.txt",
                  sep = "\t", header = TRUE)
con[,c(10:17)] <- NULL
con <- con[,c(1,3,6,9)]
con$type <- "Control_3SS"
con$source <- "Control_3SS"
con$plot <- "Control_3SS"


res <- read.table(file = "Desktop/Figures/figure_3_new/VastDB_adjusted_indexed_SHAPE.txt",
                  sep = "\t", header = TRUE)

ses <- read.table(file = "Desktop/Figures/figure_3_new/K700E_SHAPE_C3SS_Canonical_merged.txt",
                  sep = "\t", header = TRUE)
ses <- ses[,c(1,4,7:10)]
ses$type <- ses$source
ses$source <- "Sensitive"
ses$c3ss_direction <- NULL
res <- res[,c(1,3,6:8)]
res$source <- "Resistant"
res$plot <- ifelse(res$type == "C3SS", res$plot <- "Resistant_C3SS", res$plot <- "Resistant_Canonical")
colnames(res)[1] <- "gene"
ses$plot <- ifelse(ses$type == "C3SS", ses$plot <- "Sensitive_C3SS", ses$plot <- "Sensitive_Canonical")

tmp <- rbind(res, ses, con)
rm(res, ses, con)
res_c3ss <- tmp %>%
  filter(plot == "Resistant_C3SS")
res_c3ss <- median(res_c3ss$adjusted_profile, na.rm = TRUE)
ses_c3ss <- tmp %>%
  filter(plot == "Sensitive_C3SS")
ses_c3ss <- median(ses_c3ss$adjusted_profile, na.rm = TRUE)



pdf(file = "Desktop/Figures/figure_3_new/K700E_SHAPE.pdf",
    width = 9, height = 4, paper = "letter")
tmp %>%
  filter(Position >= 149, Position <= 151,
         !gene %in% c("FBXO41", "SRRT")) %>%
  mutate(plot = factor(plot, levels = c("Resistant_C3SS", "Sensitive_C3SS", "Resistant_Canonical","Sensitive_Canonical", "Control_3SS"))) %>%
  ggplot(aes(x = as.factor(Position), y = adjusted_profile, fill = plot)) +
  geom_boxplot() +
  theme_bw() +
  geom_hline(yintercept = res_c3ss, linetype = "dashed", linewidth = 1, col = "red") +
  geom_hline(yintercept = ses_c3ss, linetype = "dashed", linewidth = 1, col = "blue") +
  scale_fill_manual(values = c("Resistant_C3SS" = "#bdc9e1",
                               "Sensitive_C3SS" = "#67a9cf",
                               "Resistant_Canonical" = "#fdc086",
                               "Sensitive_Canonical" = "#fdc086",
                               "Control_3SS" = "#969696"))
dev.off()

pdf(file = "Desktop/RNASoc_Nov24/Sensitive_vs_Resistant_NAG_SHAPE.pdf",
    width = 8, height = 4, paper = "letter")
tmp %>%
  filter(Position >= 149, Position <= 151,
         !gene %in% c("FBXO41", "SRRT"),
         plot %in% c("Resistant_C3SS", "Sensitive_C3SS",
                     "Control_3SS")) %>%
  mutate(plot = factor(plot, levels = c("Resistant_C3SS", "Sensitive_C3SS", "Resistant_Canonical","Sensitive_Canonical", "Control_3SS"))) %>%
  ggplot(aes(x = as.factor(Position), y = adjusted_profile, fill = plot)) +
  geom_boxplot() +
  theme_bw() +
  geom_hline(yintercept = res_c3ss, linetype = "dashed", linewidth = 1, col = "#a6cee3") +
  geom_hline(yintercept = ses_c3ss, linetype = "dashed", linewidth = 1, col = "#fbb4b9") +
  scale_fill_manual(values = c("Resistant_C3SS" = "#a6cee3",
                               "Sensitive_C3SS" = "#fbb4b9",
                               "Resistant_Canonical" = "#67a9cf",
                               "Sensitive_Canonical" = "#f768a1",
                               "Control_3SS" = "#969696")) +
  labs(x = "Source",
       y = "Shape Reactivity") +
  theme(text = element_text(size = 18),
        axis.text.x = element_blank())
dev.off()

pdf(file = "Desktop/Figures/figure_3_new/Sensitive_vs_Resistant_NAG_SHAPE_no_legend.pdf",
    width = 8, height = 4, paper = "letter")
tmp %>%
  filter(Position >= 149, Position <= 151,
         !gene %in% c("FBXO41", "SRRT")) %>%
  mutate(plot = factor(plot, levels = c("Resistant_C3SS", "Sensitive_C3SS", "Resistant_Canonical","Sensitive_Canonical", "Control_3SS"))) %>%
  ggplot(aes(x = as.factor(Position), y = adjusted_profile, fill = plot)) +
  geom_boxplot() +
  theme_bw() +
  geom_hline(yintercept = res_c3ss, linetype = "dashed", linewidth = 1, col = "#a6cee3") +
  geom_hline(yintercept = ses_c3ss, linetype = "dashed", linewidth = 1, col = "#fbb4b9") +
  scale_fill_manual(values = c("Resistant_C3SS" = "#a6cee3",
                               "Sensitive_C3SS" = "#fbb4b9",
                               "Resistant_Canonical" = "#67a9cf",
                               "Sensitive_Canonical" = "#f768a1",
                               "Control_3SS" = "#969696")) +
  labs(x = "",
       y = "") +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        text = element_text(size = 18))
dev.off()

pdf(file = "Desktop/RNASoc_Nov24/all_shape.pdf",
    width = 8, height = 8, paper = "letter")
ggarrange(sensitive, resistant, all, ncol = 1, nrow = 3)
dev.off()


#a6cee3

test_shape_reactivity <- function(data, plot_values, position) {
  # Filter data based on plot values and position
  subset_data <- data %>%
    filter(plot %in% plot_values, Position == position) 
  
  # Separate the groups
  group1_data <- subset_data %>%
    filter(plot == plot_values[1]) %>%
    pull(adjusted_profile)
  
  group2_data <- subset_data %>%
    filter(plot == plot_values[2]) %>%
    pull(adjusted_profile)
  
  # Shapiro-Wilk test for normality
  shapiro_group1 <- shapiro.test(group1_data)
  shapiro_group2 <- shapiro.test(group2_data)
  
  # Perform Wilcoxon test
  wilcox_result <- wilcox.test(adjusted_profile ~ plot, data = subset_data)
  
  # Output the results
  list(
    Shapiro_Wilk_Group1 = shapiro_group1,
    Shapiro_Wilk_Group2 = shapiro_group2,
    Wilcoxon_Test = wilcox_result
  )
}

# Example of using the function
test_shape_reactivity(tmp, c("Resistant_C3SS", "Resistant_Canonical"), 150)

test_shape_reactivity(tmp, c("Resistant_C3SS", "Sensitive_C3SS"), 150)
test_shape_reactivity(tmp, c("Resistant_C3SS", "Sensitive_Canonical"), 150)
test_shape_reactivity(tmp, c("Resistant_C3SS", "Control_3SS"), 150)



medians <- tmp %>%
  group_by(plot) %>%
  summarise(median = median(adjusted_profile, na.rm = TRUE))

pdf(file = "Desktop/Figures/figure_3_new/shape_reactivities_boxplot.pdf",
    width = 8, height = 8, paper = "letter")
ggplot(tmp, aes(x = plot, y = -log10(adjusted_profile), group = plot)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  labs(y = "Shape Reactivity",
       x = "Source") +
  scale_y_continuous(limits = quantile(tmp))
dev.off()

ggplot(tmp, aes(x = plot, y = -log(adjusted_profile, 2))) +
  geom_quasirandom() +
  theme_bw() +
  labs(y = "Shape Reactivity",
       x = "Source")

ggplot(tmp, aes(x = adjusted_profile, y = plot, fill = plot)) +
  geom_density_ridges()


# Function to perform pairwise Wilcoxon tests
perform_wilcox_tests <- function(data) {
  # Remove NAs from the adjusted_profile column
  data <- data[!is.na(data$adjusted_profile), ]
  
  # Perform pairwise Wilcoxon tests
  result <- pairwise.wilcox.test(data$adjusted_profile, data$plot,
                                 p.adjust.method = "none") # Or adjust the method as needed
  return(result)
}

pos <- tmp %>%
  filter(Position == 150)
# Apply the function to your dataframe
perform_wilcox_tests(pos)
