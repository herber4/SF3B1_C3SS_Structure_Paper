library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)

setwd("")

all_random_ag_shape <- read.table(file = "VastDB_RES_shortMD_outsideSJs_randomAGs_bpp_and_SHAPE.txt",
                                  sep = "\t", header = TRUE)

mas <- read.table(file = "shape_vs_no_shape_resistant/SHAPE_NOSHAPE_combined_both_MDS.txt",
                  sep = "\t", header = TRUE)
ses <- mas
ses <- ses %>%
  filter(type == "C3SS",
         !Position == 150,
         !Position == 151,
         plot == "WT_C3SS_md28_to_34")
can <- mas %>%
  filter(type == "Canonical",
         Position >= 141,
         Position <= 161)
ses <- anti_join(ses, can, by = c("gene_name", "adjusted_profile", "av_sum_p","Sequence"))
# Ensure the data is sorted by gene and Position
ses <- ses %>% arrange(gene, Position)

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

tmp <- result[,c(2,11)]
tmp <- tmp %>%
  filter(Position > 10)
# Expand the dataframe
tmp$unique_id <- paste(tmp$gene_name, tmp$Position, sep = ":")

tmp <- tmp %>%
  rowwise() %>%
  mutate(new_positions = list((Position - 10):(Position + 10))) %>%
  unnest(new_positions) %>%
  ungroup()

tmp$Position <- NULL
colnames(tmp)[3] <- "Position"

c3ss <- mas %>%
  filter(type == "C3SS",
         plot == "WT_C3SS_md28_to_34")
random_ag <- merge(tmp, c3ss, by = c("gene_name", "Position"))
random_ag[,c(4,5)] <- NULL
# Add a new column 'Position' with values 1 through 60 for each gene
random_ag <- random_ag %>%
  group_by(gene_name, unique_id, plot) %>%
  mutate(new_index = row_number())

ss_ag <- mas %>%
  filter(Position >= 141, Position <= 161)
ss_ag[,c(1,3)] <- NULL
ss_ag$new_index <- ss_ag$Position-140
ss_ag$Nucleotide <- NULL
random_ag$unique_id <- NULL
random_ag$new_plot <- paste("random_ag", random_ag$plot, sep = ":")
ss_ag$new_plot <- paste("SS_AG", ss_ag$plot, sep = ":")

mas <- rbind(ss_ag, random_ag)
rm(ss_ag, tmp, result, random_ag, c3ss, ses, can)
mas$plot <- NULL
colnames(mas)[8] <- "source"

con <- read.table(file = "shape_vs_no_shape_controls/Control_SHAPE_bpp_indexed_for_line_graphs.txt",
                  sep = "\t", header = TRUE)
con[,c(3,8)] <- NULL
colnames(con)[1] <- "gene_name"
con$new_index <- con$Position
con$new_plot <- paste("Control", con$source, sep = ":")

mas <- rbind(mas, con)
rm(con)

colnames(all_random_ag_shape)[11] <- "new_plot"
all_random_ag_shape <- all_random_ag_shape[,c(3,6,9:12)]
mas <- mas[,c(10:12, 5, 2, 1)]
colnames(all_random_ag_shape)[4] <- "gene_name"

mas <- rbind(all_random_ag_shape, mas)

prof_mean <- mas %>%
  filter(new_plot == "shortMD_outsideSJs_randomAGs")
prof_mean <- mean(prof_mean$adjusted_profile, na.rm = TRUE)  
  
acc_mean <- mas %>%
  filter(new_plot == "shortMD_outsideSJs_randomAGs")
acc_mean <- mean(acc_mean$av_sum_p, na.rm = TRUE)

stats <- mas %>%
  group_by(new_plot, new_index) %>%
  summarise(av_sum_p = mean(av_sum_p, na.rm = TRUE),
            sd_sum_p = sd(av_sum_p, na.rm = TRUE),
            av_adj_prof = mean(adjusted_profile, na.rm = TRUE),
            sd_adj_prof = sd(adjusted_profile, na.rm = TRUE))


mas %>%
  filter(new_plot %in% c("random_ag:WT_C3SS_md28_to_34",
                         "SS_AG:WT_Can_md28_to_34",
                         "SS_AG:WT_C3SS_md28_to_34",
                         "shortMD_outsideSJs_randomAGs")) %>%
  ggplot(aes(x = as.factor(new_index), y = 1-av_sum_p, fill = new_plot)) +
  geom_boxplot()

#####
short md
#####
filt <- c("SS_AG:WT_Can_md28_to_34",
          "SS_AG:WT_C3SS_md28_to_34",
          "shortMD_outsideSJs_randomAGs")

values <- c("SS_AG:WT_Can_md28_to_34" = "#2171b5",
            "SS_AG:WT_C3SS_md28_to_34" = "#a6cee3",
            "shortMD_outsideSJs_randomAGs" = "black")

pdf(file = "WT_Accessibility_linegraph.pdf",
    width = 8, height = 4, paper = "letter")
ggplot(stats %>%
         filter(new_plot %in% filt), 
       aes(x = as.factor(new_index), y = 1 - av_sum_p, group = new_plot, col = new_plot)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(x = "Position", y = "Accessibility", color = "Plot") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0,1) +
  scale_color_manual(values = values) +
  geom_hline(yintercept = 0.481, linetype = "dashed", color = "black")
dev.off()
short_acc  

pdf(file = "WT_SHAPE_linegraph.pdf",
    width = 8, height = 4, paper = "letter")
ggplot(stats %>%
         filter(new_plot %in% filt), 
       aes(x = as.factor(new_index), y = av_adj_prof, group = new_plot, col = new_plot)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(x = "Position", y = "SHAPE Reactivity", color = "Plot") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0,1) +
  scale_color_manual(values = values) +
  geom_hline(yintercept = 0.504, linetype = "dashed", color = "black")
dev.off()


sig <- mas %>%
  filter(new_plot %in% c("SS_AG:WT_Can_md28_to_34",
                         "SS_AG:WT_C3SS_md28_to_34",
                         "shortMD_outsideSJs_randomAGs"))


# Load necessary library
library(dplyr)

# Initialize an empty data frame to store the results
wilcoxon_results <- data.frame(
  new_index = integer(),
  group1 = character(),
  group2 = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each unique position in `new_index`
for (pos in unique(sig$new_index)) {
  # Subset the data for the current position
  subset_data <- sig %>% filter(new_index == pos)
  
  # Get all unique pairs of groups in `new_plot`
  group_combinations <- combn(unique(subset_data$new_plot), 2, simplify = FALSE)
  
  # Perform Wilcoxon test for each pair of groups
  for (groups in group_combinations) {
    group1 <- groups[1]
    group2 <- groups[2]
    
    # Subset data for the two groups
    data_group1 <- subset_data %>% filter(new_plot == group1) %>% pull(av_sum_p)
    data_group2 <- subset_data %>% filter(new_plot == group2) %>% pull(av_sum_p)
    
    # Perform Wilcoxon rank-sum test
    test_result <- wilcox.test(data_group1, data_group2, exact = FALSE)
    
    # Append the result to the results data frame
    wilcoxon_results <- rbind(
      wilcoxon_results,
      data.frame(
        new_index = pos,
        group1 = group1,
        group2 = group2,
        p_value = test_result$p.value
      )
    )
  }
}

# View the results
wilcoxon_results

p <- wilcoxon_results %>%
  filter(p_value < .05)


wilcoxon_results <- data.frame(
  new_index = integer(),
  group1 = character(),
  group2 = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each unique position in `new_index`
for (pos in unique(sig$new_index)) {
  # Subset the data for the current position
  subset_data <- sig %>% filter(new_index == pos)
  
  # Get all unique pairs of groups in `new_plot`
  group_combinations <- combn(unique(subset_data$new_plot), 2, simplify = FALSE)
  
  # Perform Wilcoxon test for each pair of groups
  for (groups in group_combinations) {
    group1 <- groups[1]
    group2 <- groups[2]
    
    # Subset data for the two groups
    data_group1 <- subset_data %>% filter(new_plot == group1) %>% pull(adjusted_profile)
    data_group2 <- subset_data %>% filter(new_plot == group2) %>% pull(adjusted_profile)
    
    # Perform Wilcoxon rank-sum test
    test_result <- wilcox.test(data_group1, data_group2, exact = FALSE)
    
    # Append the result to the results data frame
    wilcoxon_results <- rbind(
      wilcoxon_results,
      data.frame(
        new_index = pos,
        group1 = group1,
        group2 = group2,
        p_value = test_result$p.value
      )
    )
  }
}

# View the results
wilcoxon_results

p <- wilcoxon_results %>%
  filter(p_value < .05)


























#######
visualize long md
#######

library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)

setwd("")


mas <- read.table(file = "shape_vs_no_shape_resistant/SHAPE_NOSHAPE_combined_both_MDS.txt",
                  sep = "\t", header = TRUE)
ses <- mas
ses <- ses %>%
  filter(type == "C3SS",
         !Position == 150,
         !Position == 151,
         plot == "WT_C3SS_md200_to_225")
can <- mas %>%
  filter(type == "Canonical",
         Position >= 141,
         Position <= 161)
ses <- anti_join(ses, can, by = c("gene_name", "adjusted_profile", "av_sum_p","Sequence"))
# Ensure the data is sorted by gene and Position
ses <- ses %>% arrange(gene, Position)

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

tmp <- result[,c(2,11)]
tmp <- tmp %>%
  filter(Position > 10)
# Expand the dataframe
tmp$unique_id <- paste(tmp$gene_name, tmp$Position, sep = ":")

tmp <- tmp %>%
  rowwise() %>%
  mutate(new_positions = list((Position - 10):(Position + 10))) %>%
  unnest(new_positions) %>%
  ungroup()

tmp$Position <- NULL
colnames(tmp)[3] <- "Position"

c3ss <- mas %>%
  filter(type == "C3SS",
         plot == "WT_C3SS_md200_to_225")
random_ag <- merge(tmp, c3ss, by = c("gene_name", "Position"))
random_ag[,c(4,5)] <- NULL
# Add a new column 'Position' with values 1 through 60 for each gene
random_ag <- random_ag %>%
  group_by(gene_name, unique_id, plot) %>%
  mutate(new_index = row_number())

ss_ag <- mas %>%
  filter(Position >= 141, Position <= 161)
ss_ag[,c(1,3)] <- NULL
ss_ag$new_index <- ss_ag$Position-140
ss_ag$Nucleotide <- NULL
random_ag$unique_id <- NULL
random_ag$new_plot <- paste("random_ag", random_ag$plot, sep = ":")
ss_ag$new_plot <- paste("SS_AG", ss_ag$plot, sep = ":")

mas <- rbind(ss_ag, random_ag)
rm(ss_ag, tmp, result, random_ag, c3ss, ses, can)
mas$plot <- NULL
colnames(mas)[8] <- "source"

con <- read.table(file = "shape_vs_no_shape_controls/Control_SHAPE_bpp_indexed_for_line_graphs.txt",
                  sep = "\t", header = TRUE)
con[,c(3,8)] <- NULL
colnames(con)[1] <- "gene_name"
con$new_index <- con$Position
con$new_plot <- paste("Control", con$source, sep = ":")

mas <- rbind(mas, con)
rm(con)


prof_mean <- mean(mas$adjusted_profile, na.rm = TRUE)
acc_mean <- mean(mas$av_sum_p, na.rm = TRUE)

stats <- mas %>%
  group_by(new_plot, new_index) %>%
  summarise(av_sum_p = mean(av_sum_p, na.rm = TRUE),
            sd_sum_p = sd(av_sum_p, na.rm = TRUE),
            av_adj_prof = mean(adjusted_profile, na.rm = TRUE),
            sd_adj_prof = sd(adjusted_profile, na.rm = TRUE))



filt <- c("random_ag:WT_C3SS_md200_to_225",
          "SS_AG:WT_Can_md200_to_225",
          "SS_AG:WT_C3SS_md200_to_225",
          "Control:long")

values <- c("random_ag:WT_C3SS_md200_to_225" = "#bae4b3",
            "SS_AG:WT_Can_md200_to_225" = "#67a9cf",
            "SS_AG:WT_C3SS_md200_to_225" = "#a6cee3",
            "Control:long" = "#969696")
long_acc <- ggplot(stats %>%
         filter(new_plot %in% filt), 
       aes(x = as.factor(new_index), y = 1 - av_sum_p, group = new_plot, col = new_plot)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(x = "Position", y = "Accessibility", color = "Plot") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0,1) +
  scale_color_manual(values = values) +
  geom_hline(yintercept = 0.493, linetype = "dashed", color = "black")

long_acc 

ggplot(stats %>%
         filter(new_plot %in% filt), 
       aes(x = as.factor(new_index), y = av_adj_prof, group = new_plot, col = new_plot)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(x = "Position", y = "SHAPE Reactivity", color = "Plot") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0,1) +
  scale_color_manual(values = values) +
  geom_hline(yintercept = 0.482, linetype = "dashed", color = "black")



short_acc
long_acc
