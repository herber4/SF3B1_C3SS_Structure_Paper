library(dplyr)
library(ggplot2)
{set.seed(12345)
setwd("/Users/herber4/Desktop/Figures/Figure_2/data/")

rmats <- read.table(file = "/Users/herber4/Desktop/Figures/Figure_1/data/K700E_rMATS_GC_Content_Meta_Data.txt",
                    sep = "\t", header = TRUE)
rmats$source <- "rmats_k700e"
rmats <- rmats %>%
  mutate(gene = sub(":.*", "", splice_id))
rmats <- rmats %>%
  filter(c3ss_length < 500,
         !gene %in% c("RPRD1A", "BRD9", "COMT", "UBA1", "EHMT1"))
r <- rmats %>%
  summarise(mean_length = mean(intron_length),
            median_length = median(intron_length),
            mean_gc_content = mean(GC_content),
            median_gc_content = median(GC_content),
            median_c3ss = median(c3ss_length),
            mean_c3ss = mean(c3ss_length))
r$source <- "rmats_k700e"
r$bootstrap_number <- "0"
can <- read.table(file = "Nalm6_WT_3SS_Canonical_annotations.regtools.txt",
                  sep = "\t", header = TRUE)
ref <- read.table(file = "gblock_meta_data.txt",
                  sep = "\t", header = TRUE)
ref$source <- "K700E"
ref <- ref %>%
  filter(!gene %in% c("RPRD1A", "BRD9", "COMT", "UBA1", "EHMT1"))

ref <- ref %>%
  summarise(mean_length = mean(length),
            median_length = median(length),
            mean_gc_content = mean(GC_content),
            median_gc_content = median(GC_content),
            mean_c3ss = mean(c3ss_length),
            median_c3ss = median(c3ss_length))
ref$source <- "K700E"
ref$bootstrap_number <- "0"
ref <- rbind(ref, r)
c3ss <- read.table(file = "Nalm6_WT_C3SS_annotations.regtools.txt",
                   sep = "\t", header = TRUE)
c3ss <- c3ss %>%
  filter(c3ss_length < 500,
         c3ss_length > -500)

library(dplyr)
{colnames(can)[14] <- "gc_content"
# Set the number of bootstrap iterations
num_bootstraps <- 1000

# Create an empty dataframe to store the results
bootstrap_results <- data.frame(
  bootstrap_number = integer(),
  mean_gc_content = numeric(),
  median_gc_content = numeric(),
  mean_length = numeric(),
  median_length = numeric()
)

# Perform bootstrapping
for (i in 1:num_bootstraps) {
  # Sample 256 rows with replacement
  sampled_data <- can[sample(nrow(can), 123, replace = TRUE), ]
  
  # Calculate mean and median of gc_content and length columns
  mean_gc <- mean(sampled_data$gc_content)
  median_gc <- median(sampled_data$gc_content)
  mean_length <- mean(sampled_data$intron_length)
  median_length <- median(sampled_data$intron_length)
  
  # Store results in the dataframe
  bootstrap_results <- bind_rows(
    bootstrap_results,
    data.frame(
      bootstrap_number = i,
      mean_gc_content = mean_gc,
      median_gc_content = median_gc,
      mean_length = mean_length,
      median_length = median_length
    )
  )
}

# View the first few rows of the bootstrap results
head(bootstrap_results)
bootstrap_results$source <- "Nalm_6_Can"}
can <- bootstrap_results
#write.table(bootstrap_results, file = "Nal.txt",
#            sep = "\t", row.names = F, quote = F)


colnames(c3ss)[18] <- "gc_content"
c3ss$c3ss_length <- abs(c3ss$c3ss_length)
num_bootstraps <- 1000

# Create an empty dataframe to store the results
bootstrap_results <- data.frame(
  bootstrap_number = integer(),
  mean_gc_content = numeric(),
  median_gc_content = numeric(),
  mean_length = numeric(),
  median_length = numeric(),
  mean_c3ss = numeric(),
  median_c3ss = numeric()
)

# Perform bootstrapping
for (i in 1:num_bootstraps) {
  # Sample 256 rows with replacement
  sampled_data <- c3ss[sample(nrow(c3ss), 123, replace = TRUE), ]
  
  # Calculate mean and median of gc_content and length columns
  mean_gc <- mean(sampled_data$gc_content)
  median_gc <- median(sampled_data$gc_content)
  mean_length <- mean(sampled_data$intron_length)
  median_length <- median(sampled_data$intron_length)
  mean_c3ss <- mean(sampled_data$c3ss_length)
  median_c3ss <- median(sampled_data$c3ss_length)
  # Store results in the dataframe
  bootstrap_results <- bind_rows(
    bootstrap_results,
    data.frame(
      bootstrap_number = i,
      mean_gc_content = mean_gc,
      median_gc_content = median_gc,
      mean_length = mean_length,
      median_length = median_length,
      mean_c3ss = mean_c3ss,
      median_c3ss = median_c3ss
    )
  )
}

# View the first few rows of the bootstrap results
head(bootstrap_results)
bootstrap_results$source <- "Nalm6_C3SS"
c3ss <- bootstrap_results

can$median_c3ss <- "0"
can$mean_c3ss <- "0"

int <- ref$mean_length[ref$source == "K700E"]
gc <- ref$mean_gc_content[ref$source == "K700E"]
ss <- ref$mean_c3ss[ref$source == "K700E"]
rmats_int <- ref$mean_length[ref$source == "rmats_k700e"]
rmats_gc <- ref$mean_gc_content[ref$source == "rmats_k700e"]
rmats_ss <- ref$mean_c3ss[ref$source == "rmats_k700e"]

tmp <- rbind(can, c3ss, ref)
level <- c("rmats_k700e", "Nalm6_C3SS", "Nalm_6_Can")

pdf(file = "/Users/herber4/Desktop/Figures/Figure_1/K700E_GC_Content.pdf",
    width = 8, height = 10, paper = "letter")
tmp %>%
  filter(!source == "K700E",
         !source == "rmats_k700e") %>%
  ggplot(aes(x = source, y = mean_gc_content)) +
  geom_boxplot() +
  geom_point(data = ref %>% filter(!source == "K700E"), aes(x = source, y = mean_gc_content),
             col = "black", size = 4) +
  theme_bw() +
  labs(x = "Source", y = "GC Content (%)") +
  theme(legend.position = "") +
  scale_x_discrete(limits = level)
dev.off()

ggsave(filename = "../GC_content_vs_bootstraps.pdf",
       dpi = 1000, width = 8, height = 10, device = "pdf", plot = g)

mean(c3ss$mean_gc_content >= gc)
mean(can$mean_gc_content >= gc)
mean(can$mean_gc_content >= rmats_gc)
mean(c3ss$mean_gc_content >= rmats_gc)
mean(can$mean_gc_content <= rmats_gc)
mean(c3ss$mean_gc_content <= rmats_gc)
mean(c3ss$mean_gc_content <= gc)
mean(can$mean_gc_content <= gc)
mean(rmats$GC_content <= gc)
mean(rmats$GC_content >= gc)



ggplot(tmp, aes(x = source, y = mean_length, fill = source)) +
  geom_boxplot() +
  geom_hline(yintercept = int, linetype = "dashed") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Source", y = "Intron Length (bp)")


level <- c("rmats_k700e", "Nalm_6_Can", "Nalm6_C3SS")
pdf(file = "/Users/herber4/Desktop/Figures/Figure_1/K700E_intron_length.pdf",
    width = 8, height = 10, paper = "letter")
tmp %>%
  filter(!source == "K700E",
         !source == "rmats_k700e") %>%
  ggplot(aes(x = source, y = mean_length)) +
  geom_boxplot() +
  geom_point(data = ref %>% filter(!source == "K700E"), aes(x = source, y = mean_length),
             col = "black", size = 4) +
  theme_bw() +
  labs(x = "Source", y = "Intron Length (bp)") +
  theme(legend.position = "") +
  scale_x_discrete(limits = level)
dev.off()
mean(c3ss$mean_length >= int)
mean(can$mean_length >= int)
mean(can$mean_length >= rmats_int)
mean(c3ss$mean_length >= rmats_int)
mean(can$mean_length <= rmats_int)
mean(c3ss$mean_length <= rmats_int)
mean(c3ss$mean_length <= int)
mean(can$mean_length <= int)
mean(rmats$intron_length <= int)
mean(rmats$intron_length >= int)

level <- c("rmats_k700e", "Nalm6_C3SS")
pdf(file = "/Users/herber4/Desktop/Figures/Figure_1/C3SS_length.pdf",
    width = 8, height = 10, paper = "letter")
ggplot(c3ss, aes(x = source, y = mean_c3ss)) +
  geom_boxplot() +
  theme_bw() +
  geom_point(data = ref %>% filter(!source == "K700E"), aes(x = source, y = mean_c3ss), col = "black", size = 4) +
  labs(x = "Source",
       y = "C3SS Length (bp)") +
  theme(legend.position = "none") +
  scale_x_discrete(limits = level)
dev.off()
d
mean(c3ss$mean_c3ss >= rmats_ss)
mean(c3ss$mean_c3ss <= rmats_ss)
mean(rmats$c3ss_length <= ss)
mean(rmats$c3ss_length >= ss)

var.test(rmats$c3ss_length, ref$c3ss_length)
t.test(rmats$c3ss_length, ref$c3ss_length, var.equal = F)}



library(dplyr)
library(ggplot2)
set.seed(12345)
setwd("/Users/herber4/Desktop/Figures/Figure_2/data/")

rmats <- read.table(file = "/Users/herber4/Desktop/Figures/Figure_1/data/K700E_rMATS_GC_Content_Meta_Data.txt",
                    sep = "\t", header = TRUE)
rmats$source <- "rmats_k700e"
rmats <- rmats %>%
  mutate(gene = sub(":.*", "", splice_id))
rmats <- rmats %>%
  filter(c3ss_length < 500,
         !gene %in% c("RPRD1A", "BRD9", "COMT", "UBA1", "EHMT1"))
r <- rmats %>%
  summarise(mean_length = mean(intron_length),
            median_length = median(intron_length),
            mean_gc_content = mean(GC_content),
            median_gc_content = median(GC_content),
            median_c3ss = median(c3ss_length),
            mean_c3ss = mean(c3ss_length))
r$source <- "rmats_k700e"
r$bootstrap_number <- "0"
can <- read.table(file = "Nalm6_WT_3SS_Canonical_annotations.regtools.txt",
                  sep = "\t", header = TRUE)
ref <- read.table(file = "gblock_meta_data.txt",
                  sep = "\t", header = TRUE)
ref$source <- "K700E"
ref <- ref %>%
  filter(!gene %in% c("RPRD1A", "BRD9", "COMT", "UBA1", "EHMT1"))

ref <- ref[,c(2:4, 15:16)]
rmats <- rmats[,c(4,5,8,10:11)]
colnames(ref)[1] <- "intron_length"

mas <- rbind(rmats, ref)

shapiro.test(rmats$GC_content)
shapiro.test(ref$GC_content)
wilcox.test(rmats$GC_content, ref$GC_content)

pdf(file = "/Users/herber4/Desktop/Figures/Supplementary/Supp_Figure_1/C3SS_length_boxie.pdf",
    width = 8, height = 10, paper = "letter")
ggplot(mas, aes(x = source, y = c3ss_length, fill = source)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Source", 
       y = "C3SS Length (bp)")+
  theme(legend.position = "none")
dev.off()

wilcox.test(ref$c3ss_length, rmats$c3ss_length)
pdf(file = "/Users/herber4/Desktop/Figures/Supplementary/Supp_Figure_1/C3SS_intron_length_boxie.pdf",
    width = 8, height = 10, paper = "letter")
ggplot(mas, aes(x = source, y = intron_length, fill = source)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Source", 
       y = "Intron Length (bp)") +
  theme(legend.position = "none")
dev.off()
shapiro.test(rmats$intron_length)
wilcox.test(rmats$intron_length, ref$intron_length)


pdf(file = "/Users/herber4/Desktop/Figures/Supplementary/Supp_Figure_1/C3SS_GC_Content_boxie.pdf",
    width = 8, height = 10, paper = "letter")
ggplot(mas, aes(x = source, y = GC_content, fill = source)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Source", 
       y = "GC Content")+
  theme(legend.position = "none")
dev.off()
wilcox.test(rmats$GC_content, ref$GC_content)
