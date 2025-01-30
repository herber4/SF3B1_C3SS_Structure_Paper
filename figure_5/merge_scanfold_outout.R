library(tidyr)
library(stringr)
library(dplyr)
library(ggplot2)

process_profile_files("base_pairs/vienna_rnafold/pub/scanfold/resistant/long_shape/mono_out",
                      "Resistant", "r")

process_profile_files("base_pairs/vienna_rnafold/pub/scanfold/k700e/long_shape/mono_out",
                      "K700E", "k")
process_profile_files("base_pairs/vienna_rnafold/pub/scanfold/con/long_shape/mono_out",
                      "Control", "c")

mas <- rbind(r, c, k)

mas$plot <- paste(mas$i, mas$j, sep = "_")

mas %>%
  filter(i > 25,
         j < 900) %>%
  ggplot(aes(x = plot, y = Z_score, fill = source)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90))

write.table(mas, file = "base_pairs/vienna_rnafold/pub/scanfold/long_shape_scanfold_windowed_zscores.txt",
            sep = "\t", row.names = F, quote = F)
