library(rtracklayer)
library(dplyr)
library(ggplot2)
library(tidyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(GenomicRanges)
library(ggseqlogo)
library(stringr)

setwd("Figure_1/data/")

genome <- BSgenome.Hsapiens.UCSC.hg38

ref <- read.table(file = "Figure_3/data/gblock_meta_data.txt",
                  sep = "\t", header = TRUE)
ref$source <- "K700E"
ref <- ref %>%
  filter(!gene %in% c("RPRD1A", "BRD9", "COMT", "UBA1", "EHMT1"))


gtf <- readGFF("/Users/herber4/gencode_annos/gencode.v39.annotation.gtf")
colnames(gtf)[1] <- "chr"
gtf$ss <- ifelse(gtf$strand == "-", gtf$end, gtf$start)

utrs <- gtf %>%
  filter(type == "UTR" |
           exon_number == 1 | 
           tag == "cds_start_NF" |
           tag == "cds_end_NF" |
           tag == "readthrough_transcript" |
           transcript_type == "unprocessed_pseudogene" |
           transcript_type == "lncRNA")

df <- read.table(file = "Figure_1/data/K700E_rMATS.txt",
                 sep = "\t", header = TRUE)
k7 <- df
k7 <- k7 %>%
  filter(!sample %in% c("HepG2_KD", "K562_Deletion"))
k7$color <- "Not Significant"
k7$color[k7$FDR <= .09 & k7$IncLevelDifference > .05 & k7$median_IJC_SAMPLE_2 >= 10 & k7$median_SJC_SAMPLE_1 >= 10] <- "Sig"
k7$color[k7$FDR <= .09 & k7$IncLevelDifference < -.05 & k7$median_IJC_SAMPLE_2 >= 10 & k7$median_SJC_SAMPLE_1 >= 10] <- "Sig"
tmp2 <- k7
tmp2$splice_id <- apply(tmp2[, c(3:4,6:11,5)], 1, function(row) paste(row, collapse = ":"))
tmp2$splice_id <- gsub(" ", "", tmp2$splice_id)
k7 <- tmp2
k7 <- k7 %>%
  filter(type == "A3SS")

up <- k7 %>%
  filter(IncLevelDifference > .05)
# Find splice_ids that occur in all three samples

common_splice_ids <- up %>%
  group_by(splice_id) %>%
  summarise(sample_count = n_distinct(sample)) %>%
  filter(sample_count >= 2) %>%
  pull(splice_id)

# Filter the original data for those splice_ids
up <- k7 %>%
  filter(IncLevelDifference > .05,
         splice_id %in% common_splice_ids)
down <- k7 %>%
  filter(IncLevelDifference < -.05)

# Find splice_ids that occur in all three samples
common_splice_ids <- down %>%
  group_by(splice_id) %>%
  summarise(sample_count = n_distinct(sample)) %>%
  filter(sample_count >= 2) %>%
  pull(splice_id)

# Filter the original data for those splice_ids
down <- k7 %>%
  filter(IncLevelDifference < -.05,
         splice_id %in% common_splice_ids)
up$direction <- "Positive"
down$direction <- "Negative"

mas <- rbind(up, down)
mas$c3ss_length <- ifelse(mas$strand == "-", mas$longExonEnd-mas$shortEE, mas$shortES-mas$longExonStart_0base)
rm(k7, df, down, up, tmp2, common_splice_ids)

mas <- mas %>%
  filter(c3ss_length < 500)

mas$sjc_sum <- mas$median_SJC_SAMPLE_1+mas$median_SJC_SAMPLE_2
mas$ijc_sum <- mas$median_IJC_SAMPLE_1+mas$median_IJC_SAMPLE_2

mas$long_inclusion_ratio <- mas$ijc_sum/(mas$ijc_sum+mas$sjc_sum)
mas$c3ss <- ifelse(mas$strand == "-", mas$longExonEnd, mas$longExonStart_0base)
mas$can_ss <- ifelse(mas$strand == "-", mas$shortEE, mas$shortES)
mas$intron_length <- ifelse(mas$strand == "-", mas$flankingES-mas$longExonEnd, mas$longExonStart_0base-mas$flankingEE)
mas$c3ss <- ifelse(mas$strand == "+", mas$c3ss+1, mas$c3ss)
mas$can_ss <- ifelse(mas$strand == "+", mas$can_ss+1, mas$can_ss)
mas$five_ss <- ifelse(mas$strand == "-", mas$flankingES, mas$flankingEE)
mas$five_ss <- ifelse(mas$strand == "+", mas$five_ss+1, mas$five_ss)

dedup <- mas
colnames(dedup)[3] <- "gene_name"
dedup$longExonStart_0base <- dedup$longExonStart_0base+1
dedup$shortES <- dedup$shortES+1

short <- dedup
colnames(short)[8:9] <- c("start", "end")
short <- anti_join(short, utrs, by = c("gene_name","chr", "start", "end","strand"))
colnames(short)[8:9] <- c("shortES", "shortEE")
long <- short
colnames(long)[6:7] <- c("start", "end")
long <- anti_join(long, utrs, by = c("gene_name","chr", "start", "end","strand"))

coords <- long[,c(4,5,31,33,37:40)]
coords <- unique(coords)

write.table(coords, file = "two_or_more/filter_junk/rmats_coords_for_dedup.txt",
            sep = "\t", row.names = F, quote = F)



######
i went through by hand and dedupped frauds from this rmats_coords_for_dedup.txt
######
df <- read.table(file = "two_or_more/filter_junk/rmats_coords_for_dedup.txt",
                 sep = "\t", header = TRUE)

df$bad <- ifelse(df$c3ss_length > df$intron_length, "bad", "good")
df <- df %>%
  filter(bad == "good")

mean(df$c3ss_length)
mean(df$intron_length)
coords <- df



c3ss <- GRanges(
  seqnames = coords$chr,
  ranges = IRanges(start = ifelse(coords$strand == "-", coords$c3ss, coords$c3ss-150),
                   end = ifelse(coords$strand == "-", coords$c3ss+150, coords$c3ss)),
  strand = coords$strand
)
mcols(c3ss) <- coords$splice_id
seqs <- getSeq(genome, c3ss)
mcols(seqs) <- mcols(c3ss)
names(seqs) <- c3ss@elementMetadata@listData[["X"]]

lc_gc_content <- sapply(seqs, function(seq) {
  gc_count <- sum(letterFrequency(DNAString(seq), letters = c("G", "C")))
  gc_percent <- gc_count / length(seq) * 100
  gc_percent
})


lc_gc_content <- as.data.frame(cbind(lc_gc_content, length = sapply(seqs, length),
                                     splice_id = mcols(seqs)$X))
rownames(lc_gc_content) <- NULL

meta <- merge(coords, lc_gc_content, by = "splice_id")
colnames(meta)[10] <- "GC_content"
meta$GC_content <- as.numeric(meta$GC_content)
ggplot(meta, aes(x = GC_content)) +
  geom_histogram()

write.table(meta, file = "two_or_more/filter_junk/K700E_rMATS_GC_Content_Meta_Data.txt",
            sep = "\t", row.names = F, quote = F)

c3ss <- GRanges(
  seqnames = coords$chr,
  ranges = IRanges(start = coords$c3ss-450, end = coords$c3ss+450),
  strand = coords$strand
)
mcols(c3ss) <- coords$splice_id
seqs <- getSeq(genome, c3ss)
mcols(seqs) <- mcols(c3ss)
names(seqs) <- c3ss@elementMetadata@listData[["X"]]
writeXStringSet(seqs, filepath = "two_or_more/filter_junk/rMATS_K700E_shared_C3SS900.fasta",
                format = "fasta")

ss <- GRanges(
  seqnames = coords$chr,
  ranges = IRanges(start = coords$can_ss-450, end = coords$can_ss+450),
  strand = coords$strand
)
mcols(ss) <- coords$splice_id
seqs <- getSeq(genome, ss)
mcols(seqs) <- mcols(ss)
names(seqs) <- ss@elementMetadata@listData[["X"]]

writeXStringSet(seqs, filepath = "two_or_more/filter_junk/rMATS_K700E_shared_C3SS_Canonical900.fasta",
                format = "fasta")

five <- GRanges(
  seqnames = coords$chr,
  ranges = IRanges(start = coords$five_ss-450, end = coords$five_ss+450),
  strand = coords$strand
)
mcols(five) <- coords$splice_id
seqs <- getSeq(genome, five)
mcols(seqs) <- mcols(five)
names(seqs) <- five@elementMetadata@listData[["X"]]

writeXStringSet(seqs, filepath = "two_or_more/filter_junk/rMATS_K700E_5SS_MAXENT.fasta",
                format = "fasta")




allthreeshared_tmp <- semi_join(coords, all_three_shared, by = "splice_id")
write.table(red, file = "manuscript_v1.0/Table_Two_rMATS.txt",
            sep = "\t", row.names = F, quote = F)
write.table(tmp, file = "manuscript_v1.0/Table_Two_Reduced.txt",
            sep = "\t", row.names = F, quote = F)
