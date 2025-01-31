library(rtracklayer)
library(dplyr)
library(ggplot2)
library(tidyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(GenomicRanges)
library(ggseqlogo)
library(stringr)

genome <- BSgenome.Hsapiens.UCSC.hg38

df1 <- read.table(file = "/WT1_junc_annotated.txt",
                 sep = "\t", header = TRUE)
df2 <- read.table(file = "/WT2_junc_annotated.txt",
                  sep = "\t", header = TRUE)
df3 <- read.table(file = "/WT3_junc_annotated.txt",
                  sep = "\t", header = TRUE)

df1$sample <- "WT1"
df2$sample <- "WT2"
df3$sample <- "WT3"
df <- rbind(df1, df2, df3)

df <- df %>%
  filter(!str_starts(chrom, "K") & !str_starts(chrom, "G"),
         exons_skipped == 0)
colnames(df)[1] <- "chr"
df$ss <- ifelse(df$strand == "-", df$start, df$end)
scores <- df %>%
  group_by(chr, ss, strand) %>%
  summarise(mean = mean(score),
            sd = sd(score))
scores <- scores %>%
  filter(mean > 10)
df[,c(4:5, 8:14, 18)] <- NULL
df <- unique(df)
df <- merge(df, scores, by = c("chr", "ss", "strand"))
rm(df1, df2, df3, scores)

gtf <- readGFF("/Users/herber4/gencode_annos/gencode.v39.annotation.gtf")
colnames(gtf)[1] <- "chr"
gtf$ss <- ifelse(gtf$strand == "-", gtf$end, gtf$start)

utrs <- gtf %>%
  filter(type == "UTR" |
           transcript_type == "retained_intron" |
           exon_number == 1 | 
           tag == "cds_start_NF" |
           tag == "cds_end_NF" |
           tag == "readthrough_transcript" |
           transcript_type == "unprocessed_pseudogene" |
           transcript_type == "lncRNA" |
           transcript_type == "processed_transcript")

test <- anti_join(df, utrs, by = c("chr", "strand", "ss"))

vdb <- read.table(file = "master_scripts/Nalm6_splice_junctions/VastDB_3SS_Annotations_Reduced.txt",
                  sep = "\t", header = TRUE)
vdb$ss <- ifelse(vdb$strand == "-", vdb$start, vdb$end)

test <- merge(df, vdb, by = c("chr","strand", "start", "end"))
#test[,c(8:14)] <- NULL
test <- test %>%
  filter(mean > 10)
pos <- test %>%
  filter(strand == "+") %>%
  group_by(chr, strand, start) %>%
  count(start, name = "counts")
pos <- pos %>%
  filter(counts > 1)
pos <- semi_join(test, pos, by = c("chr", "strand", "start"))

neg <- test %>%
  filter(strand == "-") %>%
  group_by(chr, strand, end) %>%
  count(end, name = "counts")
neg <- neg %>%
  filter(counts > 1)
neg <- semi_join(test, neg, by = c("chr", "strand", "end"))
test <- rbind(pos, neg)
rm(tmp, pos, neg)
test$ss.x <- NULL
colnames(test)[17] <- "ss"

zeros <- unique(test)

zeros[,c(5:8,13,15:16)] <- NULL
mas <- zeros %>%
  filter(LE_o > 0)
rm(test, gtf, utrs, vdb, df)

can <- zeros %>%
  group_by(group = ifelse(strand == "+", start, end)) %>%
  filter(mean == max(mean)) %>%
  ungroup()

can$source <- "Canonical"


#zeros <- zeros %>% group_by(group = ifelse(strand == "+", start, end))
tmp <- anti_join(zeros, can, by = "ss")

pos <- tmp %>%
  filter(strand == "+")
pos[,c(4, 7:9)] <- NULL
colnames(pos)[4:5] <- c("c3ss_mean", "c3ss_sd")
pos <- merge(can, pos, by = c("chr", "strand", "start"))
pos$c3ss_length <- abs(pos$ss.x-pos$ss.y)
colnames(pos)[15] <- "c3ss"
colnames(pos)[10] <- "ss"

min <- tmp %>%
  filter(strand == "-")
min[,c(3, 7:9)] <- NULL
colnames(min)[4:5] <- c("c3ss_mean", "c3ss_sd")
min <- merge(can, min, by = c("chr", "strand", "end"))
min$c3ss_length <- abs(min$ss.x-min$ss.y)
colnames(min)[15] <- "c3ss"
colnames(min)[10] <- "ss"

mas <- rbind(pos, min)

mas %>%
  filter(c3ss_length < 1000) %>%
  ggplot(aes(x = c3ss_length)) +
  geom_density()
rm(min, pos, tmp, can)

k7 <- read.table(file = "/gblock_meta_data.txt",
                 sep = "\t", header = TRUE)
colnames(mas)[7] <- "gene"
mas <- anti_join(mas, k7, by = "gene")
mas$LE_o <- NULL
mas$c3ss_length <- mas$ss-mas$c3ss
mas <- mas %>%
  mutate(direction = case_when(
    strand == "+" & c3ss_length > 0 ~ "upstream",
    strand == "+" & c3ss_length < 0 ~ "downstream",
    strand == "-" & c3ss_length > 0 ~ "downstream",
    strand == "-" & c3ss_length < 0 ~ "upstream"
  ))
mas$c3ss_length <- mas$c3ss_length*-1
mas$id <- paste(mas$gene, mas$chr, sep = ":")
mas$id <- paste(mas$id, mas$ss, sep = ":")
mas$id <- paste(mas$id, mas$c3ss, sep = ":")
{gr <- GRanges(seqnames = mas$chr,
              ranges = IRanges(start = mas$ss-450, end = mas$ss+450),
              strand = mas$strand)

mcols(gr) <- mas$id
seqs <- getSeq(genome, gr)
mcols(seqs) <- mcols(gr)
seqs@ranges@NAMES <- seqs@elementMetadata@listData[["X"]]
writeXStringSet(seqs, filepath = "/Nalm6_Regtools_Canonical_Junctions.fasta",
                format = "fasta")}

{gr <- GRanges(seqnames = mas$chr,
               ranges = IRanges(start = mas$c3ss-450, end = mas$c3ss+450),
               strand = mas$strand)
  
  mcols(gr) <- mas$id
  seqs <- getSeq(genome, gr)
  mcols(seqs) <- mcols(gr)
  seqs@ranges@NAMES <- seqs@elementMetadata@listData[["X"]]
  writeXStringSet(seqs, filepath = "/Nalm6_Regtools_C3SS_Junctions.fasta",
                  format = "fasta")}

gc <- GRanges(seqnames = mas$chr,
              ranges = IRanges(start = ifelse(mas$strand == "-", mas$c3ss, mas$c3ss-150),
                               end = ifelse(mas$strand == "-", mas$c3ss+150, mas$c3ss)),
              strand = mas$strand)
mcols(gc) <- mas[,c(17)]
gc_seqs <- getSeq(genome, gc)
mcols(gc_seqs) <- mcols(gc)
names(gc_seqs) <- gc_seqs@elementMetadata@listData[["X"]]
seqs <- gc_seqs
lc_gc_content <- sapply(seqs, function(seq) {
  gc_count <- sum(letterFrequency(DNAString(seq), letters = c("G", "C")))
  gc_percent <- gc_count / length(seq) * 100
  gc_percent
})

lc_gc_content <- as.data.frame(cbind(lc_gc_content, length = sapply(seqs, length),
                                     id = mcols(seqs)$X))
rownames(lc_gc_content) <- NULL
mas <- merge(mas, lc_gc_content, by = "id")
mas$length <- NULL
mas$intron_length <- abs(mas$start-mas$end)

write.table(mas, file = "/Nalm6_WT_C3SS_annotations.regtools.txt",
            sep = "\t", row.names = F, quote = F)

gc <- GRanges(seqnames = mas$chr,
              ranges = IRanges(start = mas$c3ss-100, end = mas$c3ss+100),
              strand = mas$strand)
mcols(gc) <- mas[,c(1)]
gc_seqs <- getSeq(genome, gc)
mcols(gc_seqs) <- mcols(gc)
names(gc_seqs) <- gc_seqs@elementMetadata@listData[["X"]]
seqs <- gc_seqs

pfm <- consensusMatrix(seqs, as.prob = TRUE)

ggplot() +
  geom_logo(data = pfm, method = "probability", col_scheme = "nucleotide") +
  theme_classic()

# get all 3' SS junctions

df1 <- read.table(file = "/WT1_junc_annotated.txt",
                  sep = "\t", header = TRUE)
df2 <- read.table(file = "/WT2_junc_annotated.txt",
                  sep = "\t", header = TRUE)
df3 <- read.table(file = "/WT3_junc_annotated.txt",
                  sep = "\t", header = TRUE)

df1$sample <- "WT1"
df2$sample <- "WT2"
df3$sample <- "WT3"
df <- rbind(df1, df2, df3)

df <- df %>%
  filter(!str_starts(chrom, "K") & !str_starts(chrom, "G"),
         exons_skipped == 0)
colnames(df)[1] <- "chr"
df$ss <- ifelse(df$strand == "-", df$start, df$end)
scores <- df %>%
  group_by(chr, ss, strand) %>%
  summarise(mean = mean(score),
            sd = sd(score))
scores <- scores %>%
  filter(mean > 10)
df <- df %>%
  filter(anchor == "DA",
         acceptors_skipped == 0,
         exons_skipped == 0,
         donors_skipped == 0,
         known_donor == 1,
         known_acceptor == 1,
         known_junction == 1)
df[,c(4:5, 8:14, 18)] <- NULL
df <- unique(df)
df <- merge(df, scores, by = c("chr", "ss", "strand"))
rm(df1, df2, df3, scores)

gtf <- readGFF("/Users/herber4/gencode_annos/gencode.v39.annotation.gtf")
colnames(gtf)[1] <- "chr"
gtf$ss <- ifelse(gtf$strand == "-", gtf$end, gtf$start)

utrs <- gtf %>%
  filter(type == "UTR" |
           transcript_type == "retained_intron" |
           exon_number == 1 | 
           tag == "cds_start_NF" |
           tag == "cds_end_NF" |
           tag == "readthrough_transcript" |
           transcript_type == "unprocessed_pseudogene" |
           transcript_type == "lncRNA" |
           transcript_type == "processed_transcript")

test <- anti_join(df, utrs, by = c("chr", "strand", "ss"))


c3ss <- read.table(file = "/Nalm6_WT_C3SS_annotations.regtools.txt",
                   sep = "\t", header = TRUE)
test <- anti_join(test, c3ss, by = c("strand", "chr", "start", "end"))
test$c3ss <- test$ss
test <- anti_join(test, c3ss, by = "c3ss")

tmp <- semi_join(test, gtf, by = "ss")

tmp$id <- paste(tmp$gene_names, tmp$start, sep = ":")
tmp$id <- paste(tmp$id, tmp$end, sep = ":")

gr <- GRanges(seqnames = tmp$chr,
              ranges = IRanges(start = tmp$ss-450, end = tmp$ss+450),
              strand = tmp$strand)
mcols(gr) <- tmp$id
seqs <- getSeq(genome, gr)
mcols(seqs) <- mcols(gr)
names(seqs) <- seqs@elementMetadata@listData[["X"]]

writeXStringSet(seqs, filepath = "/Nalm6_WT_3SS_regtools.fasta",
                format = "fasta")

gc <- GRanges(seqnames = tmp$chr,
              ranges = IRanges(start = ifelse(tmp$strand == "-", tmp$ss, tmp$ss-150),
                               end = ifelse(tmp$strand == "-", tmp$ss+150, tmp$ss)),
              strand = tmp$strand)
mcols(gc) <- tmp$id
gc_seqs <- getSeq(genome, gc)
mcols(gc_seqs) <- mcols(gc)
names(gc_seqs) <- gc_seqs@elementMetadata@listData[["X"]]
seqs <- gc_seqs
lc_gc_content <- sapply(seqs, function(seq) {
  gc_count <- sum(letterFrequency(DNAString(seq), letters = c("G", "C")))
  gc_percent <- gc_count / length(seq) * 100
  gc_percent
})

lc_gc_content <- as.data.frame(cbind(lc_gc_content, length = sapply(seqs, length),
                                     id = mcols(seqs)$X))
rownames(lc_gc_content) <- NULL
mas <- merge(tmp, lc_gc_content, by = "id")
mas$length <- NULL
mas$intron_length <- abs(mas$start-mas$end)

write.table(mas, file = "/Nalm6_WT_3SS_Canonical_annotations.regtools.txt",
            sep = "\t", row.names = F, quote = F)
