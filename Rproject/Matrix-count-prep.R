
# READ count matrix Merged_polyA_hisat_SuperDuper_genes.matrix
# FILTER BY LOW-EXPRESION
# CALCULATE
# PCA and dendogram
# SAVE output for WGCNA and DE analysis

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

dir <- "/Users/cigom/Documents/GitHub/conopeptides/06.Quantification/MATRIX_RSEM_dir"

f <- list.files(dir, pattern = "Merged_polyA_hisat_SuperDuper_isoforms.matrix", full.names = T)

# count_f <- list.files(path = path, pattern = "Counts.txt", full.names = T)

# MTD_f <- list.files(path = path, pattern = "METADATA.tsv", full.names = T)


library(tidyverse)

COUNTS <- read_tsv(f)

colNames <- gsub("_L[0-9]|.genes|.isoforms.results", "", basename(names(COUNTS)))

colNames[1] <- "rowid"

colnames(COUNTS) <- c(colNames)

rowNames <- COUNTS$rowid

COUNTS <- COUNTS %>% select_if(is.double) %>% as(., "matrix")

rownames(COUNTS) <- rowNames

# 1) Filter data by removing low-abundance genes ----

by_count <- 1; by_freq <- 2

keep <- rowSums(COUNTS > by_count) >= by_freq

sum(keep)/nrow(COUNTS) # N transcripts

nrow(COUNTS <- COUNTS[keep,])

COUNTS <- round(COUNTS)

# Save

file_out <- gsub(".matrix", ".filt.rds",f)

write_rds(COUNTS, file = file_out)

# exit


library(DESeq2)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = COUNTS,
  colData = .colData,
  design = ~ 1 )

dds <- estimateSizeFactors(ddsFullCountTable) 

dds <- estimateDispersions(dds)

# However for other downstream analyses – e.g. for visualization or clustering – it might be useful to work with transformed versions of the count data. Ex:

vst <- DESeq2::varianceStabilizingTransformation(dds) # vst if cols > 10

# vst <- DESeq2::vst(dds) # vst if cols > 10

ntr <- DESeq2::normTransform(dds)

DESeq2::plotPCA(ntr, intgroup = "pH")
DESeq2::plotPCA(vst, intgroup = "pH")

raw_df <- vsn::meanSdPlot(assay(dds), plot = F)

vst_df <- vsn::meanSdPlot(assay(vst), plot = F)
ntr_df <- vsn::meanSdPlot(assay(ntr), plot = F)

rbind(data.frame(py = vst_df$sd, px = vst_df$rank, col = "vst"),
  data.frame(py = ntr_df$sd, px = ntr_df$rank, col = "ntr"),
  data.frame(py = raw_df$sd, px = raw_df$rank, col = "raw")) %>%
  filter(col != "raw") %>%
  ggplot(aes(px, py, color = col)) +
  labs(x = "Ranks", y = "sd", color = "") +
  geom_line(orientation = NA, position = position_identity(), size = 2) +
  theme_bw(base_family = "GillSans", base_size = 20) +
  theme(legend.position = "top")