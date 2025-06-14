# Create a single databaset set of Diamond blastp results
# In this case  Diamond (blastp conoServer and Tox-prot)
# but it can extend to any outfmt6 result

# Due to multiple hits per sequence, BLAST output will be drive individually from the DB

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

read_outfmt6 <- function(f) {
  
  # seqid = transcript_id

  # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp
  
    outfmt6.names <- c("protein_id", "subject", "identity", "coverage", "mismatches", "gaps", "seq_start", "seq_end", "sub_start", "sub_end", "e", "bitscore", "qcovhsp","scovhsp")
  
  
  df <- read_tsv(f, col_names = F) 
  
  colnames(df) <- outfmt6.names
  
  df <- df %>% dplyr::mutate(FileName = basename(f))
  
  return(df)
  
  
}

dir <- "/Users/cigom/Documents/GitHub/conopeptides/05.Prediction/Merged_polyA_hisat_SuperDuper.transdecoder_dir/outfmt6_dir/"

pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR"

f <- list.files(dir, pattern = "outfmt6", full.names = T) 

DB <- do.call(rbind, lapply(f, read_outfmt6))

# as hits retriveved are 3, lets select the best score
# ex

DB %>% filter(protein_id %in% "Cluster-15813.148410.p2") %>% view()

# query: MKLTC-VLIVAVLILTACQFTAADDMEYPKWLRGLSTDX-SERGCWLCLGPNACCRG-SVCHD-YCPR
# target: MKLT-GVLIVAVLILTACQFTAADDMEYPKWLRGLSTD-KSERGCWLCLGPNACCRG-DVCH-SYCPR

# How to keep accurate name for blast mapping?

DB <- DB %>% 
  group_by(protein_id, FileName) %>%
  arrange(bitscore, .by_group = T) %>%
  # group_by(sampleB, y_axis) %>%
  filter(bitscore == max(bitscore)) %>%
  slice_head(n = 1)

DB <- DB %>% 
  select(protein_id, subject, FileName) %>%
  mutate(FileName = gsub("Merged_polyA_hisat_SuperDuper_longest_orfs_vs_", "",FileName)) %>%
  mutate(FileName = gsub(".diamond.blastp.outfmt6", "",FileName)) %>%
  pivot_wider(names_from = FileName, values_from = subject)


names(DB) <- c("protein_id", "uniprotkb_toxprot", "conoserver_protein")

write_rds(DB, file = file.path(pub_dir, "diamond_blastp_sources.rds"))
