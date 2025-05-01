# For best method to assembly full-length transcripts
# transcript_id, conosorter classification, signalP identification, and emmapper data
# LOAD, eggnog_mapper.emapper.annotations.rds from long-orfs OR predicted peptide step
#        These DB includes gene_id (keyid), protein_id (ID.p[0-9]+$) and annotation fields (GOs, PFAm, COGs, etc)
# LOAD ConoSorter results (from .transdecoder.pep query), including fields as hydrophobicity, PHMM/RegExp and family_class ()
# LOAD conodictor (cross-check)
# LOAD, SignalP6 results, including protein_id (ID.p[0-9]+$)  and SP[sec/SPI] value ()
# LOAD Diamond blastp conosorter (include cols: ) 
# LOAD Diamond blastp Tox-prot
# LOAD sequence from Merged_polyA_hisat_SuperDuper.fasta.transdecoder.pep
# LOAD paste0(pub_dir, "/WGCNA.tsv")
# For conotoxin DB, generate a clustering-sequence label using DECIPHER::clustering
# JOIN in the follow order:
# 

library(tidyverse)

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)


pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR"


# 1) Emapper -----

eggNOG_cols <- c("gene_id","protein_id", "Preferred_name","Description", "GO","PFAMs","BRITE", "CAZy", "COG_category")

mapper_f <- list.files(path = pub_dir, pattern = "eggnog_mapper.emapper.annotations.rds",full.names = T)

MAPPERDB <- read_rds(mapper_f) %>% select_at(vars(contains(eggNOG_cols), starts_with("KEGG"))) 


# Conopeptide annotation
# 2.1) CONOSORTER -----

cnsrtr_f <- list.files(path = pub_dir, pattern = "ConoSorter_regex_pHMM.rds", full.names = T)

SORTERDB <- read_rds(cnsrtr_f)
# 2.2 CONODICTOR ------

cndctr_f <- file.path(pub_dir, "Conodictor2.rds")
DICTORDB <- read_rds(cndctr_f)

# 3) SIGNALp -----
signalp6_f <- list.files(pub_dir, pattern = "signalp6_Merged_polyA_hisat_SuperDuper_transdecoder.rds", full.names = T)

SIGNALPDB <- read_rds(signalp6_f)

# 4) Diamond (blastp conoServer and Tox-prot) -----
blastp_f <- list.files(pub_dir, pattern = "diamond_blastp_sources.rds", full.names = T)

# file.path(pub_dir, "diamond_blastp_sources.rds")

BLASTPDB <- read_rds(blastp_f)

# 5) Fasta file
# bind gene_id, dna_seq, protein_id, pep_seq

dna <- list.files(pub_dir, "Merged_polyA_hisat_SuperDuper.fasta$", full.names = T)

pep <- list.files(pub_dir, "Merged_polyA_hisat_SuperDuper.fasta.transdecoder.pep$", full.names = T)


# IOBUF_SIZE=200002

library(Rsamtools)
fa = FaFile(pep)
indexFa(fa,as = "AAStringSet")
(param = scanFaIndex(fa))

pep <- scanFa(fa, param=param, as = "AAStringSet")

# pep <- Biostrings::readAAStringSet(pep, format = "fasta")

head(pepdf <- data.frame(protein_id = names(pep), pep_seq = c(pep)) %>% as_tibble())

dna <- Biostrings::readDNAStringSet(dna)

gene_id <- sapply(strsplit(names(dna), " "), `[`, 1)

# dnadf must filtered to peptide sequences

head(dnadf <- data.frame(gene_id = gene_id, dna_seq = c(dna)) %>% as_tibble())

seqdf <- pepdf %>% 
  mutate(gene_id = gsub(".p[0-9]+$","", protein_id)) %>%
  left_join(dnadf) 


# Bind Mapper -----

any(seqdf$protein_id %in% MAPPERDB$protein_id)

DB1 <- seqdf %>% left_join(MAPPERDB) 

# DB1 %>% filter(gene_id %in% "Cluster-15813.136460") %>% view()

SORTERDB <- SORTERDB %>% select(-Method)

SIGNALPDB <- SIGNALPDB %>% select(-Method)

# Bind Conopeptide DB ====

# (Conodictor + conoSorter + blastp + signal )

# 1 BOTH             1162
# 2 ConoSorter       3622
# 3 Conodictor       7352
# 4 Total.           12136


PEPTIDESDB <- DICTORDB %>% full_join(SORTERDB)

PEPTIDESDB <- PEPTIDESDB %>% left_join(BLASTPDB)

# Sanity check
# ===== PEPTIDESDB
sum(PEPTIDESDB$protein_id %in% SORTERDB$protein_id)/nrow(SORTERDB)
# ===== ConoDictor
sum(PEPTIDESDB$protein_id %in% DICTORDB$protein_id)/nrow(DICTORDB)


# Add cross-check column

DF1 <- DICTORDB %>% distinct(protein_id) %>% mutate(prediction_tool = "Conodictor")
DF2 <- SORTERDB %>% ungroup() %>% distinct(protein_id) %>% mutate(prediction_tool = "ConoSorter")

which_tools <- function(x) { 
  x <- x[!is.na(x)] 
  n <- length(unique(x))
  x <- unique(x)
  
  if(n > 1) {
    x <- "BOTH"
  } else
    x <- paste(x, sep = '|', collapse = '|') }

CrossCheckdf <- rbind(DF1, DF2) %>%
  group_by(protein_id) %>%
  summarise(
    across(prediction_tool, .fns = which_tools), 
    .groups = "drop_last")

CrossCheckdf %>% count(prediction_tool)
  
PEPTIDESDB <- PEPTIDESDB %>% left_join(CrossCheckdf)

# Bind to signalp

sum(PEPTIDESDB$protein_id %in% SIGNALPDB$protein_id)

PEPTIDESDB <- SIGNALPDB %>% right_join(PEPTIDESDB)

nrow(PEPTIDESDB)

PEPTIDESDB %>% count(prediction_tool, Signalp_class)

# Bind emapper to PeptideDB ======

DB1 <- DB1 %>% left_join(PEPTIDESDB)


# Sanity check (1)

# ===== SignalP

sum(DB1$protein_id %in% SIGNALPDB$protein_id)/nrow(SIGNALPDB)

# ===== PEPTIDESDB

sum(DB1$protein_id %in% PEPTIDESDB$protein_id)/nrow(PEPTIDESDB)

# ===== BlastP to Toxprot and conoserver

sum(DB1$protein_id %in% BLASTPDB$protein_id)

# Outpts =====

write_rds(DB1, file = paste0(pub_dir, "/structured_db.rds"))


# AA

pepseqs <- PEPTIDESDB %>% 
  left_join(seqdf) %>%
  drop_na(prediction_tool) %>% 
  # filter(Signalp_class == "SP") %>%
  # select(protein_id, prediction_tool, Signalp_class, Superfamily, hmm_pred_conodictor, uniprotkb_toxprot, conoserver_protein) %>%
  mutate(uniprotkb_toxprot = sapply(strsplit(uniprotkb_toxprot, "[|]"), `[`, 2)) %>%
  mutate(uniprotkb_toxprot = ifelse(is.na(uniprotkb_toxprot), "Unknown", uniprotkb_toxprot)) %>%
  mutate(conoserver_protein = ifelse(is.na(conoserver_protein), "Unknown", conoserver_protein)) %>%
  select(pep_seq, protein_id, prediction_tool, Signalp_class, uniprotkb_toxprot, conoserver_protein) %>%
  unite("protein_id", protein_id:conoserver_protein, sep = "|") %>%
  mutate(pep_seq = gsub("[*]$", "", pep_seq)) %>%
  pull(pep_seq, name = protein_id) 
  
pepseqs <- Biostrings::AAStringSet(pepseqs)

Biostrings::writeXStringSet(pepseqs, file.path(pub_dir, "conopeptides.pep"))


# DNA

seqs <- PEPTIDESDB %>% 
  left_join(seqdf) %>%
  drop_na(prediction_tool) %>% 
  # filter(Signalp_class == "SP") %>%
  mutate(uniprotkb_toxprot = sapply(strsplit(uniprotkb_toxprot, "[|]"), `[`, 2)) %>%
  mutate(uniprotkb_toxprot = ifelse(is.na(uniprotkb_toxprot), "Unknown", uniprotkb_toxprot)) %>%
  mutate(conoserver_protein = ifelse(is.na(conoserver_protein), "Unknown", conoserver_protein)) %>%
  select(dna_seq, gene_id, prediction_tool, Signalp_class, uniprotkb_toxprot, conoserver_protein) %>%
  unite("gene_id", gene_id:conoserver_protein, sep = "|") %>%
  pull(dna_seq, name = gene_id)


seqs <- DB1 %>% 
  drop_na(tab) %>% 
  filter(Signalp_class == "SP") %>%
  distinct(dna_seq, gene_id) %>%
  # mutate(pep_seq = gsub("[*]$", "", pep_seq)) %>%
  pull(dna_seq, name = gene_id) 

seqs <- Biostrings::DNAStringSet(seqs)

Biostrings::writeXStringSet(seqs, file.path(pub_dir, "conopeptides.fasta"))

# ===== Prior to save peptide DB, create a clusters

clusters <- DECIPHER::Clusterize(pepseqs,
  cutoff=0.5, # < 50% distant
  minCoverage=0.5, # > 50% coverage
  processors=NULL) # use all CPUs

barplot(sort(table(clusters)))

clustersdf <- clusters %>% as_tibble(rownames = "protein_id")

clustersdf <- clustersdf %>% mutate(cluster = paste0("seqgroup_", cluster))

clustersdf <- clustersdf %>% mutate(protein_id =  sapply(strsplit(protein_id, "[|]"), `[`, 1))


# 
# clustersdf <- data.frame(pepseqs) %>% 
#   as_tibble(rownames = "protein_id") %>% 
#   left_join(clustersdf) %>%
#   select(-protein_id) %>%
#   dplyr::rename("pep_seq" = "pepseqs")

PEPTIDESDB %>% 
  left_join(seqdf) %>%
  left_join(clustersdf) %>% 
  write_tsv(paste0(pub_dir, "/conopeptides.tsv"))

PEPTIDESDB %>% count(Signalp_class, prediction_tool)


# Exit

# library(msa)
# 
# align <- msa::msa(seqs, method = "Muscle")
# 
# .align <- msaConvert(align)$seq
# 
# DECIPHER::BrowseSeqs(DNAStringSet(.align))