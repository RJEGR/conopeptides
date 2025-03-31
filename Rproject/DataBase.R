# For best method to assembly full-length transcripts
# transcript_id, conosorter classification, signalP identification, and emmapper data
# LOAD, eggnog_mapper.emapper.annotations.rds from long-orfs OR predicted peptide step
#        These DB includes gene_id (keyid), protein_id (ID.p[0-9]+$) and annotation fields (GOs, PFAm, COGs, etc)
# LOAD ConoSorter results (from .transdecoder.pep query), including fields as hydrophobicity, PHMM/RegExp and family_class ()
# LOAD, SignalP6 results, including protein_id (ID.p[0-9]+$)  and SP[sec/SPI] value (yet running ...)
# LOAD sequence from Merged_polyA_hisat_SuperDuper.fasta.transdecoder.pep
# LOAD paste0(pub_dir, "/WGCNA.tsv")
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

  
# 2) CONOSORTER -----

cnsrtr_f <- list.files(path = pub_dir, pattern = "ConoSorter_regex_pHMM.rds", full.names = T)

SORTERDB <- read_rds(cnsrtr_f)

# 3) SIGNALp -----
signalp6_f <- list.files(pub_dir, pattern = "signalp6_Merged_polyA_hisat_SuperDuper_transdecoder.rds", full.names = T)

SIGNALPDB <- read_rds(signalp6_f)

# 4) Fasta file
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


# Bind -----

any(seqdf$protein_id %in% MAPPERDB$protein_id)

DB1 <- seqdf %>% left_join(MAPPERDB) # %>% drop_na()

# DB1 %>% count(gene_id, sort = T)

# DB1 %>% filter(gene_id %in% "Cluster-15813.136460") %>% view()

SORTERDB %>% count(Method)
SIGNALPDB  %>% count(Method)

SORTERDB <- SORTERDB %>% select(-Method)
SIGNALPDB <- SIGNALPDB %>% select(-Method)

# Scheck

sum(DB1$protein_id %in% SIGNALPDB$protein_id)


DB1 <- DB1 %>% left_join(SIGNALPDB)


any(DB1$protein_id %in% SORTERDB$protein_id)

DB1 <- DB1 %>% left_join(SORTERDB)

write_rds(DB1, file = paste0(pub_dir, "/structured_db.rds"))
