# Q: Why some identical cds rise from distinc transcript? (see Merged_polyA_hisat_SuperDuper.fasta)
# Q: Does putative conotoxin without SP signal are true non-conotoxin?

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)


library(tidyverse)
library(Rsamtools)

pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR/"

# LOAD data -----

dir <- "/Users/cigom/Documents/GitHub/conopeptides/06.Quantification/MATRIX_RSEM_dir"

subdirs <- list.files(dir, pattern = "KALLISTO_Merged_polyA_hisat_SuperDuper.fasta.transdecoder_DIR", full.names = T) 

f <- list.files(path = subdirs, pattern = "KALLISTO_Merged_polyA_hisat_SuperDuper.fasta.transdecoder.matrix", full.names = T)

kallisto_matrix <- read_rds(f)


DNA <- list.files(pub_dir, "Merged_polyA_hisat_SuperDuper.fasta$", full.names = T)

fa = FaFile(DNA)
indexFa(fa,as = "DNAStringSet")
(param = scanFaIndex(fa))

DNA <- scanFa(fa, param=param, as = "DNAStringSet")


DB <- read_tsv(paste0(pub_dir, "/conopeptides.tsv"))

DB %>% 
  count(Signalp_class, prediction_tool, pep_seq, sort = T) %>%
  ggplot(aes(y = prediction_tool, x = n, fill = Signalp_class)) + geom_col()

DB %>% 
  # if unique peptides, select filter & prediction_tool == "BOTH"
  filter(Signalp_class == "SP") %>%
  count(pep_seq, sort = T) %>% filter(n == 5) %>% head(1)

# Same peptide and cds 


query_pep <- "MFRLTAVCCFLLVIVQMDMILATSIPCTATGHPCGDFMWCCDAGDVCCDHLGPGFCMERRQCPLFGK"

DB %>% filter(grepl(query_pep, pep_seq)) %>% distinct(dna_seq) %>% pull()

query_cds <- DB %>% filter(grepl(query_pep, pep_seq)) %>% distinct(dna_seq) %>% pull()

query_cds <- structure(query_cds, names = paste0("cds_", seq(length(query_cds)))) 

# But different transcript

DB %>% 
  filter(grepl("MFRLTAVCCFLLVIVQMDMILATSIPCTATGHPCGDFMWCCDAGDVCCDHLGPGFCMERRQCPLFGK", pep_seq)) %>% view()
  # mutate(protein_id = gsub(".p[0-9]+$","", protein_id)) %>%
  distinct(protein_id) %>%
  pull() -> query

# paste(query, collapse = "|")

# Test if quantification is redundant 

kallisto_matrix[rownames(kallisto_matrix) %in% query,]


myXStringSet <- DNAStringSet(c(query_cds, as.character(DNA[names(DNA) %in% query])))


library(msa)

align <- msa::msa(myXStringSet, method = "ClustalW", order = "input")

.align <- msa::msaConvert(align)$seq

names(.align) <- msa::msaConvert(align)$nam

# data(BLOSUM62)
# msaConservationScore(align, BLOSUM62)

library(ggsci)

pat <- c("-", alphabet(myXStringSet, baseOnly=TRUE))

colors <- structure(pal_aaas()(length(pat)), names = rev(pat))

scales::show_col(colors)

DECIPHER::BrowseSeqs(DNAStringSet(.align), colWidth = 200, colors = colors, colorPatterns = T)



