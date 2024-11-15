# /Users/cigom/Documents/GitHub/conopeptides/05.Prediction
# Merge blast from conoserver 

# filter fasta by toxin, and calculate Nx
# LOAD DNA AND ORF

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)


read_outfmt6 <- function(f) {
  
  # seqid = transcript_id
  outfmt6.names <- c("transcript_id", "subject", "identity", "coverage", "mismatches", "gaps", "seq_start", "seq_end", "sub_start", "sub_end", "e", "score")
  
  
  df <- read_tsv(f, col_names = F) 
  
  colnames(df) <- outfmt6.names
  
  df <- df %>% mutate(db = basename(f))
  
  return(df)
  
  
}

library(tidyverse)

dir <- "~/Documents/GitHub/conopeptides/05.Prediction/"

dir_out <- dir

f <- list.files(dir, pattern = "blastx.outfmt6", full.names = T)

dbname <- gsub(".fa.diamond.blastx.outfmt6","",basename(f))

DB1 <- do.call(rbind, lapply(f, read_outfmt6)) %>% mutate(db = dbname) %>%
  dplyr::rename("identifier" = "subject") %>%
  mutate(identifier =  sapply(strsplit(identifier, "[|]"), `[`, 1)) %>%
  select(transcript_id, identifier, identity, e) %>%
  dplyr::rename("dna_identity" = "identity", "dna_e" = "e") 

# DB1 %>% count(identifier, sort = T)

# Load conodictor

f <- list.files(dir, pattern = "summary.csv", full.names = T)

cndictor <- read_csv(f) %>% dplyr::rename("transcript_id" = "sequence") %>%
  separate(transcript_id,into = c("transcript_id", "frame"), sep = "_frame=") 

# Load conoserverDB

dir <- "~/Documents/GitHub/conopeptides/"

conoserverDB <- read_tsv(file.path(dir, "conoserver_protein.tsv")) 

names(conoserverDB) <- paste0(gsub(" ","_", names(conoserverDB)), "_conoserver")

names(conoserverDB)[1] <- "identifier"

# Bind conoserver + conodictor

DB1 <- DB1 %>% 
  # distinct(transcript_id, identifier) %>% 
  left_join(conoserverDB) %>%
  left_join(cndictor, by = "transcript_id") 


# add sequences

dir <- "/Users/cigom/Documents/GitHub/conopeptides/04.Merge/transdecoder_dir/"

f <- list.files(dir, pattern = "Merged_clusters.fasta$", full.names = T) 

DNA <- Biostrings::readDNAStringSet(f)

keepdna <- sapply(strsplit(names(DNA), " "), `[`, 1)

names(DNA) <- keepdna

str(QUERY <- DB1 %>% distinct(transcript_id) %>% pull())

# Sanity check

sum(keepdna %in% QUERY) # must match str(QUERY) # 2526

DNA <- DNA[keepdna %in% QUERY]

width <- Biostrings::width(DNA)

DNADB <- data.frame(dna_width = width, dna_sequence = DNA, transcript_id = names(DNA)) %>% as_tibble() 


DB1 <- DNADB %>% right_join(DB1)

# ORFS

f <- list.files(dir, pattern = "Merged_clusters.fasta.transdecoder.pep", full.names = T) 

PEP <- Biostrings::readAAStringSet(f)

keeppep <- sapply(strsplit(names(PEP), " "), `[`, 1)

head(transcript_id <- sapply(strsplit(keeppep, ".p[0-9]"), `[`, 1))

names(PEP) <- keeppep

sum(keep <- transcript_id %in% QUERY) # must match str(QUERY) # 499

PEP <- PEP[keep]


# Get the longest
width <- Biostrings::width(PEP)

PEPDB <- data.frame(width = width, PEP, transcript_id = transcript_id[keep], peptide_id = names(PEP)) %>% as_tibble() 

# PEPDB %>% filter(transcript_id %in% "TRINITY_DN102516_c0_g1_i1") %>% filter(width == max(width))
PEPDB %>% count(transcript_id, sort = T)

PEPDB <- PEPDB %>%
  group_by(transcript_id) %>%
  arrange(desc(width)) %>%
  slice_head(n = 1) %>% ungroup()
  # filter(width == max(width)) %>%

# Sanity check, n == 1 per transcript

PEPDB %>% count(transcript_id, sort = T)

# BIND ORF TO BLASTP

f <- list.files(dir, pattern = ".blastp.outfmt6", full.names = T) 

dbname <- gsub(".fa.diamond.blastp.outfmt6","",basename(f))

DB2 <- do.call(rbind, lapply(f, read_outfmt6)) %>% mutate(db = dbname) %>%
  dplyr::rename("identifier" = "subject", "peptide_id" = "transcript_id") %>%
  mutate(identifier =  sapply(strsplit(identifier, "[|]"), `[`, 1))


DB2 <- PEPDB %>% right_join(DB2, by = "peptide_id")

DB2 <- DB2 %>% select(transcript_id, peptide_id, PEP, width, identifier, identity, e) %>%
  dplyr::rename("orf_width" = "width", "orf_sequences" = "PEP", "orf_id" = "peptide_id",
    "orf_identity" = "identity", "orf_e" = "e") 

DB2 %>% distinct(transcript_id) # 291

sum(DB2$transcript_id %in% QUERY) # must match only 499 with orfs

DB <- DB1 %>% left_join(DB2) 

write_tsv(DB, file = file.path(dir_out, "Merged_clusters_conoserver_and_conodictor.tsv"))


contig_width <- sort(Biostrings::width(DNA), decreasing = T)

contig_width <- sort(Biostrings::width(PEP), decreasing = T)

metrics_Nx <- function(contig_width) {
  
  # width <- contig_Nx(f)
  
  width <- contig_width
  
  assembly_seqs <- length(width)
  
  v <- seq(0.1,1, by = 0.1)
  
  Lx <- function(x) { sum(cumsum(width) < (sum(width) * x)) + 1 }
  
  l <- unlist( lapply(v, Lx))
  
  # add the number of sequences per Nx (field assembly_seqs)
  
  n_seqs <- sum(width > width[l[1]]) + 2
  
  # Cut the widths into chunks based on the Nx breakpoints 
  breakpts <- l
  
  # chunks <- cut(seq_along(width), breaks = l, labels = FALSE)
  
  chunks <- cut(rev(seq_along(width)), breaks = breakpts, labels = FALSE) 
  
  # Split the widths into chunks based on the cuts 
  chunks <- split(width, chunks) 
  
  n_seqs <- c(n_seqs,   as.vector(unlist(lapply(chunks, length))))
  
  n_frac <- n_seqs/assembly_seqs
  
  metrics_df <- data.frame(x = paste0("N", v*100), n = width[l], l = l, 
    n_seqs, n_frac, Assembly = basename(f))
  
  
  return(metrics_df)
  
}

metrics_Nx(contig_width) %>%
  mutate(x = factor(x, levels = unique(x))) %>%
  mutate(facet = "Conopeptides (346) transcriptome assembly") %>%
  ggplot(aes(x = x, y = n, group = Assembly, color = Assembly)) +
  geom_vline(xintercept = "N50", linetype="dashed", alpha=0.5) +
  ggplot2::geom_path(linewidth = 1.5, lineend = "round") +
  geom_point(shape = 21, size = 4) +
  facet_grid(~ facet) +
  labs(x = "Nx", y = "Contig length", color = "Assembly method") +
  scale_color_grey("") +
  scale_fill_grey("") +
  # scale_fill_manual("Assembly method", values = c("black", "grey89")) +
  guides(color=guide_legend(title = "", nrow = 1)) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()) 


metrics_Nx(contig_width) %>%
  mutate(x = factor(x, levels = unique(x))) %>%
  ggplot(aes(y = n_frac, x = x, group = Assembly, fill = Assembly)) +
  # ggplot2::geom_col() +
  ggplot2::geom_col(position = position_dodge2(reverse = T)) +
  geom_text(aes(label = n_seqs), vjust = -0.25) +
  labs(y = "Frac. of Scaffolds", x = "Nx", fill = "Assembly method") +
  scale_color_grey("") +
  scale_fill_grey("") +
  ylim(0,0.5) +
  guides(fill=guide_legend(title = "", ncol = 1)) +
  theme_bw(base_size = 16, base_family = "GillSans") +
  theme(legend.position = "top", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank())


