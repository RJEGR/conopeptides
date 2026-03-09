


rm(list = ls())

if(!is.null(dev.list())) dev.off()
# ============================================================================
# Complete Example: Split sequences into chunks of 500 and save as FASTA files
# ============================================================================

library(dplyr)
library(Biostrings)


# Set output directory
# 

# pub_dir <- "C://Users//cinai/OneDrive/Documentos/PUBLICATION_DIR/"
# 

pub_dir <- "C://Users//cinai/OneDrive/Documentos/PUBLICATION_DIR/"

CONOPEPDB <- read_tsv(paste0(pub_dir, "/conopeptides.tsv")) %>%  #view()
  mutate(len = nchar(pep_seq)-1) %>%
  filter(effective_length_frac > 0.5) %>%
  # To be consistent w/ RES
  filter(Signalp_class == "SP") %>%
  drop_na(prediction_tool, Superfamily) 



outdir <- "C://Users//cinai/OneDrive/Documentos/PUBLICATION_DIR/fasta_chunks"

dir.create(outdir, showWarnings = FALSE)

# ============================================================================
# STEP 1: Define the dedup_StringSet function
# ============================================================================
dedup_StringSet <- function(seqset) {
  
  
  require(Biostrings)
  
  seq_chars <- as.character(seqset)
  split_names <- split(names(seqset), seq_chars)
  unique_seqs <- Biostrings::AAStringSet(names(split_names))
  names(unique_seqs) <- sapply(split_names, paste, collapse = "|")
  unique_seqs
}

# ============================================================================
# STEP 2: Define the split_chunks function
# ============================================================================
split_chunks <- function(df, fasta_file) {
  
  cat("Writing", nrow(df), "sequences to", fasta_file, "\n")
  
  df %>%
    mutate(pep_seq = gsub("\\*", "", pep_seq)) %>%
    pull(pep_seq, name = protein_id) %>%
    Biostrings::AAStringSet() %>%
    dedup_StringSet() %>%
    Biostrings::writeXStringSet(fasta_file)
}

# ============================================================================
# STEP 3: Split data into chunks of 500 sequences
# ============================================================================
df <- CONOPEPDB

n_sequences <- 300

# Calculate number of chunks needed
n_chunks <- round(nrow(df) / n_sequences)

# Create grouping factor
group_factor <- cut(
  seq_len(nrow(df)),
  breaks = n_chunks,
  labels = FALSE
)

# Split the data frame into a list of data frames
list_of_dfs <- split(df, group_factor)

# View results
print(paste("Total sequences:", nrow(df)))
print(paste("Number of chunks:", length(list_of_dfs)))
print(paste("Sequences per chunk:", sapply(list_of_dfs, nrow)))

# ============================================================================
# STEP 4: Save each chunk as FASTA file
# ============================================================================

lapply(
  names(list_of_dfs),
  function(chunk) {
    fasta_name <- paste0("sequences_chunk_", chunk, ".fasta")
    fasta_file <- file.path(outdir, fasta_name)
    split_chunks(list_of_dfs[[chunk]], fasta_file)
  }
)

# ============================================================================
# OPTIONAL: Verify the output files
# ============================================================================
fasta_files <- list.files(outdir, pattern = ".fasta", full.names = TRUE)

print(paste("Generated", length(fasta_files), "FASTA files:"))

print(fasta_files)

# Check sequence counts in each file
for (file in fasta_files) {
  seqs <- Biostrings::readAAStringSet(file)
  cat(basename(file), ":", length(seqs), "sequences\n")
}


# READ outputs and viz
# 
# The boundaries of the mature peptide regions are identified using sequence patterns that 
# describe cleavage sites of proprotein convertases (Table 1), 
# and of two exo-peptidases that are known to act on mollusks 
# (carboxypeptidase E [CPE] and peptidylglycine alpha-amidating monooxygenase [PAM]). 
# CPE cleaves the Lys or Arg terminal residues and PAM cleaves yje C-terminal glycine and 
# amidated the C-terminus. Finally the cysteine framework of the peptide is identified as 
# well as the mature conopeptide with the most similar sequence in ConoServer.
# 

library(tidyverse)

dir <- "C://Users//cinai/Downloads/"

f <- list.files(path = dir, pattern = "csv__", full.names = T)

sel_col <- c("Name","Signal","Pre","Mature","Class","Framework")

DB <- read_csv(f) |> select(any_of(sel_col)) |> 
  distinct(Name, Mature, Class, Framework) |>
  dplyr::rename("protein_id" = Name)

nrow(DB)

DB |>
  count(Class, sort = T)

recode_to <- structure(c("Control","Shrimp", "Mollusk", "Polychaete", "Mixed"))

recode_to <- structure(recode_to, names = c("Ctrl","Cam","Lit","Pol","Mix"))

dir <- "C://Users//cinai/OneDrive/Documentos/PUBLICATION_DIR/06.Quantification/MATRIX_RSEM_dir/"

f <- "glmLRT_multiple_contrast_ctrl_and_treatments_kallisto.rds"

DEGS <- read_rds(file.path(dir, f)) %>% filter(abs(logFC) > 2 & FDR < 0.05) 

dat <- DEGS |> 
  distinct(protein_id, sam_group, sampleX) |>
  mutate(Feeding = sapply(strsplit(sampleX, "_"), `[`, 1)) |>
  dplyr::mutate(Feeding = dplyr::recode_factor(Feeding, !!!recode_to)) |> 
  distinct(protein_id, Feeding)  |>
  group_by(protein_id) |>
  summarise(
    combination = paste(sort(unique(Feeding)), collapse = ","),
    n = n_distinct(Feeding),
    .groups = "drop"
  )

dat   |> 
  left_join(distinct(CONOPEPDB, protein_id, tab, Superfamily)) |> 
  left_join(DB) |> drop_na()
  write_csv(file.path(pub_dir, "S3_relational_db_conotoxin_expressión.csv"))

read_rds(file.path(pub_dir, "S1_table_global_conotoxin_expressión.rds")) |> 
  # distinct(pep_seq, Superfamily, Feeding)  |>
  group_by(protein_id, pep_seq, Superfamily, tab) |>
  summarise(
    Expression = sum(fill),
    combination = paste(sort(unique(Feeding)), collapse = ","),
    n = n_distinct(Feeding),
    .groups = "drop"
  ) |> left_join(DB) |> left_join(dat) |>
  write_csv(file.path(pub_dir, "S3_relational_db_conotoxin_expressión.csv"))


DEGS |> 
  left_join(DB)


library(dplyr)
library(ggplot2)


# Selecting unique degs (see upset plot 89+87+58+25+2)
# 

unique_degs <- DEGS |> 
  distinct(protein_id, sam_group, sampleX) |>
  mutate(Feeding = sapply(strsplit(sampleX, "_"), `[`, 1)) |>
  dplyr::mutate(Feeding = dplyr::recode_factor(Feeding, !!!recode_to)) |> 
  distinct(protein_id, Feeding)  |>
  group_by(protein_id) |>
  summarise(
    combination = paste(sort(unique(Feeding)), collapse = ","),
    n = n_distinct(Feeding),
    .groups = "drop"
  ) |> filter(n == 1) |> distinct(protein_id)  |> pull()

# Filter for top 20 genes by absolute logFC

str(unique(DEGS$protein_id))

sum(unique_degs %in% unique(DB$protein_id))
sum(unique_degs %in% unique(CONOPEPDB$protein_id))

hist(DEGS$logFC)
hist(DEGS$logFC[DEGS$protein_id %in% unique_degs])

top20_genes <- DEGS %>%
  filter(protein_id %in% unique_degs) %>%
  left_join(DB) %>%
  filter(Class == "conotoxin") %>%
  group_by(sam_group) %>%
  arrange(desc(abs(logFC))) %>%
  slice_head(n = 5) %>%
  mutate(lfcSE = logFC/qnorm(FDR / 2)) %>%
  # Calculate error bar positions
  mutate(
    ymin = logFC - lfcSE,
    ymax = logFC + lfcSE,
    # Create factor for ordering
    # Mature = factor(Mature, levels = Mature),
    # Add significance indicator
    sig = ifelse(FDR < 0.05, "*", "")
  ) %>%
  arrange(logFC)  # Sort by log2FC for better visualization

# Create the plot
top20_genes %>%
  ggplot(aes(x = reorder(Mature, logFC), 
             y = logFC,
             fill = logFC > 0)) +
  geom_col(width = 0.7, color = "black", size = 0.3) +
  facet_grid(~ sam_group)  +
  geom_errorbar(aes(ymin = ymin, ymax = ymax),
                width = 0.2,
                color = "black",
                size = 0.4) +
  # Add significance stars
  geom_text(aes(y = ymax + 0.2, label = sig),
            size = 5,
            color = "black") +
  # Aesthetics
  scale_fill_manual(values = c("TRUE" = "red2", "FALSE" = "forestgreen"),
                    name = "Direction",
                    labels = c("TRUE" = "Up-regulated", "FALSE" = "Down-regulated")) +
  labs(
    title = "Top 20 Differentially Expressed Genes",
    x = "Gene Name",
    y = expression(Log[2] ~ "Fold Change"),
    caption = "* padj < 0.05"
  ) +
  theme_bw(base_family = "GillSans", base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  ) +
  coord_flip() -> p_top20

p_top20

# Save the plot
ggsave(p_top20, 
       filename = "top20_log2fc_genes.png",
       width = 8, height = 6, dpi = 300)
