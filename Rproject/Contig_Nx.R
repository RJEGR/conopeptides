# Ricardo Gomez-Reyes
# Visualize Nx and BUSCO completness 
# Calculate nx distribution
# Include sizes after transrating
# Try to tag conopeptides according to conoserver dataset


# Include Merge > superduper
# /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/californicus/04.Merge/MMseqs/MMap/S1_HISAT2_SAM_BAM_FILES_Merged_clusters_DIR/LACE_Merged_clusters_DIR

rm(list = ls())

if(!is.null(dev.list())) dev.off()

library(tidyverse, help, pos = 2, lib.loc = NULL)


require(Biostrings)
require(dplyr)
require(ggplot2)

# dir <- "/Users/cigom/Documents/GitHub/conopeptides/02.Assembly/"
# dir <- "/Users/cigom/Documents/GitHub/conopeptides/07.Reference/"

pepdir <- "/Users/cigom/Documents/GitHub/conopeptides/05.Prediction/ConoSorter_dir/SORTED_regex_pHMM_dir/"
pepf <- list.files(pepdir, "pep$", full.names = T)


dir <- "/Users/cigom/Documents/GitHub/conopeptides/Nx_Metrics_dir/"

f <- list.files(dir, "fasta$", full.names = T)

# f <- f[grepl("^transcripts.fasta|blastx", basename(f))]


contig_Nx <- function(f, stringSet = "DNA") {
  
  if(stringSet != "DNA") {
    stringSet <- Biostrings::readAAStringSet(f)
  } else
  
  stringSet <- Biostrings::readDNAStringSet(f)
  
  contig_width <- sort(Biostrings::width(stringSet), decreasing = T)
  
  # contig_df <- data.frame(width = contig_width, Assembly = basename(f))
  
  return(contig_width)
}


metrics_df <- function(f, stringSet = "DNA") {
  
  # width <- sort(seq(1,100), decreasing = TRUE)
  
  width <- contig_Nx(f, stringSet = stringSet)
  
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

# sum(metrics_df(f[1])$n_seqs)

pepdf <- lapply(pepf, metrics_df, stringSet = "AA")
pepdf <- do.call(rbind, pepdf)
pepdf <- pepdf %>% mutate(stringSet = "B) Orfs (Transcriptome)")


df <- lapply(f, metrics_df)

df <- do.call(rbind, df)

df <- mutate(df, x = factor(x, levels = unique(df$x)))

unique(df$Assembly)


recode_to <- c(
  "spades_hisat_superDuper.fasta",
  "spades.fasta",
  "Merged_hisat_SuperDuper.fasta",
  "MMseqs_hisat_SuperDuper.fasta",
  "MMseqs.fasta",
  "trinity_hisat_superDuper.fasta",
  "Trinity.fasta")

recode_to <- structure(
  c(
    "Spades-Lace (Hisat)", 
    "Spades (S)", 
    "Concat-Lace T-S (Hisat)",
    "MMseqs-Lace T-S (Hisat)",
    "MMseqs (T-S)",
    "Trinity-Lace (Hisat)",
    "Trinity (T)"), 
  names = recode_to)

df <- mutate(df, Assembly = dplyr::recode_factor(Assembly, !!!recode_to))

recode_to <- c(
  "SuperDuper_Spades_regex_pHMM.pep",
  "spades_regex_pHMM.pep",
  "Merged_hisat_SuperDuper_regex_pHMM.pep",
  "MMseqs_hisat_SuperDuper_regex_pHMM.pep",
  "MMseqs_regex_pHMM.pep",
  "SuperDuper_Trinity_regex_pHMM.pep",
  "Trinity_regex_pHMM.pep")

recode_to <- structure(
  c(
    "Spades-Lace (Hisat)", 
    "Spades (S)", 
    "Concat-Lace T-S (Hisat)",
    "MMseqs-Lace T-S (Hisat)",
    "MMseqs (T-S)",
    "Trinity-Lace (Hisat)",
    "Trinity (T)"), 
  names = recode_to)

pepdf <- mutate(pepdf, Assembly = dplyr::recode_factor(Assembly, !!!recode_to))

# dir <- "/Users/cigom/Documents/GitHub/conopeptides/04.Merge/transdecoder_dir/"
# # dir <- "/Users/cigom/Documents/GitHub/conopeptides/07.Reference/"
# 
# f <- list.files(dir, "fasta$", full.names = T)
# df2 <- lapply(f, metrics_df)
# df2 <- do.call(rbind, df2)

# df <- rbind(df, df2)


df <- df %>% mutate(stringSet = "DNA (Transcriptome)") %>% rbind(pepdf)

# rnsps <- mean(contig_Nx(f[[1]]))
# trnt <- mean(contig_Nx(f[[2]]))

p <- ggplot(df, aes(x = x, y = n, group = Assembly, color = Assembly)) +
  geom_vline(xintercept = "N50", linetype="dashed", alpha=0.5) +
  ggplot2::geom_path(linewidth = 1.5, lineend = "round") +
  facet_wrap(~  stringSet, scales = "free_y") +
  geom_point(shape = 21, size = 4, aes(size = n_frac)) +
  labs(x = "Nx", y = "Contig length", color = "Assembly method") +
  # ggsci::scale_color_jco() +
  ggsci::scale_color_startrek() +
  # guides(color=guide_legend(title = "", nrow = 2, ncol = 4)) +
  theme_bw(base_size = 14, base_family = "GillSans") +
  theme(legend.position = "top", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()) 

p <- p + guides(color = guide_legend(title = "", nrow = 2, ncol = 4))

p  

# annotate("text", y = rnsps, x = "N50", angle = 90, label = "label")

# p

ggsave(p, filename = 'Nx-methods.png', path = dir, width = 10, height = 5, device = png, dpi = 300)

p2 <- ggplot(df, aes(x = x, y = n_frac, group = Assembly, fill = Assembly)) +
  # ggplot2::geom_col() +
  ggplot2::geom_col(position = position_dodge2(reverse = T)) +
  labs(x = "", y = "Frac. of Scaffolds", fill = "Assembly method") +
  # scale_color_grey("") +
  # scale_fill_grey("") +
  ggsci::scale_fill_startrek() +
  guides(fill=guide_legend(title = "", ncol = 1)) +
  facet_wrap(~  stringSet, scales = "free_y") + 
  scale_x_discrete(position = "top") +
  theme_bw(base_size = 14, base_family = "GillSans") +
  theme(legend.position = "none", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()) 

df %>% 
  group_by(stringSet, x) %>% mutate(frac = n/sum(n)) %>%
  ggplot(aes(x = x, y = n_frac, group = Assembly, fill = Assembly)) +
  facet_grid(Assembly ~  stringSet, scales = "free_y") + 
  ggplot2::geom_col() +
  labs(x = "", y = "Frac. of Scaffolds", fill = "Assembly method") +
  # scale_color_grey("") +
  # scale_fill_grey("") +
  ggsci::scale_fill_startrek() +
  scale_x_discrete(position = "top") +
  theme_bw(base_size = 14, base_family = "GillSans") +
  theme(legend.position = "none", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()) 

# p2 + guides(fill = guide_legend(title = "", nrow = 2, ncol = 4))

library(patchwork)

p2 <- p / p2


ggsave(p2, filename = 'Nx-methods-2.png', path = dir, width = 7.5, height = 10, device = png, dpi = 300)

df %>% group_by(Assembly) %>% summarise(n_seqs = sum(n_seqs), n_frac = sum(n_frac)) %>% arrange(desc(n_seqs))
