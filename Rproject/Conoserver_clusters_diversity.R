# Count the profile of gene superfamilies and multi, and full annotations
# therefore plot dendogram of sequences per method
# OPTIMIZED: Using vfold_cv with stratified resampling and clustering per fold

# Based on the number of sequence clusters and the number of members within each cluster, 
# you can estimate parameters related to protein family size, conservation, and diversity. These metrics provide a quantitative measure of the evolutionary and functional landscape of your dataset.

# 
# Key Parameters and Their Interpretation
# Protein Family Size Distribution: By plotting a histogram of the number of members per cluster, 
# you can visualize the distribution of protein family sizes. 
# This helps you identify whether your dataset is dominated by a few large, highly expanded families (e.g., kinases) or 
# by many small, singleton families.
#
# Protein Family Diversity: The total number of clusters is a direct measure of the number of distinct protein families in your dataset. 
# A higher number of clusters relative to the total number of sequences indicates greater diversity.


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)
library(rsample)
library(patchwork)

my_custom_theme <- function(base_size = 10, legend_pos = "top", ...) {
  theme_bw(base_family = "GillSans", base_size = base_size) +
    theme(legend.position = legend_pos,
          strip.placement = "outside", 
          strip.background = element_rect(fill = 'gray90', color = 'white'),
          strip.text = element_text(angle = 0, size = base_size, hjust = 0), 
          axis.text = element_text(size = rel(0.7), color = "black"),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          ...
    )
}

# ============================================================================
# LOAD DATA
# ============================================================================

dir <- "/Users/rjegr/Documents/Windows/Documents/ConotoxinBenchmark/INPUTS/ConoSorter_dir/"
outdir <- "/Users/rjegr/Documents/GitHub/ConotoxinBenchmark/INPUTS/"

recode_to <- c("STRINGTIE","SPADES", "TRINITY","IDBA", "MEGAHIT", "RNABLOOM" ,"BRIDGER", "TRANSABBYS", "BINPACKER","SOAPDENOVO" ,"CSTONE", "TRANSLIG", "Baseline", "PLASS")

recode_to <- structure(c("StringTie","rnaSPAdes", "Trinity", "IDBA", "MEGAHIT", "RNA-Bloom","BRIDGER","Trans-ABySS", "BinPacker", "SOAP-denovo", "Cstone", "TransLiG", "Baseline", "PinguiN (nuclassemble)"), names = recode_to)

recode_col <- c("full","multi", "fragment","chimera", "error", "NA")
recode_col <- structure(c("Full","Multi", "Fragment", "Chimera", "Error", "Unnanotated"), names = recode_col)

outName <- "protein_1" 

f <- file.path(dir, paste0(outName, "_ConoSorter.rds"))

# DB <- read_rds(f) |> mutate(Assembler = gsub(".fa.transdecoder", "", Method)) |>
#   mutate(vfold_set = sapply(strsplit(Assembler, "_"), `[`, 1)) |> 
#   mutate(Assembler = sapply(strsplit(Assembler, "_"), `[`, 5)) |> 
#   dplyr::mutate(Assembler = dplyr::recode(Assembler, !!!recode_to))

dir <- "/Users/rjegr/Documents/Windows/Documents/ConotoxinBenchmark/INPUTS/BLAST_based_annotation_dir/"
f <- list.files(dir, full.names = T, pattern = ".rds")

annotation_results <- do.call(rbind, lapply(f, read_rds)) |>
  mutate(file_name = gsub("_into_Fold[0-9]+[0-9]+.[1|2].blast", "", file_name)) |> 
  dplyr::rename("protein_id" = qseqid, "Assembler" = file_name) |>
  mutate(vfold_set = sapply(strsplit(Assembler, "_"), `[`, 1), 
         Assembler = sapply(strsplit(Assembler, "_"), `[`, 5)) |> 
  dplyr::mutate(Assembler = dplyr::recode(Assembler, !!!recode_to)) |>
  dplyr::mutate(final_annotation = dplyr::recode_factor(final_annotation, !!!recode_col))


dir <- "/Users/rjegr/Documents/Windows/Documents/ConotoxinBenchmark/INPUTS/"

f <- list.files(path = dir, pattern = "curated_nuc_conoServerDB.rds", full.names = T)

conoServerDB <- read_rds(f) |> 
  mutate(genesuperfamily = ifelse(is.na(genesuperfamily), "Other", genesuperfamily)) |>
  mutate(genesuperfamily = gsub(" superfamily", "", genesuperfamily),
         # name2 = sapply(strsplit(name, "[|]"), `[`, 2),
         genesuperfamily = ifelse(grepl("Divergent", genesuperfamily), "Divergent", genesuperfamily)
         ) |>
  dplyr::rename("hits" = "entry_id")

.DB <- read_tsv(file.path(outdir, "benchmark_assemblers.tsv")) |>
  dplyr::rename("protein_id" = contig_name) |> ungroup() |> 
  # mutate(protein_id = gsub(".p[0-9]+$", "", protein_id)) |>
  dplyr::mutate(Assembler = dplyr::recode(Assembler, !!!recode_to)) |>
  distinct(protein_id, Assembler, vfold_set, reference_coverage, hits) |>
  left_join(annotation_results, by = c("protein_id", "Assembler", "vfold_set")) 


DB <- conoServerDB |>
  distinct(hits, name, genesuperfamily, sequence) |> 
  # right_join(DB, by = c("protein_id", "Assembler", "vfold_set")) |>
  right_join(.DB, by = c("hits"), relationship = "many-to-many") 

# ============================================================================
# CLUSTERIZE FUNCTION (from Conotoxin_protein_fam_diversity.R)
# ============================================================================

clusterize <- function(data, seq_type = "sequence") {
  
  identify_sequence_type <- function(seq_vector) {
    require(Biostrings)
    seq_string <- paste(seq_vector, collapse = "")
    
    is_dna <- tryCatch({
      DNAString(seq_string)
      TRUE
    }, error = function(e) {
      FALSE
    })
    
    if (!is_dna) {
      is_aa <- tryCatch({
        AAString(seq_string)
        TRUE
      }, error = function(e) {
        FALSE
      })
      
      if (is_aa) {
        return("Amino Acid")
      } else {
        return("Not a valid DNA or Amino Acid sequence")
      }
    } else {
      return("DNA")
    }
  }
  
  seq_subject <- data |>
    dplyr::rename("subject" = seq_type) |>
    distinct(subject) |> 
    pull(subject, name = subject) 
  
  is_dna <- identify_sequence_type(sample(seq_subject, 1))
  
  if(is_dna == "DNA") {
    cat("\n is DNA, using DNAStringSet\n ")
    seq_subject <- Biostrings::DNAStringSet(seq_subject)
  } else {
    cat("\n is AA, using AAStringSet\n ")
    seq_subject <- Biostrings::AAStringSet(seq_subject)
  }
  
  clusters <- DECIPHER::Clusterize(seq_subject,
                                   cutoff = 0.5,  
                                   minCoverage = 0.95, 
                                   processors = NULL)
  
  seq_name <- paste0(seq_type, "_clus")
  
  clustersdf <- data.frame(
    seq_type = rownames(clusters), 
    seq_name = paste0(seq_type, "_", clusters$cluster)) |> as_tibble()
  
  names(clustersdf) <- c(seq_type, seq_name)
  
  return(clustersdf)
}

# ============================================================================
# OPTIMIZED: STRATIFIED VFOLD_CV WITH BOOTSTRAPPING AND CLUSTERING
# ============================================================================

# Step 1: Filter superfamilies with n > 10
data_filtered <- conoServerDB |>
  dplyr::count(genesuperfamily, sort = TRUE) |>
  filter(n >= 20) |>
  left_join(conoServerDB, by = "genesuperfamily") |>
  distinct(genesuperfamily, sequence, proteinsequence, hits)


cat("\nFiltered data: ", nrow(data_filtered), " sequences from ", 
    n_distinct(data_filtered$genesuperfamily), " superfamilies\n")

# Step 2: Create 12-fold cross-validation with stratification by genesuperfamily
folds <- rsample::vfold_cv(data_filtered, v = 12, strata = "genesuperfamily") #, 

# Step 3: Define stratify_data function - applies to each fold
stratify_data <- function(split) {
  
  # Get analysis set (training set)
  analysis_data <- rsample::analysis(split)
  fold_id <- labels(split)$id
  
  cat("\n========================================")
  cat("\nProcessing fold: ", fold_id, "\n")
  cat("  Analysis set size: ", nrow(analysis_data), "\n")
  
  
  # Step 3.1: Bootstrap resampling to equalize n_size = 20 per superfamily
  
  
  analysis_bootstrapped <- 
    split(analysis_data, analysis_data$genesuperfamily) |>
    lapply(function(sf_data) {
      n_size <- 20
      # n_available <- nrow(sf_data)
      # FIRST: Get unique sequences per superfamily
      sf_unique <- sf_data %>% distinct(sequence, .keep_all = TRUE)
      n_available <- nrow(sf_unique)
      
      
      
      if (n_available >= n_size) {
        # Downsample if more than n_size
        sf_data |> slice_sample(n = n_size, replace = FALSE)
      } else {
        # Oversample (bootstrap) if less than n_size
        sf_data |> slice_sample(n = n_size, replace = TRUE)
      }
    }) |>
    dplyr::bind_rows() |>
    distinct(sequence, .keep_all = TRUE)  # Remove duplicates from bootstrap
  
  cat("  After bootstrap resampling: ", nrow(analysis_bootstrapped), "\n")
  
  # Step 3.2: Clusterize sequences
  cat("  Clusterizing sequences...\n")
  sequence_clusters <- clusterize(analysis_bootstrapped, seq_type = "sequence")
  
  # Step 3.3: Count sequences and clusters per superfamily
  results <- analysis_bootstrapped |>
    distinct(genesuperfamily, sequence) |>
    left_join(sequence_clusters, by = "sequence") |>
    group_by(genesuperfamily) |>
    summarise(
      seq_number = n_distinct(sequence),
      clus_number = n_distinct(sequence_clus),
      .groups = "drop"
    ) |>
    mutate(
      index = 1 - (clus_number / seq_number),
      fold_id = fold_id
    ) |>
    mutate(seq_number = as.numeric(seq_number))
  
  cat("  Summary: ", nrow(results), " superfamilies processed\n")
  
  return(results)
}

# Step 4: Apply stratify_data to all folds using lapply + bind_rows
fold_results <- dplyr::bind_rows(
  lapply(folds$splits, stratify_data)
)

cat("\n========================================\n")
cat("Total results: ", nrow(fold_results), " rows\n")
cat("Folds: ", n_distinct(fold_results$fold_id), "\n")
cat("Superfamilies: ", n_distinct(fold_results$genesuperfamily), "\n")

# ============================================================================
# SUMMARY STATISTICS: Mean and SD by Superfamily
# ============================================================================

summary_stats <- fold_results |>
  group_by(genesuperfamily) |>
  summarise(
    mean_index = mean(index, na.rm = TRUE),
    sd_index = sd(index, na.rm = TRUE),
    mean_seq = mean(seq_number, na.rm = TRUE),
    sd_seq = sd(seq_number, na.rm = TRUE),
    n_folds = n_distinct(fold_id),
    .groups = "drop"
  ) |>
  arrange(desc(mean_index))

print(summary_stats)

# ============================================================================
# VISUALIZATION: Index (Diversity) by Superfamily with Mean ± SD
# ============================================================================

fold_results <- fold_results |>
  group_by(genesuperfamily) |>
  mutate(mean_index = mean(index, na.rm = TRUE)) |>
  ungroup() |>
  mutate(genesuperfamily = fct_reorder(genesuperfamily, mean_index)) 

p1 <- fold_results |>
  ggplot(aes(y = genesuperfamily, x = index, alpha = seq_number)) +
  geom_jitter(width = 0.1, height = 0.2, shape = 21, size = 2, color = "navy") +
  stat_summary(fun = "mean", geom = "point", size = 4, color = "red", shape = 18) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "errorbarh", height = 0.3, color = "red", size = 1) +
  scale_alpha_continuous("Seq. Number", range = c(0.3, 0.9)) +
  labs(
    title = "Protein Family Diversity Index Across 12 Folds",
    y = "Gene Superfamily",
    x = "Diversity Index (1 - clus_number/seq_number)",
    caption = "Points: individual fold results | Red diamond: mean ± 1 SD"
  ) +
  my_custom_theme(base_size = 11, legend_pos = "bottom") +
  theme(panel.grid.major.x = element_line(color = "gray90", size = 0.2))

print(p1)

ggsave(p1, filename = "diversity_index_by_superfamily.png", 
       path = outdir, width = 10, height = 8, dpi = 300)

# ============================================================================
# ADDITIONAL VISUALIZATION: Boxplot with individual points
# ============================================================================

p2 <- fold_results |>
  ggplot(aes(y = genesuperfamily, x = index, fill = genesuperfamily)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(alpha = seq_number), width = 0, height = 0.1, 
              shape = 21, size = 1.5, color = "black") +
  scale_alpha_continuous("Seq. Number", range = c(0.2, 0.8)) +
  scale_fill_viridis_d(guide = "none") +
  labs(
    title = "Distribution of Diversity Index Across 12 Folds",
    y = "Gene Superfamily",
    x = "Diversity Index (1 - clus_number/seq_number)"
  ) +
  my_custom_theme(base_size = 11, legend_pos = "bottom") +
  theme(panel.grid.major.x = element_line(color = "gray90", size = 0.2))

print(p2)

ggsave(p2, filename = "diversity_index_boxplot.png", 
       path = outdir, width = 10, height = 8, dpi = 300)

# ============================================================================
# SAVE RESULTS
# ============================================================================

# write_rds(fold_results, file.path(outdir, "fold_clustering_results.rds"))
# write_csv(summary_stats, file.path(outdir, "superfamily_diversity_summary.csv"))

cat("\n✓ Analysis complete! Results saved.\n")


# ============================================================================
# PLOT WITH STEP 2
# ============================================================================

discrete_scale <- DB %>% ungroup() %>% distinct(Assembler) %>% pull()

n <- length(discrete_scale)

scale_col <- c(ggsci::pal_startrek()(7), ggsci::pal_cosmic()(n-7))

scale_col <- structure(scale_col, names = sort(discrete_scale))

f_assembler <- c("StringTie","Spades", "Trinity", "IDBA", "MEGAHIT", "RNA-bloom", "Plass")

p3 <- DB |>
  filter(Assembler %in% f_assembler) |>
  dplyr::filter(!is.na(hits)) |>
  filter(final_annotation %in% c("Multi", "Full")) |>
  mutate(final_annotation
         = ifelse(final_annotation %in% c("Multi", "Full"), "Full and Multi", final_annotation)) |>
  # dedup for not inflate, just want to count intersection between assemblers
  distinct(hits, final_annotation, Assembler, genesuperfamily) |>
  filter(genesuperfamily %in% unique(summary_stats$genesuperfamily)) |>
  mutate(genesuperfamily = factor(genesuperfamily, levels = levels(fold_results$genesuperfamily))) |>
  count(hits, final_annotation, Assembler, genesuperfamily, sort = T) |>
  group_by(final_annotation, Assembler, genesuperfamily) |> tally(sort = T) |>
  ggplot(aes(y = genesuperfamily, x = Assembler, fill = Assembler, alpha = n)) +
  facet_grid(~final_annotation) +
  # geom_col() + # position = position_fill(reverse = T)
  geom_tile(height = 1, width = 0.9 , alpha = 0.2) +
  geom_text(aes(label = n), color = "black", size = 4) +
  scale_fill_manual(values = scales::alpha(scale_col, 0.9)) +
  my_custom_theme() +
  guides(fill = "none", color = "none", alpha = "none") 



colors_ <- c("#357EBD", "#486983", "#FEA31A", "#F1E688", "gray")


# DB |>
#   drop_na(hits, genesuperfamily) |>
#   filter(genesuperfamily %in% unique(summary_stats$genesuperfamily)) |>
#   mutate(genesuperfamily = factor(genesuperfamily, levels = levels(fold_results$genesuperfamily))) |>
#   dplyr::mutate(final_annotation = dplyr::recode_factor(final_annotation, !!!recode_col)) |>
#   ggplot(aes(x = reference_coverage, y = genesuperfamily, fill=final_annotation)) +
#   # facet_wrap(~ Method) +ggplot2::stat_ecdf()
#   ggridges::geom_density_ridges_gradient(
#     jittered_points = T,
#     position = ggridges::position_points_jitter(width = 0.05, height = 0),
#     point_shape = '|', point_size = 0.5, point_alpha = 1, alpha = 0.7) +
#   scale_fill_manual(values = colors_) +
#   my_custom_theme()


p4 <- DB |>
  filter(!final_annotation %in% c("Multi", "Full")) |>
  # dedup for not inflate, just want to count intersection between assemblers
  distinct(hits, final_annotation, Assembler, genesuperfamily) |>
  filter(genesuperfamily %in% unique(summary_stats$genesuperfamily)) |>
  mutate(genesuperfamily = factor(genesuperfamily, levels = levels(fold_results$genesuperfamily))) |>
  count(hits, final_annotation, Assembler, genesuperfamily, sort = T) |>
  group_by(final_annotation, Assembler, genesuperfamily) |> tally(sort = T) |>
  ggplot(aes(y = genesuperfamily, x = n, fill = final_annotation)) +
  # geom_col(position = position_fill(reverse = T)) +
  geom_col(position = position_stack(reverse = T)) +
  scale_fill_manual(values = colors_) +
  my_custom_theme() 

p1 + p3 + p4

# Dendogram ======

heatmapdf <- .DB %>% drop_na(hits) 

# labels <- heatmapdf %>% distinct(hits, name) %>% pull(name, name = hits)

# heatmapdf |> distinct(Assembler, hits)

# table(heatmapdf$genesuperfamily)

DataViz <- .DB |> 
  drop_na(hits) |>
  filter(reference_coverage > 0.95) |>
  distinct(Assembler, hits) |>
  group_by(hits) |>
  summarise(
    combination = paste(sort(unique(Assembler)), collapse = ","),
    n_assemblers = n_distinct(Assembler),
    .groups = "drop"
  ) |> 
  right_join(distinct(conoServerDB, hits, genesuperfamily, name, sequence))


# AS COUNT
# DataViz <- heatmapdf %>%
#   filter(reference_coverage > 0.95) |>
#   distinct(Assembler, hits) |>
#   count(hits)  |> 
#   right_join(distinct(conoServerDB, hits, genesuperfamily, name, sequence))

# USE THEN DATAVIZ TO RUN DENDOGRAM AND HEATMAP

seqs <- DataViz |> 
  distinct(hits, sequence) |> pull(name = hits, sequence)
  
library(DECIPHER)
library(phangorn)

alignment <- DECIPHER::AlignSeqs( Biostrings::DNAStringSet(seqs))

# This method is time-computer consuming
phangAlign <- phangorn::phyDat(as(alignment, "matrix"), type="DNA") 
dm <- phangorn::dist.ml(phangAlign)
treeNJ <- phangorn::NJ(dm) # Note, tip order != sequence order
# plot(treeNJ)
fit <- phangorn::pml(treeNJ, data=phangAlign) # this step is fast
fitGTR <- stats::update(fit, k = 4, inv = 0.2) # this step is fast
fitGTR <- phangorn::optim.pml(fitGTR) # this step is fast

t <- fitGTR$tree

# calculate Pairwise Distances from a Phylogenetic Tree 
# disimilariy matrx using cophenetic.phylo
hclust <- hclust(dist(ape::cophenetic.phylo(t)))

# plot(t)

# plot(hclust)

# PLOT

# recode_to <- structure(c("StringTie","rnaSPAdes", "Trinity", "IDBA", "MEGAHIT", "RNA-Bloom","BRIDGER","Trans-ABySS", "BinPacker", "SOAP-denovo", "Cstone", "TransLiG", "Baseline", "PinguiN (nuclassemble)")
                       
f_assembler <- c("StringTie","rnaSPAdes", "Trinity", "IDBA", "MEGAHIT", "RNA-Bloom", "PinguiN (nuclassemble)")


p <- DataViz %>%
  mutate(
    assemblers = str_split(combination, ",")
  ) |>
  unnest(assemblers) |>
  filter(assemblers %in% f_assembler) |>
  count(assemblers, hits, genesuperfamily) |>
  mutate(hits = factor(hits, levels = hclust$labels[hclust$order])) |>
  drop_na(assemblers) |>
  # mutate(col = ifelse(is.na(Method), "Not-assembled", "Assembled")) %>%
  ggplot(aes(y = assemblers, x = hits, fill = n)) + #  fill = as.factor(Score_sf)
  geom_tile(width = 0.5, height = 0.25) +
  # facet_grid(~ genesuperfamily, space = "free", scales = "free") + # switch = "x"
  # scale_fill_viridis_c(option = "magma", direction = -1) +
  scale_x_discrete("Conopeptides (hits of contigs assembled)",position = "top") +
  # scale_y_discrete(limits = rev(recode_to), breaks = rev(recode_to)) +
  # theme_bw(base_family = "GillSans", base_size = 12) + 
  my_custom_theme(legend_pos = "none") +
  theme(strip.text = element_text(angle = 90, size = 7, hjust = 0, vjust = 0.5),
        axis.ticks.length=unit(14,"pt")) +
  labs (title = "Co-ocurrance of conotoxin identification", 
        subtitle = "Artificial fastq true set analysis",
        x = "Sequence alignment (Pairwise Distances)") 

library(ggh4x)

p <- p + 
  ggh4x::scale_x_dendrogram(hclust = hclust, position = 'bottom', labels = NULL) +
  guides(x.sec = guide_axis_manual(labels = labels, label_size = 5, angle = 90, label_family = "GillSans"))  

p
