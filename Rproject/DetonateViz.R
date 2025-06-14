
# Purpose: Determine wheter the conflict in conotoxin nomenclature are true new putative conotoxins or artifacts

# DetonateViz.R
# Doint this for the 12,136 CDS
# For transrate :
# Read contigs.csv
# Read bam_info.csv
# select columns score (further read paper to understand values https://pmc.ncbi.nlm.nih.gov/articles/PMC4971766/)
# parse to DataBase.R

# For detonate, 
# Read score.isoforms.results
# select column score effective_length, contig_impact_score and IsoPct (omit because isoforms are unmasked)
# IsoPct = (Transcript abundance / Total gene abundance) × 100
# contig_impact_score = Relative contribution of a contig to explaining RNA-Seq data (Log ratio of true contig vs. noise hypothesis; used for trimming)

# Technically, for each contig, the contig impact score is the log ratio of the probability that the contig represents a true biological transcript versus the probability that the contig is just background noise.

# A higher contig impact score means the contig is more likely to be a true and important part of the transcriptome; a negative score suggests the contig is likely noise or redundant and could be trimmed from the assembly to improve overall accuracy.

# This score enables assembly refinement by identifying and removing contigs that do not meaningfully explain the data, thereby improving assembly quality.

#  AS INDIVIDUAL VALUES DOES NOT SUPPORT ANY COMPARISON, LETS COMPARE FROM THE TRUE-SET ANALYSIS 

library(tidyverse)

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR/"

dir <- "/Users/cigom/Documents/GitHub/conopeptides/03.Coverage/Contiguity_dir/Detonate_conopeptides_dir/"

# LOAD data -----

# DB <- read_rds(paste0(pub_dir, "/structured_db.rds"))

DB <- read_tsv(paste0(pub_dir, "/conopeptides.tsv"))  %>% mutate(gene_id = gsub(".p[0-9]+$", "", protein_id))


list.files(dir, "", full.names = F)

# f <- list.files(dir, "conopeptides.isoforms.results", full.names = T)

f <- list.files(dir, "conopeptides.score.isoforms.results", full.names = T)

scoresdf <- read_tsv(f)

scoresdf <- scoresdf %>%  mutate(protein_id = sapply(strsplit(gene_id, "[|]"), `[`, 1)) %>% select(-transcript_id, -gene_id)

scoresdf <- scoresdf %>% left_join(DB) # %>% 

# IsoPct is commonly used as a threshold to remove lowly expressed or potentially misassembled isoforms from the assembly. For example, isoforms with an IsoPct of zero (i.e., not supported by any reads or fragments) are often removed, as they are likely artifacts

scoresdf %>% count(IsoPct)

# Filter artifacts: Remove isoforms with impact_score > 0.5 and TPM < 1.

scoresdf %>% select(contig_impact_score) %>% arrange(contig_impact_score)

scoresdf %>% 
  mutate(Annotation = ifelse(is.na(Conflict), "A) No redundant", "B) Redundant (conflict)")) %>%
  mutate(x = contig_impact_score, x = sign(x) * log(1+abs(x)) ) %>%
  ggplot(aes(effective_length/length, x)) + geom_point(aes(color = Annotation)) +
  ggthemes::scale_color_calc() +
  facet_grid(~ Signalp_class) +
  labs(y = "Contig impact score (pseudo-log transform)", x = "Frac. of Effective length (Detonate)") +
  theme_bw(base_family = "GillSans", base_size = 14) + theme(legend.position = "top")


scoresdf %>% 
  select(contig_impact_score, CPM) %>% 
  mutate(Bad = ifelse(contig_impact_score > 0.5 & CPM < 1, "Bad", "Good")) %>% 
  count(Bad)

scoresdf %>% ggplot(aes(contig_impact_score)) + geom_density()

# Focus only in the conflicts for potential NEW

scoresdf %>% 
  mutate(Conflict = ifelse(is.na(Conflict), "A) No redundant", "B) Redundant (conflict)")) %>%
  dplyr::count(tab, Conflict, Signalp_class)

caption <- "A higher contig impact score means the contig is more likely to be a true and important part of the transcriptome\nA negative score suggests the contig is likely noise or redundant and could be trimmed from the assembly."

scoresdf %>% 
  drop_na(Superfamily, Conflict) %>%
  filter(Signalp_class == "SP" & Score_sf != 1) %>%
  mutate(facet = ifelse(is.na(Conflict), "A) No redundant", "B) Redundant (conflict)")) %>%
  mutate(tab = ifelse(tab == "pHMM", "Novel candidates (pHMM)", "Homology predicted")) %>%
  mutate(x = contig_impact_score, x = sign(x) * log(1+abs(x)) ) %>%
  ggplot(aes(y = Superfamily, x = x, fill = after_stat(x))) +
  # geom_violin() +
  facet_grid(~ tab) +
  ggridges::geom_density_ridges_gradient(
    jittered_points = F,
    position = ggridges::position_points_jitter(width = 0.05, height = 0)) +
  geom_point(shape = 124, size = 3) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  scale_fill_viridis_c(option = "C") +
  labs(y = "Putative conotoxin with conflict assignment", x = "Contig impact score (pseudo-log transform)", caption = caption) +
  theme_bw(base_family = "GillSans", base_size = 14) + theme(legend.position = "top")

# Focus only in the conflicts

scoresdf %>% 
  drop_na(Region) %>%
  mutate(x = contig_impact_score, x = sign(x) * log(1+abs(x)) ) %>%
  mutate(facet = ifelse(is.na(Conflict), "A) No redundant", "B) Redundant (conflict)")) %>%
  ggplot(aes(y = x, x = facet, fill = after_stat(x))) +
  geom_violin() +
  geom_jitter(alpha = 0.5) +
  facet_grid( ~ Region) +
  scale_fill_viridis_c(option = "C") +
  labs(y = "Frac. of Effective length (Detonate)", x = "") +
  theme_bw(base_family = "GillSans", base_size = 14) + theme(legend.position = "none")

DB %>% 
  drop_na(Superfamily, Conflict) %>%
  # filter(Signalp_class == "SP" & Score_sf != 1) %>%
  mutate(facet = ifelse(is.na(Conflict), "", "B) Redundant (conflict)")) %>%
  mutate(x = contig_impact_score, x = sign(x) * log(1+abs(x)) ) %>% 
  ggplot(aes(x, effective_length_frac, color = Signalp_class)) + geom_point(shape = 1, size = 3) +
  facet_grid( ~ facet) +
  labs(y = "Frac. of Effective length (Detonate)", x = "Contig impact score (pseudo-log transform)") +
  theme_bw(base_family = "GillSans", base_size = 14) + theme(legend.position = "top")

DB %>% 
  # drop_na(Region) %>%
  mutate(x = contig_impact_score, x = sign(x) * log(1+abs(x)) ) %>%
  mutate(facet = ifelse(is.na(Conflict), "A) No redundant", "B) Redundant (conflict)")) %>%
  ggplot(aes(y = x, x = Signalp_class, fill = after_stat(x))) +
  geom_violin() +
  geom_jitter(alpha = 0.5) +
  # facet_grid( ~ Region) +
  scale_fill_viridis_c(option = "C") +
  labs(y = "Contig impact score (pseudo-log transform)", x = "") +
  theme_bw(base_family = "GillSans", base_size = 14) + theme(legend.position = "none")


DB %>% 
  # drop_na(Region) %>%
  mutate(x = contig_impact_score, x = sign(x) * log(1+abs(x)) ) %>%
  mutate(facet = sign(x)) %>%
  mutate(seq_len = nchar(pep_seq)-1) %>%
  ggplot(aes(seq_len * sign(x))) +
  geom_histogram() +
  facet_grid(facet ~ .) 
  scale_fill_viridis_c(option = "C") +
  labs(y = "Contig impact score (pseudo-log transform)", x = "") +
  theme_bw(base_family = "GillSans", base_size = 14) + theme(legend.position = "none")
