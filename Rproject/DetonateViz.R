# DetonateViz.R
# Read contigs.csv
# Read bam_info.csv
# select columns score (further read paper to understand values https://pmc.ncbi.nlm.nih.gov/articles/PMC4971766/)
# parse to DataBase.R

#  AS INDIVIDUAL VALUES DOES NOT SUPPORT ANY COMPARISON, LETS COMPARE FROM THE TRUE-SET ANALYSIS 

library(tidyverse)

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR/"

dir <- "/Users/cigom/Documents/GitHub/conopeptides/03.Coverage/Contiguity_dir/conopeptides_detonate_dir/"

# LOAD data -----

DB <- read_rds(paste0(pub_dir, "/structured_db.rds"))

list.files(dir, "", full.names = F)

# f <- list.files(dir, "conopeptides.isoforms.results", full.names = T)

f <- list.files(dir, "conopeptides.score.isoforms.results", full.names = T)

scoresdf <- read_tsv(f)


scoresdf <- scoresdf %>% left_join(DB) %>% drop_na(tab) %>% filter(Signalp_class == "SP")


scoresdf %>% count(IsoPct)

# Critical metric: Higher values indicate transcripts likely to be artifacts (prioritize low-impact isoforms).

# Filter artifacts: Remove isoforms with impact_score > 0.5 and TPM < 1.

scoresdf %>% select(contig_impact_score) %>% arrange(contig_impact_score)

scoresdf %>% 
  ggplot(aes(effective_length/length, contig_impact_score)) + geom_point()

scoresdf %>% 
  select(contig_impact_score, CPM) %>% 
  mutate(Bad = ifelse(contig_impact_score > 0.5 & CPM < 1, "Bad", "Good")) %>% 
  count(Bad)

scoresdf %>% ggplot(aes(contig_impact_score)) + geom_density()

scoresdf %>% 
  mutate(facet = ifelse(is.na(Conflict), "A) No redundant", "B) Redundant (conflict)")) %>%
  ggplot(aes(y = Superfamily, x = effective_length/length, fill = after_stat(x))) +
  # geom_violin() +
  facet_grid(~ facet) +
  ggridges::geom_density_ridges_gradient(
    jittered_points = T,
    position = ggridges::position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 3, point_alpha = 1, alpha = 1) +
  scale_fill_viridis_c(option = "C") +
  labs(y = "", x = "Effective length (Detonate)") +
  theme_bw(base_family = "GillSans", base_size = 14) + theme(legend.position = "none")


