# Transrate.R
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

dir <- "/Users/cigom/Documents/GitHub/conopeptides/03.Coverage/Contiguity_dir/conopeptides_to_conoserver_nucleic_transrate_dir//conopeptides/"

# LOAD data -----

DB <- read_rds(paste0(pub_dir, "/structured_db.rds"))

f <- list.files(dir, "contigs.csv", full.names = T)

transratedf <- read_csv(f) %>% dplyr::rename("protein_id" = "contig_name")

transratedf <- transratedf %>% left_join(DB) # %>% drop_na(tab) %>% filter(Signalp_class == "SP")

transratedf %>% dplyr::count(Superfamily, sort = T)

transratedf %>%
  drop_na(hits) %>%
  ggplot(aes(y = tab, x = reference_coverage, fill = after_stat(x))) +
  # geom_violin() +
  # facet_grid(~ facet) +
  ggridges::geom_density_ridges_gradient(
    jittered_points = T,
    position = ggridges::position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 3, point_alpha = 1, alpha = 1) +
  scale_fill_viridis_c(option = "C") +
  labs(y = "", x = "Reference coverage (TransRate)") +
  theme_bw(base_family = "GillSans", base_size = 12) + 
  theme(legend.position = "none", 
    # axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
    axis.title.x = element_text(size = 7)) +
  scale_x_continuous(position = "top", breaks = c(0.8, 0.9, 1))

# as there was a strong monotonic relationship between contig accuracy and TransRate contig score (detonate score), use this column to represent accuracy of the assembly contig


f <- list.files(dir, "_bam_info.csv", full.names = T)


baminfo_df <- read_csv(f) %>% dplyr::rename("gene_id" = "name")

names(baminfo_df)

baminfo_df <- baminfo_df %>% left_join(DB) %>% drop_na(tab) %>% filter(Signalp_class == "SP")

baminfo_df %>% 
  mutate(facet = ifelse(is.na(Conflict), "A) No redundant", "B) Redundant (conflict)")) %>%
  ggplot(aes(y = Superfamily, x = bases_uncovered, fill = after_stat(x))) +
  # geom_violin() +
  facet_grid(~ facet) +
  ggridges::geom_density_ridges_gradient(
    jittered_points = T,
    position = ggridges::position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 3, point_alpha = 1, alpha = 1) +
  scale_fill_viridis_c(option = "C") +
  labs(y = "", x = "TransRate contig score") +
  theme_bw(base_family = "GillSans", base_size = 14) + theme(legend.position = "none")

