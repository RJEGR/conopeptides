# 
# Dataviz EDGE.R results
# Se barplot of conopeptides from the multiple contrast ( )
# run upsetR data and plot


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)


library(tidyverse)

pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR/"

# LOAD data -----

DB <- read_rds(paste0(pub_dir, "/structured_db.rds"))

dir <- "/Users/cigom/Documents/GitHub/conopeptides/06.Quantification/MATRIX_RSEM_dir/"


RES <- read_rds(paste0(dir, "/glmLRT_multiple_contrast.rds")) %>% dplyr::rename("gene_id" = "ids")


# Match only conopeptides

CONOPEPDB <- DB %>% drop_na(tab) %>% 
  # dplyr::count(Signalp_class)
  filter(Signalp_class == "SP")

CONOPEPDB %>% distinct(Superfamily)

query_genes <- CONOPEPDB %>% distinct(gene_id) %>% pull()

RES <- RES %>% filter(gene_id %in% query_genes)

RES %>% dplyr::count(sampleA, sampleB) 


RES %>%
  ggplot(aes(PValue)) + 
  facet_wrap(sampleA ~ sampleB) +
  geom_histogram()

