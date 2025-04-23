# Conotoxin composition diversity as Phuong et al., 2016:
# N mature conotoxin
# N gene superfamilies
# N cys frameworks

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

time_levs <- c("Ctrl", "2", "4", "6")
recode_time <- structure(c("Control", "2 months", "4 months", "6 months"), names = time_levs)

Diatery_levs <- c("Ctrl","Cam", "Lit", "Pol", "Mix")
recode_Diatery <- structure(c("Control","Shrimp", "Mollusk", "Polychaete", "Mixed"), names = Diatery_levs)

# recode_Diatery <- structure(c("CT","SD", "LD", "PD", "MD"), names = Diatery_levs)


library(tidyverse)

pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR/"

# LOAD data -----

DB <- read_rds(paste0(pub_dir, "/structured_db.rds"))

dir <- "/Users/cigom/Documents/GitHub/conopeptides/06.Quantification/MATRIX_RSEM_dir/"

RES <- read_rds(paste0(dir, "/glmLRT_multiple_contrast_ctrl_and_treatments.rds")) %>%  filter(FDR < 0.05 & abs(logFC) > 2)

.count_vst <- read_rds(paste0(dir, "/counts_vst_nt_raw.rds"))$vst

count_vst <- .count_vst

.colData <- list.files(dir, pattern = "Manifest", full.names = T) 

.colData <- read_tsv(.colData) %>% 
  mutate(Diatery = ifelse(grepl("cam2v_", LIBRARY_ID), "camv", Diatery)) %>%
  mutate(design = ifelse(!Time == "Ctrl", paste0(Diatery, "_", Time), "Ctrl")) %>%
  select(LIBRARY_ID, Time, Diatery, design) %>%
  mutate_if(is.character, as.factor)


# 
# As ORFs include exclusively peptides with start and stop codon. Considere include all classes of regions

DB %>% drop_na(tab) %>% filter(Signalp_class == "SP") -> CONOPEPDB

CONOPEPDB %>% drop_na(tab) %>% filter(Signalp_class == "SP") %>% count(Region, sort = T)

# recode_peptide <- structure(c("Mature","Propeptide", "Precursor",), names = Which_regions)

# CONOPEPDB <- CONOPEPDB %>% mutate(Region = recode_factor(Region, !!!recode_peptide, .ordered = T)) 

CONOPEPDB %>% distinct(pep_seq)
CONOPEPDB %>% distinct(gene_id)
CONOPEPDB %>% count(Superfamily, sort = T)

CONOPEPDB %>% filter(!is.na(Conflict)) %>% count(tab, Superfamily, sort = T)
CONOPEPDB %>% filter(is.na(Conflict)) %>% count(tab, Superfamily, sort = T) 

dim(count_vst <- count_vst[rownames(count_vst) %in% query_genes,])

# agglomerate gene_matrix by same superfamily (not!! because gene profile going to be masked, just use for previz purpose )

barvizA <- CONOPEPDB %>% 
  count(Superfamily, tab, sort = T)

barvizB <- .count_vst %>%
  as_tibble(rownames = "gene_id") %>%
  left_join(distinct(CONOPEPDB, Superfamily, gene_id, tab)) %>%
  group_by(Superfamily, tab) %>%
  summarise_at(vars(all_of(colnames(.count_vst))), sum) %>% ungroup() 


plotdf <- barvizB %>% 
  pivot_longer(cols = all_of(colnames(.count_vst)), values_to = 'fill', names_to = "LIBRARY_ID") %>%
  left_join(.colData, by = "LIBRARY_ID") %>%
  group_by(Superfamily,tab) %>%
  summarise(Treads = sum(fill), Mean = mean(fill), sd = sd(fill)) %>%
  left_join(barvizA) %>%
  drop_na(Superfamily) %>%
  mutate(x = Treads) %>% ungroup() %>%
  arrange(x) %>% mutate(Superfamily = factor(Superfamily, levels = unique(Superfamily))) %>%
  mutate(label = paste0("(", n,")"))

p <- plotdf %>% 
  ggplot(aes(y = Superfamily, x = x)) +
  # facet_grid(~ tab, scales = "free_x", switch = "y") +
  geom_col(aes(fill = tab)) +
  scale_x_continuous(labels = scales::comma_format(scale = 0.001, suffix = "M")) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  scale_fill_grey("") +
  theme(legend.position = "top", 
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.12, "cm"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    # panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = 'grey95', color = 'white')) +
  geom_text(aes(label = label), size = 1.5,
    hjust = -0.1, vjust = 0, 
    family = "GillSans", position = position_dodge(width = 1)) 

p


