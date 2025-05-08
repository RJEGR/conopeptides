
# This is a copy from EDGERViz.R 
# Groups direction: positive logFC == sampleA & negative logFC == sampleB
# selecting groups of contrasts as Edith suggest

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

DB <- read_tsv(paste0(pub_dir, "/conopeptides.tsv"))

dir <- "/Users/cigom/Documents/GitHub/conopeptides/06.Quantification/MATRIX_RSEM_dir/"

subdir <- "KALLISTO_Merged_polyA_hisat_SuperDuper.fasta.transdecoder_DIR/"

f <- "p05_exactTest_multiple_contrast.rds"

f <- list.files(file.path(dir, subdir), f, full.names = T)

RES <- 
  read_rds(f) %>%
  dplyr::rename("protein_id" = "ids") %>%
  mutate(sign = sign(logFC)) %>%
  mutate(sampleX = ifelse(sign == 1, sampleA, sampleB))

# In sampple X we can know the independent/cumulative number of DEGs for each group

RES %>% dplyr::count(sampleX, sort = T)

RES %>% dplyr::count(sampleA, sampleB, sort = T) 

# count_vst <- read_rds(paste0(dir, "/counts_vst_nt_raw.rds"))$vst

.colData <- list.files(dir, pattern = "Manifest", full.names = T) 

.colData <- read_tsv(.colData) %>% 
  mutate(Diatery = ifelse(grepl("cam2v_", LIBRARY_ID), "camv", Diatery)) %>%
  mutate(design = ifelse(!Time == "Ctrl", paste0(Diatery, "_", Time), "Ctrl")) %>%
  select(LIBRARY_ID, Time, Diatery, design) %>%
  mutate_if(is.character, as.factor)

# Match putative conopeptides (Precursor and Pro-peptide and mature)

DB %>% 
  # filter(Signalp_class == "SP") %>%
  dplyr::count(Signalp_class, prediction_tool)

CONOPEPDB <- DB %>% 
  # filter(Signalp_class == "SP") %>%
  drop_na(prediction_tool) 

str(query_genes <- CONOPEPDB %>% distinct(protein_id) %>% pull()) # 12136 putative conopeptide genes

sum(query_genes %in% unique(RES$protein_id)) # 7624/12136 as EDGE.R remove low expressed transcripts (freq > 1 exp > 1)

RES <- RES %>% filter(protein_id %in% query_genes) 

nrow(RES %>% distinct(protein_id)) # 7624 putative conopeptides (not DEGs filtered yet)

RES %>%
  ggplot(aes(FDR)) + 
  # facet_wrap(sampleA ~ sampleB) +
  geom_histogram()

# Select multiple contrast of interest
# shrimp ----
# sampleA: Shrimp (2, 4,6 months)
# sampleB: Shrimp (4,6, months and Ctrl)


omitpattern <- c('camv|Pol|Lit|Mix|Cam_6$')


Shrimpdf <- RES %>% 
  filter(!grepl(omitpattern, sampleA)) %>%
  filter(!grepl(omitpattern, sampleB)) %>%
  mutate(sam_group = "Shrimp")
# filter(grepl("Cam_[0-9]$", sampleA)) %>%
# filter(grepl("Cam_[0-9]$|Ctrl", sampleB)) 

Shrimpdf %>%
  dplyr::count(sampleA, sampleB)


# Polychaete ----
# sampleA: Polychaete (2, 4 months)
# sampleB: Polychaete (4,6, months and Ctrl)

omitpattern <- c('camv|Cam|Lit|Mix')

polypdf <- RES %>% 
  filter(!grepl(omitpattern, sampleA)) %>%
  filter(!grepl(omitpattern, sampleB)) %>%
  mutate(sam_group = "Polychaete")
# filter(if_any(where(is.character), ~ grepl(pattern = 'Pol_[0-9]$', x = .x, ignore.case = T))) %>%
# filter(if_any(where(is.character), ~ grepl(pattern = 'Pol_[0-9]$|Ctrl', x = .x, ignore.case = T)))

polypdf %>%
  dplyr::count(sampleA, sampleB)

# Mollusk ----
# sampleA: Lit (2, 4 months)
# sampleB: Lit (4,6, months and Ctrl)


omitpattern <- c('camv|Cam|Pol|Mix')

Litdf <- RES %>% 
  filter(!grepl(omitpattern, sampleA)) %>%
  filter(!grepl(omitpattern, sampleB)) %>%
  mutate(sam_group = "Mollusk")



Litdf %>%
  dplyr::count(sampleA, sampleB)

# Mixed
# sampleA: Mix (2, 4 months)
# sampleB: Mix (4,6, months and Ctrl)


omitpattern <- c('camv|Cam|Pol|Lit')

Mixdf <- RES %>% 
  filter(!grepl(omitpattern, sampleA)) %>%
  filter(!grepl(omitpattern, sampleB)) %>%
  mutate(sam_group = "Mixed")


Mixdf %>%
  dplyr::count(sampleA, sampleB)


# Dataviz
DataViz <- rbind(
  Shrimpdf,
  Litdf,
  polypdf,
  Mixdf)


# Omit by now DEGs enriched in Ctrl (ie sampleX != "Ctrl)
DataViz <- DataViz %>% filter(sampleX != "Ctrl")


# PLOT DEGS (summary) -----

DataViz %>%
  ggplot(aes(FDR)) + 
  geom_histogram()

nrow(DataViz %>% distinct(protein_id)) # 6086 putative conopeptides presented in the selected contrast (not DEGs filtered yet)

write_rds(DataViz, file = paste0(file.path(dir, subdir), "/cds_exactTest_multiple_contrast_ctrl_and_treatments.rds"))


# Exit ------

DataViz <- polypdf


# global view of changes in expression :

DataViz %>% 
  filter(FDR < 0.05) %>%
  dplyr::count(sam_group, sampleA, sampleB, sampleX) %>%
  mutate(sampleX = ifelse(sampleX == sampleA, paste0(sampleA, " (", sampleB,")"), paste0(sampleB, " (", sampleA,")")))


# DataViz %>% dplyr::count(sam_group, sampleA, sampleB, sampleX)  %>%
#   separate(sampleX, into = c("Diatery", "Time"), sep = "_", remove = T) %>%
#   mutate(Time = ifelse(is.na(Time), "0 months", Time)) %>%
#   mutate(Time = recode_factor(Time, !!!recode_time, .ordered = T)) %>%
#   mutate(Diatery = recode_factor(Diatery, !!!recode_Diatery, .ordered = T)) %>%
#   ggplot(aes(sampleA, sampleB, fill = n)) + 
#   facet_grid(~sam_group, scales = "free_x") + geom_tile() + geom_text(aes(label =n))


# what are the number of intersected?


# .UPSETDF <- DataViz %>% 
#   mutate(DEG = sign(log2FoldChange)) %>%
#   # dplyr::mutate(sampleB = dplyr::recode_factor(sampleB, !!!recode_to)) %>%
#   group_by(transcript_id, DEG) %>%
#   summarise(across(sampleB, .fns = list), n = n())

DataViz %>% 
  filter(FDR < 0.05) %>%
  select(protein_id, sam_group, sampleA, sampleB) %>%
  pivot_longer(cols = c("sampleA", "sampleB"), names_to = "summarise_group", values_to = "summarise_col") %>%
  filter(summarise_col != "Ctrl") %>%
  dplyr::count(sam_group, summarise_group, summarise_col)


UPSETDFA <- DataViz %>% 
  filter(FDR < 0.05) %>%
  select(gene_id, sam_group, sampleA, sampleB) %>%
  pivot_longer(cols = c("sampleA", "sampleB"), names_to = "summarise_group", values_to = "summarise_col") %>%
  filter(summarise_col != "Ctrl") %>%
  # filter(logFC > 0 ) %>% # selecting only sampleA groups
  # mutate(summarise_col = sampleA) %>%
  # dplyr::mutate(sampleB = dplyr::recode_factor(sampleB, !!!recode_to)) %>%
  group_by(gene_id, sam_group, summarise_group) %>%
  summarise(across(summarise_col, .fns = list), n = n())

UPSETDFB <- DataViz %>% 
  filter(logFC < 0 ) %>% # selecting only sampleA groups
  mutate(summarise_col = sampleB) %>%
  # dplyr::mutate(sampleB = dplyr::recode_factor(sampleB, !!!recode_to)) %>%
  group_by(gene_id) %>%
  summarise(across(summarise_col, .fns = list), n = n())

UPSETDF <- rbind(UPSETDFA, UPSETDFB)

library(ggupset)

UPSETDFA %>%
  # mutate(col = ifelse(n == 1, "A", "B")) %>%
  ggplot(aes(x = summarise_col)) + # , fill = SIGN fill = as.factor(DEG))
  geom_bar(position = position_dodge(width = 1)) +
  geom_text(stat='count', aes(label = after_stat(count)), 
    position = position_dodge(width = 1), vjust = -0.5, family = "GillSans", size = 3.5) +
  scale_x_upset(order_by = "degree", reverse = F) +
  theme_combmatrix(
    combmatrix.panel.point.color.fill = "black",
    combmatrix.label.make_space = F,
    # combmatrix.panel.point.color.fill = panel.point.color.fill,
    combmatrix.panel.line.size = NA, 
    base_family = "GillSans", base_size = 16) 
# axis_combmatrix(levels = recode_to) +
# labs(x = '', y = 'Number of transcripts') +
# # scale_color_manual("", values = col) +
# scale_fill_manual("", values =  c("gray20", "grey70"), 
#   labels = c("Up-expressed","Down-expressed")) +
# guides(fill = guide_legend(title = "", nrow = 1)) 

# .... continue =====

nrow(DataViz %>% distinct(gene_id)) # 729 conopeptide DEGs (FDR < 0.05 & abs( logFC) > 2, from the filter

# write_rds(DataViz, file = paste0(dir, "/glmLRT_multiple_contrast_ctrl_and_treatments.rds"))

# filter to only up under  Diatery (ie exclude sampleX == Ctrl)

DataVizTop <- DataViz %>% filter(sampleX != "Ctrl") 



# Create a upset format



# DataViz %>%
#   ggplot(aes(FDR)) + 
#   facet_wrap(sampleA ~ sampleB) +
#   geom_histogram()

# To recode 

RES %>% 
  select(gene_id, sampleA, sampleB) %>%
  separate(sampleA, into = c("Diatery", "Time"), sep = "_", remove = T) %>%
  mutate(Time = recode_factor(Time, !!!recode_time, .ordered = T)) %>%
  mutate(Diatery = recode_factor(Diatery, !!!recode_Diatery, .ordered = T)) 
# mutate(sampleA = recode_factor(Diatery, !!!recode_Diatery, .ordered = T)) 

