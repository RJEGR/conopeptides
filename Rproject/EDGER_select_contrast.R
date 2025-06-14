
# This is a copy from EDGERViz.R 
# Groups direction: positive logFC == sampleA & negative logFC == sampleB
# selecting groups of contrasts as Edith suggest
# Because Quantification performed at CDS level, including protein_ids with identical CDS, lets to collapse DEG results based on the CDS sequence. This is posible as redundancy spread to identical CDS having identical expression patterns. 


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

# time_levs <- c("Ctrl", "2", "4", "6")
# recode_time <- structure(c("Control", "2 months", "4 months", "6 months"), names = time_levs)

# Diatery_levs <- c("Ctrl","Cam", "Lit", "Pol", "Mix")
# recode_Diatery <- structure(c("Control","Shrimp", "Mollusk", "Polychaete", "Mixed"), names = Diatery_levs)

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

DB %>%  dplyr::count(Signalp_class, prediction_tool)

CONOPEPDB <- DB %>% 
  # filter(Signalp_class == "SP") %>%
  # filter(prediction_tool == "BOTH") %>%
  drop_na(prediction_tool) 

str(query_genes <- CONOPEPDB %>% distinct(protein_id) %>% pull()) # 12136 putative conopeptide genes

sum(query_genes %in% unique(RES$protein_id)) # 7624/12136 as EDGE.R remove low expressed transcripts (freq > 1 exp > 1)

RES <- RES %>% filter(protein_id %in% query_genes) 

nrow(RES %>% distinct(protein_id)) # 7624 putative conopeptides (not DEGs filtered yet)

RES %>%
  ggplot(aes(FDR)) + 
  # facet_wrap(sampleA ~ sampleB) +
  geom_histogram()

# Previz global number of conotoxin degs by contrast group =====

levs <- c("Ctrl","Cam_2","Camv_2", "Cam_4", "Lit_2", "Lit_4","Pol_2", "Pol_4", "Mix_2", "Mix_4")

recode_to <- structure(c("Control","Shrimp", "Shrimp","Shrimp","Mollusk", "Mollusk","Polychaete", "Polychaete","Mixed","Mixed"), names = levs)

RES %>% dplyr::distinct(sampleB)

# positive logFC == sampleA & negative logFC == sampleB

DataVizdf <- RES %>%
  # If collapse DEGS to CDS Level use:
  # CONOPEPDB %>% distinct(protein_id, dna_seq, pep_seq) %>% right_join(RES) %>%  select(-protein_id) %>% distinct() %>%
  # Currently are only FDR < 0.05
  filter(FDR < 0.05 & abs(logFC) > 2 ) %>%
  mutate(facet = ifelse( sign(logFC) == 1, "up in sampleA", "up in sampleB")) %>%
  dplyr::count(sampleA, sampleB, facet, sort = T) 


DataVizdf %>%
  # dplyr::mutate(facet = dplyr::recode_factor(sampleA, !!!recode_to)) %>%
  ggplot(aes(y = sampleA, x = sampleB, fill = n)) +
  facet_grid(~ facet, scales = "free") +
  geom_tile(color = 'white', linewidth = 0.5) +
  geom_text(aes(label = n), size = 3, family = "GillSans", color = "white") +
  theme_bw(base_family = "GillSans", base_size = 10) +
  labs(subtitle = "DEGS: FDR < 0.05 & abs(logFC) > 2 & Signalp_class == SP") +
  theme(
    legend.position = "none",
    # panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    # panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    # panel.grid.major.x = element_blank(),
    # axis.text.y.right = element_text(angle = 0, hjust = 1, vjust = 0, size = 2.5),
    axis.text.y = element_text(angle = 0, size = 7),
    axis.text.x = element_text(angle = 0, size = 7),
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 0, size = 10)) -> P

P <- P +
  geom_segment(
     data = filter(DataVizdf, facet == "up in sampleA"), 
      x = 8, xend = 1, 
      y = 10, yend = 10, 
      colour = "gray7", 
      arrow = arrow(ends = "last", length = unit(0.15, "cm"))) +
  annotate("text", x = 4, y = 10.3, size = 3, label = "Contrast sence",  color = "gray7", family = "GillSans") +
  geom_segment(
    data = filter(DataVizdf, facet == "up in sampleB"), 
    x = 8, xend = 1, 
    y = 10, yend = 10, 
    colour = "gray7", 
    arrow = arrow(ends = "first", length = unit(0.15, "cm")))

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

nrow(DataViz %>% distinct(protein_id)) # 6086 putative conopeptides presented in the selected contrast (not log2FC and SIgnalP filtered yet)

write_rds(DataViz, file = paste0(file.path(dir, subdir), "/cds_exactTest_multiple_contrast_ctrl_and_treatments.rds"))

DataViz <- CONOPEPDB %>% distinct(protein_id, dna_seq, pep_seq) %>% right_join(DataViz) %>%  select(-protein_id) %>% distinct()

write_rds(DataViz, file = paste0(file.path(dir, subdir), "/cds_exactTest_multiple_contrast_ctrl_and_treatments_cds_level.rds"))


frames_df <- DataViz %>%
  filter(FDR < 0.05 & abs(logFC) > 2 ) %>%
  mutate(facet = ifelse( sign(logFC) == 1, "up in sampleA", "up in sampleB")) %>%
  dplyr::count(sampleA, sampleB, facet, sort = T) 


  
P <- P + geom_tile(data=frames_df, color="orange", fill = NA, linewidth = 1)

ggsave(P, filename = 'Marginals_degs_cds_level.png', 
  path = pub_dir, width = 7, height = 3, device = png, dpi = 600)

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

