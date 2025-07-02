
# Use logFC > 4 for more astringency
# Step0 select up in control (Note Cam_2 and Cam_4 are invert to control in sampleA/B)
# Step1 Select up putative conotoxins in time (2 and 4) vs control
# Step2 Select up putative conotoxin dfferent between time 2 and 4 per diet
# Step3 Contrast temporalidad and diet specificity using venn diagram
# Step 4 Plot a facet of Metatranscriptome Sf, overal and exclusive DEGs per diet group

c("#146179", "#09BC9F", "#FEB65F", "#C55E2D")

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)


library(tidyverse)

pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR/"

# LOAD data -----

DB <- read_tsv(paste0(pub_dir, "/conopeptides.tsv"))

# Filter transcripts by TRUE conotoxin

CONOPEPDB <- DB %>%
  filter(Signalp_class == "SP" & contig_impact_score > 0 & prediction_tool == "BOTH" & nchar(pep_seq) < 200) %>%
  mutate(uniprotkb_toxprot = sapply(strsplit(uniprotkb_toxprot, "[|]"), `[`, 2)) %>%
  mutate(uniprotkb_toxprot = ifelse(is.na(uniprotkb_toxprot), "uID", uniprotkb_toxprot)) %>%
  mutate(conoserver_protein = ifelse(is.na(conoserver_protein), "uID", conoserver_protein)) %>%
  mutate(conoserver_protein = ifelse(is.na(conoserver_protein), "uID", conoserver_protein)) %>%
  mutate(hmm_pred_conodictor = stringr::str_to_sentence(hmm_pred_conodictor)) %>%
  mutate(hmm_pred_conodictor = ifelse(is.na(hmm_pred_conodictor), "uID", hmm_pred_conodictor)) %>%
  mutate(Superfamily = ifelse(is.na(Superfamily), "uID", Superfamily)) %>%
  mutate(prediction_tool = ifelse(tab %in% "pHMM", paste0(prediction_tool,"_",tab), prediction_tool)) %>%
  mutate(sf = Superfamily) %>%
  select(protein_id, pep_seq, sf, Signalp_class, prediction_tool, hmm_pred_conodictor, Superfamily, uniprotkb_toxprot, conoserver_protein) %>%
  unite("y_axis", Signalp_class:conoserver_protein, sep = "|") 


CONOPEPDB %>% drop_na()

dir <- "/Users/cigom/Documents/GitHub/conopeptides/06.Quantification/MATRIX_RSEM_dir/"

subdir <- "KALLISTO_Merged_polyA_hisat_SuperDuper.fasta.transdecoder_DIR/"

f <- "cds_exactTest_multiple_contrast_ctrl_and_treatments.rds"

# f <- "cds_exactTest_multiple_contrast_ctrl_and_treatments_cds_level.rds"

f <- list.files(file.path(dir, subdir), f, full.names = T)

RES <- read_rds(f)

RES %>% dplyr::count(sampleX, sort = T)

# As duplicated CDS have same expression patters and so on, same gene-pairwise values in DE analysis, 
# just deduplicate redundancy omiting protein_id but using pep_seq column

CONOPEPDB %>% count(pep_seq, sort = T) # If there is not any filter, the pep_seq will be duplicated, otherwise, if filtering are stringency not duplicates will be found

# deduplicate to cds_level/pep_level using CONOPEPDB and join to RES

nrow(RES)

# number of intersected
# Only fraction low 100 % must be found between unique(RES$protein_id)/unique(CONOPEPDB$protein_id)

sum(unique(RES$protein_id) %in% unique(CONOPEPDB$protein_id)) 

nrow(RES <- RES %>% right_join(CONOPEPDB) %>% select(-protein_id) %>% distinct())

RES %>% right_join(CONOPEPDB) %>% 
  write_rds(file = file.path(pub_dir, "EDGER_EffectSize_input.rds"))

# Separate by now DEGs enriched in Ctrl (ie sampleX != "Ctrl)

Controldf <- RES %>% filter(sampleX == "Ctrl")  %>% filter(abs(logFC) > 4 & FDR < 0.05) 

Controldf %>% dplyr::count(sampleA, sampleB, sampleX, sort = T)

# The Controldf contain the 'basal' conotoxin found in empiriral/natural conditions

# In addition  filter significat DEGS

RES %>% ggplot(aes(logFC)) + geom_histogram() + facet_grid(sam_group ~.)

DataViz <- RES %>% filter(sampleX != "Ctrl") %>% filter(abs(logFC) > 4 & FDR < 0.05) 

DataViz %>% ggplot(aes(logFC)) + geom_histogram() + facet_grid(sam_group ~.)

DataViz %>% dplyr::count(sampleA, sampleB, sampleX)


# Step1: Contrasting results against control (ie. omit contrast dietA_time1 vs dietA_time2)

cols_to_check <- c("sampleA", "sampleB")

DataViz <- DataViz %>%
  filter(
    if_any(all_of(cols_to_check), ~ str_detect(.x, "Ctrl"))) 


DataViz %>% dplyr::count(sam_group, sampleA, sampleB, sampleX)

DataVizTop <- 
  DataViz %>% 
  #dplyr::mutate(sampleB = dplyr::recode_factor(sampleB, !!!recode_to)) %>%
  group_by(sam_group, sampleA, sampleB, sampleX) %>%
  arrange(logFC, .by_group = T) %>%
  group_by(sampleB, y_axis) %>%
  filter(abs(logFC) == max(abs(logFC))) %>%
  group_by(sampleX) %>%
  slice_head(n = 10)


DataVizTop %>%
  dplyr::count(sam_group, sampleA, sampleB, sampleX)

DataVizTop %>% ungroup() %>% distinct(y_axis)

# Split indistinc from unique, for labeling color

q <- DataVizTop %>% distinct(y_axis) %>% ungroup() %>% dplyr::count(y_axis, sort = T) %>% filter(n == 1) %>% pull(y_axis)

DataVizTop <- DataVizTop %>% mutate(col = ifelse(y_axis %in% q, "grey70", "grey89"))

col_vals <- c("grey70", "grey89")

col_vals <- structure(col_vals, names = col_vals)

DataVizTop %>% ungroup() %>% dplyr::count(col, sort = T) 

"log2fold-change (vst)"

x_label <- expression(log[2]~"fold-change")


gene_names <- DataVizTop %>% pull(y_axis, name = protein_id)

DataVizTop %>% 
  mutate(facet = paste0("Ctrl", "_vs_", sampleX)) %>%
  mutate(Label = y_axis, row_number = row_number(Label)) %>%
  #mutate(Label = paste0(Label, " (", protein_name, ")")) %>%
  mutate(row_number = paste(row_number, sampleX, sep = ":")) %>%
  mutate(Label = factor(paste(Label, row_number, sep = "__"),
    levels = rev(paste(Label, row_number, sep = "__")))) %>%
  mutate(star = ifelse(FDR <.001, "***", 
    ifelse(FDR <.01, "**",
      ifelse(FDR <.05, "*", "")))) %>%
  mutate(logFC = abs(logFC)) %>%
  mutate(
    # ymin = (abs(logFC) - lfcSE) * sign(logFC),
    # ymax = (abs(logFC) + lfcSE) * sign(logFC),
    y_star = logFC + (0.15)* sign(logFC)) %>% 
  ggplot(aes(y = Label, x = logFC)) + 
  # facet_grid(facet ~. , scales = "free_y", space = "free")+
  ggforce::facet_col(facet ~. , scales = "free_y", space = "free") +
  scale_y_discrete(labels = function(x) gsub("__.+$", "", x)) +
  scale_fill_manual(values = col_vals) +
  geom_col(aes(fill = col), width = 0.5, size = 0.25,
    position = position_stack(reverse = T), color = "white") +
  # geom_errorbar(aes(xmin = ymin, xmax = ymax), width = 0.05,
  # position = position_identity(), color = "black") +
  geom_text(aes(x = y_star, label=star), 
    # hjust = .7 ,
    # vjust=  .7,  
    color="black", position = position_identity(), 
    family = "GillSans", size = 1.5)+
  labs(y = NULL, x = x_label, title = "Up-expressed conotoxins under different diets") +
  theme_bw(base_family = "GillSans", base_size = 10) +
  theme(
    legend.position = "none",
    # panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    # axis.text.y.right = element_text(angle = 0, hjust = 1, vjust = 0, size = 2.5),
    axis.text.y = element_text(angle = 0, size = 3),
    axis.text.x = element_text(angle = 0),
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 1, size = 7)) -> P


P

ggsave(P, filename = 'EDGERLOG2FCTOP10.png', 
  path = pub_dir, width = 5, height = 7, device = png, dpi = 600)

# Filter Control by intersected down-expressed in all diets experiments
# Find the upexpressed degs from Ctrl frequently co-occurring in diets

ControlDEGS <- Controldf %>% 
  mutate(Time = ifelse(grepl("_2$", sampleB), "2 months", "4 months")) %>% 
  dplyr::count(Time, pep_seq, sort = T) %>% filter(n == 4) %>% 
  distinct(pep_seq) %>%
  left_join(DB) %>%
  distinct(protein_id, pep_seq)

ExperimentalDEGS <- DataViz %>% distinct(pep_seq) %>% left_join(DB) %>% distinct(protein_id, pep_seq)


queries_for_zscore_heatmap <- rbind(ControlDEGS, ExperimentalDEGS) %>% distinct(protein_id) %>% pull()

# omit -----

myXStringSet <- DataVizTop %>% 
  ungroup() %>% filter(sam_group == "Polychaete") %>% 
  distinct(y_axis, pep_seq) %>% 
  # mutate(pep_seq = gsub("[*]$","", pep_seq)) %>%
  pull(pep_seq, name = y_axis)


myXStringSet <- Biostrings::AAStringSet(c(myXStringSet))


library(msa)

align <- msa::msa(myXStringSet, method = "ClustalW", order = "input")

.align <- msa::msaConvert(align)$seq

names(.align) <- msa::msaConvert(align)$nam

# data(BLOSUM62)
# msaConservationScore(align, BLOSUM62)

library(ggsci)

pat <- c("-", alphabet(myXStringSet, baseOnly=TRUE))

# colors <- structure(pal_aaas()(length(pat)), names = rev(pat))

# scales::show_col(colors)

DECIPHER::BrowseSeqs(AAStringSet(.align), colWidth = 120)


# Ven diagram or upset -----

# Q: Those upexpressed are unique from time or expression increase by time?
# 

# Filter

# If NA in DataViz, is because not preserved in CONOPEPDB (filtered good assembled conotoxins)

DataViz <- DataViz %>% drop_na(sam_group)

DataViz %>% 
  dplyr::count(sam_group, sampleA, sampleB, sampleX)


DataViz %>% 
  select(pep_seq, sam_group, sampleX) %>%
  mutate(sampleX = ifelse(grepl("_2$", sampleX), "2 months", "4 months")) %>%
  # filter(grepl("_2$", sampleX)) %>%
  distinct() %>%
  dplyr::count(sam_group, sampleX)

# Reformat
# 
UPSETDF <- DataViz %>% 
  select(pep_seq, sam_group, sampleX) %>%
  mutate(sampleX = ifelse(grepl("_2$", sampleX), "2 months", "4 months")) %>%
  distinct() %>%
  # filter(grepl("_2$", sampleX)) %>%
  mutate(summarise_col = sam_group) %>% 
  # dplyr::mutate(sampleB = dplyr::recode_factor(sampleB, !!!recode_to)) %>%
  group_by(pep_seq, sampleX) %>%
  summarise(across(summarise_col, .fns = list), n = n())


# Plot upset
recode_to <- structure(c("Shrimp", "Mollusk", "Polychaete", "Mixed"))

DataViz %>% 
  select(pep_seq, sam_group, sampleX) %>%
  mutate(sampleX = ifelse(grepl("_2$", sampleX), "2 months", "4 months")) %>%
  distinct() %>%
  mutate(summarise_col = sam_group) %>% 
  count(sampleX, summarise_col) %>%
  mutate(summarise_col = factor(summarise_col, levels = rev(recode_to))) %>%
  ggplot(aes(y = summarise_col, x = n)) + 
  geom_col(position = position_stack(reverse = T), fill = "black") +
  facet_grid(~ sampleX, scales = "free_x", space = "free_x") +
  # geom_text(aes(label = summarise_col), color = "white")
  theme_bw(base_family = "GillSans", base_size = 7) +
  labs(x = "Set size (Number of putative conotoxins)", y = "Diet") +
  scale_fill_manual("", values = c("#F3E0F7","#63589F")) +
  theme(
    legend.position = "top",
    # panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    # axis.text.y.right = element_text(angle = 0, hjust = 1, vjust = 0, size = 2.5),
    axis.text.y = element_text(angle = 0, size = 10),
    axis.text.x = element_text(angle = 0),
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 1, size = 10)) -> P

# P

ggsave(P, filename = 'Intersections_setSize.png', 
  path = pub_dir, width = 3, height = 2, device = png, dpi = 600)


library(ggupset)

my_font <- "GillSans"

P <- UPSETDF %>%
  # filter(sampleX == "4 months") %>%
  ggplot(aes(x = summarise_col)) +
  facet_grid(sampleX ~ ., scales = "free_y") +
  geom_bar(fill = "black") +
  scale_y_reverse("Number of putative conotoxins") +
  geom_text(stat='count', aes(label = after_stat(count)), 
    position = position_dodge(width = 1.2), vjust = 1, family = "GillSans", size = 2) +
  scale_x_upset(order_by = "degree", reverse = F, position = "top") +
  labs(x = "Degree of intersections") +
  theme_bw() + 
  theme_combmatrix(
    # combmatrix.label.text = element_text(color = "blue", size=10),
    combmatrix.panel.point.color.fill = "black",
    combmatrix.label.make_space = F,
    # combmatrix.panel.point.color.fill = panel.point.color.fill,
    combmatrix.panel.point.size = 0.15,
    combmatrix.panel.line.size = 0.15,
    base_family = "GillSans", base_size = 12,
    strip.background = element_rect(fill = 'gray90', color = 'white'),
    strip.text = element_text(color = "black", size = 10, family = "GillSans",hjust = 0.5),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_text(family = my_font),
    axis.title.y = element_text(family = my_font),
    axis.text.x = element_text(family = my_font),
    axis.text.y = element_text(family = my_font)) +
  axis_combmatrix(levels = recode_to) 

# P +   scale_y_reverse("Number of conotoxins", breaks = seq(0,40, by = 10)) 

ggsave(P, filename = 'Intersections.png', 
  path = pub_dir, width = 4, height = 4, device = png, dpi = 600)

# 


# Facets ======
# Filter conotoxin 

# CONOPEPDB <- DB %>%
  # drop_na(Superfamily) %>%
  # filter(Signalp_class == "SP" & contig_impact_score > 0 & prediction_tool == "BOTH") %>%
  # select(pep_seq, hmm_pred_conodictor, Superfamily, prediction_tool, tab) %>%
  # mutate(hmm_pred_conodictor = stringr::str_to_sentence(hmm_pred_conodictor)) %>%
  # mutate(hmm_pred_conodictor = ifelse(is.na(hmm_pred_conodictor), "uID", hmm_pred_conodictor)) %>%
  # mutate(Superfamily = ifelse(is.na(Superfamily), "uID", Superfamily)) %>%
  # mutate(prediction_tool = ifelse(tab %in% "pHMM", paste0(prediction_tool,"_",tab), prediction_tool)) %>%
  # unite("y_axis", hmm_pred_conodictor:tab, sep = "|") %>%
  # mutate(y_axis = Superfamily)

# Global metatrascriptome
# Here we can see number of allelic variants for different conotoxins families

CONOPEPDB <- CONOPEPDB %>% mutate(y_axis = sf)

Globaldf <- CONOPEPDB %>%
  dplyr::count(y_axis, sort = T) %>% drop_na() %>%
  mutate(y_axis = factor(y_axis, levels = rev(unique(y_axis)))) %>%
  mutate(facet = "A) Global diversity")

Globaldf %>%
  ggplot(aes(y = y_axis, x = n)) + geom_col() +
  labs(x = "Number of transcripts")


# 1) Split global from exclusive DEGs

# Add or omit sampleX if want to facet 2 from 4 moths ()

overall_degs_df <- UPSETDF %>% 
  ungroup() %>%
  unnest(summarise_col) %>%
  left_join(CONOPEPDB) %>% 
  dplyr::count(sampleX, summarise_col, y_axis, sort = T) %>% drop_na() %>%
  mutate(facet = "A) Overall DEGs")
  
unique_degs_df <- 
  UPSETDF %>% 
  ungroup() %>%
  filter(n == 1) %>% 
  unnest(summarise_col) %>%
  # Fix redundancy <----
  left_join(CONOPEPDB) %>%
  dplyr::count(sampleX, summarise_col, y_axis, sort = T) %>% drop_na() %>%
  mutate(facet = "B) Exclusive DEGs")

unique_degs_df %>%
  rbind(overall_degs_df) %>%
  mutate(y_axis = factor(y_axis, levels = levels(Globaldf$y_axis))) %>%
  mutate(summarise_col = factor(summarise_col, levels = recode_to)) %>%
  mutate(label = scales::comma(n)) %>%
  ggplot(aes(y = y_axis, x = summarise_col, fill = n)) + 
  # facet_grid(~ sampleX + summarise_col) +
  ggh4x::facet_nested(~ facet +sampleX, nest_line = T) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  geom_text(aes(label= label), hjust= 1, vjust = 0.5, size = 3, family = "GillSans", color = "white") +
  theme_bw(base_family = "GillSans", base_size = 10) +
  labs(x = "", y = "Superfamily") +
  theme(
    legend.position = "none",
    # panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    # axis.text.y.right = element_text(angle = 0, hjust = 1, vjust = 0, size = 2.5),
    axis.text.y = element_text(angle = 0, size = 12),
    axis.text.x = element_text(angle = 0),
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 1, size = 12)) -> P


ggsave(P, filename = 'Superfamilies_by_degs.png', 
  path = pub_dir, width = 10, height = 10, dpi = 500, device = png)

# By N reads (z-score) ====
# in addition to number of transcripts, summarise number of reads per family (or maybe zscore?)
# Caution!!!
# agglomerate gene_matrix by same superfamily going to mask allelic variation 

dir <- "/Users/cigom/Documents/GitHub/conopeptides/06.Quantification/MATRIX_RSEM_dir/"

subdir <- "KALLISTO_Merged_polyA_hisat_SuperDuper.fasta.transdecoder_DIR/"

f <- "KALLISTO_Merged_polyA_hisat_SuperDuper.fasta.transdecoder.matrix"

f <- list.files(file.path(dir, subdir), f, full.names = T)

dim(gene_matrix <- round(readRDS(f)))


.colData <- list.files(dir, pattern = "Manifest", full.names = T) 

recode_to <- structure(c("Control","Shrimp", "Mollusk", "Polychaete", "Mixed"))

recode_to <- structure(recode_to, names = c("Ctrl","Cam","Lit","Pol","Mix"))

recode_time <- structure(c("","2 month", "4 month"))
recode_time <- structure(recode_time, names = c("ctrl","2","4"))


.colData <- read_tsv(.colData) %>% 
  mutate(LIBRARY_ID = ifelse(grepl("cam2v_", LIBRARY_ID), NA, LIBRARY_ID)) %>%
  mutate(LIBRARY_ID = ifelse(grepl("cam6", LIBRARY_ID), NA, LIBRARY_ID)) %>%
  # mutate(LIBRARY_ID = ifelse(grepl("ctrl", LIBRARY_ID), NA, LIBRARY_ID)) %>%
  select(LIBRARY_ID, Time, Diatery) %>% drop_na(LIBRARY_ID) %>%
  dplyr::mutate(Diatery = dplyr::recode_factor(Diatery, !!!recode_to)) %>%
  dplyr::mutate(Time = dplyr::recode_factor(Time, !!!recode_time)) %>%
  mutate_if(is.character, as.factor)

overall_degs_df <- UPSETDF %>% 
  ungroup() %>%
  unnest(summarise_col)

# Matrix to zscore
# matrix to tibble
# filter degs
# pivot longer
# summarise by sample group
# join to sf (CONOPEPDB)
# 

z_scores <- function(x) {(x-mean(x))/sd(x)}

sum(keepRows <- rownames(gene_matrix) %in% queries_for_zscore_heatmap)

sum(keepCols <- colnames(gene_matrix) %in% .colData$LIBRARY_ID)

str(HeatmapViz <- t(apply(gene_matrix[keepRows,keepCols], 1, z_scores)))

# Normalize to vst---

HeatmapViz[is.na(HeatmapViz)] <- 0

h <- heatmap(HeatmapViz, keep.dendro = TRUE )

HeatmapViz <- HeatmapViz %>%
  as_tibble(rownames = "protein_id") %>%
  left_join(distinct(CONOPEPDB, protein_id, pep_seq, y_axis)) %>%
  pivot_longer(cols = colnames(HeatmapViz), names_to = "LIBRARY_ID", values_to = "n") %>%
  left_join(.colData) %>%
  # mutate(Time = ifelse(grepl("2", Time), "2 months", "4 months")) %>%
  filter(!is.na(n))

HeatmapViz<- HeatmapViz %>% filter(y_axis == "I1")

lo = floor(min(HeatmapViz$n))
up = ceiling(max(HeatmapViz$n))
mid = (lo + up)/2

P <- HeatmapViz %>%
  drop_na(pep_seq) %>%
  # sample_n(50) %>%
  # group_by(Diatery) %>% mutate(n = z_scores(n))
  ggplot(aes(y = pep_seq, x = Diatery, fill = n)) + 
  ggh4x::facet_nested(y_axis ~ Time, scales = "free",space = "free", nest_line = T,switch = "y") +
  geom_raster() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
    na.value = "white", midpoint = mid, limit = c(lo, up), 
    breaks = seq(lo, up, by = 3),
    name = NULL) +
  labs(y = "", x = "") +
  theme_classic(base_family = "GillSans", base_size = 12) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 7),
    # axis.text.y = element_blank(), 
    # axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1,size = 12),
    strip.text.y.left = element_text(
      angle = 0, hjust = 1,
      size = 5),
    strip.background = element_rect(colour = "transparent", fill = "transparent", size = 1)
  )

P

ggsave(P, filename = 'Conotoxins_degs_under_diet_treatment.png', 
  path = pub_dir, width = 5, height = 10, device = png, dpi = 600)


# Overal by superfam ====

# correlate n transcripts vs expresion



DF1 <- gene_matrix %>%
  as_tibble(rownames = "protein_id") %>%
  filter(protein_id %in% queries_for_zscore_heatmap) %>%
  left_join(distinct(CONOPEPDB, protein_id, pep_seq, y_axis)) %>%
  select(-protein_id, -pep_seq) %>%
  group_by(y_axis) %>%
  summarise_at(vars(all_of(colnames(gene_matrix))), sum) %>% 
  # mutate_at(vars(all_of(colnames(gene_matrix))), z_scores)
  pivot_longer(cols = colnames(gene_matrix), names_to = "LIBRARY_ID", values_to = "n") %>%
  right_join(.colData) %>%
  mutate(facet = "A) Overall DEGs") %>%
  # group_by(Diatery) %>% mutate(n = z_scores(n))
  group_by(y_axis) %>% mutate(n = z_scores(n))
  

unique_degs_df <- 
  UPSETDF %>% 
  ungroup() %>%
  filter(n == 1) %>% 
  unnest(summarise_col)

DF2 <- gene_matrix %>%
  as_tibble(rownames = "protein_id") %>%
  filter(protein_id %in% queries_for_zscore_heatmap) %>%
  left_join(distinct(CONOPEPDB, protein_id, pep_seq, y_axis)) %>%
  right_join(distinct(unique_degs_df, pep_seq)) %>%
  select(-protein_id, -pep_seq) %>%
  group_by(y_axis) %>%
  summarise_at(vars(all_of(colnames(gene_matrix))), sum) %>% 
  pivot_longer(cols = colnames(gene_matrix), names_to = "LIBRARY_ID", values_to = "n") %>%
  right_join(.colData) %>%
  filter(Time != "Ctrl") %>%
  mutate(facet = "B) Exclusive DEGs") %>%
  # group_by(Diatery) %>% mutate(n = z_scores(n))
  group_by(y_axis) %>% mutate(n = z_scores(n))

HeatmapViz <- DF1 %>%
  rbind(DF2) %>%
  drop_na(y_axis) %>%
  mutate(label = scales::comma(n)) %>%
  mutate(y_axis = factor(y_axis, levels = levels(Globaldf$y_axis))) 

lo = floor(min(HeatmapViz$n))
up = ceiling(max(HeatmapViz$n))
mid = (lo + up)/2


HeatmapViz %>%
  mutate(Time = factor(Time, levels = c("Ctrl", "2 month", "4 month"))) %>%
  drop_na(Diatery) %>%
  ggplot(aes(y = y_axis, x = Diatery, fill = n)) + 
  # facet_grid(~ sampleX + summarise_col) +
  ggh4x::facet_nested(~ facet + Time, nest_line = T, scales = "free_x", space = "free_x") +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  # geom_raster() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
    na.value = "white", midpoint = 0, limit = c(lo, up), 
    breaks = c(-2, 0, 3), # seq(lo, up, by = 3),
    name = NULL) +
  # geom_text(aes(label= label), hjust= 1, vjust = 0.5, size = 3, family = "GillSans", color = "white") +
  theme_bw(base_family = "GillSans", base_size = 10) +
  labs(x = "", y = "Superfamily") +
  labs(x = "", y = "Superfamily") +
  theme(
    legend.position = "top",
    # panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    # axis.text.y.right = element_text(angle = 0, hjust = 1, vjust = 0, size = 2.5),
    axis.text.y = element_text(angle = 0, size = 12),
    axis.text.x = element_text(angle = 0),
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 1, size = 12)) -> P

P  <- P + guides(
  fill = guide_colorbar(barwidth = unit(1.5, "in"),
    barheight = unit(0.1, "in"), label.position = "bottom",
    alignd = 0.5,
    title = "Row Z-score",
    title.position  = "top",
    title.theme = element_text(size = 10, family = "GillSans", hjust = 1),
    ticks.colour = "black", ticks.linewidth = 0.35,
    frame.colour = "black", frame.linewidth = 0.35,
    label.theme = element_text(size = 10, family = "GillSans")
  ))


ggsave(P, filename = 'Superfamilies_by_degs_reads.png', 
  path = pub_dir, width = 10, height = 10, dpi = 500, device = png)



# Test if allelic ----

# split by sf (y_axis)
# write in a list of vectors
# apply msa::msa
# format to tydy format
# unlist object of lists
# plot individual
# OR
# Calculate ConsensusSequence score per sequences group, and plot in a single plot all the sf


# myXStringSet <- distinct(CONOPEPDB, protein_id, pep_seq, y_axis) %>%
#   right_join(distinct(overall_degs_df, pep_seq)) %>%
#   drop_na() %>%
#   mutate(seqname = paste(y_axis, protein_id, sep = "|")) %>% 
#   distinct(y_axis, seqname, pep_seq) %>%
#   pull(pep_seq, name = y_axis)
  
myXStringSet <- Globaldf %>% filter(y_axis == "MTFLLLLVSV") %>% 
  left_join(CONOPEPDB) %>%  
  mutate(y_axis = paste0(y_axis, "|", protein_id)) %>% pull(pep_seq, name = y_axis)


myXStringSet <- Biostrings::AAStringSet(c(myXStringSet))

library(msa)

align <- msa::msa(myXStringSet, method = "ClustalW", order = "input")

.align <- msa::msaConvert(align)$seq
# 
names(.align) <- msa::msaConvert(align)$nam

library(ggsci)

pat <- c("-", alphabet(myXStringSet, baseOnly=TRUE))
# 
DECIPHER::BrowseSeqs(AAStringSet(.align), colWidth = Inf)


# Method 1
# conMat <- consensusMatrix(align)

data(BLOSUM62)
msa::msaConservationScore(align, BLOSUM62)

calculate_entropy <- function(alignment) {
  
  
  require(stringr)
  require(dplyr)
  require(tidyr)
  
  # Read the alignment file (assuming FASTA format)
  # alignment <- Biostrings::readAAStringSet(alignment_file)
  #alignment <- read.fasta(alignment_file, as.string = TRUE) #Alternative using ape package
  
  alignment <- AAStringSet(alignment)
  
  # Convert to a matrix where rows are sequences and columns are positions
  alignment_matrix <- str_split(as.character(alignment), "", simplify = TRUE)
  
  # Get the number of positions
  num_positions <- ncol(alignment_matrix)
  
  # Initialize a list to store entropy values for each position
  entropy_values <- numeric(num_positions)
  
  # Iterate through each position and calculate entropy
  for (i in 1:num_positions) {
    # Get the column (position) from the matrix
    position_data <- alignment_matrix[, i]

    # Calculate frequencies
    frequencies <- table(position_data) / length(position_data)
    
    
    # Recalculate to zero gaps from the entropy
    # position_data <- position_data[!grepl("-", position_data)]
    # position_data <- position_data[!grepl("-", names(frequencies))]
    
    # Calculate entropy using Shannon entropy formula
    entropy <- -sum(frequencies * log2(frequencies), na.rm = TRUE)
    
    # entropy <- -sum(frequencies * log2(frequencies + 1e-10)) # Add a small value to avoid log(0)
    
    max_H <- log2(length(frequencies))
    
    # If normalize (range 0 to 1)
    entropy_values[i] <- 1 - (entropy / max_H)
    
    # Store the entropy value
    entropy_values[i] <- entropy
    
    # of if want to use the max frequent residue
    # entropy_values[i] <- max(frequencies)
  }
  
  # Return the entropy values
  return(entropy_values)
}


plot(calculate_entropy(align))

ggseqlogo::ggseqlogo(.align)

library(ggseqlogo)
# 
# LOGO <- geom_logo(.align, p = F, method = "probability")
# 
# lo = floor(min(LOGO$y))
# up = ceiling(max(LOGO$y))
# mid = (lo + up)/2
# 
# LOGO %>% group_by(position) %>% summarise(sum(y))
# 
# LOGO %>% 
#   # mutate()
#   ggplot(aes(y = y, x = position)) + 
#   geom_po
