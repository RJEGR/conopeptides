# Conotoxin composition diversity as Phuong et al., 2016:
# Open conotoxins.tsv
# Filter peptides with signalp and predicted as conotoxin by at least one tool
# filter(Signalp_class == "SP") %>% drop_na(prediction_tool) 
# Additionally, filter conopeptide genes with low expression levels (cpm < 1)
# N mature conotoxin
# N gene superfamilies
# N cys frameworks
# Split by treatment

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

# DB <- read_rds(paste0(pub_dir, "/structured_db.rds"))

DB %>% mutate(len = nchar(pep_seq)) %>% 
  # count(prediction_tool, len) %>%
  # drop_na(prediction_tool) %>%
  ggplot() + stat_ecdf(aes(len, color = prediction_tool))
  # geom_col(position = position_stack(), aes(y = n/90497, x = len, fill = prediction_tool)) 

CONOPEPDB <- read_tsv(paste0(pub_dir, "/conopeptides.tsv"))

CONOPEPDB %>% count(Signalp_class, prediction_tool)

CONOPEPDB <- CONOPEPDB %>% filter(Signalp_class == "SP") %>% drop_na(prediction_tool) 


dir <- "/Users/cigom/Documents/GitHub/conopeptides/06.Quantification/MATRIX_RSEM_dir/"

# EDGER RES: 1678 putative conopeptides presented in the contrast selected (not DEGs filtered yet)
RES.P <- read_rds(paste0(dir, "/glmLRT_multiple_contrast_ctrl_and_treatments.rds")) %>%  
  filter(FDR < 0.05 & abs(logFC) > 2)

RES.P %>% distinct(gene_id) # 1,224 putative pep are degs

.count_vst <- read_rds(paste0(dir, "/counts_vst_nt_raw.rds"))$vst

.colData <- list.files(dir, pattern = "Manifest", full.names = T) 

.colData <- read_tsv(.colData) %>% 
  mutate(Diatery = ifelse(grepl("cam2v_", LIBRARY_ID), "camv", Diatery)) %>%
  mutate(design = ifelse(!Time == "Ctrl", paste0(Diatery, "_", Time), "Ctrl")) %>%
  select(LIBRARY_ID, Time, Diatery, design) %>%
  mutate_if(is.character, as.factor)

RES.P %>% dplyr::count(sam_group, sampleA, sampleB, sampleX)

# 
# As ORFs include exclusively peptides with start and stop codon. Considere include all classes of regions

CONOPEPDB %>% mutate(len = nchar(pep_seq)) %>% ggplot() + geom_histogram(aes(len)) + facet_grid(prediction_tool ~ .)

# recode_peptide <- structure(c("Mature","Propeptide", "Precursor",), names = Which_regions)

# CONOPEPDB <- CONOPEPDB %>% mutate(Region = recode_factor(Region, !!!recode_peptide, .ordered = T)) 

CONOPEPDB %>% distinct(pep_seq)
CONOPEPDB %>% distinct(gene_id)
CONOPEPDB %>% count(Superfamily, sort = T)

# 0) choose if DEGs or total N of putative conopeptide genes ========

# here is where the subset of genes must be filtered, for example by exclusively DEGS or DEGs w/ assignment by BOTH predited tool

str(query_genes <- CONOPEPDB %>% distinct(gene_id) %>% pull()) # total 3514
  
str(query_degs <- RES.P %>% distinct(gene_id) %>% pull()) # 1224 degs


# query_degs <- CONOPEPDB %>%
#   filter(prediction_tool == "BOTH") %>%
#   distinct(gene_id) %>%
#   left_join(distinct(RES.P, gene_id)) %>%
#   distinct(gene_id) %>% pull()


str(query_degs <- RES %>% distinct(gene_id) %>% pull()) # 1224 degs

count_vst <- .count_vst

# Using degs as they not include low expression levels

dim(count_vst <- count_vst[rownames(count_vst) %in% query_degs,])


head(sort(rowSums(count_vst)))

z_scores <- function(x) {(x-mean(x))/sd(x)}

count_vst <- apply(count_vst, 1, z_scores)

h <- heatmap(t(count_vst), , keep.dendro = T)

hc_samples <- as.hclust(h$Colv)
plot(hc_samples)
hc_order <- hc_samples$labels[h$colInd]

hc_genes <- as.hclust(h$Rowv)
# plot(hc_genes)
order_genes <- hc_genes$labels[h$rowInd]

# Plot 1-----

# The short form for "Unidentified" is ID.

gene_names <- CONOPEPDB %>% 
  mutate(uniprotkb_toxprot = sapply(strsplit(uniprotkb_toxprot, "[|]"), `[`, 2)) %>%
  mutate(uniprotkb_toxprot = ifelse(is.na(uniprotkb_toxprot), "uID", uniprotkb_toxprot)) %>%
  mutate(conoserver_protein = ifelse(is.na(conoserver_protein), "uID", conoserver_protein)) %>%
  mutate(conoserver_protein = ifelse(is.na(conoserver_protein), "uID", conoserver_protein)) %>%
  mutate(hmm_pred_conodictor = stringr::str_to_sentence(hmm_pred_conodictor)) %>%
  mutate(hmm_pred_conodictor = ifelse(is.na(hmm_pred_conodictor), "uID", hmm_pred_conodictor)) %>%
  mutate(Superfamily = ifelse(is.na(Superfamily), "uID", Superfamily)) %>%
  select(gene_id, prediction_tool, hmm_pred_conodictor, Superfamily, uniprotkb_toxprot, conoserver_protein) %>%
  unite("gene_name", prediction_tool:conoserver_protein, sep = "|") %>%
  distinct(gene_id, gene_name) %>% 
  pull(gene_name, name = gene_id)

# gene_names <- CONOPEPDB %>%  distinct(Superfamily) %>% pull(Superfamily, name = Superfamily)

sum(order_genes %in% names(gene_names))

gene_names <- gene_names[match(order_genes, names(gene_names))]

identical(names(gene_names), order_genes) # TRUE

# head(gene_names <- gsub("(Mature)", "M", gene_names))
# head(gene_names <- gsub("(Pro-region)", "P", gene_names))
# head(gene_names <- gsub("(Signal)", "S", gene_names))


HeatmapViz <- count_vst %>% t() %>%
  as_tibble(rownames = 'gene_id') %>%
  pivot_longer(cols = rownames(count_vst), values_to = 'fill', names_to = "LIBRARY_ID") %>%
  left_join(.colData, by = "LIBRARY_ID") 


lo = floor(min(HeatmapViz$fill))
up = ceiling(max(HeatmapViz$fill))
mid = (lo + up)/2

HeatmapViz %>% distinct(gene_id)

HeatmapViz %>%
  filter(abs(fill) > 1) %>%
  filter(Diatery != "camv") %>% filter(Time != "6") %>%
  mutate(Time = recode_factor(Time, !!!recode_time, .ordered = T)) %>%
  mutate(Diatery = recode_factor(Diatery, !!!recode_Diatery, .ordered = T)) %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = hc_order)) %>%
  mutate(gene_id = factor(gene_id, levels = order_genes)) %>%
  ggplot(aes(x = Time, y = gene_id, fill = fill)) + 
  facet_grid(~ Diatery, scales = "free_x", space = "free_x") +
  # geom_tile(color = 'white', linewidth = 0.2) +
  geom_raster() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
    na.value = "white", midpoint = mid, limit = c(lo, up), 
    breaks = seq(lo, up, by = 3),
    name = NULL) +
  labs(x = "", y = "") +
  scale_y_discrete(position = 'left', labels = gene_names) +
  # scale_x_discrete(position = "bottom") 
  # ggh4x::scale_y_dendrogram(hclust = hc_genes, position = "left", labels = NULL) +
  # guides( y.sec = ggh4x::guide_axis_manual(labels = gene_names, label_size = 5)) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top",
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 0, size = 8),
    axis.line.y = element_line(color = 'white'),
    axis.text.y = element_text(hjust = 1, size = 7),
    # axis.ticks.y = element_blank(),
    axis.ticks.length = unit(2.5, "pt")) 

# Exit or continue

# 1) Create the column
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


