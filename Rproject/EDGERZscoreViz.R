# Verify direction of logFC sign
# Heatmap of zscore for conopeptides 
# PCA only for conopeptides
# Diversity 

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

# RES <- read_rds(paste0(dir, "/glmLRT_multiple_contrast.rds")) %>% 
#   filter(FDR < 0.05 & abs(logFC) > 2) %>%
#   dplyr::rename("gene_id" = "ids") 

# glmLRT_multiple_contrast_ctrl_and_treatments.rds came from EDGERViz.R from subset of groups

RES <- read_rds(paste0(dir, "/glmLRT_multiple_contrast_ctrl_and_treatments.rds")) %>%  filter(FDR < 0.05 & abs(logFC) > 2)

count_vst <- read_rds(paste0(dir, "/counts_vst_nt_raw.rds"))$vst

.colData <- list.files(dir, pattern = "Manifest", full.names = T) 

.colData <- read_tsv(.colData) %>% 
  mutate(Diatery = ifelse(grepl("cam2v_", LIBRARY_ID), "camv", Diatery)) %>%
  mutate(design = ifelse(!Time == "Ctrl", paste0(Diatery, "_", Time), "Ctrl")) %>%
  select(LIBRARY_ID, Time, Diatery, design) %>%
  mutate_if(is.character, as.factor)



# Match only conopeptides
Which_regions <- c("(Mature)","(Mature)-(Pro-region)", "(Mature)-(Pro-region)-(Signal)")

CONOPEPDB <- DB %>% drop_na(tab) %>% 
  # dplyr::count(Signalp_class)
  filter(Signalp_class == "SP") %>%
  filter(Region %in% Which_regions)

# recode peptides
recode_peptide <- structure(c("Mature","Propeptide", "Precursor"), names = c("(Mature)", "(Mature)-(Pro-region)", "(Mature)-(Pro-region)-(Signal)"))

CONOPEPDB <- CONOPEPDB %>% 
  mutate(Region = recode_factor(Region, !!!recode_peptide, .ordered = T)) 

CONOPEPDB %>% distinct(Region, Superfamily)


RES %>% distinct(gene_id)

str(query_genes <- RES %>% distinct(gene_id) %>% left_join(CONOPEPDB) %>% distinct(gene_id) %>% pull())

CONOPEPDB %>% 
  select(gene_id, Region, Superfamily,tab, conoserver_protein, uniprotkb_toxprot) %>%
  right_join(RES) %>%
  write_tsv(file = paste0(pub_dir, "glmLRT_multiple_contrast_ctrl_and_treatments_annot_peptides.tsv"))
  
dim(count_vst <- count_vst[rownames(count_vst) %in% query_genes,])

# agglomerate gene_matrix by same superfamily (not!! because gene profile going to be masked, just use for previz purpose )

CONOPEPDB %>% count(Superfamily, sort = T)

# nrow(.count <- count_vst %>% 
#   as_tibble(rownames = "gene_id") %>%
#   left_join(distinct(CONOPEPDB, Superfamily, gene_id)) %>% 
#   group_by(Superfamily) %>%
#   summarise_at(vars(all_of(colnames(count_vst))), sum) %>% ungroup() )

# count_vst <- .count %>% elect_if(is.double) %>% as("matrix")

# rownames(count_vst) <- .count$Superfamily

z_scores <- function(x) {(x-mean(x))/sd(x)}

str(HeatmapViz <- t(apply(count_vst, 1, z_scores)))

HeatmapViz[is.na(HeatmapViz)] <- 0

# prep sanity check for elimin duplications in CONOPEPDB

HeatmapViz <- HeatmapViz %>% 
  as_tibble(rownames = 'gene_id') %>%
  pivot_longer(cols = colnames(HeatmapViz), values_to = 'fill', names_to = "LIBRARY_ID") %>%
  left_join(.colData, by = "LIBRARY_ID") %>% 
  left_join(distinct(CONOPEPDB, gene_id, Superfamily), by = "gene_id") 


h <- heatmap(count_vst, col = cm.colors(12), keep.dendro = T)

hc_samples <- as.hclust(h$Colv)
plot(hc_samples)
hc_order <- hc_samples$labels[h$colInd]

hc_genes <- as.hclust(h$Rowv)
# plot(hc_genes)
order_genes <- hc_genes$labels[h$rowInd]

# re-label gene_names
gene_names <- CONOPEPDB %>% 
  # mutate(gene_id = Superfamily) %>%
  # mutate(uniprot = gsub("_HUMAN", "", uniprot)) %>%
  # mutate(uniprot = paste0(uniprot, " (", protein_name,")")) %>% 
  distinct(gene_id, Superfamily) %>% 
  pull(Superfamily, name = gene_id)

# gene_names <- CONOPEPDB %>%  distinct(Superfamily) %>% pull(Superfamily, name = Superfamily)

sum(order_genes %in% names(gene_names))

gene_names <- gene_names[match(order_genes, names(gene_names))]

identical(names(gene_names), order_genes) # TRUE

# head(gene_names <- gsub("(Mature)", "M", gene_names))
# head(gene_names <- gsub("(Pro-region)", "P", gene_names))
# head(gene_names <- gsub("(Signal)", "S", gene_names))

lo = floor(min(HeatmapViz$fill))
up = ceiling(max(HeatmapViz$fill))
mid = (lo + up)/2

# hc_order2 <- DataViz %>% drop_na(CONTRASTE_D) %>% group_by(CONTRASTE_D, LIBRARY_ID) %>% summarise(fill = sum(fill)) %>% pull(LIBRARY_ID)
# hc_order2 <- sort(gene_names)


# DataViz %>% view()

HeatmapViz %>%
  filter(Diatery != "camv") %>%
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
  # scale_y_discrete(position = 'left', labels = gene_names) +
  scale_x_discrete(position = "bottom") +
  ggh4x::scale_y_dendrogram(hclust = hc_genes, position = "left", labels = NULL) +
  guides( y.sec = ggh4x::guide_axis_manual(labels = gene_names, label_size = 5)) +
  theme_classic(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top",
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 0, size = 8),
    axis.line.y = element_line(color = 'white'),
    axis.text.y = element_text(hjust = 1, size = 7),
    # axis.ticks.y = element_blank(),
    axis.ticks.length = unit(2.5, "pt")) -> p

# p


p <- p +  theme(
  # axis.ticks.x = element_blank(), 
  axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
  # axis.text.x = element_blank(), 
  axis.line.x = element_blank()
)

p <- p + theme(panel.spacing.x = unit(0, "mm"))


p <- p + guides(
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


ggsave(p, filename = 'Row_Zscore_heatmap_short.png', 
  path = pub_dir, width = 5, height = 5, dpi = 500, device = png)

# Diversity plot -----


SummarizedExperiment::assay(dds) %>% 
  as_tibble(rownames = 'transcript_id') %>%
  pivot_longer(cols = colnames(datExpr), values_to = 'fill', names_to = "LIBRARY_ID") %>%
  left_join(datTraits, by = "LIBRARY_ID") %>% 
  left_join(distinct(DB,transcript_id, uniprot), by = "transcript_id") %>%
  separate(uniprot,into = c("uniprot", "genus"), sep = "_") %>%
  mutate(transcript_id = factor(transcript_id, levels = rev(order_genes))) %>%
  drop_na(CONTRASTE_D) %>%
  ggplot(aes(x = transcript_id, y = fill, color = CONTRASTE_D)) +
  geom_boxplot(outlier.shape=NA)+
  # geom_point(position=position_jitterdodge(), alpha = 0.5) +
  scale_x_discrete(labels = gene_names) +
  labs(y = "vst", x = "") +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(
    # axis.ticks.x = element_blank(), 
    axis.text.x = element_text(angle = -90, hjust = 0, size = 10),
    # axis.text.x = element_blank(), 
    # axis.line.x = element_blank()
  ) + coord_flip() -> p


ggsave(p, filename = 'Oncogenes_boxplot.png', 
  path = dir, width = 10, height = 7, dpi = 500, device = png)

