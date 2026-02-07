# Figures for pub
# LOAD count matrix, conopeptidesDB (astringent filter) and DEGS,
# Transform count matrix to vst
# Contrast deg-subset vs raw matrix and calculate PCA
# USE deg-subset to calculate zscore and plot heatmap
# 

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

extrafont::loadfonts(device = "win")

base_text_fam <- "Gill Sans MT"

my_custom_theme <- function(base_size = 14, legend_pos = "top", ...) {
  base_size = 14
  theme_bw(base_family = "Gill Sans MT", base_size = base_size) +
    theme(legend.position = legend_pos,
          strip.placement = "outside", 
          strip.background = element_rect(fill = 'gray90', color = 'white'),
          strip.text = element_text(angle = 0, size = base_size, hjust = 0), 
          axis.text = element_text(size = rel(0.7), color = "black"),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          ...
    )
}

pub_dir <- "C://Users//cinai/OneDrive/Documentos/PUBLICATION_DIR"

# dir <- "/Users/cigom/Documents/GitHub/conopeptides/06.Quantification/MATRIX_RSEM_dir/"

dir <- "C://Users//cinai/OneDrive/Documentos/PUBLICATION_DIR/06.Quantification/MATRIX_RSEM_dir/"



library(tidyverse)


recode_to <- structure(c("Control","Shrimp", "Mollusk", "Polychaete", "Mixed"))
 
# col_values <- c("#282419", "#C43726", "#F1DCBD", "#449A6D", "#366B4C")

diet_col <- c("gray20","#146179", "#09BC9F", "#FEB65F", "#C55E2D")

diet_col <- structure(diet_col, names = recode_to)

recode_to <- structure(recode_to, names = c("Ctrl","Cam","Lit","Pol","Mix"))

recode_time <- structure(c("","2 months", "4 months"))
recode_time <- structure(recode_time, names = c("Ctrl","2","4"))

# data 1 ----

subdir <- "KALLISTO_Merged_polyA_hisat_SuperDuper.fasta.transdecoder_DIR/"

f <- "KALLISTO_Merged_polyA_hisat_SuperDuper.fasta.transdecoder.matrix"

f <- list.files(file.path(dir, subdir), f, full.names = T)

dim(datExpr <- round(readRDS(f)))


# data 2 ----

Manifest <- list.files(dir, pattern = "Manifest", full.names = T)  %>%
  read_tsv() %>% 
  mutate(LIBRARY_ID = ifelse(grepl("cam2v_", LIBRARY_ID), NA, LIBRARY_ID)) %>%
  mutate(LIBRARY_ID = ifelse(grepl("cam6", LIBRARY_ID), NA, LIBRARY_ID)) %>%
  # mutate(LIBRARY_ID = ifelse(grepl("ctrl", LIBRARY_ID), NA, LIBRARY_ID)) %>%
  select(LIBRARY_ID, Time, Diatery) %>% drop_na(LIBRARY_ID) %>%
  dplyr::mutate(Diatery = dplyr::recode_factor(Diatery, !!!recode_to)) %>%
  dplyr::mutate(Time = dplyr::recode_factor(Time, !!!recode_time)) %>%
  mutate(Design = paste(Time, Diatery, sep = "_")) %>%
  mutate(LIBRARY_ID = gsub("_|-",".", LIBRARY_ID)) %>%
  mutate_if(is.character, as.factor)

# data 2.5


CONOPEPDB <- read_tsv(paste0(pub_dir, "/conopeptides.tsv")) %>%  #view()
  mutate(len = nchar(pep_seq)-1) %>%
  filter(effective_length_frac > 0.5) %>%
  # To be consistent w/ RES
  filter(Signalp_class == "SP") %>%
  drop_na(prediction_tool, Superfamily) 

# data 3 ----

f <- "glmLRT_multiple_contrast_ctrl_and_treatments_kallisto.rds"

DEGS <- read_rds(file.path(dir, f)) %>% filter(abs(logFC) > 4 & FDR < 0.05) 

# try individual Ctrl from Diet at timeX

alldf <- DEGS %>% dplyr::count(sampleA, sampleB, sampleX, sort = T)

DEGS <- DEGS %>% left_join(CONOPEPDB) %>% drop_na(Superfamily)

DEGS %>% dplyr::count(sampleA, sampleB, sampleX, sort = T) %>% right_join(alldf)

ControlDEGs <- DEGS %>%
  filter(sampleX == "Ctrl") %>%
  distinct(protein_id, Superfamily, sampleX)

# from DEGS, select only UPexpressed for every Diet at timeX

cols_to_check <- c("sampleA", "sampleB")

DEGS <- DEGS %>%
  filter(sampleX != "Ctrl") %>%
  filter(
    if_any(all_of(cols_to_check), ~ str_detect(.x, "Ctrl"))) 


DataViz <- rbind(ControlDEGs %>% dplyr::count(sampleX, Superfamily, sort = T),
                 DEGS %>% dplyr::count(sampleX, Superfamily, sort = T)) 

yaxis_level <- DataViz %>% group_by(Superfamily) %>% 
  tally(n,sort = T)  %>%  pull(Superfamily)


DataViz %>%
  mutate(Superfamily = factor(Superfamily, levels = rev(yaxis_level))) %>%
  ggplot(aes(sampleX, Superfamily, label = n)) + geom_text() + my_custom_theme()

DEGS %>% 
  mutate(Superfamily = factor(Superfamily, levels = rev(yaxis_level))) %>%
  ggplot(aes(y = Superfamily, x = abs(logFC))) + 
  facet_grid(~ sam_group) +
  # geom_jitter(position = position_jitter(0.1), shape = 1) +
  stat_summary(fun = "mean", geom = "line", color="red") +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red") +
  labs(x = expression(Log[2] ~FC), y = "Gene Superfamily") +
  my_custom_theme()

# data 4

CONOPEPDB <- DEGS %>% distinct(protein_id, Superfamily)%>% dplyr::rename("yaxis" = "Superfamily")

yaxis_level <- CONOPEPDB %>% dplyr::count(yaxis, sort = T) %>% pull(yaxis)


rbind(ControlDEGs %>% dplyr::count(sampleX, Superfamily, sort = T),
      DEGS %>% dplyr::count(sampleX, Superfamily, sort = T)) %>%
  mutate(Superfamily = factor(Superfamily, levels = rev(yaxis_level))) %>%
  ggplot(aes(sampleX, Superfamily, label = n)) + geom_text()

# CONOPEPDB <- read_tsv(paste0(pub_dir, "/conopeptides.tsv")) %>%
#   filter(Signalp_class == "SP" & contig_impact_score > 0 & nchar(pep_seq) < 200) %>% drop_na(Superfamily) %>%
#   distinct(protein_id, Superfamily) %>% dplyr::rename("yaxis" = "Superfamily")
# 
# yaxis_level <- CONOPEPDB %>% dplyr::count(yaxis, sort = T) %>% pull(yaxis)

# CONOPEPDB %>% dplyr::count(yaxis) %>% view()


# remove samples from the plot
colnames(datExpr) <- gsub("_|-",".", colnames(datExpr))

keep <-  colnames(datExpr) %in% levels(Manifest$LIBRARY_ID) 

sum(keepRows <- rownames(datExpr) %in% CONOPEPDB$protein_id)

sum(keepCols <- colnames(datExpr) %in% levels(Manifest$LIBRARY_ID) )


dim(datExpr <- datExpr[keepRows,keepCols])



# NOT test using only degs for PCA

# PCA -----



#  blind=TRUE should be used for comparing samples in a manner unbiased by prior information on samples, 
# for example to perform sample QA (quality assurance). 

dfPCA <- DESeq2::varianceStabilizingTransformation(datExpr, blind = T)

# datExpr <- assay(dds)


PCA = prcomp(t(dfPCA), center = T, scale. = F)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
# percentVar <- round(PCA$sdev/sum(PCA$sdev)*100,1)

sd_ratio <- sqrt(percentVar[2] / percentVar[1])

PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

k <- 4

PCAdf %>% 
  dist(method = "euclidean") %>% 
  hclust() %>% 
  cutree(., k) %>% 
  as_tibble(rownames = 'LIBRARY_ID') %>% 
  mutate(cluster = paste0('C', value)) %>% 
  dplyr::select(-value) -> hclust_res


PCAdf %>%
  mutate(LIBRARY_ID = rownames(.)) %>%
  left_join(Manifest) %>% 
  left_join(hclust_res) %>%
  mutate(col = Diatery, label = ifelse(Time != "Ctrl", as.character(Time), "Control")) %>%
  ggplot(., aes(PC1, PC2, label = label)) +
  ggforce::geom_mark_ellipse(aes(label = label, group = label)) +
  # coord_fixed(ratio = sd_ratio) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  scale_color_manual(values = diet_col) +
  geom_point(size = 5, alpha = 1, aes(color = col), 
    shape = 21, stroke = 1.5, fill = "white") +
  ylim(-50, 50) +
  xlim(-50, 50) +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  my_custom_theme() +
  guides(color = guide_legend(title = "", byrow = T)) +
  theme(plot.title = element_text(hjust = -0.5), 
    legend.position = 'top',
    legend.spacing.x = unit(0.1, 'cm'),
    legend.text = element_text(size = 8),
    panel.grid = element_blank()
    # legend.spacing.y = unit(-1, 'mm')
  ) -> plot_pca

plot_pca

# Heatmap of degs -----
# Summarise first to Superfamilies category

queries_for_zscore_heatmap <- DEGS %>% 
  distinct(protein_id, pep_seq) %>% 
  pull(pep_seq, name = protein_id)

z_scores <- function(x) {(x-mean(x))/sd(x)}

sum(keepRows <- rownames(datExpr) %in% names(queries_for_zscore_heatmap))

sum(keepCols <- colnames(datExpr) %in% levels(Manifest$LIBRARY_ID))

# Note: here we must to secure Ctrl levels must be down. When transform to vst the pattern in Ctrl bias 
# Therore omiting this step datExpr1 <- transform_vst(datExpr[keepRows,keepCols])


HeatmapViz <- datExpr[keepRows,keepCols] %>%
  as_tibble(rownames = "protein_id") %>%
  right_join(CONOPEPDB) %>%
  group_by(yaxis) %>%
  summarise_at(vars(all_of(colnames(datExpr))), sum) %>%
  data.frame(row.names = .$yaxis) %>%
  select(-yaxis) %>% as.matrix()

str(HeatmapViz <- t(apply(HeatmapViz, 1, z_scores)))


HeatmapViz <- HeatmapViz %>% 
  # mutate_at(vars(all_of(colnames(datExpr))), z_scores) %>% # <- here is an issue as all control are down-epressed not red color m
  # data.frame(row.names = .$yaxis) %>%
  # select_at(vars(all_of(colnames(datExpr)))) %>%
  as_tibble(rownames = "yaxis") %>%
  pivot_longer(cols = colnames(HeatmapViz), names_to = "LIBRARY_ID", values_to = "n") %>%
  left_join(Manifest) %>% 
  mutate(yaxis = factor(yaxis, levels = rev(yaxis_level))) %>%
  # mutate(facet = "Global conotoxins expression")
  mutate(facet = "Differential up-expressed conotoxins under diet")
  


lo = floor(min(HeatmapViz$n))
up = ceiling(max(HeatmapViz$n))
mid = (lo + up)/2

plot_heatmap <- HeatmapViz %>%
  ggplot(aes(y = yaxis, x = Diatery, fill = n)) + 
  facet_grid(facet ~ Time, scales = "free",space = "free") +
  # geom_raster() +
  geom_tile(color = 'white', linewidth = 0.2) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
    na.value = "white", midpoint = 0, limit = c(lo, up), 
    breaks = c(lo, 0, up),
    name = NULL) +
  labs(y = "", x = "") +
  my_custom_theme() +
  theme(
    # legend.position = "top",
    # panel.grid.minor.y = element_blank(),
    # panel.grid.major.y = element_blank(),
    # panel.grid.minor.x = element_blank(),
    # panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1,size = 7),
    strip.text = element_text(
      angle = 0, hjust = 0,
      size = 10),
    strip.background = element_rect(colour = "transparent", fill = "transparent", size = 1)
  )

plot_heatmap <- plot_heatmap + guides(
  fill = guide_colorbar(
    barwidth = unit(2, "in"),
    barheight = unit(0.1, "in"), 
    label.position = "bottom",
    label.hjust = 0.5,
    title = "Row Z-score",
    title.position  = "top", title.hjust = 0,
    title.theme = element_text(size = 10, family = base_text_fam, hjust = 1),
    ticks.colour = "black", ticks.linewidth = 0.35,
    frame.colour = "black", frame.linewidth = 0.35,
    label.theme = element_text(size = 10, family = base_text_fam)
  ))


plot_heatmap

ggsave(plot_heatmap, filename = 'PUB_heatmap3.png', 
       path = pub_dir, width = 4, height = 5, device = png, dpi = 600)

# PLOT Number of intersected and unique

UPSETDF <- DEGS %>% 
  separate(sampleX, into = c("Time", "Diatery"), sep = "_")
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



ggsave(plot_pca, filename = 'PUB_PCA.png', 
  path = pub_dir, width = 4, height = 4, device = png, dpi = 600)

