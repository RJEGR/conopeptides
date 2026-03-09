
# Read .filt.rds files
# Calculate PCA
# Plot 
# Calculate samples-dendogram
# Plot

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


pub_dir <- "C://Users//cinai/OneDrive/Documentos/PUBLICATION_DIR/"
  
# dir <- "/Users/cigom/Documents/GitHub/conopeptides/06.Quantification/MATRIX_RSEM_dir"

# f <- list.files(dir, pattern = "Merged_polyA_hisat_SuperDuper_isoforms.filt.rds", full.names = T)

# f <- list.files(dir, pattern = "Merged_polyA_hisat_SuperDuper.fasta.transdecoder_genes.filt.rds", full.names = T)
# f <- list.files(dir, pattern = "Merged_polyA_hisat_SuperDuper.fasta.transdecoder_genes.rds", full.names = T)

dir <- "C://Users//cinai/OneDrive/Documentos/PUBLICATION_DIR/06.Quantification/MATRIX_RSEM_dir/KALLISTO_Merged_polyA_hisat_SuperDuper.fasta.transdecoder_DIR/"

f <- list.files(dir, pattern = "KALLISTO_Merged_polyA_hisat_SuperDuper.fasta.transdecoder.matrix", full.names = T)

library(DESeq2)
library(tidyverse)

dim(COUNT <- read_rds(f))

COUNT <- round(COUNT)

#
#
#

time_levs <- c("Ctrl", "2", "4", "6")

recode_time <- structure(c("Control", "2 months", "4 months", "6 months"), names = time_levs)

Diatery_levs <- c("Ctrl","Cam", "Lit", "Pol", "Mix")

recode_Diatery <- structure(c("Control","Shrimp", "Mollusk", "Polychaete", "Mixed"), names = Diatery_levs)



dir <- "C://Users//cinai/OneDrive/Documentos/PUBLICATION_DIR/03.Coverage/"

Manifest <- list.files(dir, pattern = "Manifest.tsv", full.names = T)

Manifest <- read_tsv(Manifest) %>% distinct() %>%
  mutate(LIBRARY_ID = ifelse(grepl("cam2v_", LIBRARY_ID), NA, LIBRARY_ID)) %>%
  mutate(LIBRARY_ID = ifelse(grepl("cam6", LIBRARY_ID), NA, LIBRARY_ID)) %>%
  drop_na() %>%
  dplyr::mutate(Feeding = dplyr::recode_factor(Feeding, !!!recode_Diatery)) %>%
  dplyr::mutate(Time = dplyr::recode_factor(Time, !!!recode_time))

# colData <- data.frame(LIBRARY_ID = factor(Manifest$LIBRARY_ID)) %>% 
#   filter(Diatery != "camv") %>% filter(Time != "6") 

# Manifest <- as.data.frame(Manifest)
# rownames(Manifest) <- Manifest$LIBRARY_ID

# dds <- DESeqDataSetFromMatrix(COUNT,
  # Manifest,
  # design = ~ LIBRARY_ID)


COUNT <- DESeq2::varianceStabilizingTransformation(COUNT)

# datExpr <- assay(dds)


# 1) PCA ------


PCA = prcomp(t(COUNT), center = T, scale. = F)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
# percentVar <- round(PCA$sdev/sum(PCA$sdev)*100,1)

sd_ratio <- sqrt(percentVar[2] / percentVar[1])

PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

PCAvar <- data.frame(
  Eigenvalues = PCA$sdev,
  percentVar = percentVar,
  Varcum = cumsum(percentVar))

PCAvar %>%
  mutate(Dim = row_number()) %>%
  ggplot(., aes(y = percentVar, x = as.factor(Dim))) +
  geom_col(position = position_dodge2()) +
  geom_path(aes(y = Varcum), group = 1) +
  geom_point(aes(y = Varcum), size = 0.5) +
  # labs(x = "Component Number", y = "Eigenvalue") +
  labs(x = "Principal component", y = "Fraction variance explained (%)") +
  scale_fill_grey("") +
  scale_color_grey("") +
  guides(color=guide_legend(nrow = 1)) +
  theme_classic(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()) 

# ggsave(p, filename = 'Eigenevalues.png', path = dir, width = 4, height = 3, device = png, dpi = 300)

k <- 4

PCAdf %>% 
  dist(method = "euclidean") %>% 
  hclust() %>% 
  cutree(., k) %>% 
  as_tibble(rownames = 'LIBRARY_ID') %>% 
  mutate(cluster = paste0('C', value)) %>% 
  dplyr::select(-value) -> hclust_res


# scales::show_col(see::pizza_colors())

# col_values <- c("gray20","#ad6aea","#ffd700","#63b8ff", "#ee2c2c")

col_values <- c("#282419", "#C43726",  "#449A6D", "#366B4C") #"#F1DCBD",


recode_time <- c(  `Ctrl` = "Control",
  `2`= "2 months", `4` = "4 months",
  `6` = "6 months")



PCAdf %>%
  mutate(LIBRARY_ID = rownames(.)) %>%
  right_join(Manifest) %>% 
  # left_join(hclust_res) %>%
  # dplyr::mutate(Time = dplyr::recode_factor(Time, !!!recode_time)) %>%
  mutate(col = Time, label = Feeding) %>%
  ggplot(., aes(PC1, PC2, label = label)) +
  # coord_fixed(ratio = sd_ratio) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  # scale_color_grey() +
  scale_color_manual(values = col_values) +
  geom_point(size = 5, alpha = 1, aes(color = col), 
    shape = 21, stroke = 1.5, fill = "white") +
  # geom_text( family = "GillSans", mapping = aes(label = label), size = 5) +
  ggrepel::geom_text_repel(family = base_text_fam, mapping = aes(label = label), size = 5) +
  # ylim(-200, 200) + xlim(-400, 400) +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  my_custom_theme() +
  guides(color = guide_legend(title = "", byrow = T)) -> p

p

ggsave(p, filename = 'SAMPLE_PCA_FOR_PUB.png', path = pub_dir, width = 5, height = 5, device = png, dpi = 800)


# 2) as heatmap ----

sample_cor = cor(datExpr, method='pearson', use='pairwise.complete.obs')

sample_dist = dist(sample_cor, method='euclidean')

hc_samples = hclust(sample_dist, method='complete')

hc_order <- hc_samples$labels[hc_samples$order]

# sample_cor <- sample_cor[match(hc_order, rownames(sample_cor)), ]

# sample_cor <-sample_cor[match(hc_order, colnames(sample_cor)), ]

sample_cor %>%
  as_tibble(rownames = 'LIBRARY_ID') %>%
  pivot_longer(cols = colnames(sample_cor), values_to = 'cor') %>%
  left_join(Manifest) %>%
  # mutate(Sample_group = factor(Sample_group, levels = unique(Sample_group))) %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = hc_order)) %>%
  mutate(name = factor(name, levels = hc_order)) -> sample_cor

scale_x <- structure(Manifest$Sample_group, names = Manifest$LIBRARY_ID)

library(ggh4x)

# recode_to <- c(  `Control` = "(A) Control",
#   `METASTASIS` = "Metastasis", `NO METASTASIS`= "No metastasis",
#   `Grado I` = "(B) Stage I", `Grado II` = "(C) Stage II", `Grado III` = "(D) Stage III",
#   `Indiferenciado` = NA)

sample_cor %>%
  ggplot(aes(x = LIBRARY_ID, y = name, color = cor, fill = cor)) + 
  geom_tile() + 
  ggsci::scale_color_material(name = "", "blue-grey") +
  ggsci::scale_fill_material(name = "", "blue-grey") +
  ggh4x::scale_y_dendrogram(hclust = hc_samples, position = "left", labels = NULL) +
  ggh4x::scale_x_dendrogram(hclust = hc_samples, position = "top", labels = NULL) +
  # scale_x_discrete(labels = scale_x) +
  guides(
  # y.sec = guide_axis_manual(labels = hc_order, label_size = 3.5)
  x.sec = guide_axis_manual(labels = hc_order, label_size = 3.5)
    ) +
  labs(x = '', y = '') +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(legend.position = 'bottom',
    axis.text.x = element_text(angle = 90,
      hjust = 1, vjust = 0.5, size = 10)
  ) -> p

p <- p +  guides(
  colour = guide_colorbar(barwidth = unit(2, "in"),
    barheight = unit(0.05, "in"), label.position = "bottom",
    alignd = 0.5,
    ticks.colour = "black", ticks.linewidth = 0.5,
    frame.colour = "black", frame.linewidth = 0.5,
    label.theme = element_text(family = "GillSans", size = 7))) 

p <- p + theme(
  # strip.background = element_rect(fill = 'white', color = 'white'),
  panel.border = element_blank(),
  plot.background = element_rect(fill='transparent', color = 'transparent'),
  plot.margin = unit(c(0,0,0,0), "pt"),
  panel.grid.minor = element_blank(),
  # axis.ticks.x = element_blank(),
  panel.grid.major = element_blank())

p

