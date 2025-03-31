# Supplementary material for pub
# LOAD data.frame from Assembly Nx metrics for DNA and orf (transdecoder.Orfs) [Nx_metrix.R] for every assembly method
# LOAD BUSCO completeness for every assembly method (BUSCO.R)
# LOAD superfm_df data.frame from ConoSorter_viz.R 

# PLOT lineplot of Nx facet DNA from ORF
# PLOT BUSCO scores by complete (single + duplicate) ~ assembly method (trinity, spades, concat), coloring by assembly step (raw vs full-length)
# PLOT summary of known/novel Superfamily classes detected by methods. Including (Mature)-(Pro-region)-(Signal) class 

# as mmseq method is only a deduplicate method, omit from the figures

rm(list = ls())

if(!is.null(dev.list())) dev.off()

library(tidyverse, help, pos = 2, lib.loc = NULL)

pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR"

Nxmetrix_df <- read_rds(paste0(pub_dir, "/Contig_nx.rds")) %>% filter(!grepl("Mmseq", Assembly))

# superfm_df %>% distinct(Method)
Nxmetrix_df %>% distinct(Assembly)

# Scale colors

library("scales")

# show_col(ggsci::pal_d3("category20")(20))
# show_col(ggsci::pal_gsea("default")(5))

# scale_color <- c("#D60C00FF", "#FF8888FF", "#8C564BFF", "#C49C94FF", "black", "gray40" )

recode_to <- unique(Nxmetrix_df$Assembly) %>% levels()

scale_color <- c("#D60C00FF", "#FF8888FF", "black", "gray40", "#2600D1FF","#8787FFFF")

scale_color <- structure(scale_color, names = recode_to)


# Plot Nx =====

# segment.curvature = 1 increases right-hand curvature, negative values would increase left-hand curvature, 0 makes straight lines
# segment.ncp = 3 gives 3 control points for the curve
# segment.angle = 20 skews the curve towards the start, values greater than 90 would skew toward the end

recode_to <- structure(c("A) DNA contigs", "B) Complete ORFs"), names = c("DNA contigs", "Complete ORFs"))

.Nxmetrix_df <- Nxmetrix_df %>% 
  mutate(stringSet = dplyr::recode_factor(stringSet, !!!recode_to))
  
p1 <- 
  .Nxmetrix_df %>% 
  mutate(stringSet = dplyr::recode_factor(stringSet, !!!recode_to)) %>%
  ggplot(aes(x = x, y = n, group = Assembly, color = Assembly)) +
  # geom_vline(xintercept = "N50", linetype="dashed", alpha=0.5) +
  ggplot2::geom_path(linewidth = 1, lineend = "round") +
  facet_wrap(~  stringSet, scales = "free_y") +
  # geom_point(shape = 21, size = 4, aes(size = n_frac)) +
  labs(x = "Nx statistics", y = "Sequence length", color = "Assembly method") +
  # ggsci::scale_color_igv() +
  scale_color_manual(values = scale_color) +
  scale_y_continuous(labels = scales::comma_format()) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = 'white', color = 'white', size = 12),
    strip.text = element_text(color = "black",hjust = 0)) + 
  guides(color = guide_legend(title = "", nrow = 2, ncol = 4))


p1 <- p1 +
  geom_point(data = filter(.Nxmetrix_df, x == "N50"), shape = 21, size = 2) +
  ggrepel::geom_text_repel(data = filter(.Nxmetrix_df, x == "N50"),
    aes(label= paste0(Assembly, " (", scales::comma(n), ")")),
    nudge_x      = 1,
    # nudge_y = 1,
    direction    = "y",
    hjust        = 0,
    vjust = -15,
    segment.curvature = 0.1,
    # segment.ncp = 3,
    # segment.angle = 20,
    size = 3, family = "GillSans") 

p1

ggsave(p1, filename = 'Nx-methods.png', path = pub_dir, width = 9, height = 3.5, device = png, dpi = 300)

# ConoSorter =====

# count the number of sequences per method and level (nuc and pep)
# 
# recode_to <- unique(superfm_df$Method) %>% levels()

assembly_lev <- Nxmetrix_df %>%
  filter(stringSet %in% "DNA contigs") %>%
  group_by(Assembly, stringSet) %>%
  summarise(n = sum(n_seqs)) %>% arrange(desc(n)) %>% 
  pull(Assembly) %>% as.character()

# recode_to <- structure(c("DNA contigs", "Complete ORFs"), )

pleft <- Nxmetrix_df %>%
  group_by(Assembly, stringSet) %>%
  summarise(n = sum(n_seqs)) %>%
  mutate(Assembly = factor(Assembly, levels = rev(assembly_lev))) %>%
  mutate(topfacet = "A) Transcriptome") %>%
  mutate(stringSet = factor(stringSet, levels = c("DNA contigs", "Complete ORFs"))) %>%
  ggplot(aes(y = Assembly, x = n/100000, fill = Assembly)) +
  geom_col(fill = "black") + 
  labs(y = "Method") +
  scale_fill_manual(values = scale_color) +
  # facet_wrap(~  stringSet, scales = "free_x") +
  ggh4x::facet_nested(~ topfacet + stringSet, scales = "free_x") +
  geom_text(aes(label= scales::comma(n)), hjust= 1.05, vjust = 0.5, size = 2.5, family = "GillSans", color = "white") +
  # labs( =) +
  scale_x_continuous( "Number of sequences", breaks = c(0,1, 2,4),
    labels = scales::comma_format(suffix = "K")) +
  # scale_x_continuous( "Number of sequences", labels = scales::scientific_format()) +
  theme_bw(base_size = 10, base_family = "GillSans") +
  theme(axis.text.x = element_text(size = 7),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 0))

pleft

# filter conoSorter only to Score_df == 3 and mature
#  matches for only the pro and mature regions (as described in Lavergne et al. (2013) and Barghi et al., 2016). 

superfm_df <- read_rds(paste0(pub_dir, "/superfm_df.rds")) %>% ungroup() %>% filter(!grepl("Mmseq", Method))

superfm_df %>% count(Region)

prodf <- superfm_df %>% filter(Region == "(Mature)-(Pro-region)") %>% mutate(Region = "Propeptide")
Maturedf <- superfm_df %>% filter(Region == "(Mature)") %>% mutate(Region = "Mature")
precursordf <- superfm_df %>% filter(Score_sf == 3) %>% mutate(Region = "Precursor")

superfm_df <- Maturedf %>% rbind(prodf, precursordf)

# Count number of transcripts with several orfs 

superfm_viz <- superfm_df %>%
  # if use gene level
  select(-gene) %>% distinct() %>%
  count(Method, Region, tab) %>% ungroup() %>% dplyr::rename("n_seqs" = "n")

# Count number of unique superfamilies (sf) annotated

superfm_viz <- superfm_df %>%
  ungroup() %>% distinct(Method, tab, Score_sf, Superfamily, Region) %>%
  count(Method, Region, tab) %>% ungroup() %>% dplyr::rename("n_sf" = "n") %>%
  left_join(superfm_viz)

# example
#  Trinity (T) method, assembly 10,978 conopetides (i.e 10,288 genes with 10,978 orfs) ranging 8 (Mature) superfamilies
#
superfm_df %>% filter(Method == "Trinity (T)" & tab == "Regex" & Region == "Mature") %>% count(Superfamily)

superfm_viz %>%
 group_by(Method) %>% summarise(n = sum(n_sf)) %>% arrange(desc(n)) %>% pull(Method) -> method_levels
  

recode_to <- structure(c("A) Precursor", "B) Propeptide", "C) Mature"), names = c("Precursor", "Propeptide", "Mature"))

pright <- superfm_viz %>%
  mutate(Region = dplyr::recode_factor(Region, !!!recode_to)) %>%
  mutate(label = paste0(scales::comma(n_seqs), " (", scales::comma(n_sf),")")) %>%
  group_by(Method) %>% mutate(frac = n_seqs/sum(n_seqs)) %>%
  mutate(Method = factor(Method, levels = rev(assembly_lev))) %>%
  # mutate(Region = factor(Region, levels = reg_lev)) %>%
  ggplot(aes(y = Method, x = n_seqs)) + 
  geom_col(position = position_dodge2(reverse = T), fill = "black") +
  facet_grid(~ Region + tab, scales = "free_x") +
  scale_fill_manual(values = scale_color) +
  ggh4x::facet_nested(~ Region + tab, scales = "free_x") +
  labs(y = "") +
  scale_x_continuous( "Number of ORFs annotated (Number of superfamily classes)", labels = scales::comma_format()) +
  theme_bw(base_size = 10, base_family = "GillSans") +
  geom_text(aes(label= label), hjust= 1.2, vjust = 0.5, size = 2.5, family = "GillSans", color = "white") +
  theme(
    axis.text.x = element_text(size = 7),
    panel.grid.minor = element_blank(), 
    axis.text.y = element_blank(), 
    # axis.title.x = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 0))



library(patchwork)


psave <- pleft + plot_spacer() + pright + plot_layout(widths = c(1, -0.05, 2))

ggsave(psave, filename = 'ConoSorter_.png', path = pub_dir, width = 10, height = 2, device = png, dpi = 700)


heatmapdf <- superfm_viz %>%
  mutate(label = paste0(scales::comma(n_seqs), " (", scales::comma(n_sf),")")) %>%
  group_by(Region,tab) %>%
  mutate(frac = n_seqs/max(n_seqs)) %>%
  mutate(Method = factor(Method, levels = rev(assembly_lev))) 

lo = floor(min(heatmapdf$frac))
up = ceiling(max(heatmapdf$frac))
mid = (lo + up)/2

heatmapdf %>% group_by(Method) %>% summarise(n_sf = sum(n_sf), n_seqs = sum(n_seqs))

pheatmap <- heatmapdf %>%
  mutate(topfacet = "B) Conopeptide annotation") %>%
  mutate(tab = factor(tab, levels = c("Regex", "pHMM"))) %>%
  mutate(Region = factor(Region, levels = c("Precursor", "Propeptide","Mature"))) %>%
  ggplot(aes(y = Method, x = Region, fill = frac)) +
  ggh4x::facet_nested(~ topfacet+tab, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 10, base_family = "GillSans") +
  scale_fill_gradient2(low = "white", high = "black", mid = "gray40",
    na.value = "white", midpoint = mid, limit = c(lo, up),
    name = NULL) +
  labs(x = "Superfamily classes", y = "") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 7),
    panel.grid.minor = element_blank(), 
    axis.text.y = element_blank(), 
    # axis.title.x = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 0)) +
  geom_tile(color = "white",
    lwd = 0.5,
    linetype = 1) +
  geom_text(aes(label= label), hjust= 0.5, vjust = 0.5, size = 2, family = "GillSans", color = "white") 

psave <- pleft + plot_spacer() + pheatmap + plot_layout(widths = c(1.5, -0.01, 1.5))

ggsave(psave, filename = 'ConoSorter_heatmap.png', path = pub_dir, width = 7.2, height = 2.2, device = png, dpi = 700)

# BUSCO =====

col <- c("#ED4647", "#EFE252", "#3A93E5", "#5BB5E7")

names <- c("Complete", "Duplicated", "Fragmented", "Missing")

#problem reading empty records

col <- structure(col, names = rev(names))

BUSCOdf <- read_tsv(
  paste0(pub_dir, "/busco_summaries.tsv")) %>%
  filter(!grepl("Mmseq", Method))


names <- c("Complete", "Duplicated", "Fragmented", "Missing")

labels <- c("Complete (C) and single-copy (S)",
  "Complete (C) and duplicated (D)",
  "Fragmented (F)  ",
  "Missing (M)")


labels <- structure(labels, names = names)

col <- structure(col, names = rev(names))

my_sp_lev <- c("Mollusca","Metazoa","Eukaryota","Mammalia","Bacteria")

figure <- BUSCOdf %>%
  # filter(Status != "Missing") %>%
  group_by(Method, odb) %>% mutate(percentage = n/sum(n)) %>%
  mutate(facet = "Completeness") %>%
  filter(!odb %in% c("Bacteria", "Mammalia")) %>%
  mutate(Status = factor(Status, levels = rev(names))) %>%
  mutate(odb = factor(odb, levels = my_sp_lev)) %>%
  ggplot(aes(x = Method, y = percentage, fill = Status)) +
  facet_grid(odb~ facet) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.75) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(y = "% BUSCOs", x = "Assembly method") +
  coord_flip() +
  scale_fill_manual("", values = col, labels = rev(labels)) +
  guides(fill=guide_legend(nrow = 4)) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "bottom", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank())

figure

BUSCOdf %>%
  count(Method)

BUSCOdf %>%
  filter(Status != "Missing") %>%
  filter(!odb %in% c("Bacteria", "Mammalia")) %>%
  mutate(facet = ifelse(grepl("Full-length", Method),"Post-assembly step", "Only assembly step")) %>%
  mutate(group = "Full-length T-S") %>%
  mutate(group = ifelse(grepl("Trinity", Method),"Trinity (T)", group)) %>%
  mutate(group = ifelse(grepl("Spades", Method),"Spades (S)", group)) %>%
  mutate(label = ifelse(Method %in% "Full-length T-S (Mmseq)", "Mmseq", NA)) %>%
  mutate(label = ifelse(Method %in% "Full-length T-S (Lace)","Lace", label)) %>% 
  # mutate(group2 = ifelse(grepl("Lace", Method),"Lace", group2)) %>%
  # mutate(group = factor(group, levels = rev(c("Trinity", "Spades", "Lace", "Mmseq")))) %>%
  filter(!Status %in% c("Missing")) %>%
  mutate(Status = dplyr::recode_factor(Status, !!!labels)) %>%
  mutate(odb = factor(odb, levels = my_sp_lev)) %>%
  group_by(Method, odb) %>% mutate(percentage = n/sum(n)) %>%
  ggplot(aes(x = percentage, y = group, fill = facet)) +
  facet_grid(odb~Status, scales = "free_x") +
  # ggh4x::facet_nested(odb~Status, scales = "free_x", switch = "x", ) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_text(aes(label= label, group = group), 
    position = position_dodge2(width = 0.9),
    hjust= 1.2,
    vjust = 0.5,
    size = 1.5, family = "GillSans", color = "white") +
  scale_x_continuous(labels = scales::percent_format(scale = 100)) +
  labs(x = "% Completeness (BUSCOs)", y = "Assembly method") +
  scale_fill_manual("", values = c("gray90","grey20")) +
  guides(fill=guide_legend(nrow = 1)) +
  theme_bw(base_size = 7, base_family = "GillSans") +
  theme(legend.position = "top", 
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.12, "cm"),
    panel.grid.minor = element_blank(), 
    strip.background = element_rect(fill = 'white', color = 'white'), 
    strip.placement = "outside",
    strip.text = element_text(color = "black",hjust = 0, size = 5),
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()) -> psave

psave

ggsave(psave, filename = 'BUSCO_C_D_F.png', path = pub_dir, width = 4, height = 2.5, device = png, dpi = 700)
