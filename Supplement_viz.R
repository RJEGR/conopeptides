# Supplementary material for pub
# LOAD data.frame from Assembly Nx metrics for DNA and orf (transdecoder.Orfs) [Nx_metrix.R] for every assembly method
# LOAD BUSCO completeness for every assembly method (BUSCO.R)
# LOAD superfm_df data.frame from ConoSorter_viz.R 

# PLOT lineplot of Nx facet DNA from ORF
# PLOT BUSCO scores by complete (single + duplicate) ~ assembly method (trinity, spades, concat), coloring by assembly step (raw vs full-length)
# PLOT summary of known/novel Superfamily classes detected by methods. Including (Mature)-(Pro-region)-(Signal) class 


rm(list = ls())

if(!is.null(dev.list())) dev.off()

library(tidyverse, help, pos = 2, lib.loc = NULL)

pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR"

superfm_df <- read_rds(paste0(pub_dir, "/superfm_df.rds")) %>% ungroup()

Nxmetrix_df <- read_rds(paste0(pub_dir, "/Contig_nx.rds"))

superfm_df %>% distinct(Method)

assembly_methods <- c(
  "Trinity (T)", 
  "Trinity-Lace (Hisat)", 
  "Spades (S)", 
  "Spades-Lace (Hisat)", 
  "Concat-Lace T-S (Hisat-polA)", 
  "MMseqs (T-S)"
  # "MMseqs-Lace T-S (Hisat)",
  # "Concat-Lace T-S (Hisat-polA-transdecoder)"
  )


recode_to <- c(
  "Trinity (T)", 
  "Full-length Trinity", 
  "Spades (S)", 
  "Full-length Spades", 
  "Full-length T-S (Lace)", 
  "Full-length T-S (Mmseq)"
  # "Full-length T-S (Mmseq-lace)",
  # "Concat-Lace T-S (Hisat-polA-transdecoder)"
  )

recode_to <- structure(
  recode_to,
  names = assembly_methods)

# Scale colors

library("scales")

show_col(ggsci::pal_igv("default")(6))

scale_color <- structure(ggsci::pal_igv("default")(6), names = recode_to)

# Recode datasets ----

Nxmetrix_df <- Nxmetrix_df %>% 
  mutate(Assembly = dplyr::recode_factor(Assembly, !!!recode_to)) 

superfm_df <- 
  superfm_df %>%
  # omit Concat-Lace T-S (Hisat-polA-transdecoder) cause here I intent to compare DNA assembly 
  filter(Method %in% assembly_methods) %>%
  mutate(Method = dplyr::recode_factor(Method, !!!recode_to))


# Plot Nx

p1 <- ggplot(Nxmetrix_df, aes(x = x, y = n, group = Assembly, color = Assembly)) +
  geom_vline(xintercept = "N50", linetype="dashed", alpha=0.5) +
  ggplot2::geom_path(linewidth = 1, lineend = "round") +
  facet_wrap(~  stringSet, scales = "free_y") +
  # geom_point(shape = 21, size = 4, aes(size = n_frac)) +
  labs(x = "Nx", y = "Contig length", color = "Assembly method") +
  # ggsci::scale_color_igv() +
  scale_color_manual(values = scale_color) +
  scale_y_continuous(labels = scales::comma_format()) +
  theme_bw(base_size = 14, base_family = "GillSans") +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 0)) + 
  guides(color = guide_legend(title = "", nrow = 2, ncol = 4))

# ggsave(p1, filename = 'Nx-methods.png', path = pub_dir, width = 10, height = 5, device = png, dpi = 300)

# plot only the number of nx

Nxmetrix_df %>%
  filter(x %in% "N90") %>%
  ggplot(aes(y = Assembly, x = n_seqs)) +
  geom_col() +
  facet_wrap(~  stringSet, scales = "free_y") 

# Nxmetrix_df %>%
#   filter(stringSet %in% "B) Complete ORFs (Transcriptome)") %>%
#   group_by(Assembly) %>% summarise(n = sum(n_seqs)) %>% arrange(desc(n)) %>%
  

# add busco

# filter conoSorter only to Score_df == 3 and mature
#  matches for only the pro and mature regions (as described in Lavergne et al. (2013) and Barghi et al., 2016). 

superfm_df %>% count(Region)

prodf <- superfm_df %>% filter(Region == "(Mature)-(Pro-region)") %>% mutate(Region = "B) Pro-region-mature")
Maturedf <- superfm_df %>% filter(Region == "(Mature)") %>% mutate(Region = "C) Mature")
precursordf <- superfm_df %>% filter(Score_sf == 3) %>% mutate(Region = "A) Precursor")

superfm_df <- Maturedf %>% rbind(prodf, precursordf)

# Count number of transcripts with several orfs 

superfm_viz <- superfm_df %>%
  # if use gene level
  select(-transcript) %>% distinct() %>%
  count(Method, Region, tab) %>% ungroup() %>% dplyr::rename("n_transcripts" = "n")

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
  
p2 <- superfm_viz %>%
  mutate(label = paste0(scales::comma(n_transcripts), " (", scales::comma(n_sf),")")) %>%
  group_by(Method) %>% mutate(frac = n_transcripts/sum(n_transcripts)) %>%
  # mutate(Method = factor(Method, levels = method_levels)) %>%
  # mutate(Region = factor(Region, levels = reg_lev)) %>%
  ggplot(aes(y = Method, x = n_transcripts)) + 
  geom_col(position = position_dodge2(reverse = T),fill = "black") +
  facet_grid(~ Region + tab, scales = "free_x") +
  scale_fill_manual(values = scale_color) +
  ggh4x::facet_nested(~ Region + tab, scales = "free_x") +
  scale_x_continuous( "Transcripts annotated (Number of conopetide classes)", labels = scales::comma_format()) +
  theme_bw(base_size = 14, base_family = "GillSans") +
  geom_text(aes(label= label), hjust= 1.2, vjust = 0.5, size = 3.5, family = "GillSans", color = "white") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 0))

p2

ggsave(p2, filename = 'ConoSorter.png', path = pub_dir, width = 10, height = 3, device = png, dpi = 300)



heatmapdf <- superfm_viz %>%
  mutate(label = paste0(scales::comma(n_transcripts), " (", scales::comma(n_sf),")")) %>%
  group_by(Region) %>% mutate(frac = n_transcripts/sum(n_transcripts)) 
  
lo = floor(min(heatmapdf$frac))
up = ceiling(max(heatmapdf$frac))
mid = (lo + up)/2


p22 <- heatmapdf %>%
  ggplot(aes(y = Method, x = tab, fill = frac)) +
  facet_grid(~ Region) +
  theme_bw(base_size = 14, base_family = "GillSans") +
  geom_text(aes(label= label), hjust= 1.2, vjust = 0.5, size = 3.5, family = "GillSans", color = "white") +
  scale_fill_gradient2(low = "white", high = "blue", mid = "gray",
    na.value = "white", midpoint = 0.05, limit = c(0, 0.25),
    name = NULL) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 0)) +
  geom_tile(color = "white",
    lwd = 0.5,
    linetype = 1) +
  geom_text(aes(label= label), hjust= 0.5, vjust = 0.5, size = 3.5, family = "GillSans", color = "white") 


ggsave(p22, filename = 'ConoSorter_heat.png', path = pub_dir, width = 7, height = 3, device = png, dpi = 300)




library(patchwork)

p1/p2

