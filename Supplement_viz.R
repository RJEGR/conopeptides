# Supplementary material for pub
# LOAD data.frame from Assembly Nx metrics for DNA and orf (transdecoder.Orfs) [Nx_metrix.R] for every assembly method
# LOAD BUSCO completeness for every assembly method (BUSCO.R)
# LOAD superfm_df data.frame from ConoSorter_viz.R 

# PLOT lineplot of Nx facet DNA from ORF
# PLOT BUSCO scores by complete (single + duplicate) ~ assembly method (trinity, spades, concat), coloring by assembly step (raw vs full-length)
# PLOT summary of known/novel Superfamily classes detected by methods. Including (Mature)-(Pro-region)-(Signal) class 

pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR"

superfm_df <- read_rds(paste0(pub_dir, "/superfm_df.rds")) %>% ungroup()

superfm_viz <- superfm_df %>%
  # if use gene level
  # select(-transcript) %>% distinct() %>%
  count(Method, Region, tab) %>% ungroup() %>% dplyr::rename("n_transcripts" = "n")

superfm_viz <- superfm_df %>%
  ungroup() %>% distinct(Method, tab, Score_sf, Superfamily, Region) %>%
  count(Method, Region, tab) %>% ungroup() %>% dplyr::rename("n_sf" = "n") %>%
  left_join(superfm_viz)

superfm_viz %>%
  filter(tab == "Regex") %>%
  mutate(label = paste0(scales::comma(n_transcripts), " (", scales::comma(n_sf),")")) %>%
  group_by(Method) %>% mutate(frac = n_transcripts/sum(n_transcripts)) %>%
  mutate(Method = factor(Method, levels = method_levels)) %>%
  mutate(Region = factor(Region, levels = reg_lev)) %>%
  ggplot(aes(y = Method, x = n_transcripts)) + 
  geom_col(position = position_dodge2(reverse = T),fill = "black") +
  facet_grid(~ Region, scales = "free_x") +
  scale_x_continuous( "Fraction (Transcript annotated/Assembled transcripts)", labels = scales::comma_format()) +
  theme_bw(base_size = 14, base_family = "GillSans") +
  geom_text(aes(label= label), hjust= 1.2, vjust = 0.5, size = 5, family = "GillSans", color = "white")
