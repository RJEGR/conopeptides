# UPSETPLOT
# Count diversity of sig. up-degs per diet-group
# Contrast against global diversity
# Include a sub-group of conotoxins per gene superfamily
# Include the upset version of kallisto cds pipeline
# describe publication


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

# dir <- "C://Users//cinai/OneDrive/Documentos/PUBLICATION_DIR/06.Quantification/MATRIX_RSEM_dir/"

# f <- "glmLRT_multiple_contrast_ctrl_and_treatments_kallisto.rds"

dir <- "C://Users//cinai/OneDrive/Documentos/PUBLICATION_DIR/03.Coverage/"

Manifest <- list.files(dir, pattern = "Manifest.tsv", full.names = T)

Manifest <- read_tsv(Manifest) %>% distinct() %>%
  mutate(LIBRARY_ID = ifelse(grepl("cam2v_", LIBRARY_ID), NA, LIBRARY_ID)) %>%
  mutate(LIBRARY_ID = ifelse(grepl("cam6", LIBRARY_ID), NA, LIBRARY_ID)) %>%
  drop_na()


pub_dir <- "C://Users//cinai/OneDrive/Documentos/PUBLICATION_DIR/"

CONOPEPDB <- read_tsv(paste0(pub_dir, "/conopeptides.tsv")) %>%  #view()
  mutate(len = nchar(pep_seq)-1) %>%
  filter(effective_length_frac > 0.5) %>%
  # To be consistent w/ RES
  filter(Signalp_class == "SP") %>%
  drop_na(prediction_tool, Superfamily) 


# 1. LOAD DEGS
# 

# recode_to <- structure(c("Control","Shrimp", "Mollusk", "Polychaete", "Mixed"))

dir <- "C://Users//cinai/OneDrive/Documentos/PUBLICATION_DIR/06.Quantification/MATRIX_RSEM_dir/"

f <- "glmLRT_multiple_contrast_ctrl_and_treatments_kallisto.rds"

DEGS <- read_rds(file.path(dir, f)) %>% filter(abs(logFC) > 2 & FDR < 0.05) 

DEGS <- DEGS %>% left_join(CONOPEPDB) %>% drop_na(Superfamily)

DEGS |> distinct(pep_seq, Superfamily, sam_group, sampleX) |> count(sampleX)

recode_time <- structure(c("0 months","2 months", "4 months"))

recode_time <- structure(recode_time, names = c("Ctrl","2","4"))

recode_to <- structure(c("Control","Shrimp", "Mollusk", "Polychaete", "Mixed"))

recode_to <- structure(recode_to, names = c("Ctrl","Cam","Lit","Pol","Mix"))


dat <- DEGS |> 
  distinct(pep_seq, Superfamily, sam_group, tab, sampleX) |>
  # mutate(Time = sapply(strsplit(sampleX, "_"), `[`, 2)) |>
  # mutate(Time = ifelse(is.na(Time), "Ctrl", Time)) |>
  # dplyr::mutate(Time = dplyr::recode_factor(Time, !!!recode_time)) |>
  # rename("Feeding" = sam_group) 
  mutate(Feeding = sapply(strsplit(sampleX, "_"), `[`, 1)) |>
  dplyr::mutate(Feeding = dplyr::recode_factor(Feeding, !!!recode_to))

dat |> count(Feeding)

dat |> 
  distinct(pep_seq, Superfamily, tab, Feeding)  |>
  group_by(pep_seq, tab, Superfamily) |>
  summarise(
    combination = paste(sort(unique(Feeding)), collapse = ","),
    n = n_distinct(Feeding),
    .groups = "drop"
  ) |>
  write_csv(file.path(pub_dir, "S2_table_degs_conotoxin_expressión.csv"))


# 2. LOAD GLOBA

# dat <- read_rds(file.path(pub_dir, "S1_table_global_conotoxin_expressión.rds")) |> ungroup()

# Create combinations of assemblers for each hit
# 

hits_combinations <- dat |>
  distinct(pep_seq, Superfamily, Feeding)  |>
  group_by(pep_seq) |>
  summarise(
    combination = paste(sort(unique(Feeding)), collapse = ","),
    n = n_distinct(Feeding),
    .groups = "drop"
  ) |>
  count(combination, name = "count")
  


# Reorder by count (descending)
hits_combinations <- hits_combinations |>
  mutate(combination = fct_reorder(combination, count, .desc = TRUE))

# Get assembler counts
assembler_counts <- dat |>
  distinct(pep_seq, Superfamily, Feeding)  |>
  count(Feeding, name = "hits_count") |>
  arrange(hits_count) |>
  mutate(Feeding = fct_reorder(Feeding, hits_count, .desc = TRUE))


# Create intersection points data (assembler x combination)
points_data <- hits_combinations |>
  mutate(
    Feeding = str_split(combination, ",")
  ) |>
  unnest(Feeding) |>
  mutate(
    Feeding = factor(Feeding, levels = levels(assembler_counts$Feeding)),
    combination = factor(combination, levels = levels(hits_combinations$combination))
  ) |>
  filter(Feeding != "")

combination_lev <- rev(levels(hits_combinations$combination))

# 1. Main bar chart (hits per combination)
bar_chart <- hits_combinations |>
  mutate(label = paste0(" (", count,")")) |>
  mutate(col = ifelse( grepl(",", combination), "Intersect", "Unique")) |>
  ggplot(aes(y = combination, x = count)) +
  # facet_grid(~ Time) +
  geom_col(width = 0.6, aes(fill = col)) +
  geom_text(aes(label = label),
            vjust = 0.5, hjust = -0.15, size= 2.5,
            color="black",
            # position=position_dodge(0.5),
            family =  "Gill Sans MT") +
  # scale_x_discrete(position = "top", drop = FALSE) + 
  scale_y_discrete(limits = combination_lev, labels = NULL) +
  scale_x_continuous(position = "top", limits = c(0, 750), breaks = c(0,300,600), labels = NULL) +
  my_custom_theme() +
  scale_fill_manual(values = c("gray90", "black")) +
  theme(
    legend.position = "none",
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x.top = element_text(angle = 0, vjust = -0.5, hjust = 1)
  ) +
  coord_cartesian(expand = FALSE) +
  labs(x = element_blank(), y = element_blank())

# 2. Intersection point chart (shows which assemblers in each combination)
# 

point_chart <- points_data |>
  mutate(col = ifelse( grepl(",", combination), "Intersect", "Unique")) |>
  mutate(facets = "Intersections") |>
  ggplot(aes(y = combination, x = Feeding)) +
  geom_line(aes(group = combination, colour = col), size = 1) +
  geom_point(aes(colour = col), size = 3.5) +
  # facet_grid(facets ~ ., switch = "y") +
  my_custom_theme() +
  ggstats::geom_stripped_cols() +
  scale_x_discrete(position = "top", drop = FALSE) + 
  scale_y_discrete(limits = combination_lev) +
  scale_color_manual(values = c("gray90", "black")) +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x.top = element_text(angle = -45, vjust = -0.5, hjust = 1)
  ) +
  # coord_cartesian(expand = FALSE) +
  labs(x = element_blank(), y = element_blank()) 




# 3. Side bar chart (total hits per assembler)

assembler_bars <- assembler_counts |>
  ggplot(aes(y = hits_count, x = Feeding)) +
  geom_col(width = 0.6, fill = 'gray90') +
  # geom_segment(aes(x = Assembler, y = 0, yend = -hits_count), size = 12) +
  geom_text(aes(label = Feeding), size = 1.5, hjust = 1.5, vjust = 0.5, angle = -90, color = "white") +
  scale_y_reverse() +
  # ggstats::geom_stripped_cols() +
  my_custom_theme() +
  theme(
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
  ) +
  coord_cartesian(expand = FALSE) +
  labs(x = element_blank(), y = element_blank())

library(patchwork)

design <- "#AC
           #B#"


PSAVE <- wrap_plots(C = bar_chart, 
                    B = assembler_bars, 
                    A = point_chart, design = design) +
  plot_layout(heights  = c(0.65,0.15,1))

PSAVE

ggsave(PSAVE, filename = 'UPSET_FOR_PUB_DEGS.png', path = pub_dir, width = 5, height = 7, device = png, dpi = 800)

#
#
#

