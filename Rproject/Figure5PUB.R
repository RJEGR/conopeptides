# 
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

dir <- "C://Users//cinai/OneDrive/Documentos/PUBLICATION_DIR/06.Quantification/MATRIX_RSEM_dir/"


f <- "glmLRT_multiple_contrast_ctrl_and_treatments_kallisto.rds"


dir <- "C://Users//cinai/OneDrive/Documentos/PUBLICATION_DIR/03.Coverage/"

Manifest <- list.files(dir, pattern = "Manifest.tsv", full.names = T)

Manifest <- read_tsv(Manifest) %>% distinct() %>%
  mutate(LIBRARY_ID = ifelse(grepl("cam2v_", LIBRARY_ID), NA, LIBRARY_ID)) %>%
  mutate(LIBRARY_ID = ifelse(grepl("cam6", LIBRARY_ID), NA, LIBRARY_ID)) %>%
  drop_na()