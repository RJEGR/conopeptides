# Read EMAPPER
# Wrangling NOGs
# Save GOs by gene_id
# Count NOGs
# Find Disulfide isomerase enzymes


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)


library(tidyverse)

pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR"

dir <- "~/Documents/GitHub/conopeptides/05.Prediction//EMAPPER_DIR/"

f <- list.files(dir, pattern = "eggnog_mapper.emapper.annotations", full.names = T)

MAPPER_DB <- read_tsv(f, comment = "##") %>% dplyr::rename("protein_id" = "#query")

MAPPER_DB <- MAPPER_DB %>% mutate(gene_id=gsub(".p[0-9]+$", "", protein_id)) %>% select(protein_id, gene_id) %>% left_join(MAPPER_DB)

MAPPER_DB %>% count(protein_id, sort = T) # 28,264 proteins

MAPPER_DB %>% count(gene_id, sort = T) # 24,255 transcripts

# MAPPER_DB %>% filter(gene_id %in% "Cluster-15813.95759") %>% select(protein_id, gene_id)

# Bind to the COG_category

NOG.col <- c('J@Translation, ribosomal structure and biogenesis','A@RNA processing and modification','K@Transcription','L@Replication, recombination and repair','B@Chromatin structure and dynamics','D@Cell cycle control, cell division, chromosome partitioning','Y@Nuclear structure','V@Defense mechanisms','T@Signal transduction mechanisms','M@Cell wall/membrane/envelope biogenesis','N@Cell motility','Z@Cytoskeleton','W@Extracellular structures','U@Intracellular trafficking, secretion, and vesicular transport','O@Posttranslational modification, protein turnover, chaperones','X@Mobilome: prophages, transposons','C@Energy production and conversion','G@Carbohydrate transport and metabolism','E@Amino acid transport and metabolism','F@Nucleotide transport and metabolism','H@Coenzyme transport and metabolism','I@Lipid transport and metabolism','P@Inorganic ion transport and metabolism','Q@Secondary metabolites biosynthesis, transport and catabolism','R@General function prediction only','S@Function unknown')

into <- c("EggNM.COG_category", "COG_name")

NOG.col <- data.frame(NOG.col) %>% separate(NOG.col, sep = "@", into = into)

# Category

MAPPER_DB %>% count(COG_category,sort = T) %>% tail() # other categories?

# MAPPER_DB <- MAPPER_DB %>% 
#   left_join(NOG.col, by = c("COG_category" = "EggNM.COG_category")) %>%
#   mutate(COG_category = paste0(COG_category, ", ", COG_name))



plotdf <- MAPPER_DB %>%
  filter(COG_category != "-") %>%
  right_join(NOG.col, by = c("COG_category" = "EggNM.COG_category")) %>%
  mutate(COG_category = paste0(COG_category, ", ", COG_name)) %>%
  dplyr::count(COG_category, sort = T) %>%
  mutate(frac = n/sum(n)) %>%
  arrange(desc(frac)) %>%
  mutate(COG_category = factor(COG_category, levels = unique(COG_category))) %>%
  mutate(facet = "A) Transcriptome")

MAPPER_DB %>%
  filter(if_any(where(is.character), ~ grepl(pattern = 'Thioredoxin', x = .x, ignore.case = T))) %>%
  filter(COG_category != "-") %>%
  right_join(NOG.col, by = c("COG_category" = "EggNM.COG_category")) %>%
  mutate(COG_category = paste0(COG_category, ", ", COG_name)) %>%
  dplyr::count(COG_category, sort = T) %>%
  mutate(frac = n/sum(n), facet = "B) Thioredoxin (PDI family)") 
  # rbind(plotdf)

data_text <- plotdf %>% 
  # group_by(COG_name, facet) %>% summarise(n = sum(n)) %>%
  mutate(n = paste0("(", n, ")")) %>%
  mutate(n = ifelse(grepl("^O,", COG_category), paste0(n, " Thioredoxin (PDI family)"), n))

plotdf %>% 
  ggplot(aes(y = COG_category, x = frac, fill = facet)) +
  labs(y = "Nested Orthologous Gene Group (NOGS)", x = "Enrichment ratio (N transcripts/Total transcripts)") +
  facet_grid(~ facet, scales = "free_x", space = "free_x", switch = "y") +
  geom_col(fill = "black") +
  scale_x_continuous(labels = scales::percent_format()) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  xlim(0,0.4) +
  theme(legend.position = "top", 
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    # panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = 'grey95', color = 'white')) +
  geom_text(data = data_text, 
    aes(label = n), size = 5,
    hjust = -0.1, vjust = 0, 
    family = "GillSans", position = position_dodge(width = 1)) -> p

p

ggsave(p, filename = 'NOGS.png', path = pub_dir, width = 8, height = 7, device = png, dpi = 300)

# Save output for Database generation

eggNOG_cols <- c("protein_id", "Preferred_name","Description", "GO","PFAMs","BRITE", "CAZy", "COG_category")

MAPPER_DB <- MAPPER_DB %>% select_at(vars(contains(eggNOG_cols), starts_with("KEGG"))) 

# Thioredoxin (protein disulfide isomerase family)

MAPPER_DB %>% 
  filter(if_any(where(is.character), ~ grepl(pattern = 'Thioredoxin', x = .x, ignore.case = T))) %>%
  

