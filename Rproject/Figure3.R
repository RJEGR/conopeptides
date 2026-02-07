
# Estimate Duct and venom vulb of C californicus transcriptome annotation of expressed conotoxins
# 
# Calculate :
# 
# Use CWE (Method 1) if you believe that "more isoforms = more functional importance."
# 
# Use Shannon (Method 2) if you want to identify gene groups that are undergoing intense Alternative Splicing regulation.
# 
# Use Simple Sum if you only care about total transcriptional output, regardless of how it's packaged.
# 

# LOAD COUNT and CONOSERVER DB
# AGGLOMERATE 
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

time_levs <- c("Ctrl", "2", "4", "6")

recode_time <- structure(c("Control", "2 months", "4 months", "6 months"), names = time_levs)

Diatery_levs <- c("Ctrl","Cam", "Lit", "Pol", "Mix")

recode_Diatery <- structure(c("Control","Shrimp", "Mollusk", "Polychaete", "Mixed"), names = Diatery_levs)


library(tidyverse)

pub_dir <- "C://Users//cinai/OneDrive/Documentos/PUBLICATION_DIR/"

# DB <- read_rds(paste0(pub_dir, "/structured_db.rds"))
# 
# CONOPEPDB %>% 
#   filter(effective_length_frac > 0.5) %>%
#   count(Superfamily,sort = T)

# CONOPEPDB %>% ggplot(aes(effective_length_frac )) + geom_histogram()

CONOPEPDB <- read_tsv(paste0(pub_dir, "/conopeptides.tsv")) %>%  #view()
  mutate(len = nchar(pep_seq)-1) %>%
  filter(effective_length_frac > 0.5) %>%
  # To be consistent w/ RES
  filter(Signalp_class == "SP") %>%
  drop_na(prediction_tool, Superfamily) 

CONOPEPDB %>% count(Superfamily, sort = T) %>% View

dir <- "C://Users//cinai/OneDrive/Documentos/PUBLICATION_DIR/06.Quantification/MATRIX_RSEM_dir/KALLISTO_Merged_polyA_hisat_SuperDuper.fasta.transdecoder_DIR/"

f <- list.files(dir, pattern = "KALLISTO_Merged_polyA_hisat_SuperDuper.fasta.transdecoder.matrix", full.names = T)

COUNTS <- read_rds(f)

keep <- rownames(COUNTS) %in% CONOPEPDB$protein_id

sum(keep)

dim(COUNTS <- COUNTS[keep,])

# normalize
# 

# COUNTS <- edgeR::cpm(COUNTS)

# file 3
# 
# 
# 
dir <- "C://Users//cinai/OneDrive/Documentos/PUBLICATION_DIR/03.Coverage/"

Manifest <- list.files(dir, pattern = "Manifest.tsv", full.names = T)

Manifest <- read_tsv(Manifest) %>% distinct() %>%
  mutate(LIBRARY_ID = ifelse(grepl("cam2v_", LIBRARY_ID), NA, LIBRARY_ID)) %>%
  mutate(LIBRARY_ID = ifelse(grepl("cam6", LIBRARY_ID), NA, LIBRARY_ID)) %>%
  drop_na()

barvizA <- CONOPEPDB %>% 
  count(Superfamily, tab, sort = T)

barvizB <- COUNTS %>%
  as_tibble(rownames = "protein_id") %>%
  left_join(distinct(CONOPEPDB, Superfamily, protein_id, tab)) %>%
  group_by(Superfamily, tab) %>%
  # summarise_at(vars(all_of(colnames(COUNTS))), sum) %>% ungroup() %>%
  pivot_longer(cols = all_of(colnames(COUNTS)), values_to = 'fill', names_to = "LIBRARY_ID") %>%
  right_join(Manifest, by = "LIBRARY_ID") %>% filter(fill > 0)

barvizB %>% 
  count(Superfamily, tab, sort = T)

barvizB %>% group_by(Feeding) %>% tally(fill)

# calculate the diversity-abundanxe index (shannon-based)
# 
# 
# 3. Calcular el Índice de Diversidad y el Peso Combinado
shannon_results <- barvizB %>%
  # Agrupamos por gen para que los cálculos sean relativos a cada gen
  group_by(Superfamily,tab,Feeding ) %>% # Feeding
  mutate(
    # Calcular la proporción (pi) de cada isoforma dentro de su gen
    total_gene_abundance = sum(fill),
    p_i = fill / total_gene_abundance
  ) %>%
  # ggplot(aes(p_i, fill)) + geom_point(aes())
  # Filtramos abundancias p_i para evitar errores con el logaritmo
  filter(p_i > 0.01) %>%
  summarise(
    # 1. Calculamos el Índice de Shannon (H)
    # H = -sum(pi * ln(pi))
    shannon_index = -sum(p_i * log(p_i)),
    
    # 2. Guardamos la abundancia total del gen
    total_abundance = first(total_gene_abundance),
    
    # 3. Calculamos el Weighted Index (Diversidad x Abundancia)
    weight = shannon_index * total_abundance,
    
    # 4. Contamos cuántas isoformas contribuyeron (opcional)
    n = n()
  ) %>%
  # Ordenar por el peso final para ver los genes más "relevantes"
  arrange(desc(weight)) 

threshold_value <- quantile(shannon_results$weight, probs = 0.5)

shannon_results <- shannon_results %>%
  ungroup() %>%
  mutate(facet = "B) Conotoxin expression") %>% 
  # filter(total_abundance > 100) %>% 
  # --- NEW STEP: Create a grouping variable for the facet ---
  mutate(range_group = ifelse(weight > threshold_value, "High Abundance", "Low Abundance")) %>%
  # Optional: Set levels so "Low" appears on the left and "High" on the right
  mutate(range_group = factor(range_group, levels = c("Low Abundance", "High Abundance"))) 

shannon_results %>% 
  count(range_group)

sf_levels <- shannon_results %>%
  group_by(Superfamily) %>% summarise(weight = sum(weight)) %>%
  arrange(desc(weight)) %>%
  distinct(Superfamily) %>% pull()
  
text_df <- shannon_results %>%
  group_by(Superfamily) %>% 
  summarise(n = sum(n), weight = sum(weight)) %>%
  mutate(label = paste0("(", n,")"))
  
p <- shannon_results %>% 
  mutate(Superfamily = factor(Superfamily, levels = rev(sf_levels))) %>%
  ggplot(aes(y = Superfamily, x = weight   )) +
  facet_grid(~ facet, scales = "free_x", space = "free") +
  geom_col(aes(fill = tab)) +
  scale_x_continuous("Weight (Total Expression x N isoforms)",labels = scales::comma_format(scale = 1E-6, suffix = "M"), limits = c(0,2E6)) +
  ylab("Gene Superfamily") +
  theme_bw(base_size = 12, base_family = base_text_fam) +
  scale_fill_grey("") +
  geom_text(data = text_df, aes(label = label), size = 2.5,
            hjust = -0.1, vjust = 0.25,
            family = base_text_fam, position = position_dodge(width = 1)) +
  my_custom_theme() +
  theme(panel.spacing.x = unit(1, "cm"))

p

ggsave(p, filename = 'Superfamilies.png', path = pub_dir, width = 5, height = 6, device = png, dpi = 500)


