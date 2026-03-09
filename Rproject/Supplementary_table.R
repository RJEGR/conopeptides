



rm(list = ls())

if(!is.null(dev.list())) dev.off()

library(tidyverse)

outdir <- "C://Users//cinai/OneDrive/Documentos/PUBLICATION_DIR/Sup_tables_dir"

dir.create(outdir, showWarnings = FALSE)

# LOAD transrate
# 

f <- "//wsl.localhost/Debian/home/ricardo/contigs_transrate_conopeptides.csv"

read_csv(f) |> 
  ggplot(aes(score/tpm)) + geom_histogram()

# [ INFO] 2025-04-23 12:55:20 : TRANSRATE OPTIMAL SCORE      0.0241
# [ INFO] 2025-04-23 12:55:20 : TRANSRATE OPTIMAL CUTOFF     0.0431
# [ INFO] 2025-04-23 12:55:20 : good contigs                   1507 of 1725
# [ INFO] 2025-04-23 12:55:20 : p good contigs                 0.87

# Transrate also calculates an expression-weighted assembly score by multiplied by its relative expression before being included in the assembly score

transratedf <- read_csv(f) |> 
  # filter(score >  0.0431) |>
  # ggplot(aes(score/tpm)) + geom_histogram()
  dplyr::rename("gene_id" = contig_name) |>
  # mutate(protein_id = paste0(protein_id, ".p"))
  select(gene_id, score, p_bases_covered)


dir <- "C://Users//cinai/Downloads/"

f <- list.files(path = dir, pattern = "csv__", full.names = T)

sel_col <- c("Name","Signal","Pre","Mature","Class","Framework")

DB <- read_csv(f) |> select(any_of(sel_col)) |> 
  distinct(Name, Signal, Pre, Mature, Class, Framework) |>
  dplyr::rename("protein_id" = Name)

nrow(DB)


DB |>
  count(Class, sort = T)

pub_dir <- "C://Users//cinai/OneDrive/Documentos/PUBLICATION_DIR/"

read_tsv(paste0(pub_dir, "/conopeptides.tsv")) |> count(tab, Region)

DB <- read_tsv(paste0(pub_dir, "/conopeptides.tsv")) |>  #view()
  filter(effective_length_frac > 0.5) |>
  # To be consistent w/ RES
  filter(Signalp_class == "SP") |>
  drop_na(prediction_tool, Superfamily) |>
  # count(tab, Region)
  mutate(gene_id = gsub(".p[0-9]+$","", protein_id)) |> left_join(transratedf) |>
  distinct(protein_id, dna_seq, effective_length_frac, score, p_bases_covered,
           pep_seq, Signal_peptide_probs, 
           uniprotkb_toxprot, conoserver_protein, Superfamily, tab) |>
  dplyr::rename("primary_dna_seq" = dna_seq, "primary_pep_seq" = pep_seq) |>
  left_join(DB) 


# DB |> ggplot(aes(effective_length_frac, p_bases_covered)) + geom_point()

COUNT_DB <- read_rds(file.path(pub_dir, "S1_table_global_conotoxin_expressión.rds")) |> 
  group_by(protein_id, Feeding) |>
  summarise(fill = round(sum(fill))) |> ungroup() |> 
  # group_by(Feeding) |> tally(fill) |> tally(n)
  pivot_wider(values_from = fill, names_from = Feeding, values_fill = 0, names_glue = "{Feeding}_cpm")


read_rds(file.path(pub_dir, "S1_table_global_conotoxin_expressión.rds")) |> 
  group_by(protein_id) |>
  summarise(fill = round(sum(fill))) |>
  left_join(read_tsv(paste0(pub_dir, "/conopeptides.tsv"))) |>
  # ggplot(aes(log10(contig_impact_score), effective_length_frac, color = Region, alpha = log10(fill), size = 10)) + geom_point() + 
  ggplot(aes(log10(contig_impact_score), fill = Region)) + geom_density() +
  ggthemes::scale_fill_pander() + theme_bw(base_size = 16)
  group_by(Region) |> summarise(fill = sum(fill), n = n_distinct(protein_id)) 

# COUNT_DB |> 
#   filter(protein_id %in% c("Cluster-15813.116417.p1", "Cluster-15813.36708.p1", "Cluster-15813.37436.p1"))

recode_to <- structure(c("Control","Shrimp", "Mollusk", "Polychaete", "Mixed"))

recode_to <- structure(recode_to, names = c("Ctrl","Cam","Lit","Pol","Mix"))

diet_col <- c("gray20","#146179", "#09BC9F", "#FEB65F", "#C55E2D")

diet_col <- structure(diet_col, names = recode_to)


dir <- "C://Users//cinai/OneDrive/Documentos/PUBLICATION_DIR/06.Quantification/MATRIX_RSEM_dir/"

f <- "glmLRT_multiple_contrast_ctrl_and_treatments_kallisto.rds"

DEGS <- read_rds(file.path(dir, f)) |> filter(abs(logFC) > 2 & FDR < 0.05) 

DEGS |> filter(sam_group == "Shrimp") |> left_join(COUNT_DB) |> view()

# Cluster-15813.57809.p1

Split_overlaps <- DEGS |> 
  distinct(protein_id, sam_group, sampleX) |>
  mutate(Feeding = sapply(strsplit(sampleX, "_"), `[`, 1)) |>
  dplyr::mutate(Feeding = dplyr::recode_factor(Feeding, !!!recode_to)) |> 
  distinct(protein_id, Feeding)  |>
  group_by(protein_id) |>
  summarise(
    Diet_groups = paste(sort(unique(Feeding)), collapse = ","),
    n = n_distinct(Feeding),
    # logFC = mean(abs(logFC)),
    .groups = "drop"
  ) 

Split_overlaps <- Split_overlaps |>
  # select(-n) |>
  left_join(COUNT_DB) |>
  left_join(DB)


excel_sheets_v2 <- Split_overlaps %>%
  filter(n == 1) |> select(-n) |>
  group_by(Diet_groups) %>%
  nest() %>%
  deframe() %>%
  map(~.x %>% select(-Superfamily))

file_out <- paste0(pub_dir, "/Supplementary_database_1.xlsx")

writexl::write_xlsx(excel_sheets_v2, path = file_out)

###
###
###

excel_sheets_v2 <- Split_overlaps %>%
  filter(n != 1) |> select(-n) |>
  group_by(Diet_groups) %>%
  nest() %>%
  deframe() %>%
  map(~.x %>% select(-Superfamily))

file_out <- paste0(pub_dir, "/Supplementary_database_2.xlsx")

writexl::write_xlsx(excel_sheets_v2, path = file_out)

file_out <- paste0(pub_dir, "/Supplementary_database.csv")

write_csv(Split_overlaps, file = file_out)



conotoxin_specialization <- Split_overlaps |>
  filter(n == 1) |>
  select(-n) |>
  left_join(COUNT_DB) |>
  left_join(DB)


Split_overlaps |> count(n, Class) |> pivot_wider(names_from = n, values_from = nn)

Split_overlaps |>  filter(n == 1) |> 
  count(Diet_groups, Superfamily, Class) |> 
  pivot_wider(names_from = Diet_groups, values_from = n, values_fill = 0) |> 
  view()

qgenes <- unique(conotoxin_specialization$protein_id)

# DEGS |>
#   left_join(DB) |>
#   filter(protein_id %in% qgenes) |>
#   mutate(col = ifelse(protein_id %in% qgenes, "unique", "overlap")) |>
#   ggplot(aes(y = logFC, x = logCPM, alpha = -log10(FDR))) +
#   facet_grid(~ sam_group) +
#   # scale_color_manual(values = c("black","red")) +
#   geom_point(aes(color = Class))
#   
#   


DEGS |> 
  mutate(Feeding = sapply(strsplit(sampleX, "_"), `[`, 1)) |>
  dplyr::mutate(Feeding = dplyr::recode_factor(Feeding, !!!recode_to)) |> 
  # filter(protein_id %in% qgenes) |>
  select(FDR, logFC, Feeding) |>
  mutate(logFC = abs(logFC)) |>
  pivot_longer(-Feeding) |>
  ggplot(aes(value, fill = Feeding)) +
  facet_grid( ~ name, scales = "free") +
  geom_histogram(position = position_stack()) +
  scale_color_manual("",values = diet_col) +
  scale_fill_manual("",values = diet_col)  +
  labs(x = element_blank(), y = element_blank()) +
  my_custom_theme()


DEGS |> 
  left_join(DB) |>
  # filter(protein_id %in% qgenes) |>
  count(Class, sort = T)

