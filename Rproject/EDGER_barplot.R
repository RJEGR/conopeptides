
# Step1 Select up putative conotoxins in time (2 and 4) vs control
# Step2 Select up putative conotoxin dfferent between time 2 and 4 per diet
# Step3 Contrast temporalidad and diet specificity using venn diagram

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)


library(tidyverse)

pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR/"

# LOAD data -----

DB <- read_tsv(paste0(pub_dir, "/conopeptides.tsv"))


CONOPEPDB <- DB %>%
  filter(Signalp_class == "SP") %>%
  mutate(uniprotkb_toxprot = sapply(strsplit(uniprotkb_toxprot, "[|]"), `[`, 2)) %>%
  mutate(uniprotkb_toxprot = ifelse(is.na(uniprotkb_toxprot), "uID", uniprotkb_toxprot)) %>%
  mutate(conoserver_protein = ifelse(is.na(conoserver_protein), "uID", conoserver_protein)) %>%
  mutate(conoserver_protein = ifelse(is.na(conoserver_protein), "uID", conoserver_protein)) %>%
  mutate(hmm_pred_conodictor = stringr::str_to_sentence(hmm_pred_conodictor)) %>%
  mutate(hmm_pred_conodictor = ifelse(is.na(hmm_pred_conodictor), "uID", hmm_pred_conodictor)) %>%
  mutate(Superfamily = ifelse(is.na(Superfamily), "uID", Superfamily)) %>%
  mutate(prediction_tool = ifelse(tab %in% "pHMM", paste0(prediction_tool,"_",tab), prediction_tool)) %>%
  select(protein_id, pep_seq, cluster, Signalp_class, prediction_tool, hmm_pred_conodictor, Superfamily, uniprotkb_toxprot, conoserver_protein) %>%
  unite("y_axis", cluster:conoserver_protein, sep = "|") 

dir <- "/Users/cigom/Documents/GitHub/conopeptides/06.Quantification/MATRIX_RSEM_dir/"

subdir <- "KALLISTO_Merged_polyA_hisat_SuperDuper.fasta.transdecoder_DIR/"

f <- "cds_exactTest_multiple_contrast_ctrl_and_treatments.rds"

# f <- "cds_exactTest_multiple_contrast_ctrl_and_treatments_cds_level.rds"

f <- list.files(file.path(dir, subdir), f, full.names = T)

RES <- read_rds(f)

# As the interest is to known which conotoxin over-expressed due to diet, omit Ctrl group in sampleX column
# In addition  filter significat DEGS

DataViz <- RES %>% filter(sampleX != "Ctrl") %>% filter(abs(logFC) > 2 & FDR < 0.05) %>% right_join(CONOPEPDB)

# As duplicated CDS have same expression patters and so on, same gene-pairwise values in DE analysis, just deduplicate redundancy omiting protein_id but using pep_seq column

# deduplicate to cds_level/pep_level

nrow(DataViz)

nrow(DataViz <- DataViz %>% select(-protein_id) %>% distinct())

# Step1: Contrasting results against control (ie. omit contrast dietA_time1 vs dietA_time2)

# DataViz <- DataViz %>% filter(if_any(where(is.character), ~ grepl(pattern = 'Ctrl', x = .x, ignore.case = T)))

cols_to_check <- c("sampleA", "sampleB")

DataViz <- DataViz %>%
  filter(
    if_any(all_of(cols_to_check), ~ str_detect(.x, "Ctrl"))) 

# DataViz <- DataViz %>% mutate(sampleA = ifelse(sampleA == "Ctrl", sampleB, sampleA))

DataViz %>% dplyr::count(sam_group, sampleA, sampleB, sampleX)

DataVizTop <- 
  DataViz %>% 
  #dplyr::mutate(sampleB = dplyr::recode_factor(sampleB, !!!recode_to)) %>%
  group_by(sam_group, sampleA, sampleB, sampleX) %>%
  arrange(logFC, .by_group = T) %>%
  group_by(sampleB, y_axis) %>%
  filter(abs(logFC) == max(abs(logFC))) %>%
  group_by(sampleX) %>%
  slice_head(n = 10)


DataVizTop %>%
  dplyr::count(sam_group, sampleA, sampleB, sampleX)

DataVizTop %>% ungroup() %>% distinct(y_axis)

# Split indistinc from unique, for labeling color

q <- DataVizTop %>% distinct(y_axis) %>% ungroup() %>% dplyr::count(y_axis, sort = T) %>% filter(n == 1) %>% pull(y_axis)

DataVizTop <- DataVizTop %>% mutate(col = ifelse(y_axis %in% q, "grey70", "grey89"))

col_vals <- c("grey70", "grey89")

col_vals <- structure(col_vals, names = col_vals)

DataVizTop %>% ungroup() %>% dplyr::count(col, sort = T) 

"log2fold-change (vst)"

x_label <- expression(log[2]~"fold-change")


gene_names <- DataVizTop %>% pull(y_axis, name = protein_id)

DataVizTop %>% 
  mutate(facet = paste0("Ctrl", "_vs_", sampleX)) %>%
  mutate(Label = y_axis, row_number = row_number(Label)) %>%
  #mutate(Label = paste0(Label, " (", protein_name, ")")) %>%
  mutate(row_number = paste(row_number, sampleX, sep = ":")) %>%
  mutate(Label = factor(paste(Label, row_number, sep = "__"),
    levels = rev(paste(Label, row_number, sep = "__")))) %>%
  mutate(star = ifelse(FDR <.001, "***", 
    ifelse(FDR <.01, "**",
      ifelse(FDR <.05, "*", "")))) %>%
  mutate(logFC = abs(logFC)) %>%
  mutate(
    # ymin = (abs(logFC) - lfcSE) * sign(logFC),
    # ymax = (abs(logFC) + lfcSE) * sign(logFC),
    y_star = logFC + (0.15)* sign(logFC)) %>% 
  ggplot(aes(y = Label, x = logFC)) + 
  # facet_grid(facet ~. , scales = "free_y", space = "free")+
  ggforce::facet_col(facet ~. , scales = "free_y", space = "free") +
  scale_y_discrete(labels = function(x) gsub("__.+$", "", x)) +
  scale_fill_manual(values = col_vals) +
  geom_col(aes(fill = col), width = 0.5, size = 0.25,
    position = position_stack(reverse = T), color = "white") +
  # geom_errorbar(aes(xmin = ymin, xmax = ymax), width = 0.05,
  # position = position_identity(), color = "black") +
  geom_text(aes(x = y_star, label=star), 
    # hjust = .7 ,
    # vjust=  .7,  
    color="black", position = position_identity(), 
    family = "GillSans", size = 1.5)+
  labs(y = NULL, x = x_label, title = "Up-expressed conotoxins under different ") +
  theme_bw(base_family = "GillSans", base_size = 10) +
  theme(
    legend.position = "none",
    # panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    # axis.text.y.right = element_text(angle = 0, hjust = 1, vjust = 0, size = 2.5),
    axis.text.y = element_text(angle = 0, size = 3),
    axis.text.x = element_text(angle = 0),
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 1, size = 7)) -> P


P


ggsave(P, filename = 'EDGERLOG2FCTOP10.png', 
  path = pub_dir, width = 5, height = 7, device = png, dpi = 600)


# 

myXStringSet <- DataVizTop %>% 
  ungroup() %>% filter(sam_group == "Polychaete") %>% 
  distinct(y_axis, pep_seq) %>% 
  # mutate(pep_seq = gsub("[*]$","", pep_seq)) %>%
  pull(pep_seq, name = y_axis)


myXStringSet <- Biostrings::AAStringSet(c(myXStringSet))


library(msa)

align <- msa::msa(myXStringSet, method = "ClustalW", order = "input")

.align <- msa::msaConvert(align)$seq

names(.align) <- msa::msaConvert(align)$nam

# data(BLOSUM62)
# msaConservationScore(align, BLOSUM62)

library(ggsci)

pat <- c("-", alphabet(myXStringSet, baseOnly=TRUE))

# colors <- structure(pal_aaas()(length(pat)), names = rev(pat))

# scales::show_col(colors)

DECIPHER::BrowseSeqs(AAStringSet(.align), colWidth = 120)


# Ven diagram or upset -----

# Q: Those upexpressed are unique from time or expression increase by time?
# 

# Filter

# If NA in DataViz, is because not preserved in CONOPEPDB

DataViz <- DataViz %>% drop_na(sam_group)

DataViz %>% 
  dplyr::count(sam_group, sampleA, sampleB, sampleX)


DataViz %>% 
  select(pep_seq, sam_group, sampleX) %>%
  mutate(sampleX = ifelse(grepl("_2$", sampleX), "2 months", "4 months")) %>%
  # filter(grepl("_2$", sampleX)) %>%
  distinct() %>%
  dplyr::count(sam_group, sampleX)

# Reformat
# 
UPSETDFA <- DataViz %>% 
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


library(ggupset)

my_font <- "GillSans"

P <- UPSETDFA %>%
  # filter(n > 1) %>%
  ggplot(aes(x = summarise_col)) +
  facet_grid(sampleX ~ ., scales = "free_y") +
  geom_bar(fill = "black") +
  scale_y_reverse("Number of conotoxins") +
  geom_text(stat='count', aes(label = after_stat(count)), 
    position = position_dodge(width = 1), vjust = 1.2, family = "GillSans", size = 2.5) +
  scale_x_upset(order_by = "degree", reverse = F, position = "top") +
  labs(x = "Degree of intersections") +
  theme_bw() +
  theme_combmatrix(
    combmatrix.panel.point.color.fill = "black",
    combmatrix.label.make_space = F,
    # combmatrix.panel.point.color.fill = panel.point.color.fill,
    combmatrix.panel.point.size = 0.25,
    combmatrix.panel.line.size = 0.25,
    base_family = "GillSans", base_size = 14,
    strip.background = element_rect(fill = 'gray90', color = 'white'),
    strip.text = element_text(color = "black", size = 12, family = "GillSans",hjust = 0.5),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_text(family = my_font),
    axis.title.y = element_text(family = my_font),
    axis.text.x = element_text(family = my_font),
    axis.text.y = element_text(family = my_font)) +
  axis_combmatrix(levels = recode_to) 

ggsave(P, filename = 'Intersections.png', 
  path = pub_dir, width = 5, height = 5, device = png, dpi = 600)

