# read output from conoSorter
# Wrangling Regex files and summarise according to
# Number of orfs predicted (normalize to the total N of transcripts per transcriptome)
# Number of known vs novel conopeptide candidates
# Hydrophobicity distribution
# Count classes of superFamilies (score_sf)

# Reruning peptide mode (/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/californicus/05.Prediction/TransDecoder_dir/test_m50_dir > conoSorter of complete orfs predicted)

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

dir <- "/Users/cigom/Documents/GitHub/conopeptides/05.Prediction/ConoSorter_dir"

pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR"

# cDNA is sorted initially using ConoSorter [31], 
# which translates raw cDNA sequences into six reading frames and 
# extracts sequences from the first start codon in each read to the first subsequent stop codon. 
# The results generated two files, the Regex.tab file containing unambiguously identified amino acid sequences and 
# unclassified amino acid sequences considered to be novel peptides, respectively

pHHM_f <- list.files(dir, pattern = "_pHMM.tab", full.names = T) # _pHMM.tab and _Regex.tab

Regex_f <- list.files(dir, pattern = "_Regex.tab", full.names = T) # _pHMM.tab and _Regex.tab

read_regex <- function(f, Hydrophobicity_val = 60, pwidth_val = 50) {
  
  DF <- read_delim(f, delim = "|", col_names = T) %>% mutate(Method = basename(f))
  
  DF <- DF %>%
    dplyr::rename("Hydrophobicity"="% Hydrophobicity (Signal)", "transcript" = "Read Name") %>%
    mutate(Hydrophobicity = gsub("%", "", Hydrophobicity), Hydrophobicity = as.double(Hydrophobicity)) %>%
    dplyr::rename("Protein_width"="# A.A", "Cys_number" = "# Cysteine(s)") %>%
    dplyr::rename("Score_sf"="Score Superfamily", "Score_class" = "Score Class") %>%
    dplyr::rename("seq"="Protein Sequence")
  
  
  DF <- DF %>% mutate(Method = gsub("_Regex.tab", "", Method))
  
  
  # Filter step as Borghie et al.
  
  DF %>% 
    filter(Hydrophobicity > Hydrophobicity_val) %>%
    filter(Protein_width >= pwidth_val)
  
  
}

write_fasta <- function(f) {
  
  outName <- gsub("_Regex.tab|_pHMM.tab", "", basename(f))
  
  select_headers <- c( "transcript", "Conopeptide", "Hydrophobicity","Signal","Pro","Mature")
  
  Regex <- read_delim(f, delim = "|", col_names = T) %>%
    dplyr::rename("Hydrophobicity"="% Hydrophobicity (Signal)", "transcript" = "Read Name") %>%
    mutate(Hydrophobicity = gsub("%", "", Hydrophobicity), Hydrophobicity = as.double(Hydrophobicity)) %>%
    dplyr::rename("Score_sf"="Score Superfamily", "Score_class" = "Score Class") %>%
    dplyr::rename("pseq"="Protein Sequence", "Protein_width"="# A.A") %>%
    dplyr::rename("Signal"="Superfamily (Signal)", "Pro" = "Superfamily (Pro-region)",  "Mature" =  "Superfamily (Mature)") %>%
    # Any filter step? 
    # filter(Protein_width < 10) %>%
    select(starts_with(c(select_headers, "pseq"))) %>%
    mutate(Hydrophobicity = round(Hydrophobicity, digits = 0)) %>%
    # mutate(Score_sf = gsub("[^0-9.-]", "", Score_sf)) %>%
    distinct() %>%
    unite("header", transcript:Mature, sep = "|") %>%
    pull(pseq, name = header)
  
  # Search pHMM.tab file for 
  
  # f2 <-  file.path(dirname(f), paste0(outName, "_pHMM.tab"))
  
  f2 <- list.files(path = dirname(f), pattern = paste0(outName, "_pHMM.tab"), full.names = T)
  
  # if(!is.empty(f2)) 
  
  pHMM <- read_delim(f2, delim = "|", col_names = T) %>%
    dplyr::rename("Hydrophobicity"="% Hydrophobicity (Signal)", "transcript" = "Read Name") %>%
    mutate(Hydrophobicity = gsub("%", "", Hydrophobicity), Hydrophobicity = as.double(Hydrophobicity)) %>%
    dplyr::rename("pseq"="Protein Sequence", "Protein_width"="# A.A") %>%
    dplyr::rename("Signal"="Superfamily (Signal)", "Pro" = "Superfamily (Pro-region)",  "Mature" =  "Superfamily (Mature)") %>%
    select(starts_with(c(select_headers, "pseq"))) %>%
    distinct() %>%
    unite("header", transcript:Mature, sep = "|") %>%
    pull(pseq, name = header)
    
  
  OUT <- c(Regex, pHMM)
  
  # else
  # OUT <- Regex

  # Any sanity check?
  
  any(names(Regex) %in% names(pHMM)) # MUST BE FALSE
  
  
  # Out file
  
  # outName <- gsub("_Regex.tab|_pHMM.tab", "", basename(f))
  
  outFile <- file.path(dirname(f), paste0(outName, "_regex_pHMM.pep"))
  
  Biostrings::writeXStringSet(Biostrings::AAStringSet(OUT), file = outFile)
  
  
}

# lapply(Regex_f[2], write_fasta)

# read_regex(Regex_f[2])

DF <- lapply(Regex_f, read_regex)

DF <- do.call(rbind,DF)

DF %>% dplyr::count(Method)

recode_to <- c("Trinity",  
  "SuperDuper_Trinity",
  "spades", 
  "SuperDuper_Spades",
  "MMseqs",
  "MMseqs_hisat_SuperDuper",
  "Merged_hisat_SuperDuper",
  "Merged_polyA_hisat_SuperDuper",
  "Merged_polyA_hisat_SuperDuper.fasta.transdecoder")

recode_to <- structure(
  c("Trinity (T)", 
    "Trinity-Lace (Hisat)",
    "Spades (S)", 
    "Spades-Lace (Hisat)", 
    "MMseqs (T-S)",
    "MMseqs-Lace T-S (Hisat)",
    "Concat-Lace T-S (Hisat)",
    "Concat-Lace T-S (Hisat-polA)",
    "Concat-Lace T-S (Hisat-polA-transdecoder)"), 
  names = recode_to)

DF <- mutate(DF, Method = dplyr::recode_factor(Method, !!!recode_to))


# Number of orfs predicted (normalize to the total N of transcripts per transcriptome) -----

dir <- "/Users/cigom/Documents/GitHub/conopeptides/Nx_Metrics_dir/"

seqs_f <- list.files(dir, "fasta$", full.names = T)

contig_len <- function(f) {
  
  require(tidyverse)
  
  recode_to <- c("Trinity.fasta",  
    "trinity_hisat_superDuper.fasta",
    "spades.fasta", 
    "spades_hisat_superDuper.fasta",
    "MMseqs.fasta",
    "MMseqs_hisat_SuperDuper.fasta",
    "Merged_hisat_SuperDuper.fasta",
    "Merged_polyA_hisat_SuperDuper.fasta")
  
  recode_to <- structure(
    c("Trinity (T)", 
      "Trinity-Lace (Hisat)",
      "Spades (S)", 
      "Spades-Lace (Hisat)", 
      "MMseqs (T-S)",
      "MMseqs-Lace T-S (Hisat)",
      "Concat-Lace T-S (Hisat)",
      "Concat-Lace T-S (Hisat-polA)"), 
    names = recode_to)
  
  DNA <- Biostrings::readDNAStringSet(f)
  
  contig_len <- length(Biostrings::readDNAStringSet(f))
  
  out <- data.frame(contig_len, Method = basename(f))
  
  
  out <- mutate(out, Method = dplyr::recode_factor(Method, !!!recode_to))
  
  return(out)
}

LEN_DF <- lapply(seqs_f, contig_len)

LEN_DF <- do.call(rbind,LEN_DF)

# DF %>% 
#   count(Method) %>% left_join(LEN_DF) %>% 
#   mutate(contig_len_frac = n/contig_len) %>% 
#   arrange(contig_len_frac) %>%
#   distinct(Method) %>% pull() %>% as.character() -> method_levels

str(gene <- DF %>% pull(transcript))

str(gene <- gsub("_[0-9]+_[0-9]+$","", gene))


DF %>% 
  # if using gene level instead of orf level:
  mutate(gene = gene) %>% distinct(Method, gene) %>%
  count(Method) %>% left_join(LEN_DF) %>% 
  mutate(contig_len = ifelse(is.na(contig_len), 90497, contig_len)) %>%
  mutate(contig_len_frac = n/contig_len) %>% 
  arrange(contig_len_frac) %>%
  mutate(Method = factor(Method, levels = unique(Method))) %>%
  ggplot(aes(y = Method, x = contig_len_frac)) +
  geom_col(fill = "black") + 
  geom_text(aes(label= scales::comma(n)), hjust= 1.2, vjust = 0.5, size = 7, family = "GillSans", color = "white") +
  # labs(x =) +
  scale_x_continuous( "Fraction of Orfs (Transcript with Orf/Assembled transcripts)", labels = scales::percent_format()) +
  theme_bw(base_size = 16, base_family = "GillSans") 

# if only orfs w/ score > 0

DF %>% 
  # if using gene level instead of orf level:
  mutate(gene = gene) %>% 
  filter(Score_sf > 0) %>%
  distinct(Method, gene) %>%
  count(Method) %>% left_join(LEN_DF) %>% 
  mutate(contig_len = ifelse(is.na(contig_len), 90497, contig_len)) %>%
  mutate(contig_len_frac = n/contig_len) %>% 
  arrange(contig_len_frac) %>%
  mutate(Method = factor(Method, levels = unique(Method))) %>%
  ggplot(aes(y = Method, x = contig_len_frac)) +
  geom_col(fill = "black") + 
  geom_text(aes(label= scales::comma(n)), hjust= 1.2, vjust = 0.5, size = 7, family = "GillSans", color = "white") +
  # labs(x =) +
  scale_x_continuous( "Fraction of Orfs (Transcript with Orf/Assembled transcripts)", labels = scales::percent_format()) +
  theme_bw(base_size = 16, base_family = "GillSans") 

# DF %>% ggplot(aes(Hydrophobicity, color = Method)) + ggplot2::stat_ecdf()
# DF %>% ggplot(aes(Hydrophobicity)) + geom_histogram() + facet_grid(~ Method)

# DF %>% ggplot(aes(Hydrophobicity)) + geom_histogram()

# Number of known vs novel conopeptide candidates -----

# DF %>% 
#   mutate(Conopeptide = ifelse(Conopeptide == "-", "Novel", "Known")) %>%
#   filter(Conopeptide == "Novel") %>%
#   count(Method, Conopeptide) %>%
#   left_join(LEN_DF) %>%
#   mutate(contig_len_frac = n/contig_len) %>% 
#   arrange(contig_len_frac) %>%
#   mutate(Method = factor(Method, levels = unique(Method))) %>%
#   # filter(Conopeptide %in% "Known") %>%
#   ggplot(aes(y = Method, x = contig_len_frac, group = Conopeptide)) + 
#   geom_col() + 
#   geom_text(aes(label= scales::comma(n)), hjust= 1.2, vjust = 0.5, size = 7, family = "GillSans", color = "white") +
#   facet_grid(~ Conopeptide, scales = "free_x") +
#   theme_bw(base_size = 14, base_family = "GillSans")


# Count classes of superFamilies (score_sf) ----

# Score classification consist if sequence presents a single (score: 1), either (score: 2) or all (score: 3) of the conopeptides regions:
# 0. The ER signal  
# 1. the pro-peptide 
# 2. The mature 

DF %>% 
  select(contains(c("Method", "Score"))) %>%
  filter(grepl("CONFLICT", Score_sf)) %>% count(Method, Score_sf)

# Omit \(!!CONFLICT!!)

DF <- DF %>% 
  mutate(Score_sf = gsub("[^0-9.-]", "", Score_sf)) %>%
  mutate(Score_class = gsub("[^0-9.-]", "", Score_class))

# in which method score classificatio increased to 3?

# DF %>% 
#   select(contains(c("Method", "Score"))) %>%
#   # mutate_if(is.character, as.double) %>%
#   mutate_at(vars(matches(c("Score"))), function(x) as.double(x)) %>%
#   pivot_longer(cols = c("Score_sf", "Score_class"), names_to = "Type", values_to = "Score") %>%
#   count(Method, Type, Score) %>%
#   drop_na(Score) %>%
#   # Must normalize to the total number of transcripts 
#   group_by(Method, Type) %>% mutate(frac = n/sum(n)) %>%
#   filter(Score == 3) %>%
#   ggplot(aes(y = Method, x = frac, fill = as.character(Score))) + 
#   geom_col() + facet_grid(~ Type, scales = "free_x") +
#   theme_bw(base_size = 14, base_family = "GillSans")

# scores per superfamilies -----


DF %>% 
  filter(Score_sf == 3) %>%
  select(contains(c("Method", "Score_sf","Superfamily ("))) %>%
  count(Method, Score_sf) %>%
  group_by(Method) %>% 
  ungroup() %>%
  arrange(n) %>%
  distinct(Method) %>% pull() %>% as.character() -> method_levels


DF %>% 
  # filter(Score_sf > 0) %>%
  select(contains(c("Method", "Score_sf","Superfamily ("))) %>%
  count(Method, Score_sf) %>%
  group_by(Method) %>% 
  mutate(frac = n/sum(n)) %>%
  mutate(Method = factor(Method, levels = method_levels)) %>%
  ggplot(aes(y = Method, x = n, fill = Score_sf)) + 
  geom_col() +
  geom_text(aes(label= scales::comma(n)), hjust= 1.2, vjust = 0.5, size = 7, family = "GillSans", color = "white") +
  facet_grid(~ Score_sf, scales = "free_x") +
  scale_x_continuous( "Fraction (Superfamily class / )", labels = scales::comma_format()) +
  theme_bw(base_size = 14, base_family = "GillSans")


# Which types of superFamilies per methods?

which_cols <- c("Superfamily (Signal)", "Superfamily (Pro-region)", "Superfamily (Mature)")

paste_col <- function(x) { 
  x <- x[!is.na(x)] 
  x <- unique(sort(x))
  x <- paste(x, sep = '-', collapse = '-')
}


# DF %>%
#   select(contains(c("transcript","Method", "Score_sf","Superfamily ("))) %>%
#   filter(Score_sf > 0) %>%
#   # unite("Region", all_of(which_cols), sep = "-")
#   mutate(row_number = row_number()) %>%
#   pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
#   filter(Superfamily != "-") %>%
#   mutate(Region = gsub("Superfamily ", "", Region)) %>%
#   filter(Method == "Concat-Lace T-S (Hisat)" & Region == "(Mature)") %>%
#   distinct(transcript)
#   distinct(Superfamily)

DF %>%
  select(contains(c("transcript","Method", "Score_sf","Superfamily ("))) %>%
  filter(Score_sf > 0) %>%
  # unite("Region", all_of(which_cols), sep = "-")
  mutate(row_number = row_number()) %>%
  pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
  filter(Superfamily != "-") %>%
  mutate(Region = gsub("Superfamily ", "", Region)) %>%
  filter(Method == "Concat-Lace T-S (Hisat)" & Score_sf == 3) %>%
  distinct(Region, Superfamily)


superfm_df <- DF %>% 
  select(contains(c("transcript","Method", "Score_sf","Superfamily ("))) %>%
  filter(Score_sf > 0) %>% 
  pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
  filter(Superfamily != "-") %>%
  mutate(Region = gsub("Superfamily ", "", Region)) %>%
  # if count number of gene (transcrits ids in the assemblies) instead of orfs candidates:
  # mutate(transcript =  gsub("_[0-9]+_[0-9]+$","", transcript)) %>% distinct() %>%
  group_by(Method, transcript) %>%
  summarise(across(Region,.fns = paste_col), Score_sf = n()) 

superfm_df %>%
  filter(Method == "Trinity-Lace (Hisat)" & Region == "(Mature)-(Pro-region)-(Signal)")

superfm_viz <- superfm_df %>%
  count(Method, Region) %>% ungroup() %>% dplyr::rename("n_transcripts" = "n")

# Include number of unique genes or superfamilies
# gsub("_[0-9]+_[0-9]+$","", gene))

superfm_df <- DF %>% 
  select(contains(c("transcript","Method", "Score_sf","Superfamily ("))) %>%
  filter(Score_sf > 0) %>% 
  pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
  filter(Superfamily != "-") %>%
  mutate(Region = gsub("Superfamily ", "", Region)) %>%
  # if count number of gene (transcrits ids in the assemblies) instead of orfs candidates:
  # mutate(transcript =  gsub("_[0-9]+_[0-9]+$","", transcript)) %>% distinct() %>%
  group_by(Method, transcript) %>%
  summarise(across(Superfamily, .fns = paste_col), Score_sf = n()) %>% 
  left_join(superfm_df) 



superfm_df %>% ungroup() %>% count(Region)

superfm_df %>%
  filter(Region != "(Mature)-(Signal)") %>%
  ungroup() %>%
  count(Method, Superfamily, Score_sf) %>%
  ggplot(aes(y = Superfamily, x = Method, fill = n)) +
  facet_grid(Score_sf ~., scales = "free_y", space = "free_y") +
  geom_tile(color = "white",
    lwd = 0.5,
    linetype = 1) +
  theme_classic(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top",
    strip.background = element_rect(fill = 'gray68', color = 'white'),
    strip.text = element_text(color = "black",hjust = 0, size = 7),
    axis.line.y = element_line(color = 'white'),
    axis.text.y = element_text(hjust = 1),
    axis.ticks.length = unit(5, "pt")) +
  guides(
    fill = guide_colorbar(barwidth = unit(3, "in"),
      barheight = unit(0.1, "in"), label.position = "bottom",
      alignd = 0.5,
      title = "Number of peptides",
      title.position  = "top",
      title.theme = element_text(size = 10, family = "GillSans", hjust = 1),
      ticks.colour = "black", ticks.linewidth = 0.35,
      frame.colour = "black", frame.linewidth = 0.35,
      label.theme = element_text(size = 10, family = "GillSans"))) +
  theme(axis.ticks.x = element_blank(), 
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
    # axis.text.x = element_blank(), 
    axis.line.x = element_blank())
  # geom_point(aes(size = n), shape = 21)
  

superfm_viz <- superfm_df %>%
  ungroup() %>% distinct(Method, Score_sf, Superfamily, Region) %>%
  count(Method, Region) %>% ungroup() %>% dplyr::rename("n_sf" = "n") %>%
  left_join(superfm_viz)
  

reg_lev <- c("(Mature)-(Pro-region)-(Signal)","(Mature)-(Pro-region)", "(Pro-region)-(Signal)","(Mature)-(Signal)", "(Mature)","(Pro-region)","(Signal)")

superfm_df %>% filter(Method == "Trinity-Lace (Hisat)" & Region == "(Mature)") %>% distinct(Superfamily)

superfm_viz %>%
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



# DF %>% count(Method, `Superfamily (Signal)`)

DF %>% 
  select(contains(c("Method", "Class ("))) %>%
  pivot_longer(-Method, names_to = "Class", values_to = "value") %>%
  filter(value != "-") %>%
  count(Method, Class) %>%
  group_by(Method) %>% mutate(frac = n/sum(n)) %>%
  ggplot(aes(y = Class, x = frac, fill = Method)) + 
  geom_col(position = position_dodge2(reverse = T)) +
  theme_bw(base_size = 14, base_family = "GillSans")



# pHMMM *Novel candidates -----

read_pHMM <- function(f,  Hydrophobicity_val = 60, pwidth_val = 50, eval = 0.05) {
  
  DF <- read_delim(f, delim = "|", col_names = T) %>% mutate(Method = basename(f))
  
  DF <- DF %>%
    dplyr::rename("Hydrophobicity"="% Hydrophobicity (Signal)", "transcript" = "Read Name") %>%
    mutate(Hydrophobicity = gsub("%", "", Hydrophobicity), Hydrophobicity = as.double(Hydrophobicity)) %>%
    dplyr::rename("Protein_width"="# A.A", "Cys_number" = "# Cysteine(s)") %>%
    dplyr::rename("E_value"="E-value Superfamily")
  
  # protein_id <- DF %>% pull(ID)
  # protein_id <- gsub("_[0-9]+_[0-9]+$","", protein_id)
  
  DF <- DF %>% mutate(Method = gsub("_pHMM.tab", "", Method)) 
  
  DF <- DF %>% 
    filter(Hydrophobicity > Hydrophobicity_val) %>%
    filter(Protein_width >= pwidth_val) %>%
    filter(as.numeric(E_value) < eval)

  
  DF1 <- DF %>% select(contains(c("transcript","Method","Superfamily ("))) %>%
    # filter(Score_sf > 0) %>% # Score_sf > 0 == Superfamily != "-"
    pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
    filter(Superfamily != "-") %>%
    mutate(Region = gsub("Superfamily ", "", Region)) %>%
    group_by(Method, transcript) %>%
    summarise(across(Region,.fns = paste_col), Score_sf = n()) %>% ungroup()
  
  
  OUT <- DF %>% 
    select(contains(c("transcript","Method", "Score_sf","Superfamily ("))) %>%
    pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
    filter(Superfamily != "-") %>%
    mutate(Region = gsub("Superfamily ", "", Region)) %>%
    group_by(Method, transcript) %>%
    summarise(across(Superfamily, .fns = paste_col), Score_sf = n()) %>% 
    left_join(DF1) %>%
    mutate(tab = "pHMM")
  
  return(OUT)
  
}

read_pHMM <- function(f) {
  
  DF <- read_delim(f, delim = "|", col_names = T) %>% mutate(Method = basename(f))
  
  DF <- DF %>%
    dplyr::rename("Hydrophobicity"="% Hydrophobicity (Signal)", "transcript" = "Read Name") %>%
    mutate(Hydrophobicity = gsub("%", "", Hydrophobicity), Hydrophobicity = as.double(Hydrophobicity)) %>%
    dplyr::rename("Protein_width"="# A.A", "Cys_number" = "# Cysteine(s)") # %>%
    # dplyr::rename("Score_sf"="Score Superfamily", "Score_class" = "Score Class") 
  
  
  DF %>% mutate(Method = gsub("_pHMM.tab", "", Method))

  
  
}

pHHMdf <- lapply(pHHM_f, read_pHMM)

pHHMdf <- do.call(rbind,pHHMdf)

pHHMdf <- mutate(pHHMdf, Method = dplyr::recode_factor(Method, !!!recode_to))

pHHMdf %>% dplyr::count(Method)

# pHHMdf <- pHHMdf %>%  filter(Hydrophobicity > 60) %>% filter(Protein_width >= 50) 
# 
# pHHMdf %>%
#   dplyr::count(Method) 

# pHHMdf_ <- pHHMdf %>% 
#   select(contains(c("transcript","Method", "Superfamily ("))) %>%
#   pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
#   filter(Superfamily != "-") %>%
#   mutate(Region = gsub("Superfamily ", "", Region)) %>%
#   # if count number of gene (transcrits ids in the assemblies) instead of orfs candidates:
#   # mutate(transcript =  gsub("_[0-9]+_[0-9]+$","", transcript)) %>% distinct() %>%
#   group_by(Method, transcript) %>%
#   summarise(across(Region,.fns = paste_col), Score_sf = n()) 

# pHHMdf_ <- pHHMdf %>% 
#   select(contains(c("transcript","Method","Superfamily ("))) %>%
#   pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
#   filter(Superfamily != "-") %>%
#   mutate(Region = gsub("Superfamily ", "", Region)) %>%
#   # if count number of gene (transcrits ids in the assemblies) instead of orfs candidates:
#   # mutate(transcript =  gsub("_[0-9]+_[0-9]+$","", transcript)) %>% distinct() %>%
#   group_by(Method, transcript) %>%
#   summarise(across(Superfamily, .fns = paste_col), Score_sf = n()) %>% 
#   left_join(pHHMdf_) %>%
#   select(names(superfm_df)) %>%
#   mutate(tab = "pHHM")


# write outp -----
superfm_df %>% mutate(tab = "Regex") %>% rbind(pHHMdf) %>%
  mutate(gene =  gsub("_[0-9]+_[0-9]+$","", transcript)) %>%
  write_rds(file = paste0(pub_dir, "/superfm_df.rds"))

# AMINOACID LEN

DF %>% select(transcript, Protein_width, Method) %>% mutate(tab = "Regex") %>%
  rbind(pHHMdf %>% select(transcript, Protein_width, Method) %>% mutate(tab = "pHHM")) %>%
  distinct() %>%
  # ggplot(aes(Protein_width)) + 
  # stat_ecdf(aes(color = tab)) + 
  # geom_violin() +
  ggplot(aes(y = Method, x = Protein_width,fill = after_stat(x))) +
  facet_grid(~ tab, scales = "free_x") +
  ggridges::geom_density_ridges_gradient(
    jittered_points = T,
    position = ggridges::position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 3, point_alpha = 1, alpha = 1) +
  scale_fill_viridis_c(option = "C") +
  labs(y = "", x = "Blast bit score") +
  theme_bw(base_family = "GillSans", base_size = 14) + theme(legend.position = "none")
    

DF %>% select(transcript, Protein_width, Method) %>% mutate(tab = "Regex") %>%
  rbind(pHHMdf %>% select(transcript, Protein_width, Method) %>% mutate(tab = "pHHM")) %>%
  distinct() %>%
  ggplot(aes(y = Method, x = Protein_width)) +
  geom_violin() +
  facet_grid(~ tab, scales = "free_x") 

# count unique genes


# Number of transcripts spread in orfs?
str(gene <- DF %>% pull(transcript))
str(gene <- gsub("_[0-9]+_[0-9]+$","", gene))
# str(gene <- sapply(strsplit(gene, "_"), `[`, 1))

DF %>% 
  select(Method, transcript) %>% 
  mutate(gene= gene) %>% 
  # filter(Method %in% "Concat-Lace T-S (Hisat-polA-transdecoder)")
  group_by(Method) %>% summarise(n_orfs = n(), n_genes = length(unique(gene))) %>%
  arrange(desc(n_genes)) %>%
  left_join(LEN_DF) 




# Number of transcripts spread in orfs?
str(gene <- pHHMdf %>% pull(transcript))
str(gene <- gsub("_[0-9]+_[0-9]+$","", gene))

pHHMdf %>% 
  select(Method, transcript) %>% 
  mutate(gene= gene) %>% 
  # filter(Method %in% "Concat-Lace T-S (Hisat-polA-transdecoder)")
  group_by(Method) %>% summarise(n_orfs = n(), n_genes = length(unique(gene))) %>%
  arrange(desc(n_genes))

DF %>% 
  select(transcript, Protein_width, Method) %>% mutate(tab = "Regex") %>%
  rbind(pHHMdf %>% select(transcript, Protein_width, Method) %>% mutate(tab = "pHHM")) 
