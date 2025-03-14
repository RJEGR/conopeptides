# read output from conoSorter
# Wrangling Regex files and summarise according to
# Number of orfs predicted (normalize to the total N of transcripts per transcriptome)
# Number of known vs novel conopeptide candidates
# Hydrophobicity distribution
# Count classes of superFamilies (score_sf)

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

dir <- "/Users/cigom/Documents/GitHub/conopeptides/05.Prediction/ConoSorter_dir"

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

lapply(Regex_f, write_fasta)

# read_regex(Regex_f[3])

DF <- lapply(Regex_f, read_regex)

DF <- do.call(rbind,DF)

DF %>% dplyr::count(Method)

recode_to <- c("Trinity",  
  "SuperDuper_Trinity",
  "spades", 
  "SuperDuper_Spades",
  "MMseqs",
  "MMseqs_hisat_SuperDuper",
  "Merged_hisat_SuperDuper")

recode_to <- structure(
  c("Trinity (T)", 
    "Trinity-Lace (Hisat)",
    "Spades (S)", 
    "Spades-Lace (Hisat)", 
    "MMseqs (T-S)",
    "MMseqs-Lace T-S (Hisat)",
    "Concat-Lace T-S (Hisat)"), 
  names = recode_to)

DF <- mutate(DF, Method = dplyr::recode_factor(Method, !!!recode_to))


# Number of orfs predicted (normalize to the total N of transcripts per transcriptome)

dir <- "/Users/cigom/Documents/GitHub/conopeptides/Nx_Metrics_dir/"

seqs_f <- list.files(dir, "fasta", full.names = T)

contig_len <- function(f) {
  
  require(tidyverse)
  
  recode_to <- c("Trinity.fasta",  
    "trinity_hisat_superDuper.fasta",
    "spades.fasta", 
    "spades_hisat_superDuper.fasta",
    "MMseqs.fasta",
    "MMseqs_hisat_SuperDuper.fasta",
    "Merged_hisat_SuperDuper.fasta")
  
  recode_to <- structure(
    c("Trinity (T)", 
      "Trinity-Lace (Hisat)",
      "Spades (S)", 
      "Spades-Lace (Hisat)", 
      "MMseqs (T-S)",
      "MMseqs-Lace T-S (Hisat)",
      "Concat-Lace T-S (Hisat)"), 
    names = recode_to)
  
  DNA <- Biostrings::readDNAStringSet(f)
  
  contig_len <- length(Biostrings::readDNAStringSet(f))
  
  out <- data.frame(contig_len, Method = basename(f))
  
  
  out <- mutate(out, Method = dplyr::recode_factor(Method, !!!recode_to))
  
  return(out)
}

LEN_DF <- lapply(seqs_f, contig_len)

LEN_DF <- do.call(rbind,LEN_DF)

DF %>% 
  count(Method) %>% left_join(LEN_DF) %>% 
  mutate(contig_len_frac = n/contig_len) %>% 
  arrange(contig_len_frac) %>%
  distinct(Method) %>% pull() %>% as.character() -> method_levels

DF %>% 
  count(Method) %>% left_join(LEN_DF) %>% 
  mutate(contig_len_frac = n/contig_len) %>% 
  arrange(contig_len_frac) %>%
  mutate(Method = factor(Method, levels = unique(Method))) %>%
  ggplot(aes(y = Method, x = contig_len_frac)) +
  geom_col(fill = "black") + 
  geom_text(aes(label= scales::comma(n)), hjust= 1.2, vjust = 0.5, size = 7, family = "GillSans", color = "white") +
  # labs(x =) +
  scale_x_continuous( "Fraction of Orfs (Transcript with Orf/Assembled transcripts)", labels = scales::percent_format()) +
  theme_bw(base_size = 16, base_family = "GillSans") 


DF %>% ggplot(aes(Hydrophobicity, color = Method)) + ggplot2::stat_ecdf()
DF %>% ggplot(aes(Hydrophobicity)) + geom_histogram() + facet_grid(~ Method)

# DF %>% ggplot(aes(Hydrophobicity)) + geom_histogram()


DF %>% distinct(Conopeptide)
DF %>% distinct(transcript)

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
  filter(Score_sf > 0) %>%
  select(contains(c("Method", "Score_sf","Superfamily ("))) %>%
  count(Method, Score_sf) %>%
  group_by(Method) %>% 
  mutate(frac = n/sum(n)) %>%
  mutate(Method = factor(Method, levels = method_levels)) %>%
  ggplot(aes(y = Method, x = frac, fill = Score_sf)) + 
  geom_col() +
  geom_text(aes(label= scales::comma(n)), hjust= 1.2, vjust = 0.5, size = 7, family = "GillSans", color = "white") +
  facet_grid(~ Score_sf, scales = "free_x") +
  scale_x_continuous( "Fraction (Superfamily class / )", labels = scales::percent_format()) +
  theme_bw(base_size = 14, base_family = "GillSans")


# Which types of superFamilies per methods?

which_cols <- c("Superfamily (Signal)", "Superfamily (Pro-region)", "Superfamily (Mature)")

paste_col <- function(x) { 
  x <- x[!is.na(x)] 
  x <- unique(sort(x))
  x <- paste(x, sep = '-', collapse = '-')
}


# Are unique nunbers of transcriots? Yes
DF %>%
  select(contains(c("transcript","Method", "Score_sf","Superfamily ("))) %>%
  filter(Score_sf > 0) %>%
  # unite("Region", all_of(which_cols), sep = "-")
  mutate(row_number = row_number()) %>%
  pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
  filter(Superfamily != "-") %>%
  mutate(Region = gsub("Superfamily ", "", Region)) %>%
  filter(Method == "Concat-Lace T-S (Hisat)" & Score_sf == 3) %>%
  distinct(transcript)

DF %>%
  select(contains(c("transcript","Method", "Score_sf","Superfamily ("))) %>%
  filter(Score_sf > 0) %>%
  # unite("Region", all_of(which_cols), sep = "-")
  mutate(row_number = row_number()) %>%
  pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
  filter(Superfamily != "-") %>%
  mutate(Region = gsub("Superfamily ", "", Region)) %>%
  filter(Method == "Concat-Lace T-S (Hisat)" & Region == "(Mature)") %>%
  distinct(transcript)
  distinct(Superfamily)

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
  # unite("Region", all_of(which_cols), sep = "-")
  mutate(row_number = row_number()) %>%
  pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
  filter(Superfamily != "-") %>%
  mutate(Region = gsub("Superfamily ", "", Region)) %>%
  group_by(Method, row_number) %>%
  summarise(across(Region, .fns = paste_col), n = n()) %>%
  count(Method, Region) %>% ungroup() %>% dplyr::rename("n_transcripts" = "n")

superfm_df <- DF %>% 
  select(contains(c("Method", "Score_sf","Superfamily ("))) %>%
  filter(Score_sf > 0) %>% 
  distinct() %>%
  mutate(row_number = row_number()) %>%
  pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
  filter(Superfamily != "-") %>%
  mutate(Region = gsub("Superfamily ", "", Region)) %>%
  group_by(Method, row_number) %>%
  summarise(across(Region, .fns = paste_col), n = n()) %>%
  count(Method, Region) %>% ungroup() %>% dplyr::rename("n_sf" = "n") %>%
  left_join(superfm_df)

reg_lev <- c("(Mature)-(Pro-region)-(Signal)","(Mature)-(Pro-region)", "(Pro-region)-(Signal)","(Mature)-(Signal)", "(Mature)","(Pro-region)","(Signal)")

superfm_df %>%
  mutate(label = paste0(scales::comma(n_transcripts), " (", scales::comma(n_sf),")")) %>%
  group_by(Method) %>% mutate(frac = n_transcripts/sum(n_transcripts)) %>%
  mutate(Method = factor(Method, levels = method_levels)) %>%
  mutate(Region = factor(Region, levels = reg_lev)) %>%
  ggplot(aes(y = Method, x = frac)) + 
  geom_col(position = position_dodge2(reverse = T),fill = "black") +
  facet_grid(~ Region, scales = "free_x") +
  scale_x_continuous( "Fraction (Transcript annotated/Assembled transcripts)", labels = scales::percent_format()) +
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

# pHHMdf %>% dplyr::count(Method)

pHHMdf %>%
  filter(Hydrophobicity > 60) %>%
  filter(Protein_width >= 50) %>%
  dplyr::count(Method) %>% view()



pHHMdf <-read_pHMM(pHHM_f[1])

Regexdf <-read_regex(Regex_f[1], Hydrophobicity_val = 0, pwidth_val = 0)

pHHMdf %>% distinct(transcript) %>% left_join(Regexdf)

Regexdf %>% distinct(transcript)

Regexdf %>%
  ggplot(aes(Protein_width, color = Method)) + ggplot2::stat_ecdf()

f <- Regex_f[1]

read_delim(f, delim = "|", col_names = T)

