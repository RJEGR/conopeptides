# Read conosorter peptide search
# 


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

dir <- "/Users/cigom/Documents/GitHub/conopeptides/05.Prediction/ConoSorter_dir/transdecoder.predict.conosorter_dir/"

# dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR/complete_and_partial_cds_dir/Conosorter_dir/" <- omit

pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR"

pattern <- "Merged_polyA_hisat_SuperDuper.fasta.transdecoder_pHMM.tab"

pHHM_f <- list.files(dir, pattern = pattern, full.names = T) # _pHMM.tab and _Regex.tab

pattern <- "Merged_polyA_hisat_SuperDuper.fasta.transdecoder_Regex.tab"

Regex_f <- list.files(dir, pattern = pattern, full.names = T) # _pHMM.tab and _Regex.tab


which_cols <- c("Superfamily (Signal)", "Superfamily (Pro-region)", "Superfamily (Mature)")

paste_col <- function(x) { 
  x <- x[!is.na(x)] 
  x <- unique(sort(x))
  x <- paste(x, sep = '_', collapse = '_')
}

read_regex <- function(f, Hydrophobicity_val = 60, pwidth_val = 50) {
  
  DF <- read_delim(f, delim = "|", col_names = T) %>% mutate(Method = basename(f))
  
  DF <- DF %>%
    dplyr::rename("Hydrophobicity"="% Hydrophobicity (Signal)", "ID" = "Read Name") %>%
    mutate(Hydrophobicity = gsub("%", "", Hydrophobicity), Hydrophobicity = as.double(Hydrophobicity)) %>%
    dplyr::rename("Protein_width"="# A.A", "Cys_number" = "# Cysteine(s)") %>%
    dplyr::rename("Score_sf"="Score Superfamily", "Score_class" = "Score Class") %>%
    dplyr::rename("seq"="Protein Sequence") %>%
    mutate(Conflict = Score_sf)
  
  protein_id <- DF %>% pull(ID)
  protein_id <- gsub("_[0-9]+_[0-9]+$","", protein_id)
  
  DF <- cbind(data.frame(protein_id), DF) %>% mutate(Method = gsub("_Regex.tab", "", Method)) %>% as_tibble()
  
  Conflictdf <- DF %>% filter(grepl("CONFLICT", Conflict)) %>% distinct(protein_id, Conflict)
  
  # Filter step as Borghie et al.
  
  DF <- DF %>% 
    filter(Hydrophobicity > Hydrophobicity_val) %>%
    filter(Protein_width >= pwidth_val) %>%
    # Omit conflicts
    mutate(Score_sf = gsub("[^0-9.-]", "", Score_sf)) %>%
    mutate(Score_class = gsub("[^0-9.-]", "", Score_class))
  
  
  DF1 <- DF %>% select(contains(c("protein_id","Method", "Score_sf","Cys_number","Superfamily ("))) %>%
    # filter(Score_sf > 0) %>% # Score_sf > 0 == Superfamily != "-"
    pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
    filter(Superfamily != "-") %>%
    mutate(Region = gsub("Superfamily ", "", Region)) %>%
    group_by(Method, protein_id, Cys_number) %>%
    summarise(across(Region,.fns = paste_col), Score_sf = n()) %>% ungroup()
  
  
  OUT <- DF %>% 
    select(contains(c("protein_id","Method", "Score_sf","Superfamily ("))) %>%
    # filter(Score_sf > 0) %>% 
    pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
    filter(Superfamily != "-") %>%
    mutate(Region = gsub("Superfamily ", "", Region)) %>%
    mutate(Superfamily =  gsub("\\(.*?\\)", "", Superfamily)) %>% 
    group_by(Method, protein_id) %>%
    summarise(across(Superfamily, .fns = paste_col), Score_sf = n()) %>% 
    # arrange(desc(Score_sf)) %>%
    left_join(DF1) %>%
    mutate(tab = "Regex")
  
  OUT <- OUT %>% left_join(Conflictdf)
  
  return(OUT)
  
  
}

read_pHMM <- function(f,  Hydrophobicity_val = 60, pwidth_val = 50, eval = 0.05) {
  
  DF <- read_delim(f, delim = "|", col_names = T) %>% mutate(Method = basename(f))
  
  # DF %>% select(contains(c("E-value", "Score","Biais"))) # how to include ??
  
  DF <- DF %>%
    dplyr::rename("Hydrophobicity"="% Hydrophobicity (Signal)", "ID" = "Read Name") %>%
    mutate(Hydrophobicity = gsub("%", "", Hydrophobicity), Hydrophobicity = as.double(Hydrophobicity)) %>%
    dplyr::rename("Protein_width"="# A.A", "Cys_number" = "# Cysteine(s)") %>%
    dplyr::rename("E_value"="E-value Superfamily") %>%
    mutate(Conflict = E_value)
  
  protein_id <- DF %>% pull(ID)
  protein_id <- gsub("_[0-9]+_[0-9]+$","", protein_id)
  
  DF <- cbind(data.frame(protein_id), DF) %>% mutate(Method = gsub("_pHMM.tab", "", Method)) %>% as_tibble()
  
  Conflictdf <- DF %>% filter(grepl("CONFLICT", Conflict)) %>% distinct(protein_id, Conflict)
  
  DF <- DF %>% 
    filter(Hydrophobicity > Hydrophobicity_val) %>%
    filter(Protein_width >= pwidth_val) %>%
    filter(as.numeric(E_value) < eval)

  # DF %>% 
  #   select(contains(c("protein_id","Method", "Score_sf","Superfamily ("))) %>%
  #   pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
  #   filter(Superfamily != "-") %>%
  #   mutate(Region = gsub("Superfamily ", "", Region)) %>%
  #   group_by(Method, protein_id) %>%
  #   summarise(across(Region,.fns = paste_col), Score_sf = n()) %>%
  #   mutate(tab = "Regex")
  # 
  
  DF1 <- DF %>% select(contains(c("protein_id","Method", "Cys_number","Superfamily ("))) %>%
    # filter(Score_sf > 0) %>% # Score_sf > 0 == Superfamily != "-"
    pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
    filter(Superfamily != "-") %>%
    mutate(Region = gsub("Superfamily ", "", Region)) %>%
    group_by(Method, protein_id, Cys_number) %>%
    summarise(across(Region,.fns = paste_col), Score_sf = n()) %>% ungroup()
  
  
  OUT <- DF %>% 
    select(contains(c("protein_id","Method", "Score_sf","Superfamily ("))) %>%
    pivot_longer(cols = all_of(which_cols), names_to = "Region", values_to = "Superfamily") %>%
    filter(Superfamily != "-") %>%
    mutate(Region = gsub("Superfamily ", "", Region)) %>%
    mutate(Superfamily =  gsub("\\(.*?\\)", "", Superfamily)) %>% 
    group_by(Method, protein_id) %>%
    summarise(across(Superfamily, .fns = paste_col), Score_sf = n()) %>% 
    left_join(DF1) %>%
    mutate(tab = "pHMM")
  
  OUT <- OUT %>% left_join(Conflictdf)
  
  return(OUT)
  
}

DB <- read_pHMM(pHHM_f) %>% rbind(read_regex(Regex_f))

DB %>% count(Superfamily, tab, Conflict) %>% view()

outName <- gsub("_Regex.tab|_pHMM.tab", "", basename(Regex_f))

outFile <- file.path(pub_dir, paste0(outName, "_ConoSorter_regex_pHMM.rds"))

DB %>% write_rds(file = outFile)

# read_delim(pHHM_f, delim = "|", col_names = T) %>% 
#   select(contains(c("E-value", "Score","Biais"))) %>%
#   pivot_longer(cols = all_of(names(.)), names_to = "Region", values_to = "value") %>%
#   mutate(Region = gsub("Superfamily ", "", Region)) %>%
#   filter(value != "-") %>%
#   filter(Region =="E-value Superfamily") %>%
#   mutate(facet = "E-value") %>%
#   mutate(facet = ifelse(grepl("Score", Region), "Score", facet)) %>%
#   mutate(facet = ifelse(grepl("Biais", Region), "Biais", facet)) %>%
#   ggplot(aes(as.numeric(value))) + geom_histogram() + facet_wrap(~ facet, scales = "free")


# read_delim(pHHM_f, delim = "|", col_names = T) %>% 
#   dplyr::rename("E_value"="E-value Superfamily") %>%
#   select(contains(c("Read Name","E-value"))) %>%
#   pivot_longer(-`Read Name`, names_to = "Region", values_to = "Superfamily") %>%
#   mutate(Region = gsub("E-value", "", Region))
