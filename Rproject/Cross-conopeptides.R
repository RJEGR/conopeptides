
# Cross-check evaluation of conopeptides identification
# LOAD COnodictor 
# Output and  include in Database.R
# 

# Load conodictor
pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR"

dir <- "/Users/cigom/Documents/GitHub/conopeptides/05.Prediction/conodictor_dir/"

subdirs <- list.files(dir, pattern = "_dir")

# f <- list.files(dir, pattern = "summary.csv", full.names = T)

# Conosorter -----

cnsrtr_f <- list.files(path = pub_dir, pattern = "ConoSorter_regex_pHMM.rds", full.names = T)

SORTERDB <- read_rds(cnsrtr_f)


# conodictor ----

read_conodictor <- function(path) {
  
  subdir <- path[1]
  
  # basedir <- sapply(strsplit(subdir, "_"), `[`, 1)
  
  basedir <- basename(subdir)
  
  f <- file.path(dir, subdir, "summary.csv")
  
  OUT <- read_csv(f) # %>% mutate(Method = basedir)
  
  names(OUT) <- paste0(names(OUT), "_conodictor")
  
  return(OUT)
}



DICTORDB <- read_conodictor(subdirs[2]) %>% dplyr::rename("protein_id" = "sequence_conodictor")

write_rds(DICTORDB, file.path(pub_dir, "Conodictor2.rds"))

# Crosscheck -----
# Equivalence between DICTOR = SORTER
# hmm_pred = Superfamily
# pssm_pred = Region
# definitive_pred = hmm_pred

DICTORDB %>% count(hmm_pred_conodictor, pssm_pred_conodictor, definitive_pred_conodictor)


DICTORDB %>% distinct(protein_id) # DICTOR found 8,514 from Merged_polyA_hisat_SuperDuper.fasta.transdecoder_dir

SORTERDB %>% distinct(protein_id) # SORTER found 4,243 Merged_polyA_hisat_SuperDuper.fasta.transdecoder_dir

any(SORTERDB$protein_id %in% DICTORDB$protein_id)

which_tools <- function(x) { 
  x <- x[!is.na(x)] 
  n <- length(unique(x))
  x <- unique(x)
  
  if(n > 1) {
    x <- "BOTH"
  } else
    x <- paste(x, sep = '|', collapse = '|') }



DICTORDB %>% right_join(SORTERDB) %>% count(protein_id, sort = F)

DF1 <- DICTORDB %>% distinct(protein_id) %>% mutate(predicted = "Conodictor")
DF2 <- SORTERDB %>% ungroup() %>% distinct(protein_id) %>% mutate(predicted = "ConoSorter")

DF <- rbind(DF1, DF2) %>%
  group_by(protein_id) %>%
  summarise(
    across(predicted, .fns = which_tools), 
    .groups = "drop_last")

DF %>% count(predicted)

DF %>% right_join(SORTERDB) %>% count(predicted, Conflict)

DF %>% right_join(DICTORDB) %>% count(predicted, pssm_pred_conodictor)
