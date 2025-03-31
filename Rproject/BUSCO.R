# 

rm(list = ls())

if(!is.null(dev.list())) dev.off()

library(tidyverse, help, pos = 2, lib.loc = NULL)

pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR"

dir <- "/Users/cigom/Documents/GitHub/conopeptides/03.Coverage/BUSCO_summaries/full_table_busco_dir"

f <- list.files(dir, "tsv", full.names = T)

f[1]


read_full_table <- function(f) {
  
  g <- gsub(".tsv", "",basename(f))
  g <- gsub("full_table_", "", g)
  g <- gsub("_odb[0-9]+$", "", g)
  g <- stringr::str_to_title(g)
  
  cat("\n", g, "\n")

  last_pos <- length(strsplit(g, "_")[[1]])
  
  str(odb <- sapply(strsplit(g, "_"), `[`, last_pos))
  
  Meth <- gsub(paste0("_", odb), "", g)
  
  # df <- read.delim(f, comment.char = "#", sep = "\t", header = F, na.strings = c(" ", "NA"))
  
  df <- read.delim(f, skip = 4, sep = "\t", header = T, na.strings = c(" ", "NA"))
  
  # df <- read_tsv(f, comment = "#",  col_names = F, na = " ", quoted_na = " ")
  
  names(df) <- c("Busco_id",	"Status",	"Sequence",	"Score",	"Length")
  
  df <- df %>% as_tibble() %>% mutate(Method = Meth, odb = stringr::str_to_title(odb))
  

  
  return(df)
}


df <- lapply(f, read_full_table)

head(df <- do.call(rbind, df))


# Recode

  
recode_to <- c(
  "Trinity",  
  "Trinity_hisat_superduper",
  "Rnaspades", 
  "Spades_hisat_superduper",
  "Mmseqs",
  "Merged_polya_hisat_superduper",
  "Merged_hisat_superduper",
  "Mmseqs_hisat_superduper")

recode_to <- structure(
  c(
    "Trinity (T)", 
    "Full-length Trinity (Lace)", 
    "Spades (S)", 
    "Full-length Spades (Lace)", 
    "Full-length T-S (Mmseq)", 
    "Full-length T-S (Lace)",
    NA,NA
  ), 
  names = recode_to)

df <- mutate(df, Method = dplyr::recode_factor(Method, !!!recode_to))


# Count Total BUSCO groups searched

df %>% distinct(Busco_id, odb) %>% count(odb)

df %>%
  # filter()
  count(Method)

# Count BUSCO notation for files

df %>% 
  drop_na(Method) %>%
  count(Method, odb, Status) %>%
  write_tsv(file = paste0(pub_dir, "/busco_summaries.tsv"))

# 

df %>%
  drop_na(Method) %>%
  filter(Status %in% "Complete") %>%
  ggplot(aes(Score, Length)) + 
  geom_point() +
  facet_grid(Status ~ Method)

df %>%
  drop_na(Method) %>%
  filter(Status %in% "Complete") %>%
  ggplot(aes(y = Method, x = Score)) + 
  geom_violin()
  