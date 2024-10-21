# Blast DataViz
# Identify which peptides from conoserver database was identify by assembly
# Intersected and distinct
# Bind Subjects to conopeptides main groups (search how)

library(tidyverse)

read_outfmt6 <- function(f) {
  
  # seqid = transcript_id
  outfmt6.names <- c("transcript_id", "subject", "identity", "coverage", "mismatches", "gaps", "seq_start", "seq_end", "sub_start", "sub_end", "e", "score")
  
  
  df <- read_tsv(f, col_names = F) 
  
  colnames(df) <- outfmt6.names
  
  df <- df %>% mutate(db = basename(f))
  
  return(df)
  
  
}


dir <- "~/Documents/GitHub/conopeptides/03.Coverage/"


f <- list.files(dir, pattern = "outfmt6", full.names = T) 


df <- do.call(rbind, lapply(f, read_outfmt6))

df %>% 
  filter(identity > 90) %>%
  distinct(subject, identity)

df %>% ggplot(aes(x = identity, y = db)) + ggridges::geom_density_ridges()

df %>% distinct(subject, db) %>%
  group_by(subject) %>%
  summarise(
    across(db, .fns = list), n = n(), .groups = "drop_last") %>%
  mutate(db = ifelse(n == 2, "Both", db)) %>%
  mutate(db = unlist(db)) %>%
  count(db) %>%
  ggplot(aes(x = n, y = db)) +
  geom_col()
  