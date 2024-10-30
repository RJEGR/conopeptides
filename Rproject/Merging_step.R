

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

read_outfmt6 <- function(f) {
  
  # seqid = transcript_id
  outfmt6.names <- c("transcript_id", "subject", "identity", "coverage", "mismatches", "gaps", "seq_start", "seq_end", "sub_start", "sub_end", "e", "score")
  
  
  df <- read_tsv(f, col_names = F) 
  
  colnames(df) <- outfmt6.names
  
  df <- df %>% mutate(db = basename(f))
  
  return(df)
  
  
}

dir <- "~/Documents/GitHub/conopeptides/04.Merge/"


f <- list.files(dir, pattern = "outfmt6", full.names = T) 

# 

df <- do.call(rbind, lapply(f, read_outfmt6))

# recode_to <- c() 

# structure(recode_to, names = c("Merged.0.98", "Merged.0.99", ))

# df %>%. dplyr::mutate(pH = dplyr::recode_factor(pH, !!!recode_to))

df <- df %>% mutate(db = gsub("vs_conoserver_protein.fa.diamond.blastx.outfmt6","", db))

# df %>% ggplot(aes(identity,-log10(e), color = db)) + geom_point()


# df %>% 
#   count(db, round(identity)) %>%
#   group_by(db) %>% mutate(frac = n / sum(n)) %>%
#   ggplot(aes(x = `round(identity)`, y = frac)) +
#   geom_col(aes(fill = db)) +
#   facet_grid(~ db) +
#   # ggplot2::stat_ecdf(linewidth = 1, alpha = 0.5) +
#   ggthemes::scale_fill_calc() + theme_bw()

df %>% ggplot(aes(identity, color = db, fill = db)) +
  geom_histogram(position = position_identity()) +
  facet_grid(~ db) +
  # ggplot2::stat_ecdf(linewidth = 1, alpha = 0.5) +
  ggthemes::scale_fill_calc() + ggthemes::scale_color_calc()  + theme_bw()

df %>% ggplot(aes(identity, y = db, fill = db)) + 
  ggridges::geom_density_ridges() +
  ggthemes::scale_fill_calc() + theme_bw()

# Split by conoserver_protein.tsv


dir <- "~/Documents/GitHub/conopeptides/"

conoserverDB <- read_tsv(file.path(dir, "conoserver_protein.tsv")) %>%
  rename("subject" = "identifier") 

df %>%
  mutate(subject =  sapply(strsplit(subject, "[|]"), `[`, 1)) %>%
  left_join(conoserverDB) %>%
  drop_na(`gene superfamily`) %>%
  ggplot(aes(y = `gene superfamily`, x = identity, fill = db)) +
  facet_grid(~ db) +
  ggridges::geom_density_ridges(
    jittered_points = T,
    position = ggridges::position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 1, point_alpha = 1, alpha = 0.7) +
  ggthemes::scale_fill_calc() + theme_bw()

# 

p <- df %>%
  mutate(subject =  sapply(strsplit(subject, "[|]"), `[`, 1)) %>%
  left_join(conoserverDB) %>%
  drop_na(`gene superfamily`) %>%
  mutate(facet = `gene superfamily`) %>%
  count(db, facet) %>%
  group_by(db) %>%
  arrange(desc(n)) %>% 
  mutate(facet = factor(facet, levels=unique(facet))) %>%
  ggplot(aes(y = db, x = n, fill = db)) + 
  facet_grid(facet ~., scales = "free_y", space = "free_y") +
  geom_col(position = position_dodge2(width = 0.5)) +
  ggthemes::scale_fill_calc() + theme_bw(base_size = 7, base_family = "GillSans") +
  theme(strip.background.y = element_rect(fill = 'grey89', color = 'white'),
    strip.text.y = element_text(angle = 0, size = 10, hjust = 0))


ggsave(p, filename = 'conoserver_protein_by_assembly.png', path = dir, 
  width = 5.5, height = 12, device = png, dpi = 300)
  
