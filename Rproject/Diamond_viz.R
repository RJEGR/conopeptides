# plot bit score from outfmt6 from DIAMOND blastx
# Files: dir <- "/Users/cigom/Documents/GitHub/conopeptides/Nx_Metrics_dir/"
# Reference: short- REFERENCE (PRECURSOR) conotoxin sequences ()

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

# dir <- "~/Documents/GitHub/conopeptides/04.Merge/"

# dir <- "~/Documents/GitHub/conopeptides/03.Coverage/"

dir <- "/Users/cigom/Documents/GitHub/conopeptides/Nx_Metrics_dir/outfmt6_dir/"

f <- list.files(dir, pattern = "outfmt6", full.names = T) 

f <- f[grepl("curated", f)]

# 

df <- do.call(rbind, lapply(f, read_outfmt6))

df <- df %>% mutate(db = gsub("_vs_conoserver_protein_curated.fa.diamond.blastx.outfmt6","", db))

df %>% dplyr::count(db)

# df <- df %>% mutate(db = gsub("_vs_conoserver_protein.fa.diamond.blastx.outfmt6","", db))
# df <- df %>% mutate(db = gsub("_vs_conoserver_protein.diamond.blastx.outfmt6","", db))


recode_to <- c("Trinity",  
  "trinity_hisat_superDuper",
  "spades", 
  "spades_hisat_superDuper",
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

df <- mutate(df, db = dplyr::recode_factor(db, !!!recode_to))


df %>% dplyr::count(db)


# df %>% ggplot(aes(identity,-log10(e), color = db)) + geom_point()


# df %>% 
#   count(db, round(identity)) %>%
#   group_by(db) %>% mutate(frac = n / sum(n)) %>%
#   ggplot(aes(x = `round(identity)`, y = frac)) +
#   geom_col(aes(fill = db)) +
#   facet_grid(~ db) +
#   # ggplot2::stat_ecdf(linewidth = 1, alpha = 0.5) +
#   ggthemes::scale_fill_calc() + theme_bw()

df %>% ggplot(aes(-log10(e), color = db, fill = db)) +
  geom_histogram(position = position_identity(), alpha = 0.5) +
  facet_grid(db ~.) +
  # ggplot2::stat_ecdf(linewidth = 1, alpha = 0.7) +
  ggthemes::scale_fill_calc() + ggthemes::scale_color_calc()  + theme_bw()

df %>% ggplot(aes(mismatches, y = db, fill = db)) + 
  ggridges::geom_density_ridges() +
  ggthemes::scale_fill_calc() + theme_bw()

# df %>% 
#   mutate(facet = ifelse(grepl("MMseq", db), "MMseqs", db)) %>%
#   mutate(facet = ifelse(grepl("Trinity", db), "Trinity", db)) %>%
#   mutate(facet = ifelse(grepl("Spades", db), "Spades", db)) %>%
#   ggplot(aes(x = score, y = facet, color = db)) +
#   geom_violin(position = position_fill()) + labs(y = "", x = "Blast bit score") +
#   ggthemes::scale_fill_calc() + ggthemes::scale_color_calc()  + 
#   theme_bw(base_family = "GillSans", base_size = 12)

# Bit score profile ----

df %>% ggplot(aes(x = score, y = db, color = db, fill = db)) +
  geom_violin() + labs(y = "", x = "Blast bit score") +
  ggthemes::scale_fill_calc() + ggthemes::scale_color_calc()  + theme_bw(base_family = "GillSans", base_size = 14)

p1 <- df %>% 
  # ggplot(aes(score, y = db, fill = db)) + 
  ggplot(aes(y = db, x = score,fill = after_stat(x))) +
  # geom_violin() +
  ggridges::geom_density_ridges_gradient(
    jittered_points = T,
    position = ggridges::position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 3, point_alpha = 1, alpha = 1) +
  scale_fill_viridis_c(option = "C") +
  labs(y = "", x = "Blast bit score") +
  theme_bw(base_family = "GillSans", base_size = 14) + theme(legend.position = "none")


p2 <- df %>% 
  distinct(subject, db) %>%
  dplyr::count(db) %>%
  ggplot(aes(x = n, y = db)) +
  geom_col(fill = "black") + labs(y = "", x = "Unique partially conopeptide annotated") +
  geom_text(aes(label= scales::comma(n)), hjust= 1.2, vjust = 0.5, size = 7, family = "GillSans", color = "white") +
  # ggthemes::scale_fill_calc() + ggthemes::scale_color_calc()  + 
  theme_bw(base_family = "GillSans", base_size = 14)

p2 <- p2 + theme(
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.y = element_blank(),
  axis.line.y = element_blank()) 

library(patchwork)

p1 + p2

# df %>% 
#   mutate(label = ifelse(e > 1, "Random", ))


df %>% ggplot(aes(x = score, y = -log10(e), color = db, fill = db)) +
  geom_point()
  ggthemes::scale_fill_calc() + ggthemes::scale_color_calc()  + theme_bw()


# Split by conoserver_protein.tsv


dir <- "~/Documents/GitHub/conopeptides/"

conoserverDB <- read_tsv(file.path(dir, "conoserver_protein.tsv")) %>%
  dplyr::rename("subject" = "identifier") 



conoserverDB %>% dplyr::count(`gene superfamily`, sort = T)



# Any of the Eliger et al., 2011?
QUERY <- conoserverDB %>% filter(grepl("Ca1.[0-9]", name)) %>% distinct(subject) %>% pull()


df %>% 
  mutate(subject =  sapply(strsplit(subject, "[|]"), `[`, 1)) %>% 
  filter(subject %in% QUERY)

df %>%
  mutate(subject =  sapply(strsplit(subject, "[|]"), `[`, 1)) %>%
  left_join(conoserverDB) %>%
  distinct(db, `protein type`, subject) %>%
  dplyr::count(db, `protein type`)

df %>%
  mutate(subject =  sapply(strsplit(subject, "[|]"), `[`, 1)) %>%
  left_join(conoserverDB) %>% 
  # drop_na(organism) %>%
  distinct(db, organism, subject) %>%
  dplyr::count(db, organism, sort = T) 


df %>%
  mutate(subject =  sapply(strsplit(subject, "[|]"), `[`, 1)) %>%
  left_join(conoserverDB) %>%
  filter(`protein type` %in% "Conus californicus")


# Looking for a obiquitinious superfamily found by assembly method

df %>%
  mutate(subject =  sapply(strsplit(subject, "[|]"), `[`, 1)) %>%
  left_join(conoserverDB) %>% 
  filter(grepl("T superfamily", `gene superfamily`)) %>% 
  distinct(db, subject, `protein type`) %>% 
  dplyr::count(subject, sort = T)
  dplyr::count(db, `protein type`, sort = T)

df %>% 
  mutate(subject =  sapply(strsplit(subject, "[|]"), `[`, 1)) %>%
  filter(subject %in% "P00659")
  
  
df %>%
  mutate(subject =  sapply(strsplit(subject, "[|]"), `[`, 1)) %>%
  left_join(conoserverDB) %>%
  drop_na(`gene superfamily`) %>%
  filter(grepl("T superfamily", `gene superfamily`)) %>%
  ggplot(aes(score, color = db)) +
  ggplot2::stat_ecdf(linewidth = 1, alpha = 0.5)


df %>%
  mutate(subject =  sapply(strsplit(subject, "[|]"), `[`, 1)) %>%
  left_join(conoserverDB) %>%
  drop_na(`gene superfamily`) %>%
  filter(grepl("T superfamily", `gene superfamily`)) %>% 
  ggplot(aes(y = `db`, x = score, fill = db)) +
  facet_grid(~ `gene superfamily`) +
  geom_violin()
  
df %>%
  mutate(subject =  sapply(strsplit(subject, "[|]"), `[`, 1)) %>%
  left_join(conoserverDB) %>%
  drop_na(`gene superfamily`) %>%
  filter(grepl("T superfamily", `gene superfamily`)) %>% 
  ggplot(aes(y = `db`, x = score, fill = db)) +
  facet_grid(~ `gene superfamily`) +
  ggridges::geom_density_ridges(
    jittered_points = T,
    position = ggridges::position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 1, point_alpha = 1, alpha = 0.7) +
  ggthemes::scale_fill_calc() + theme_bw()


df %>%
  mutate(subject =  sapply(strsplit(subject, "[|]"), `[`, 1)) %>%
  left_join(conoserverDB) %>%
  drop_na(`gene superfamily`) %>%
  ggplot(aes(y = `gene superfamily`, x = score, fill = db)) +
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

df %>%
  mutate(subject =  sapply(strsplit(subject, "[|]"), `[`, 1)) %>%
  left_join(conoserverDB) %>%
  mutate(facet = `gene superfamily`) %>%
  drop_na(facet) %>%
  distinct(db, subject, facet) %>%
  count(db, facet) %>%
  group_by(db) %>%
  arrange(desc(n)) %>% 
  mutate(facet = factor(facet, levels=unique(facet))) %>%
  ggplot(aes(x = db, y = facet, fill = n)) +
  geom_tile(color = 'white', linewidth = 0.7, width = 1) +
  theme_classic(base_family = "GillSans", base_size = 12) +
  theme(legend.position = 'top', 
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
  
