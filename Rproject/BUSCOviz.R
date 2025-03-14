
library(tidyverse)


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

dir <- "/Users/cigom/Documents/GitHub/conopeptides/03.Coverage/BUSCO_summaries/"

f <- list.files(dir, pattern = "BUSCO.tsv", full.names = T)


read_tsv(f) %>% count(my_species)
  
df <- read_tsv(f) %>% 
  mutate(my_species = gsub("_hisat_[s,S]uperDuper_", ".SuperDuper_", my_species)) %>% 
  separate(my_species, into = c("Method", "my_species"), sep = "_") %>%
  mutate(my_species = stringr::str_to_title(my_species)) %>%
  mutate(Method = stringr::str_to_title(Method))


df <- df %>% filter(my_species != "Mammalia")

df %>% count(Method)

recode_to <- c("Trinity",  
  "Trinity.superduper",
  "Rnaspades", 
  "Spades.superduper",
  "Merged.superduper",
  "Mmseqs",
  "Mmseqs.superduper")

recode_to <- structure(
  c("Trinity (T)", 
    "Trinity-Lace (Hisat)",
    "Spades (S)", 
    "Spades-Lace (Hisat)", 
    "Concat-Lace T-S (Hisat)",
    "MMseqs (T-S)",
    "MMseqs-Lace T-S (Hisat)"), 
  names = recode_to)

df <- mutate(df, Method = dplyr::recode_factor(Method, !!!recode_to, .ordered = T))

col <- c("#ED4647", "#EFE252", "#3A93E5", "#5BB5E7")

#names <- c("Complete", "Duplicated", "Fragmented", "Missing")

names <- c("S", "D", "F", "M")

labels <- c("Complete (C) and single-copy (S)",
  "Complete (C) and duplicated (D)",
  "Fragmented (F)  ",
  "Missing (M)")


labels <- structure(labels, names = names)

col <- structure(col, names = rev(names))

my_sp_lev <- c("Mollusca","Metazoa","Eukaryota","Mammalia","Bacteria")

figure <- df %>%
  mutate(facet = "Completeness") %>%
  # filter(my_species != "Bacteria") %>%
  mutate(category = factor(category, levels = rev(names))) %>%
  mutate(my_species = factor(my_species, levels = my_sp_lev)) %>%
  ggplot(aes(x = Method, y = my_percentage, fill = category)) +
  facet_grid(my_species~ facet) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.75) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(y = "% BUSCOs", x = "Assembly method") +
  coord_flip() +
  scale_fill_manual("", values = col, labels = rev(labels)) +
  guides(fill=guide_legend(nrow = 4)) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "bottom", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank())

# figure <- figure + facet_grid(~ category, scales = "free_x")

ggsave(figure, filename = 'BUSCO-methods.png', path = dir, width = 5, height = 7, device = png, dpi = 300)


Method_levs <- df %>% filter(category %in% c("S") & my_species == "Eukaryota") %>% arrange(my_percentage) %>% pull(Method) %>% as.character()


df %>%
  mutate(facet = "Completeness") %>%
  filter(category %in% c("F","D","S")) %>%
  mutate(category = dplyr::recode_factor(category, !!!labels)) %>%
  mutate(my_species = factor(my_species, levels = my_sp_lev)) %>%
  mutate(Method = factor(Method, levels = Method_levs)) %>%
  ggplot(aes(x = Method, y = my_percentage, fill = category)) +
  facet_grid(my_species~ category, scales = "free_x") +
  geom_col(position = position_stack(reverse = TRUE), width = 0.75) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(y = "% BUSCOs", x = "Assembly method") +
  coord_flip() +
  scale_fill_manual("", values = structure(col, names = rev(labels))) +
  guides(fill=guide_legend(nrow = 4)) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "none", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank())




df %>%
  mutate(facet = ifelse(grepl("Lace", Method),"Post-assembly step (Lace)", "Only assembly step")) %>%
  mutate(group ="Concat-Lace T-S (Hisat)") %>%
  mutate(group = ifelse(grepl("Trinity", Method),"Trinity", group)) %>%
  mutate(group = ifelse(grepl("Spades", Method),"Spades", group)) %>%
  mutate(group = ifelse(grepl("MMseqs", Method),"MMseqs", group)) %>%
  mutate(group = factor(group, levels = c("Trinity", "Spades", "MMseqs", "Concat-Lace T-S (Hisat)"))) %>%
  filter(category %in% c("F","D","S")) %>%
  mutate(category = dplyr::recode_factor(category, !!!labels)) %>%
  mutate(my_species = factor(my_species, levels = my_sp_lev)) %>%
  # mutate(Method = factor(Method, levels = Method_levs)) %>%
  ggplot(aes(x = group, y = my_values, fill = facet)) +
  facet_grid(my_species~category, scales = "free_x") +
  geom_col(position = position_dodge2(width = 1, preserve = "single"), width = 1) +
  # scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(y = "% Completeness (BUSCOs)", x = "Assembly method") +
  # scale_fill_manual("", values = c("black", "blue")) +
  guides(fill=guide_legend(nrow = 1)) +
  theme_bw(base_size = 14, base_family = "GillSans") +
  theme(legend.position = "top", 
    strip.background = element_rect(fill = 'grey89', color = 'white'), 
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line.x = element_blank(),
    axis.line.y = element_blank())



df %>%
  filter(category %in% c("F","D","S")) %>%
  mutate(category = dplyr::recode_factor(category, !!!labels)) %>%
  select(-my_values) %>%
  pivot_wider(names_from = category, values_from = my_percentage) %>%
  view()
