
# Read blastn results (optional)
# Read transrate scores (to evals accuracy of assembly methods)
# Read detonate scores (optional)
#


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

dir <- "~/Downloads/CONOSERVERDB/"

filepattern <- "conoserver_nucleic_californicus_length_100.fa"

f <- list.files(dir, pattern = filepattern, full.names = T)

seqs <- Biostrings::readDNAStringSet(f)

dir <- "/Users/cigom/Documents/GitHub/conopeptides/Artificial_rnaseq_dir/randomize_paired_evals_dir/"

subdirs <- list.files(dir, pattern = "_transrate_dir")

read_transrate_scores <- function(path) {
  
  # Path where 
  
  subdir <- path[1]
  
  basedir <- sapply(strsplit(subdir, "_"), `[`, 1)
  
  # read_csv(list.files(subdir1, "contigs.csv", full.names = T))
  
  f <- file.path(dir, subdir, basedir, "contigs.csv")
  
  read_csv(f) %>% mutate(Method = basedir)
  
}

transratedf <- lapply(subdirs, read_transrate_scores)

transratedf <- do.call(rbind,transratedf)

transratedf %>% count(Method)

# Transrate -----

# as there was a strong monotonic relationship between contig accuracy and TransRate contig score, use this column to represent accuracy of the assembly contig


recode_to <- c("SuperDuper", "Spades", "Trinity")

transratedf <- transratedf %>%
  mutate(Method = factor(Method, levels = rev(recode_to)))

# plot 1

p1 <- transratedf %>%
  drop_na(hits) %>%
  ggplot(aes(y = Method, x = reference_coverage, fill = after_stat(x))) +
  # geom_violin() +
  # facet_grid(~ facet) +
  ggridges::geom_density_ridges_gradient(
    jittered_points = T,
    position = ggridges::position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 3, point_alpha = 1, alpha = 1) +
  scale_fill_viridis_c(option = "C") +
  labs(y = "", x = "Reference coverage (TransRate)") +
  theme_bw(base_family = "GillSans", base_size = 12) + 
  theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
    axis.title.x = element_text(size = 7)) +
  scale_x_continuous(position = "top", breaks = c(0.8, 0.9, 1))

p1

# plot 2

# transratedf %>% ggplot(aes(reference_coverage, cpg_ratio)) + geom_point()

barpdf <- transratedf %>% 
  # Label those not found in the reference (Ho: chimera?)
  mutate(col = ifelse(is.na(hits), "No-hit", "Hit")) %>%
  # filter(reference_coverage)
  count(Method, col) %>%
  mutate(label = paste0("(", n,")"))

p2 <- barpdf %>%
  ggplot(aes(x = n/123, y = Method, fill = col)) +
  geom_col() + 
  theme_bw(base_family = "GillSans", base_size = 12) +
  labs(y = "", x = "Accuracy (N contigs/ N hits to reference)") +
  ggthemes::scale_fill_calc() +
  scale_x_continuous(position = "top") +
  theme(legend.position = "bottom", axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 7))

p2 <- p2 + 
  geom_text(data = filter(barpdf, col == "Hit"),
    aes(label= label), hjust= 1.05, vjust = 0.5, size = 2.5, family = "GillSans", color = "white")

# Plot 3


heatmapdf <- transratedf %>% drop_na(hits) 

heatmapdf <- clusters %>% as_tibble(rownames = "hits") %>% left_join(heatmapdf)

p3 <- heatmapdf %>%
  # mutate(col = ifelse(is.na(hits), "No-hit", "Hit")) %>%
  ggplot(aes(y = Method, x = hits, fill = reference_coverage)) +
  geom_tile(color = 'white', width = 0.5, height = 0.25) +
  facet_grid(~ cluster, space = "free", scales = "free", switch = "x") +
  scale_fill_viridis_c(option = "magma", direction = -1) +
  scale_x_discrete("Conopeptides (hits of contigs assembled)",position = "top") +
  theme_bw(base_family = "GillSans", base_size = 12) + 
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 7),
    axis.title.x = element_text(size = 7),
    strip.background = element_rect(fill = 'white', color = 'white', size = 10)) +
  labs (title = "Artificial fastq true set analysis") 

# p3

library(patchwork)

# plot_spacer() + pright + plot_layout(widths = c(1, -0.05, 2))

psave <- p3 + p2 + plot_spacer() + p1 + plot_layout(widths = c(0.5, 0.2, -0.02 ,0.2))

pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR"

ggsave(psave, filename = 'randomize_paired_evals.png', path = pub_dir, width = 10, height = 3, device = png, dpi = 700)


paste_col <- function(x) { 
  x <- x[!is.na(x)] 
  x <- unique(sort(x))
  x <- paste(x, sep = '', collapse = '||')
}

hitsdf <- transratedf %>%
  drop_na(hits) %>%
  group_by(hits) %>%
  summarise(across(Method,.fns = paste_col))

# Why some conotoxin are assembled than others?
# Measure in the reference
# GC content
# Sequence length
# Entropy
# Complexity metrics
# Seq homology

str(query <- transratedf %>% drop_na(hits) %>% distinct(hits) %>% pull())

str(Names <- sapply(strsplit(names(seqs), "[|]"), `[`, 1))

sum(keep <- Names %in% query)

# seqs <- seqs[keep]

DB <- data.frame(seqs, widthref = Biostrings::width(seqs)) %>% 
  as_tibble(rownames = "id") %>% 
  separate(id, into = c("hits","Desc", "sp"), sep = "[|]") %>%
  left_join(hitsdf)
  
DB %>% distinct(seqs) %>% pull()

# Seq homology analysis

# load the DECIPHER library in R
library(DECIPHER)

# cluster the sequences
clusters <- Clusterize(seqs,
  cutoff=0.5, # < 50% distant
  minCoverage=0.5, # > 50% coverage
  processors=NULL) # use all CPUs

# view the cluster numbers
table(clusters)

# build a maximum likelihood phylogenetic tree =====

# align coding sequences
aligned <- AlignTranslation(seqs,
  type="DNAStringSet") # choose AA or DNA

aligned <- seqs

BrowseSeqs(aligned, highlight=0)

# optimize the tree (hard to run)

tree <- DECIPHER::TreeLine(aligned,
  method="ML",
  model=MODELS, # choose a model or test all
  showPlot=TRUE,
  processors=NULL) # use all CPUs
