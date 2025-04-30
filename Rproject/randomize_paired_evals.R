
# Read blastn results (optional)
# Read transrate scores (to evals accuracy of assembly methods)
# Read detonate scores (optional)
#


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR"

dir <- "~/Downloads/CONOSERVERDB/"

filepattern <- "conoserver_nucleic_californicus_length_100.fa"

f <- list.files(dir, pattern = filepattern, full.names = T)

seqs <- Biostrings::readDNAStringSet(f)

reff <- "conoserver_nucleic.fa"
reff <- list.files(dir, pattern = reff, full.names = T)
refseq <- Biostrings::readDNAStringSet(reff)
  
str(hits <- sapply(strsplit(names(refseq), "[|]"), `[`, 1))
str(Desc <- sapply(strsplit(names(refseq), "[|]"), `[`, 2))
str(sp <- sapply(strsplit(names(refseq), "[|]"), `[`, 3))

refseqdf <- data.frame(hits, Desc, sp) %>% as_tibble()

names(refseq) <- hits

sum(keep <- hits %in% names(seqs)) # match to 123 records
refseq <- refseq[keep]


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

# Include classification of toxin or prediction (Transdecoder > ConoSorter)
recode_to <- structure(c("Precursor", "Incomplete", "Incomplete", "Incomplete"), 
  names = c("(Mature)_(Pro-region)_(Signal)", 
    "(Mature)", 
    "(Pro-region)_(Signal)","(Signal)"))


subdirs <- list.files(dir, pattern = "ConoSorter_dir", full.names = T)

subdirs <- list.files(subdirs, pattern = "_Regex.tab", full.names = T)

read_conosort <- function(path) {
  
  paste_col <- function(x) { 
    x <- x[!is.na(x)] 
    x <- unique(sort(x))
    x <- paste(x, sep = '_', collapse = '_')
  }
  
  which_cols <- c("Superfamily (Signal)", "Superfamily (Pro-region)", "Superfamily (Mature)")
  
  
  f <- path[1]
  
  basedir <- sapply(strsplit(basename(f), "[.]"), `[`, 1)
  
  basedir <- gsub("_longest_orfs_Regex", "", basedir)
  
  DF <- read_delim(f, delim = "|", col_names = T) %>% mutate(Method = basedir) %>%
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
    # filter(Hydrophobicity > Hydrophobicity_val) %>%
    # filter(Protein_width >= pwidth_val) %>%
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

conoSortdf <- do.call(rbind, lapply(subdirs, read_conosort)) %>% ungroup()

conoSortdf <- conoSortdf %>% 
  mutate(contig_name = gsub(".p[0-9]+$","", protein_id)) %>%
  mutate(Region = dplyr::recode_factor(Region, !!!recode_to))
  
conoSortdf %>% count(tab, Region)

conoSortdf %>% count(Method)

conoSortdf %>% distinct(contig_name)

any(transratedf$contig_name %in% conoSortdf$contig_name)

transratedf <- transratedf %>% left_join(conoSortdf)

transratedf %>% count(tab, Region)

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
  labs(y = "Method", x = "Accuracy (N contigs/ N hits to reference)") +
  ggthemes::scale_fill_calc() +
  scale_x_continuous(position = "top") +
  theme(legend.position = "bottom", 
    # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 7))

p2 <- p2 + 
  geom_text(data = filter(barpdf, col == "Hit"),
    aes(label= label), hjust= 1.05, vjust = 0.5, size = 2.5, family = "GillSans", color = "white")

p2 <- p2 + guides(fill=guide_legend(title = ""))

# Plot 3

# load the DECIPHER library in R
# library(DECIPHER)

# cluster the sequences
clusters <- DECIPHER::Clusterize(refseq,
  cutoff=0.5, # < 50% distant
  minCoverage=0.5, # > 50% coverage
  processors=NULL) # use all CPUs


alignment <- DECIPHER::AlignSeqs(refseq)

# This method is time-computer consuming
phangAlign <- phangorn::phyDat(as(alignment, "matrix"), type="DNA") 
dm <- phangorn::dist.ml(phangAlign)
treeNJ <- phangorn::NJ(dm) # Note, tip order != sequence order
plot(treeNJ)
fit <- phangorn::pml(treeNJ, data=phangAlign) # this step is fast
fitGTR <- stats::update(fit, k = 4, inv = 0.2) # this step is fast
fitGTR <- phangorn::optim.pml(fitGTR) # this step is fast

t <- fitGTR$tree

# calculate Pairwise Distances from a Phylogenetic Tree 
# disimilariy matrx using cophenetic.phylo
hclust <- hclust(dist(ape::cophenetic.phylo(t)))

plot(t)
plot(hclust)
# view the cluster numbers
table(clusters)


heatmapdf <- transratedf %>% drop_na(hits) 

heatmapdf <- clusters %>% as_tibble(rownames = "hits") %>% left_join(heatmapdf)

heatmapdf <- heatmapdf %>% left_join(refseqdf, by = "hits")

labels <- heatmapdf %>% distinct(hits, Desc) %>% pull(Desc, name = hits)



p3 <- heatmapdf %>%
  # mutate(col = ifelse(is.na(Method), "Not-assembled", "Assembled")) %>%
  ggplot(aes(y = Method, x = hits)) + #  fill = as.factor(Score_sf)
  geom_tile(fill = "black", color = 'white', width = 0.5, height = 0.25) +
  # facet_grid(~ cluster, space = "free", scales = "free", switch = "x") +
  # scale_fill_viridis_c(option = "magma", direction = -1) +
  # scale_x_discrete("Conopeptides (hits of contigs assembled)",position = "top") +
  scale_y_discrete(limits = rev(recode_to), breaks = rev(recode_to)) +
  theme_bw(base_family = "GillSans", base_size = 12) + 
  theme(legend.position = "bottom", 
    # axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    # axis.text.x = element_text(angle = 90, size = 5), 
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    # panel.grid.minor.x = element_blank(),
    # panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 7),
    axis.title.x = element_text(size = 7),
    strip.background = element_rect(fill = 'white', color = 'white', size = 10)) +
  labs (title = "Artificial fastq true set analysis", x = "Sequence alignment (Pairwise Distances)") 

library(ggh4x)

p3 <- p3 + 
  ggh4x::scale_x_dendrogram(hclust = hclust, position = 'bottom', labels = NULL) +
  guides(x.sec = guide_axis_manual(labels = labels, label_size = 5, angle = 70, label_family = "GillSans"))  

ggsave(p3, filename = 'randomize_paired_evals_heatmap.png', path = pub_dir, width = 10, height = 5, device = png, dpi = 700)

# plot 4: barplot of toxins


conoSortdf %>% count(Region)

bardf <- transratedf %>% 
  mutate(hits = ifelse(is.na(hits), "No-hit", "Hit")) %>%
  mutate(protein_id = ifelse(is.na(protein_id), "No-orf", "orf")) %>%
  count(Method, hits, protein_id, Region) %>%
  drop_na(Region)

bardf <- conoSortdf %>% 
  filter(Method %in% "conoserver_nucleic_californicus_legth") %>%
  mutate(hits = "Hit") %>%
  mutate(protein_id = "orf") %>%
  count(Method, hits, protein_id, Region) %>%
  rbind(bardf)

# transratedf %>% drop_na(Superfamily) %>% count(Method, Superfamily)
bardf %>%
  pivot_wider(names_from = hits, values_from = n) %>%
  # mutate(`No-hit` = ifelse(is.na(`No-hit`), "", `No-hit`)) %>%
  mutate(label = ifelse(!is.na(`No-hit`), paste(Hit, " (",`No-hit`,")"), Hit)) %>%
  ggplot(aes(y = Method, x = Region, fill = Hit)) +
  scale_y_discrete(limits = rev(recode_to), breaks = rev(recode_to)) +
  geom_tile(color = "white",
    lwd = 0.5,
    linetype = 1) +
  geom_text(aes(label= label), hjust= 0.5, vjust = 0.5, size = 5, family = "GillSans", color = "white") +
  theme_bw(base_size = 16, base_family = "GillSans") 

p4 <- bardf %>%
  pivot_wider(names_from = hits, values_from = n) %>%
  # mutate(`No-hit` = ifelse(is.na(`No-hit`), "", `No-hit`)) %>%
  mutate(label = ifelse(!is.na(`No-hit`), paste(Hit, " (",`No-hit`,")"), Hit)) %>%
  ggplot(aes(y = Method, x = Hit, fill = Region)) +
  # facet_grid(~ Region, space = "free_x", scales = "free_x") +
  geom_col() +
  labs(y = "") +
  scale_y_discrete(limits = rev(recode_to), breaks = rev(recode_to)) +
  geom_text(aes(label= label), hjust= 1.05, vjust = 0.5, size = 2.5, family = "GillSans", color = "white") +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "bottom",
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    strip.background = element_rect(fill = 'white', color = 'white', size = 10),
    axis.title.x = element_text(size = 7)) +
  ggthemes::scale_fill_calc() 

p4 <- p4 + scale_x_continuous("N contigs with assignment",position = "top") +  
  guides(fill=guide_legend(title = ""))

p4

library(patchwork)

# plot_spacer() + pright + plot_layout(widths = c(1, -0.05, 2))

psave <- p2 + p1 + plot_spacer() + p4 + plot_layout(widths = c(0.2, 0.2, -0.02 ,0.2))

psave

ggsave(psave, filename = 'randomize_paired_evals.png', path = pub_dir, width = 7.5, height = 2.2, device = png, dpi = 700)


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
