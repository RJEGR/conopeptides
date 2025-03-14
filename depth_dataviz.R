

library(tidyverse)

# Normalize position to the sequence length from fasta file 

dir <- "/Users/cigom/STATS/"


f <- list.files(dir, full.names = T, pattern = ".depth.txt")

df <- read_tsv(f, col_names = F)

names(df) <- c("Chr", "Position", "Depth")

max(df$Position)

# normalize length to width----

dir <- "~/Downloads/CONOSERVERDB/"

f <- list.files(dir, pattern = "conoserver_nucleic.fa.gz", full.names = T)

Labels <- names(Biostrings::readDNAStringSet(f))


str(sp <- sapply(strsplit(Labels, "[|]"), `[`, 3))
str(Descr <- sapply(strsplit(Labels, "[|]"), `[`, 2))
str(Chr <- sapply(strsplit(Labels, "[|]"), `[`, 1))

Descr <- gsub("Conus californicus conotoxin |Conus californicus isolate","",Descr)

f2 <- list.files(dir, pattern = "_californicus.fa", full.names = T)

seqs <- Biostrings::readDNAStringSet(f2)

keep <- Chr %in% names(seqs)

head(WIDTH <- data.frame(Chr = names(seqs), width = Biostrings::width(seqs), sp = sp[keep], Descr = Descr[keep]))

df <- df %>% 
  left_join(WIDTH)


# Order position based on alignment

library(msa)

msa_align <- msa::msa(seqs, method = "ClustalW") # Muscle perform better for 

out <- as_tibble(ggmsa::tidy_msa(DNAStringSet(msa_align))) %>% 
  dplyr::rename("Chr"="name", "pos_alignment"="position") %>%
  filter(character != "-") %>%
  group_by(Chr) %>%
  mutate(ID = paste(Chr, row_number(), sep = "_"))

out %>%
  ggplot(aes(x = pos_alignment, y = Chr)) +
  geom_tile(aes(fill = character)) + # size = 0.1, width = 0.95, color = "black"
  # geom_text(aes(label = character), vjust = 0.5, hjust = 0.5, size= 2.5, family =  "GillSans") +
  theme_classic(base_size = 7, base_family = "GillSans") +
  # scale_fill_manual("", values = c("white", "#cd201f", "#FFFC00","#00b489","#31759b")) +
  theme(
    legend.position = "bottom", 
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_text(hjust = 1),
    axis.ticks.length = unit(5, "pt"))

# Reframe Depth in the alignment

df %>%   
  group_by(Chr) %>%
  mutate(ID = paste(Chr, row_number(), sep = "_")) %>%
  left_join(out, by = "ID") %>%
  ggplot(aes(x = pos_alignment, y = Descr)) +
  # geom_tile(aes(fill = Depth)) + # size = 0.1, width = 0.95, color = "black"
  geom_point(aes(color = Depth), shape = 15) +
  # geom_text(aes(label = character), vjust = 0.5, hjust = 0.5, size= 2.5, family =  "GillSans") +
  theme_classic(base_size = 10, base_family = "GillSans") +
  # scale_fill_manual("", values = c("white", "#cd201f", "#FFFC00","#00b489","#31759b")) +
  scale_fill_viridis_c() + scale_color_viridis_c() +
  scale_x_continuous(position = "top") +
  theme(
    legend.position = "bottom", 
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_text(hjust = 1),
    axis.ticks.length = unit(5, "pt"))


df %>% 
  # filter(Chr %in% QUERY) %>%
  # group_by(Chr) %>% mutate(Depth = Depth/sum(Depth)) %>%
  # group_by(Chr) %>%
  mutate(Position = Position/width) %>%
  ggplot(aes(y = Descr, x = Position, color = Depth))  +
  geom_point(shape = 15)+
  scale_x_continuous(position = "top") +
  scale_fill_viridis_c() + theme_classic(base_family = "GillSans", base_size = 10)


####

# Install the necessary packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("GenomicAlignments", quietly = TRUE)) {
  BiocManager::install("GenomicAlignments")
}
if (!requireNamespace("Rsamtools", quietly = TRUE)) {
  BiocManager::install("Rsamtools")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

# Load the required libraries
library(GenomicAlignments)
library(Rsamtools)
library(ggplot2)

# Function to read BAM file and generate coverage table
generate_coverage_table <- function(bam_file) {
  # Read BAM file
  param <- ScanBamParam(what = c("rname", "pos", "qwidth"))
  bam <- readGAlignments(bam_file, param = param)
  
  # Calculate coverage
  coverage <- IRanges::coverage(bam)
  # coverage <- GenomicRanges::coverage(bam, )
  
  # Create a data frame from coverage
  coverage_df <- as.data.frame(coverage) # as.data.frame(as(coverage, "IRanges"))
  
  # Rename columns for clarity
  colnames(coverage_df) <- c("sequence_name", "start", "end", "coverage")
  
  return(coverage_df)
}

# Function to plot coverage using ggplot2
plot_coverage <- function(coverage_df) {
  ggplot(coverage_df, aes(x = start, y = coverage)) +
    geom_line() +
    labs(title = "Coverage Plot", x = "Position", y = "Coverage") +
    theme_minimal()
}

# Example usage
bam_dir <- "~/Documents/GitHub/conopeptides/07.Reference/S1_BOWTIE2_BAM_FILES_conoserver_nucleic_californicus_DIR/"  # Replace with the path to 

bam_file <- list.files(bam_dir, full.names = T)

# your BAM file
coverage_df <- generate_coverage_table(bam_file[2])
coverage_plot <- plot_coverage(coverage_df)

# Save the plot to a file
ggsave("coverage_plot.png", coverage_plot)

# Display the plot
print(coverage_plot)

