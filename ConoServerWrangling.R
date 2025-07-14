# Conotoxin server processing
# 
rm(list = ls())

if(!is.null(dev.list())) dev.off()

dir <- "~/Downloads/CONOSERVERDB/"

f1 <- list.files(dir, pattern = "conoserver_protein.fa$", full.names = T)

library(Biostrings)
library(tidyverse)

seq1 <- Biostrings::readAAStringSet(f1)

sum(keep <- nchar(seq1) > 1) # because Entry found which does not contain a sequence: P05971|Vi001|Conus.

seq1[!keep]

seq1 <- seq1[keep]

length(seq1) # 8522


# Deduplicate sequences

length(unique(as.character(seq1)))/length(seq1) # 7298 ( 0.8563717)

str(Names <- sapply(strsplit(names(seq1), "[|]"), `[`, 1))

# Names are unique? true (ie. no duplicates)
sum(duplicated(Names))

# write a apply to save as vector the unique identifiers

# SEQS <- sort(as.character(seq1)) # <-- this would be errors in next steps, so omit

# names(SEQS) <- Names

SEQS <- structure(as.character(seq1), names = names(seq1))

str(names(SEQS) <- sapply(strsplit(names(SEQS), "[|]"), `[`, 1))

# Sanity check
# visualice vector position 1 matches between SEQS names and sequence in seq1
SEQS[1]
seq1[1]

# Identify duplicated elements
SEQS <- sort(SEQS)

is_duplicated <- duplicated(SEQS) | duplicated(SEQS, fromLast = TRUE)

# Mark duplicated elements
marked_vec <- ifelse(is_duplicated, "duplicated seqs", "unigenes")

# Print marked vector
table(marked_vec) # must match  1224 duplicated records removed w/ seqkit

# Split duplicated elements into groups
duplicated_groups <- split(names(SEQS)[is_duplicated], SEQS[is_duplicated])

# Split the string by comma

head(duplicated_groups <- lapply(duplicated_groups, function(x) paste0(x, collapse = "|")))

head(duplicated_groups <- unlist(duplicated_groups))

duplicated_seqs <- names(duplicated_groups)

duplicated_ids <- paste0("dedup_seq_",seq(1,length(duplicated_groups)))

# duplicated_ids <- duplicated_groups

UNIGENES <- structure(duplicated_seqs, names = duplicated_ids)

head(UNIGENES <- c(SEQS[!is_duplicated], UNIGENES))

length(UNIGENES) # Close to match grep -c "^>" conoserver_protein_dedup.fa (7299)

length(UNIGENES) == length(unique(as.character(seq1))) # Must be true

# OUt file

outName <- gsub(".fa*", "", basename(f1))

outFile <- file.path(dir, paste0(outName, "_curated.fa"))

Biostrings::writeXStringSet(Biostrings::AAStringSet(UNIGENES), file = outFile)

# Filter only precursors ----

# Subset only precursor?

table(sapply(strsplit(names(seq1), "[|]"), `[`, 4))

# Patent Precursor Synthetic Wild type 
# 1292      2930      1243      3057 

sum(keep <- grepl("Precursor", names(seq1))) # Only 2930 Precursor sequences

# The size distribution does not seem different between Precursor and mature
# hist(nchar(seq1[!keep]))
hist(nchar(seq1[keep]))

# seq1[keep]

# seq1[!keep]

str(SEQS <- unique(sort(as.character(seq1[keep])))) # 2824 unique precursor sequences

sum(keep <- UNIGENES %in% unique(SEQS))

outName <- gsub(".fa*", "", basename(f1))

outFile <- file.path(dir, paste0(outName, "_precursor_curated.fa"))

Biostrings::writeXStringSet(Biostrings::AAStringSet(UNIGENES[keep]), 
  file = outFile)


# Write data.frame of duplicates 

head(data.frame(duplicated_groups))

head(data.frame(structure(duplicated_seqs, names = duplicated_ids)))

head(data.frame(SEQS[!is_duplicated]))

# library(ShortRead)
# 
# ShortRead::writeFasta(seq1, file = outFile)

# Nucleotide level ----

f2 <- list.files(dir, pattern = "conoserver_nucleic.fa$", full.names = T)

seqs <- Biostrings::readDNAStringSet(f2)

# Count size 
mean(nchar(seqs))
sd(nchar(seqs))
min(nchar(seqs))
max(nchar(seqs))

# Not work if there is not duplicates

DEDUP <- function(seqs) {
  
  # L <- length(unique(seqs))
  
  # if(L != length(seqs))
    # if not true, omit the function

  Names <- sapply(strsplit(names(seqs), "[|]"), `[`, 1)
  
  names(seqs) <- Names
  
  # SEQS <- sort(as.character(seqs))
  
  
  # Identify duplicated elements
  is_duplicated <- duplicated(SEQS) | duplicated(SEQS, fromLast = TRUE)
  
  # Mark duplicated elements
  marked_vec <- ifelse(is_duplicated, "duplicated seqs", "unigenes")
  
  # Split duplicated elements into groups
  duplicated_groups <- split(names(SEQS)[is_duplicated], SEQS[is_duplicated])
  
  if(!is.null(duplicated_groups)) {
    
    # Split the string by comma
    
    head(duplicated_groups <- lapply(duplicated_groups, function(x) paste0(x, collapse = ";")))
    
    head(duplicated_groups <- unlist(duplicated_groups))
    
    duplicated_seqs <- names(duplicated_groups)
    
    duplicated_ids <- paste0("dedup_seq_",seq(1,length(duplicated_groups)))
    
    UNIGENES <- structure(duplicated_seqs, names = duplicated_ids)
    
    head(UNIGENES <- c(SEQS[!is_duplicated], UNIGENES))
    
    length(UNIGENES)
    
    length(UNIGENES) == length(unique(as.character(seqs))) # Must be true
    
    UNIGENES <- Biostrings::DNAStringSet(UNIGENES)
    
  } else
    
    UNIGENES <- Biostrings::DNAStringSet(SEQS)
  
  return(UNIGENES)
 
}

table(sapply(strsplit(names(seqs), "[|]"), `[`, 3))


# Generate randomize sequences from size 100, keep only sequences > 150 =----

sum(keep <- nchar(seqs) >= 150)

seqssized <- seqs[keep]


table(sapply(strsplit(names(seqssized), "[|]"), `[`, 3))

str(Names <- sapply(strsplit(names(seqssized), "[|]"), `[`, 1))

names(seqssized) <- Names

# and dedup
# seqssized <- DEDUP(seqssized)

outName <- gsub(".fa", "", basename(f2))

outFile <- file.path(dir, paste0(outName, "_seq_length_150.fa"))

Biostrings::writeXStringSet(seqssized, file = outFile)


# only californiconus ====

keep <- grepl('Conus californicus', names(seqs))

seqs <- seqs[keep]

sum(keep <- nchar(seqs) >= 100)

seqs <- seqs[keep]

# seqs2 <- seqs[!keep]

str(Names <- sapply(strsplit(names(seqs), "[|]"), `[`, 1))

names(seqs) <- Names

unique(seqs)

outName <- gsub(".fa", "", basename(f2))

outFile <- file.path(dir, paste0(outName, "_californicus_length_100.fa"))

Biostrings::writeXStringSet(seqs, file = outFile)

table(sapply(strsplit(names(seqs), "[|]"), `[`, 2))



# MSA ====

library(msa)

msa_align <- msa::msa(seqs, method = "Muscle") # Muscle perform better for 

out <- as_tibble(ggmsa::tidy_msa(DNAStringSet(msa_align)))

.align <- msaConvert(msa_align)$seq

# .align <- gsub("U","T", msaConvert(msa_align)$seq)

library(bioseq)

bioseq_set <- bioseq::dna(.align)

length(table(bioseq::seq_cluster(bioseq_set, threshold = 0.05)))

bioseq::seq_consensus(bioseq_set)

# table(SEQ_CLUSTERS <- bioseq::seq_cluster(bioseq_set, threshold = 0.95, method = "single"))


#

out <- ggmsa::tidy_msa(alignment) %>% as_tibble() %>% dplyr::rename("Name"="name") 
# separate(name, into = c("Name", "KnownRNA"), sep = " ")

# see::social_colors()
# scales::show_col(see::social_colors())

ps <- out %>%
  ggplot(aes(x = position, y = name)) +
  geom_tile(aes(fill = character)) + # size = 0.1, width = 0.95, color = "black"
  # geom_text(aes(label = character), vjust = 0.5, hjust = 0.5, size= 2.5, family =  "GillSans") +
  theme_classic(base_size = 7, base_family = "GillSans") +
  # scale_fill_manual("", values = c("white", "#cd201f", "#FFFC00","#00b489","#31759b")) +
  theme(
    legend.position = "none", 
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_text(hjust = 1),
    axis.ticks.length = unit(5, "pt"))





library(ggseqlogo)

msa_align <- msa::msa(seqs, method = "ClustalW") 

.align <- msaConvert(msa_align)$seq


ggseqlogo(.align, seq_type = "dna",  method = "probability")

LOGO <- geom_logo(.align, p = F, seq_type = "dna", method = "probability")

lo = floor(min(LOGO$y))
up = ceiling(max(LOGO$y))
mid = (lo + up)/2

# LOGO %>% group_by(position, order) %>% summarise(sum(y))

LOGO %>% 
  # mutate()
  ggplot(aes(y = y, x = position, fill = "group", group = "group_by")) + 
  geom_col() 
  # facet_grid(seq_group ~.) +
  # scale_fill_gradient2(
  # low = "white", high = "gray5", mid = "gray50",
  # # low = "blue", high = "red", mid = "white", 
  # na.value = "white", midpoint = mid, limit = c(lo, up),
  # name = NULL) 


# Exit

names(seqs[width(seqs) > 200])

str(Names <- sapply(strsplit(names(seqs), "[|]"), `[`, 2))
table(Names)



library(tidyverse)

dir <- "~/Downloads/"
# f <- list.files(path = dir, pattern = "PRJNA526781 _SraRunInfo", full.names = T)

f <- list.files(path = dir, 
  pattern = "SupplementaryTablesv5-forsubmission.csv", full.names = T)

df <- read_csv(f[1])

df %>% dplyr::count(`WoRMS genus`) %>% 
  ggplot(aes(y = `WoRMS genus`, x = n)) +
  geom_col()
