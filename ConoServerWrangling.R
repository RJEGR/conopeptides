# Conotoxin server processing

dir <- "~/Downloads/"

f1 <- list.files(dir, pattern = "conoserver_protein", full.names = T)
f2 <- list.files(dir, pattern = "conoserver_nucleic", full.names = T)

library(Biostrings)

seq1 <- Biostrings::readAAStringSet(f1)

seq2 <- Biostrings::readDNAStringSet(f2)

sum(keep <- nchar(seq1) > 1) # because Entry found which does not contain a sequence: P05971|Vi001|Conus.

seq1[!keep]

seq1 <- seq1[keep]

# Subset only precursor?

table(sapply(strsplit(names(seq1), "[|]"), `[`, 4))

# Patent Precursor Synthetic Wild type 
# 1292      2930      1243      3057 

sum(keep <- grepl("Precursor", names(seq1))) # Only 2930 Precursor sequences

# The size distribution does not seem different between Precursor and mature
hist(nchar(seq1[!keep]))
hist(nchar(seq1[keep]))

seq1[keep]
seq1[!keep]




# OUt file

outName <- gsub(".fa.gz", "", basename(f1))
outFile <- file.path(dir, paste0(outName, ".curated.fa"))

Biostrings::writeXStringSet(seq1, file = outFile)

library(ShortRead)

ShortRead::writeFasta(seq1, file = outFile)
