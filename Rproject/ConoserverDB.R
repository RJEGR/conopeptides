# Parse Conoserver db
# Note: the format of the headers in the FASTA file is:
  
  # for proteins: conoserver identifier|name|organism|protein type|toxin class|gene superfamily|cysteine framework|pharmacological family|evidence
# for nucleic acids: conoserver identifier|name|organism

# But also grep https://www.conoserver.org/?page=about_conotoxins&bpage=cononames

library(Biostrings)
library(tidyverse)

dir <- "~/Documents/GitHub/conopeptides/"

f <- list.files(dir, pattern = "conoserver_protein.fa", full.names = T)

seqs <- Biostrings::readAAStringSet(f[1])

DF <- as_tibble(data.frame(seqs, row.names = names(seqs)), rownames = "ID")

cols <- "identifier|name|organism|protein type|toxin class|gene superfamily|cysteine framework|pharmacological family|evidence"

cols <- strsplit(cols,"[|]")[[1]]

DF <- DF %>% separate(ID, into = cols, sep = "[|]")

DF %>% count(name)
DF %>% count(`gene superfamily`)

DF %>% count(`cysteine framework`)
DF %>% count(`evidence`)

write_tsv(DF, file = file.path(dir, "conoserver_protein.tsv"))

# library(xml2)
# 
# xml2::rad
# read_xml(f[2])
