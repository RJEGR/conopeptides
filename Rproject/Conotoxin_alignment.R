# transcriptome conopeptides wranglin

dir <- "~/Documents/GitHub/conopeptides/05.Prediction/"


DB <- read_tsv(file.path(dir, "Merged_clusters_conoserver_and_conodictor.tsv"))

DB %>% count(cysteine_framework_conoserver, sort = T)

DB %>% count(gene_superfamily_conoserver, sort = T)

DB %>% filter(is.na(gene_superfamily_conoserver)) %>% view()

DB %>% 
  filter(grepl("M superfamily|O1 superfamily|O3 superfamily", gene_superfamily_conoserver)) %>% 
  mutate(CDS = ifelse(!is.na(orf_sequences), "True", "False")) %>%
  mutate(CDS = factor(CDS, levels = c("True", "False"))) %>%
  count(CDS, gene_superfamily_conoserver, sort = T) %>% view()

  
DB %>% filter(grepl("M superfamily|O1 superfamily|O3 superfamily", gene_superfamily_conoserver)) %>%
  mutate(CDS = ifelse(!is.na(orf_sequences), "True", "False")) %>%
  mutate(CDS = factor(CDS, levels = c("True", "False"))) %>%
  ggplot(aes(y= dna_identity, x = dna_width, col = CDS)) + 
  geom_point(shape = 21) + facet_grid(~ gene_superfamily_conoserver) +
  scale_color_manual(values = c("red", "gray30")) +
  labs(y = "% Identity (Conoserver blastx)", x = "DNA sequences width") +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()) 

seqs <- DB %>% 
  filter(grepl("M superfamily", gene_superfamily_conoserver)) %>%
  pull(dna_sequence, name = transcript_id) 
  
seqs <- Biostrings::DNAStringSet(seqs)

# DB %>% filter(!is.na(orf_sequences)) %>% count(gene_superfamily_conoserver) 

library(msa)

align <- msa::msa(seqs, method = "Muscle")

.align <- msaConvert(align)$seq

DECIPHER::BrowseSeqs(DNAStringSet(.align))


# Fasta de 48 peptidos (O1 superfamily, O3 superfamily, M superfamily)
# remove chr *

seqs <- DB %>% 
  filter(grepl("M superfamily|O1 superfamily|O3 superfamily", gene_superfamily_conoserver)) %>%
  drop_na(orf_sequences) %>%
  mutate(orf_sequences = gsub("[*]$", "", orf_sequences)) %>%
  mutate(transcript_id = paste0(transcript_id, "|", gene_superfamily_conoserver)) %>%
  pull(orf_sequences, name = transcript_id) 

seqs <- Biostrings::AAStringSet(seqs)


Biostrings::writeXStringSet(seqs, file.path(dir, "Merged_clusters_O1_O3_M_superfamily.cds"))
