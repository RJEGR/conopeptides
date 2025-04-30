# REad signal P and plot patterns per assembly method
# This code compare performance of assembly methods to find signalP
# Instead of regex_pHMM_prediction_results derived from DNA analysis, lets use from predicted (transdecoder) peptide results
# Waiting for signalp6_Merged_polyA_hisat_SuperDuper.fasta.transdecoder.pep

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

dir <- "~/Documents/GitHub/conopeptides/05.Prediction/SIGNALP6_OUT/SIGNALP_OUT_DIR/"


f <- list.files(dir, pattern = "txt", full.names = T)

library(tidyverse)

read_signalp_res <- function(f) {

  DF <- read_tsv(f, skip = 1) 
  
  keepCols <- names(DF)[1:4]
  
  Method <- gsub("_output_dir_prediction_results.txt", "", basename(f))
  
  DF %>% 
    select(all_of(keepCols)) %>% mutate(Method = Method) %>% 
    dplyr::rename("Signal_peptide_probs" = "SP(Sec/SPI)", "protein_id" = "# ID") %>%
    dplyr::rename("Signalp_class" = "Prediction", "Other_probs" = "OTHER")
  

}


DF <- lapply(f, read_signalp_res)

DF <- do.call(rbind,DF)

DF %>% distinct(protein_id)

DF <- DF %>% mutate(protein_id = sapply(strsplit(protein_id, " "), `[`, 1))

DF %>% dplyr::count(Signalp_class, Method)

DF %>% 
  filter(Signal_peptide_probs > 0.1) %>%
  dplyr::count(Method)

DF %>% 
  filter(Signalp_class == "SP") %>%
  ggplot(aes(Signal_peptide_probs, color = Method, fill = Method)) + 
  # geom_density()
  stat_ecdf(linewidth = 1, alpha = 0.5)


# save 
pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR"


write_rds(DF, file = paste0(pub_dir, "/signalp6_Merged_polyA_hisat_SuperDuper_transdecoder.rds"))
