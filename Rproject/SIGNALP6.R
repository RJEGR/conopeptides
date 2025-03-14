# REad signal P and plot patterns per assembly method


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

dir <- "~/Documents/GitHub/conopeptides/05.Prediction/SIGNALP6_OUT/"

f <- list.files(dir, pattern = "txt", full.names = T)

library(tidyverse)

read_predict_res <- function(f) {

  DF <- read_tsv(f, skip = 1) 
  
  keepCols <- names(DF)[1:4]
  
  Method <- gsub("_output_dir_prediction_results.txt", "", basename(f))
  
  DF %>% select(all_of(keepCols)) %>% mutate(Method = Method) %>% dplyr::rename("Signal_peptide_probs" = "SP(Sec/SPI)")
  

}

DF <- lapply(f, read_predict_res)

DF <- do.call(rbind,DF)

DF %>% dplyr::count(Prediction, Method)

DF %>% 
  filter(Signal_peptide_probs > 0.1) %>%
  dplyr::count(Method)

DF %>% ggplot(aes(y = Method, x = Signal_peptide_probs, fill = Method)) + geom_boxplot()
