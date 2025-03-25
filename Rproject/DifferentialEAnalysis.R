# 
# Run differential expression analysis
# Focus only on conotoxins!


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)


library(tidyverse)

pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR/"

DB <- read_rds(paste0(pub_dir, "/structured_db.rds"))

# gene2GO <- read_rds(paste0(pub_dir, "/gene2GO.rds"))

dir <- "/Users/cigom/Documents/GitHub/conopeptides/06.Quantification/MATRIX_RSEM_dir/"

.colData <- list.files(dir, pattern = "Manifest", full.names = T) 

.colData <- read_tsv(.colData) %>% 
  mutate(datTraits = paste0(Diatery, "-", Time))

table(.colData$LIBRARY_ID, .colData$Time)

f <- list.files(dir, pattern = "Merged_polyA_hisat_SuperDuper_isoforms.filt.rds", full.names = T)

dim(datExpr <- readRDS(f))

# x) Run Multiple Contrast Comparison =====

library(DESeq2)

CONTRAST <- .colData %>% dplyr::select(starts_with("CONTRAST")) %>% names()

run_contrast_DE <- function(COUNTS, colData, CONTRAST = NULL, ref = NULL) {
  
  get_res <- function(dds, contrast, alpha_cutoff = 0.1) {
    
    sA <- contrast[1]
    sB <- contrast[2]
    
    contrast <- as.character(DESeq2::design(dds))[2]
    
    keepA <- as.data.frame(colData(dds))[,contrast] == sA
    keepB <- as.data.frame(colData(dds))[,contrast] == sB
    
    contrast <- c(contrast, sA, sB)
    
    res = results(dds, contrast, alpha = alpha_cutoff)
    

    if(!sum(keepA) == 1) {
      baseMeanA <- rowMeans(DESeq2::counts(dds,normalized=TRUE)[,keepA])
    } else
      baseMeanA <- DESeq2::counts(dds,normalized=TRUE)[,keepA]
    
    if(!sum(keepB) == 1) {
      baseMeanA <- rowMeans(DESeq2::counts(dds,normalized=TRUE)[,keepB])
    } else
      baseMeanA <- DESeq2::counts(dds,normalized=TRUE)[,keepB]
    
    # baseMeanA <- rowMeans(DESeq2::counts(dds,normalized=TRUE)[,keepA])
    # baseMeanB <- rowMeans(DESeq2::counts(dds,normalized=TRUE)[,keepB])
    
    res %>%
      as.data.frame(.) %>%
      cbind(baseMeanA, baseMeanB, .) %>%
      cbind(sampleA = sA, sampleB = sB, .) %>%
      as_tibble(rownames = "Name") %>%
      mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
      mutate_at(vars(!matches("Name|sample|pvalue|padj")),
        round ,digits = 2)
  }
  
  # CONTRAST: Column in colData with character vector of Design
  
  names(colData)[1] <- "LIBRARY_ID"
  
  names(colData)[names(colData) %in% CONTRAST] <- "Design"
  
  colData <- colData %>% drop_na(Design)
  
  # any(colnames(COUNTS) %in% colData$LIBRARY_ID) # sanity check
  
  colData <- mutate_if(colData, is.character, as.factor)
  
  keep <- colnames(COUNTS) %in% colData$LIBRARY_ID 
  
  COUNTS <- COUNTS[,keep]
  
  colData <- colData %>% mutate(Design = relevel(Design, ref = "Control"))
  
  require(DESeq2)
  
  ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData = COUNTS,
    colData = colData,
    design = ~ Design )
  
  dds <- estimateSizeFactors(ddsFullCountTable) 
  
  dds <- estimateDispersions(dds)
  
  dds <- nbinomWaldTest(dds)
  
  # return(dds)
  
  contrast <- levels(colData(dds)$Design)
  
  res <- get_res(dds, contrast)
  
  return(res)
}

# run_contrast_DE(datExpr, .colData, CONTRAST = CONTRAST[1])

any(colnames(datExpr) == .colData$LIBRARY_ID)

out <- list()

for (j in 1:length(CONTRAST)) {
  
  i <- j
  
  cat("\nRunning ",CONTRAST[i], "\n")
  
  
  out [[i]] <- run_contrast_DE(datExpr, .colData, CONTRAST = CONTRAST[i]) %>% mutate(CONTRAST = CONTRAST[i])
}

do.call(rbind, out) -> RES

# Run design Time ~ Diatery 

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = datExpr,
  colData = .colData,
  design = Time ~ Diatery )

dds <- estimateSizeFactors(ddsFullCountTable) 

dds <- estimateDispersions(dds)

dds <- nbinomWaldTest(dds)

# return(dds)

contrast <- levels(colData(dds)$Design)

res <- get_res(dds, contrast)