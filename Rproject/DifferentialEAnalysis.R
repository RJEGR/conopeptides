# 
# Run differential expression analysis
# Subset DB to only on conotoxins


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)


library(tidyverse)

pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR/"

DB <- read_rds(paste0(pub_dir, "/structured_db.rds"))

dir <- "/Users/cigom/Documents/GitHub/conopeptides/06.Quantification/MATRIX_RSEM_dir/"

.colData <- list.files(dir, pattern = "Manifest", full.names = T) 

.colData <- read_tsv(.colData) %>% 
  mutate(design = paste0(Diatery, "_", Time)) %>%
  mutate_if(is.character, as.factor)

table(.colData$Diatery, .colData$Time)

f <- list.files(dir, pattern = "Merged_polyA_hisat_SuperDuper_isoforms.filt.rds", full.names = T)

dim(datExpr <- readRDS(f))

# If Filter only conopeptides
# This step is useful to run DE analysis without replicates
# reducing one or more explanatory factors from the linear model matrix

f <- list.files(dir, pattern = "Merged_polyA_hisat_SuperDuper_isoforms.rds", full.names = T)

dim(datExpr <- readRDS(f))

DBSBST <- DB %>% 
  select(gene_id, tab, Signalp_class, PFAMs, Description) %>%
  filter(!is.na(tab)) %>% filter(Signalp_class == "SP")

sum(unique(DBSBST$gene_id) %in% rownames(datExpr))

queryids <- DBSBST %>% distinct(gene_id) %>% pull()
  
keep <- rownames(datExpr) %in% queryids

sum(keep)/length(queryids) # N transcripts

nrow(datExpr <- datExpr[keep,])

datExpr <- round(datExpr)

# multiple contrast Diatery~time ----

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = datExpr,
  colData = .colData,
  design =  ~ design )

dds <- estimateSizeFactors(ddsFullCountTable) 

dds <- estimateDispersions(dds)

dds <- nbinomWaldTest(dds)

keep <- names(rowData(dds))

keep <- grepl(pattern = paste0("^","design"), keep)

contrasts <- names(rowData(dds))[keep]

contrasts <- gsub("_vs_", "_", contrasts)

contrasts <- strsplit(contrasts, "_")

split_res <- function(dds, contrast) {
  
  Level <- paste(contrast, collapse = "_")
  
  cat(Level,"\n")
  
  out <- results(dds, contrast = contrast )
  
  cbind(data.frame(contrast = Level), out)
}

res <- lapply(strsplit(contrasts, "_"), function(x) split_res(dds, contrast = x))

do.call(rbind, res) -> res

# as DESEQ2 works with replicates, lets use edgeR

# Using EDgeR and calculating dispersion value
# Section What to do if you have no replicates:
# Simply pick a reasonable dispersion value, based on your experience with similar data, and use that for exactTest or glmFit. 
# Typical values for the common BCV (square-rootdispersion) for datasets arising from well-controlled experiments are
# 0.4 for human data, 0.1 for data on genetically identical model organisms or 0.01 for technical replicates.

# from https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

# any(DB$gene_id %in% rownames(datExpr))

# DB %>% drop_na(tab) %>% view()

# 2) Run Multiple Contrast Comparison =====

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
      baseMeanB <- rowMeans(DESeq2::counts(dds,normalized=TRUE)[,keepB])
    } else
      baseMeanB <- DESeq2::counts(dds,normalized=TRUE)[,keepB]
    
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


# CONTRAST <- c("Time", "Diatery")

for (j in 1:length(CONTRAST)) {
  
  i <- j
  
  cat("\nRunning ",CONTRAST[i], "\n")
  
  
  out [[i]] <- run_contrast_DE(datExpr, .colData, CONTRAST = CONTRAST[i]) %>% mutate(CONTRAST = CONTRAST[i])
}

do.call(rbind, out) -> RES

write_rds(RES, file = paste0(pub_dir, "/multipleContrastDeseq2.rds"))

RES %>% 
  drop_na(padj) %>%
  filter(padj < 0.05) %>%
  ggplot(aes(padj)) + 
  facet_grid(sampleA ~ CONTRAST) +
  geom_histogram()

# Run design ~ global profile under Diatery or Time -----

run_global_de <- function(datExpr, colData, Design) {
  
  
  names(colData)[names(colData) %in% Design] <- "Design"
  
  colData <- mutate_if(colData, is.character, as.factor)
  
  colData <- colData %>% mutate(Design = relevel(Design, ref = "Ctrl"))
  
  ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData = datExpr,
    colData = colData,
    design =  ~ Design )
  
  dds <- estimateSizeFactors(ddsFullCountTable) 
  
  dds <- estimateDispersions(dds)
  
  dds <- nbinomWaldTest(dds)
  
  # contrast = c('factorName','numeratorLevel','denominatorLevel'),
  
  contrasts <- names(rowData(dds))[grepl(paste0("^",Design), names(rowData(dds)))]
  
  contrasts <- gsub("_vs_", "_", contrasts)
  
  split_res <- function(dds, contrast) {
    
    Level <- paste(contrast, collapse = "_")
    
    cat(Level,"\n")
    
    out <- results(dds, contrast = contrast )
    
    cbind(data.frame(contrast = Level), out)
  }
  
  res <- lapply(strsplit(contrasts, "_"), function(x) split_res(dds, contrast = x))
  
  do.call(rbind, res) -> res
  
  return(res)
}

res_time <- run_global_de(datExpr, .colData, Design = "Time")

which_contrasts <- c("Diatery", "Time")

RES <- lapply(which_contrasts, function(x) run_global_de(datExpr, .colData, Design = x))

RES <- do.call(rbind, RES) %>% as_tibble() 

write_rds(RES, file = paste0(pub_dir, "/Deseq2GlobalContrast.rds"))

res %>%
  as_tibble() %>%
  drop_na(padj) %>%
  filter(padj < 0.05) %>%
  ggplot(aes(padj)) + 
  facet_grid(~ contrast) +
  geom_histogram()

# Run design Time ~ 

# Run design Diatery ~ 