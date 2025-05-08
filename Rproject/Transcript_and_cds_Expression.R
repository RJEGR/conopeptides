

# Correlate expression and prevalence for both, isoform and cds expression levels
# Include methods RSEM, two pseudoalignments: SALMON and KAllisto tested with cds expression levels: which resulted in best read coverage?
# Calculate PCA or dendogram of samples and compare
# As bowtie/RSEM is not sensitivit to cds, run kallisto/salmon in cds and include here


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR/"

# Library	Raw PE (x2)	Poly-A trim
Total_reads <- c(
  "cam2_E_CKDL230016188-1A_H75JVDSX7"=	14406214,
  "cam2v_E_CKDL230016197-1A_H75JVDSX7" = 9317842,
  "cam4_E_CKDL230016192-1A_H75JVDSX7" = 16625616,
  "cam6_E_CKDL230016196-1A_H75JVDSX7"	= 7369370,
  "ctrl_E_CKDL230016187-1A_H75JVDSX73" =	13505191,
"lit2_E_CKDL230016189-1A_H75JVDSX7_L3"	=  14290690,
"lit4_E_CKDL230016193-1A_H75JVDSX7_L3"	= 18788783,
"mix2_E_CKDL230016191-1A_H75JVDSX7_L3" =	 	12077409,
"mix4_E_CKDL230016195-1A_H75JVDSX7_L3" = 20992813,
"pol2_E_CKDL230016190-1A_H75JVDSX7_L3" = 	15087939,
"pol4_E_CKDL230016194-1A_H75JVDSX7_L3" =	17987538)


combine_matrix <- function(dir) {
  
  # To solve >>> the Reduce section
  
  subdirs <- list.files(dir, pattern  = "transcripts_quant", full.names = T)
  
  subdirs <- list.files(subdirs, pattern  = ".tsv$|quant.sf$", full.names = T)
  
  read_matrix_kallisto_salmon <- function(f) {
    
    OUT <- read_tsv(f)
    
    keep_cols <- names(OUT) %in% c("target_id","Name", "est_counts", "NumReads") # TPM or tpm
    
    keep_cols <- names(OUT)[keep_cols]
    
    OUT <- OUT %>% select(all_of(keep_cols))
    
    sampleName <- sapply(strsplit(dirname(f), "/"), `[`, 10)
    
    sampleName <- gsub("_transcripts_quant","", sampleName)
    
    sampleName <- gsub("_L[0-9]$","", sampleName)
    
    names(OUT) <- c("gene_id", sampleName)
    
    return(OUT)
  }
  
  data_list <- lapply(subdirs, read_matrix_kallisto_salmon)
  
  # Function to bind columns and set row names
  bind_columns_with_rownames <- function(data_list) {
    
    # Function to check if the first column matches across all data frames
    check_first_column_match <- function(data_list) {
      # Extract the first column of all data frames
      first_columns <- lapply(data_list, function(df) df[, 1]$gene_id)
      
      # Check if all first columns are identical
      
      # identical(first_columns[[1]], first_columns[[11]])
      
      all_match <- Reduce(function(x,y) identical(first_columns[[1]], y), first_columns)
      
      return(all_match)
    }
    
    # Apply the function to the list
    result <- check_first_column_match(data_list)
    
    # Print the result
    if (result) {
      print("The first column matches across all data frames.")
    } else {
      print("The first column does not match across all data frames.")
    }
    
    # Extract the first data frame's first column as row names
    rownames_final <- data_list[[1]]$gene_id
    
    # Use cbind to combine all data frames, dropping the first column in each
    combined <- do.call(cbind, lapply(data_list, function(df) df[, -1, drop = FALSE]))
    
    # Set the row names of the final data frame
    rownames(combined) <- rownames_final
    
    return(combined)
  }
  
  # Apply the function to the list OUT
  OUT <- bind_columns_with_rownames(data_list)
  
  return(OUT)
}

# LOAD data -----

DB <- read_rds(paste0(pub_dir, "/structured_db.rds"))

dir <- "/Users/cigom/Documents/GitHub/conopeptides/06.Quantification/MATRIX_RSEM_dir"

subdirs <- list.files(dir, pattern = "SALMON_Merged_polyA_hisat_SuperDuper.fasta.transdecoder_DIR", full.names = T)

salmon_matrix <- combine_matrix(subdirs)

subdirs <- list.files(dir, pattern = "KALLISTO_Merged_polyA_hisat_SuperDuper.fasta.transdecoder_DIR", full.names = T) 

kallisto_matrix <-combine_matrix(subdirs)

file_out <- file.path(subdirs, "KALLISTO_Merged_polyA_hisat_SuperDuper.fasta.transdecoder.matrix")

write_rds(as(kallisto_matrix, "matrix"), file = file_out)

coveragedf_salmon <- colSums(salmon_matrix) %>% as_tibble(rownames = "LIBRARY_ID") %>% dplyr::rename("salmon_cds_level" = "value")

coveragedf_kallisto <- colSums(kallisto_matrix) %>% as_tibble(rownames = "LIBRARY_ID") %>% dplyr::rename("kallisto_cds_level" = "value")

# 2) 

f <- list.files(dir, pattern = "Merged_polyA_hisat_SuperDuper_isoforms.rds", full.names = T)

dim(COUNT <- read_rds(f))

coveragedf_isoforms <- colSums(COUNT) %>% as_tibble(rownames = "LIBRARY_ID") %>% dplyr::rename("rsem_isoform_level" = "value")

summarize_expression <- function(count) {
  
  query <- rownames(count)
  
  expr <- rowSums(count)
  
  prev <- apply(count, 1, function(x) sum(x > 0))
  
  identical(query, names(expr))
  identical(query, names(prev))
  
  head(OUT <- data.frame(query, prev, expr))
  
  as_tibble(OUT)
  
}

summarize_isoform <- summarize_expression(COUNT)

names(summarize_isoform) <- c("gene_id", paste(names(summarize_isoform)[2:3], "_rsem_isoform_level",sep = ""))

salmon_summarize <- summarize_expression(salmon_matrix)

names(salmon_summarize) <- c("protein_id", paste(names(salmon_summarize)[2:3], "_salmon_cds_level",sep = ""))

kallisto_summarize <- summarize_expression(kallisto_matrix)

names(kallisto_summarize) <- c("protein_id", paste(names(kallisto_summarize)[2:3], "_kallisto_cds_level",sep = ""))


# f <- list.files(dir, pattern = "Merged_polyA_hisat_SuperDuper.fasta.transdecoder_genes.filt.rds", full.names = T)

f <- list.files(dir, pattern = "Merged_polyA_hisat_SuperDuper.fasta.transdecoder_genes.rds", full.names = T)

dim(COUNT <- read_rds(f))

coveragedf <- colSums(COUNT) %>% as_tibble(rownames = "LIBRARY_ID") %>% 
  dplyr::rename("rsem_cds_level" = "value") %>%
  left_join(coveragedf_isoforms) %>%
  left_join(coveragedf_salmon) %>%
  left_join(coveragedf_kallisto)




data.frame(coveragedf[1], coveragedf[-1]/Total_reads) %>% 
  pivot_longer(-LIBRARY_ID) %>%
  ggplot() + 
  facet_grid(~ name) + #  scales = "free_x
  geom_col(aes(y = LIBRARY_ID, x = value, fill = name))

coveragedf %>%
  # pivot_longer(-LIBRARY_ID) %>%
  select_if(is.numeric) %>%
  rstatix::cor_mat(method = "spearman") %>%
  rstatix::pull_lower_triangle() %>%
  rstatix::cor_plot(label = TRUE)


summarize_cds <- summarize_expression(COUNT)

names(summarize_cds)<- c("protein_id", paste(names(summarize_cds)[2:3], "_rsem_cds_level",sep = ""))


summarize_cds %>%
  # mutate(protein_id = query) %>%
  mutate(gene_id = gsub(".p[0-9]+$","", protein_id)) %>%
  left_join(summarize_isoform) %>%
  left_join(salmon_summarize) %>%
  left_join(kallisto_summarize) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  select_at(vars(starts_with("expr_"))) %>%
  # select_if(is.numeric) %>%
    # mutate_if(is.character, as.factor) %>%
    # mutate_if(is.factor, as.numeric) %>%
  rstatix::cor_mat(method = "spearman") %>%
  # rstatix::cor_gather()
  # rstatix::cor_reorder() %>%
  rstatix::pull_lower_triangle() %>%
  rstatix::cor_plot(label = TRUE)

df2 %>%
  mutate(protein_id = query) %>%
  mutate(query = gsub(".p[0-9]+$","", query)) %>%
  left_join(df1) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  select_if(is.numeric) %>%
  ggplot(aes(log10(expr_isoform_level), log10(expr_cds_level))) + ggdensity::geom_hdr()


