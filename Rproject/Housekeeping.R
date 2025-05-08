# Find stable gene expression
# according to Eisenberg and Levanon

# Load necessary libraries

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

set.seed(280325)

library(tidyverse)

dir <- "/Users/cigom/Documents/GitHub/conopeptides/06.Quantification/MATRIX_RSEM_dir/"


# LOAD And transform

subdir <- "KALLISTO_Merged_polyA_hisat_SuperDuper.fasta.transdecoder_DIR/"

f <- "KALLISTO_Merged_polyA_hisat_SuperDuper.fasta.transdecoder.matrix"

f <- list.files(file.path(dir, subdir), f, full.names = T)

datExpr <- round( read_rds(f))

vst <- DESeq2::vst(datExpr) # vst if cols > 10 and varianceStabilizingTransformation if cols < 10


# Or read
rds <- read_rds(file_out <- paste0(dir, "/counts_vst_nt_raw.rds"))

vst <- rds$vst

datExpr <- rds$raw

# Function to filter genes based on the given criteria
filter_genes <- function(expression_matrix, expr_lev = 5) {
  
  # Calculating high expression leves
  
  # z_scores <- function(x) {(x-mean(x))/sd(x)}
  
  # str(gene_z <- apply(expression_matrix, 1, z_scores))

  # Calculate mean and standard deviation for each gene
  
  gene_means <- rowMeans(expression_matrix)
  
  gene_sds <- apply(expression_matrix, 1, sd)
  
  # Criterion I: Expression observed in all samples
  
  prevalence <- apply(expression_matrix, 1, function(x) sum(x > 0))
  
  # table(prevalence)
  
  keep <- prevalence ==  ncol(expression_matrix)
  
  nrow(gene_expression <- expression_matrix[keep,])
  
  # Criterion II: Low variance over tissues or developmental stages
  low_variance_genes <- gene_sds < 1
  
  # Criterion III: No exceptional expression in any samples (this is equivalent to z-score)
  no_exceptional_expression <- apply(expression_matrix, 1, function(x) {all(abs(x - mean(x)) < 2)})
  
  # Criterion IV: Medium to high expression level
  medium_high_expression <- gene_means > expr_lev
  
  
  # Combine all criteria
  filtered_genes <- low_variance_genes & no_exceptional_expression & medium_high_expression
  
  # Return the filtered gene expression matrix
  return(expression_matrix[filtered_genes, ])
  
}

# using vst to treat with parametric features as mean, sd, 

filtered_gene_expression <- filter_genes(vst) 

cv_values <- apply(filtered_gene_expression, 1, function(x) {sd(x) / mean(x)})

# input both, filtered and raw count to calculte common dispersion

boostrap_common_dispersion <- function(raw_count, count_matrix, sample_size = 0.3, cv_thres = 0.1) {
  
  sample_size <- round(nrow(count_matrix) * sample_size) # using X % of the size
  
  resampled_indices <- sample(c(0:nrow(count_matrix)), sample_size, replace = T)
  
  resampled_matrix <- count_matrix[resampled_indices, ]
  
  # Filter matrix and CV values
  
  cv_values <- apply(resampled_matrix, 1, function(x) {sd(x) / mean(x)})
  
  plot(hist(cv_values))
  
  housekeeping <- names(cv_values[cv_values < cv_thres])
  
  # str(housekeeping <- rownames(resampled_matrix))
  
  require(edgeR)
  
  count_matrix <- round(count_matrix)
  
  g <- colnames(count_matrix)
  
  y <- DGEList(counts=raw_count, group=g)
  
  y1 <- y
  
  y1$samples$group <- 1
  
  y0 <- estimateDisp(y1[housekeeping,], trend="none", tagwise=FALSE)
  
  common_dispersion <- y0$common.dispersion
  
  return(common_dispersion)
  
  
}


# Boostraping 1000 rounds the calculation of dispersion using 30% of the stable genes per resampling

boostrap_dispersion <- replicate(1000, boostrap_common_dispersion(raw_count = datExpr, count_matrix = filtered_gene_expression ))

# if(!is.null(dev.list())) dev.off()

# write_rds(boostrap_dispersion, file = paste0(dir,"/boostrap_dispersion.rds"))

write_rds(boostrap_dispersion, file = paste0(dir,"/cds_kallisto_boostrap_dispersion.rds"))

hist(boostrap_dispersion)

mean(boostrap_dispersion)


# screen annotation data

data.frame(stable_genes) %>% 
  as_tibble(rownames = "gene_id") %>%
  left_join(DB) %>% distinct(Description) %>% view()

# 

# Test
# Assume 'gene_expression' is your matrix of log2(TPM) values with rows as genes and columns as tissue types or developmental stages
# For demonstration, let's create a sample gene_expression matrix
set.seed(123)

dim(gene_expression <- matrix(runif(1000, 0, 10), nrow = 100, ncol = 10))

create_vector <- function(n, mean, sd) {
  # Generate a vector of values from a normal distribution
  values <- rnorm(n, mean = mean, sd = sd)
  return(values)
}

# Example usage
set.seed(123) # For reproducibility
n <- 10 # Number of values
mean_value <- 5 # Desired mean
sd_value <- 1 # Desired standard deviation

# Create the vector
dim(true_vector_values <- replicate(12, create_vector(n, mean_value, sd_value), simplify = F))

dim(outliers_vector_values <- replicate(12, create_vector(n, 5, 2),  simplify = F))

gene_expression <- rbind(gene_expression, do.call(rbind, true_vector_values), do.call(rbind, outliers_vector_values))

dim(gene_expression)

rownames(gene_expression) <- paste0("Gene", 1:nrow(gene_expression))
colnames(gene_expression) <- paste0("Tissue", 1:10)

