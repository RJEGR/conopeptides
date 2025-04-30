
# Ricardo Gomez-Reyes
# To do:
# READ COUNT MATRIX
# READ Manifest (metadata)
# Transform Count Matrix to z-score
# Previs w/ heatmap
# CONVERt wide to longer data.frame
# Plot ggplot using geom_tile or geom_raster

# Number of genes and samples
num_genes <- 100
num_samples <- 10

# Create a matrix of random gene expression values
# Values are randomly drawn from a normal distribution
gene_expression_matrix <- matrix(
  data = rnorm(num_genes * num_samples, mean = 10, sd = 2), # mean=10, sd=2
  nrow = num_genes,
  ncol = num_samples
)

# Add row and column names
rownames(gene_expression_matrix) <- paste0("Gene_", 1:num_genes)
colnames(gene_expression_matrix) <- paste0("Sample_", 1:num_samples)

# View the first few rows of the matrix
print(head(gene_expression_matrix))

# filter count matrix
query_ids <- sample(rownames(gene_expression_matrix), 20) # replace with true set

keep <- rownames(gene_expression_matrix) %in% query_ids

gene_expression_matrix <- gene_expression_matrix[keep,]

# Number of samples
num_samples <- 10

# Create three factors:
# 1. Experimental Group (e.g., "Control" or "Treatment")
# 2. Batch (e.g., "Batch_1", "Batch_2")
# 3. Time Point (e.g., "T1", "T2", "T3")

# Factor 1: Experimental Group
experimental_group <- factor(rep(c("Control", "Treatment"), each = num_samples / 2))

# Factor 2: Batch
batch <- factor(rep(c("Batch_1", "Batch_2"), times = num_samples / 2))

# Factor 3: Time Point
time_point <- factor(rep(c("T1", "T2", "T3", "T1"), length.out = num_samples))

# Combine factors into a data frame for better organization
sample_factors <- data.frame(
  Sample = paste0("Sample_", 1:num_samples),
  Experimental_Group = experimental_group,
  Batch = batch,
  Time_Point = time_point
)

# View the factors
print(sample_factors)


z_scores <- function(x) {(x-mean(x))/sd(x)}

gene_expression_matrix <- apply(gene_expression_matrix, 1, z_scores)



h <- heatmap(t(gene_expression_matrix), , keep.dendro = T)

hc_samples <- as.hclust(h$Colv)
plot(hc_samples)
hc_order <- hc_samples$labels[h$colInd]

hc_genes <- as.hclust(h$Rowv)
plot(hc_genes)
order_genes <- hc_genes$labels[h$rowInd]
plot(hc_genes)


#

LongerDF <- gene_expression_matrix %>% 
  as_tibble(rownames = 'Sample') %>%
  pivot_longer(cols = colnames(gene_expression_matrix), values_to = 'fill', names_to = "Genes") %>%
  left_join(sample_factors, by = "Sample") %>% 
  mutate(Genes = factor(Genes, levels = order_genes)) %>%
  mutate(Sample = factor(Sample, levels = hc_order)) 



lo = floor(min(LongerDF$fill))
up = ceiling(max(LongerDF$fill))
mid = (lo + up)/2

LongerDF %>%
  ggplot(aes(x = Sample, y = Genes, fill = fill)) +  
  geom_tile(color = 'white', linewidth = 0.7, width = 1) +
  facet_grid(~ Experimental_Group+Batch+Time_Point, scales = "free_x") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
    na.value = "white", midpoint = mid, limit = c(lo, up),
    name = NULL)

