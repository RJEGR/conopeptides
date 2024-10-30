
library(tidyverse)


dir <- "~/Documents/GitHub/conopeptides/03.Coverage/"


f <- list.files(dir, pattern = "counts.txt", full.names = T) 

Manifest <- list.files(dir, pattern = "Manifest.tsv", full.names = T)

Manifest <- read_tsv(Manifest) %>% distinct()

cols <- read_tsv(f[1], skip = 1) %>% select(contains(".sorted.bam")) %>% names()

readR <- function(f) {
  
  df <- read_tsv(f, skip = 1) 
  
  m <- df %>% select(contains(c("Chr",".sorted.bam"))) # %>% as("matrix")
  
  colNames <- colnames(m)
  
  colNames <- basename(colNames)
  
  colnames(m) <- gsub(".sorted.bam", "", colNames)
  
  data <- m %>% select_if(is.double) %>% as("matrix")
  
  rownames(data) <- m$Chr
  
  # 
  data <- DESeq2::vst(round(data))
  
  return(data)
  
}


DF <- lapply(f, readR)

DF <- do.call(rbind, DF)

DF <- readR(f[1])

m <- DF %>% 
  as_tibble(rownames = "ID") %>%
  filter(ID %in% query) %>%
  select_if(is.double) %>% as("matrix")



z_scores <- function(x) {(x-mean(x))/sd(x)}

datExpr <- t(apply(m, 1, z_scores))

heatmap(datExpr)

PCA <- function(f) {
  
  df <- read_tsv(f, skip = 1) 
  
  m <- df %>% select(contains(".sorted.bam")) %>% as("matrix")
  
  # m <- df %>% select_if(is.double) %>% as("matrix")
  
  colNames <- colnames(m)
  
  colNames <- basename(colNames)
  
  colnames(m) <- gsub(".sorted.bam", "", colNames)
  
  data <- DESeq2::vst(round(m))
  
  PCA <- prcomp(t(data), center = T, scale. = FALSE)
  
  PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], Method = basename(f))
  
  percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
  
  # percentVar <- round(PCA$sdev/sum(PCA$sdev)*100,1)
  
  PCAvar <- data.frame(
    Eigenvalues = PCA$sdev,
    percentVar = percentVar,
    Varcum = cumsum(percentVar),
    Method = basename(f))
  
  return(list(PCAdf, PCAvar))
}

PCAdf <- PCA(f[1])

# percentVar <- PCAdf[[2]]$percentVar


PCAdf <- lapply(f, PCA)

percentVar1 <- paste0(PCAdf[[2]][[2]]$percentVar[1], "% ,",PCAdf[[1]][[2]]$percentVar[1])

percentVar1 <- paste0("PC1, VarExp: ", percentVar1, "%")


percentVar2 <- paste0(PCAdf[[2]][[2]]$percentVar[2], "% ,",PCAdf[[1]][[2]]$percentVar[2])
percentVar2 <- paste0("PC2, VarExp: ", percentVar2, "%")


PCAdf_ <- rbind(PCAdf[[1]][[1]], PCAdf[[2]][[1]])

# rownames(PCAdf_) <- sapply(strsplit(rownames(PCAdf_), "_E"), `[`, 1)


recode_to <- c("Trinity_counts.txt", "Rnaspades_counts.txt")
recode_to <- structure(c("Trinity", "Rnaspades"), names = recode_to)

PCAdf_ %>%
  mutate(Method = dplyr::recode_factor(Method, !!!recode_to)) %>%
  mutate(LIBRARY_ID = rownames(.)) %>% 
  mutate(LIBRARY_ID = sapply(strsplit(LIBRARY_ID, "_L"), `[`, 1)) %>% 
  mutate(Group = sapply(strsplit(LIBRARY_ID, "_E"), `[`, 1)) %>% 
  left_join(Manifest, by = "LIBRARY_ID") %>%
  # mutate(g = substr(LIBRARY_ID, 1,2)) %>%
  ggplot(., aes(PC1, PC2, color = Time)) +
  theme_bw(base_size = 14, base_family = "GillSans") +
  # xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  # ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  xlab(percentVar1) + ylab(percentVar2) +
  # scale_color_grey("") +
  facet_grid(~ Method) +
  # ylim(-100,100) + xlim(-100,100) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  geom_point(size = 7, alpha = 0.7) +
  ggrepel::geom_text_repel( family = "GillSans", mapping = aes(label = Feeding), size = 3.5, color = "black") +
  # geom_text( family = "GillSans", mapping = aes(label = LIBRARY_ID), size = 3.5, color = "black") +
  theme(legend.position = "top") -> p

p

ggsave(p, filename = '03.Quantification.PCA.png', path = dir, width = 5, height = 4, device = png, dpi = 300)


# PCAdf[[2]]
rbind(PCAdf[[2]][[2]],PCAdf[[1]][[2]]) %>%
  mutate(Method = dplyr::recode_factor(Method, !!!recode_to)) %>%
  group_by(Method) %>%
  mutate(Dim = row_number()) %>%
  ggplot(., aes(y = percentVar, x = as.factor(Dim), fill = Method, color = Method)) +
  geom_col(position = position_dodge2()) + # fill = "gray78", color = "gray78"
  geom_line(aes(y = Varcum, group = Method)) +
  geom_point(aes(y = Varcum)) +
  # labs(x = "Component Number", y = "Eigenvalue") +
  labs(x = "Principal component", y = "Fraction variance explained (%)") +
  scale_fill_grey("") +
  scale_color_grey("") +
  guides(color=guide_legend(nrow = 1)) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()) 


# dendogram against first assembly

hc_cor_mat <- function(f) {
  
  df <- read_tsv(f, skip = 1) 
  
  m <- df %>% select(contains(".sorted.bam")) %>% as("matrix")
  
  colNames <- colnames(m)
  
  colNames <- basename(colNames)
  
  colnames(m) <- gsub(".sorted.bam", "", colNames)
  
  DATA <- DESeq2::vst(round(m))
  
  sample_cor = cor(DATA, method='pearson', use='pairwise.complete.obs')
  sample_dist = dist(sample_cor, method='euclidean')
  hc = hclust(sample_dist, method='complete')
  
  return(hc)
}

hc <- lapply(f, hc_cor_mat)

library(dendextend)

# Create two dendrograms
hc[[1]]$labels


hc[[2]]$labels <- gsub("S2_RSEM_CALCULATION_FILES/","",hc[[2]]$labels)

hc[[2]]$labels <- sapply(strsplit(hc[[2]]$labels, "_CK"), `[`, 1)

hc[[1]]$labels <- gsub("S2_RSEM_CALCULATION_FILES/","",hc[[1]]$labels)
hc[[1]]$labels <- sapply(strsplit(hc[[1]]$labels, "_CK"), `[`, 1)


identical(sort(hc[[1]]$labels), sort(hc[[2]]$labels))

dend1 <- as.dendrogram(hc[[1]])
dend2 <- as.dendrogram(hc[[2]])


# Create a list to hold dendrograms
dend_list <- dendlist("TRINITY" = dend1, "RNASPADES" = dend2)


# Align and plot two dendrograms side by side
dend_list %>%
  untangle(method = "step1side") %>% # Find the best alignment layout
  tanglegram(
    highlight_distinct_edges = FALSE,
    common_subtrees_color_lines = TRUE,
    common_subtrees_color_branches = TRUE) 


# Compute alignment quality. Lower value = good alignment quality
dendlist(dend1, dend2) %>%
  untangle(method = "step1side") %>%
  entanglement() # Alignment quality

# Baker’s Gamma Index (see baker’s paper from 1974) is a measure of association (similarity) between two trees of Hierarchical clustering (dendrograms). It is defined as the rank correlation between the stages at which pairs of objects combine in each of the two tree

cor_bakers_gamma(dend1, dend2)

