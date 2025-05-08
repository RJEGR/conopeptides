# To do: 
# Run differential expression analysis w/o replicates
# reducing one or more explanatory factors from the linear model matrixw w/o replicates
# LOAD raw matrix expression
# LOAD Structured DB 
# FILTER expression matrix to conopeptides 'filter(!is.na(tab)) %>% filter(Signalp_class == "SP")'
# From initial raw matrix expression transform to vst
# from the vst transformed raw matrix find housekeeping gene
# calculate dispersion value of housekeeping vector 'y0 <- estimateDisp(y1[housekeeping,], trend="none", tagwise=FALSE)'
# Run EDGE using 'y0$common.dispersion' value and conopeptide gene matrix


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)


library(tidyverse)

pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR/"

# LOAD data -----

DB <- read_rds(paste0(pub_dir, "/structured_db.rds"))

dir <- "/Users/cigom/Documents/GitHub/conopeptides/06.Quantification/MATRIX_RSEM_dir/"

# disp_value <- read_rds(paste0(dir,"/boostrap_dispersion.rds"))

# f <- list.files(dir, pattern = "Merged_polyA_hisat_SuperDuper_isoforms.rds", full.names = T)

subdir <- "KALLISTO_Merged_polyA_hisat_SuperDuper.fasta.transdecoder_DIR/"

f <- "KALLISTO_Merged_polyA_hisat_SuperDuper.fasta.transdecoder.matrix"

f <- list.files(file.path(dir, subdir), f, full.names = T)


dim(datExpr <- readRDS(f))

datExpr <- round(datExpr)

disp_value <- read_rds(paste0(dir,"/cds_kallisto_boostrap_dispersion.rds"))


.colData <- list.files(dir, pattern = "Manifest", full.names = T) 

.colData <- read_tsv(.colData) %>% 
  mutate(Diatery = ifelse(grepl("cam2v_", LIBRARY_ID), "camv", Diatery)) %>%
  mutate(design = ifelse(!Time == "Ctrl", paste0(Diatery, "_", Time), "Ctrl")) %>%
  select(LIBRARY_ID, Time, Diatery, design) %>%
  mutate_if(is.character, as.factor)

table(.colData$Diatery, .colData$Time)

# (omit) Transform raw matrix to VST and save ----

# library(DESeq2)

# ddsFullCountTable <- DESeqDataSetFromMatrix(
#   countData = datExpr,
#   colData = .colData,
#   design = ~ 1 )

# dds <- estimateSizeFactors(ddsFullCountTable) 

# dds <- estimateDispersions(dds)

vst <- DESeq2::vst(datExpr) # vst if cols > 10 and varianceStabilizingTransformation if cols < 10

ntr <- DESeq2::normTransform(dds)

file_out <- paste0(dir, "/counts_vst_nt_raw.rds")
  
write_rds(list("vst" = vst, "normTransform" = ntr, "raw" = datExpr), file = file_out)

raw_df <- vsn::meanSdPlot(datExpr, plot = F)
vst_df <- vsn::meanSdPlot(vst, plot = F)
ntr_df <- vsn::meanSdPlot(assay(ntr), plot = F)

rbind(
  data.frame(py = vst_df$sd, px = vst_df$rank, col = "vst"),
  # data.frame(py = ntr_df$sd, px = ntr_df$rank, col = "ntr"),
  data.frame(py = raw_df$sd, px = raw_df$rank, col = "raw")) %>%
  filter(col != "ntr") %>%
  # ggplot(aes(y = col, x = py, fill = after_stat(x))) +
  # ggridges::geom_density_ridges_gradient(
  #   jittered_points = T,
  #   position = ggridges::position_points_jitter(width = 0.05, height = 0),
  #   point_shape = '|', point_size = 3, point_alpha = 1, alpha = 1) +
  # scale_fill_viridis_c(option = "C")
  ggplot(aes(px, py, color = col)) +
  labs(x = "Ranks", y = "sd", color = "") +
  geom_line(orientation = NA, position = position_identity(), size = 2) +
  theme_bw(base_family = "GillSans", base_size = 20) +
  theme(legend.position = "top")

# creates edgeR object -----

library(edgeR)
# library(limma)

g <- levels(.colData$design)
  
# y <- DGEList(counts=datExpr, group=g)

# design <- model.matrix(~g)


# Calculate dispersion ======

# Tnsert this into the full data object and proceed:
# y$common.dispersion <- mean(disp_value)
# fit <- glmFit(y, design)
# lrt <- glmLRT(fit)
# 
# tTags <- topTags(lrt, n = NULL)
# 
# result_table <- tTags$table

# if bayes

# 
# fit <- lmFit(count, design)
# 
# contrast <- paste(levels(fl), collapse = '-')
# 
# cont.matrix <- makeContrasts(contrasts = contrast, levels = levels(fl))
# 
# fit2 <- contrasts.fit(fit, cont.matrix)
# 
# fit2 <- eBayes(fit2)
# 
# tT <- topTable(fit2, adjust="fdr", sort.by="B", 
#   number = Inf) %>% 
#   as_tibble(rownames = "ids") %>%
#   mutate_at(vars(!matches("ids|P.Value|adj.P.Val")), 
#     round, digits = 2)



# Run DE analysis ------

create_pairs <- function(vec) {
  combn(vec, 2, simplify = FALSE)
}


combinations_vector <- unlist(lapply( create_pairs(g), paste, collapse = "-"))

# Print the combinations vector
print(combinations_vector)

lapply(combinations_vector, function(x) unlist(strsplit(x, "-")))

# herence is where we apply the loop [] or lapply
i <- 1

# gstr <- unlist(strsplit(combinations_vector[i], "-"))
# 
# count <- datExpr
# 
# colData <- .colData

# continue here

run_edgeR <- function(count, gstr, colData, disp_value = NULL) {
  
  sam_gstr <- structure(levels(colData$design), names = levels(colData$LIBRARY_ID))
  
  sam_gstr <- sam_gstr[sam_gstr %in% gstr]
  
  fl <- as.factor(sam_gstr)
  
  design <- model.matrix(~fl)
  
  # filter samples from the contrast
  
  keep_cols <- colnames(count) %in% names(sam_gstr)

  # count <- round(count[,keep_cols])
  
  by_count <- 1; by_freq <- 1
  
  keep_genes <- rowSums(count > by_count) >= by_freq
  
  sum(keep_genes)
  
  dim(count <- count[keep_genes,keep_cols ])
  
  colnames(count) <- sam_gstr
  
  count <- as(count, 'matrix')
  
  count <- round(count)

  # Performing EdgerR
  
  DGE <- DGEList(counts = count, group = fl)
  
  DGE$common.dispersion <- disp_value[1]
  
  # fit <- glmFit(DGE, design)
  # 
  # lrt <- glmLRT(fit)
  # 
  # tTags <- topTags(lrt, n = NULL)
  
  
  # After adding disp_value run as below instead of glm
  DGE = estimateTagwiseDisp(DGE)
  
  exact <- exactTest(DGE, pair = levels(fl), dispersion = disp_value[1])
  
  tTags <- topTags(exact, n = NULL)

  result_table <- tTags$table
  
  sampleA <- levels(fl)[1]
  sampleB <- levels(fl)[2] 
  
  result_table <- data.frame(sampleA,
    sampleB, 
    result_table)
  
  reorder_cols <- c( "ids" ,"sampleA", "sampleB", "logFC",  "PValue", "FDR", "logCPM")
  
  result_table %>%
    as_tibble(rownames = 'ids') %>%
    mutate(logFC = -1 * logFC) %>% # reset logfc so it's A/B instead of B/A to be consistent with DESeq2
    mutate_at(vars(!matches("ids|sample|PValue|FDR")),
      round ,digits = 2) %>%
    select_at(vars(all_of(reorder_cols)))
  
}


# DF <- run_edgeR(datExpr, gstr, colData, disp_value = mean(disp_value))

comb_list <- strsplit(combinations_vector, "-")

OUT <- lapply(comb_list, 
  function(x) run_edgeR(datExpr, gstr = x, .colData, disp_value = mean(disp_value)))

OUT <- do.call(rbind, OUT)

OUT <- OUT %>% filter(PValue < 0.05)

# write_rds(OUT, file = paste0(dir, "/glmLRT_multiple_contrast_ctrl_and_treatments.rds"))
write_rds(OUT, file = paste0(file.path(dir, subdir), "/p05_exactTest_multiple_contrast.rds"))

# pre 

