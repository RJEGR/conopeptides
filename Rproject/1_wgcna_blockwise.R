# WGCNA


library(WGCNA)
library(flashClust)
library(tidyverse)


rm(list = ls());

if(!is.null(dev.list())) dev.off()

dir <- "/Users/cigom/Documents/GitHub/conopeptides/06.Quantification/MATRIX_RSEM_dir"

f <- list.files(dir, pattern = "Merged_polyA_hisat_SuperDuper_isoforms.filt.rds", full.names = T)

datExpr <- read_rds(f)

colData <- data.frame(LIBRARY_ID = factor(colnames(datExpr)))

library(DESeq2)

dds <- DESeqDataSetFromMatrix(datExpr,
  colData,
  design = ~ LIBRARY_ID)

dds <- DESeq2::varianceStabilizingTransformation(dds)

datExpr <- assay(dds)

datExpr <- t(datExpr) # log2(count+1) #

str(datExpr)

cat("\n:::::\n")

gsg = goodSamplesGenes(datExpr, verbose = 3)

gsg$allOK

if (!gsg$allOK) {
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
  datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
}

file_out <- gsub(".rds", ".wgcna.input.rds",f)

# write_rds(datExpr, file = file_out)

dim(datExpr <- read_rds(file_out))

max_power <- 30

powers = c(c(1:10), seq(from = 10, to = max_power, by=1)) 

allowWGCNAThreads()

cor_method =  "cor" # by default WGCNA::cor(method =  'pearson') is used, "bicor"

corOptionsList = list(use ='p') # maxPOutliers = 0.05, blocksize = 20000

sft <- pickSoftThreshold(datExpr, 
  powerVector = powers, 
  corFnc = cor_method,
  corOptions = corOptionsList,
  verbose = 5, 
  networkType = "unsigned")

file_out <- gsub(".rds", paste0(dir, '/wgcna_Soft.',cor_method, '.rds'),f)

# readr::write_rds(sft, file = file_out)

# Continue here tomorrow

sft <- readr::read_rds(file_out)

soft_values <- abs(sign(sft$fitIndices[,3])*sft$fitIndices[,2])

soft_values <- round(soft_values, digits = 2)

power_pct <- soft_values[which.max(soft_values)]

softPower <- sft$fitIndices[,1][which(soft_values >= power_pct)]

meanK <- sft$fitIndices[softPower,5]

hist(sft$fitIndices[,5])

softPower <- min(softPower)

cat("\nsoftPower value", softPower, '\n')


title1 = 'Scale Free Topology Model Fit,signed R^2'
title2 = 'Mean Connectivity'

caption = paste0("Lowest power for which the scale free topology index reaches the ", power_pct*100, " %")

library(tidyverse)

sft$fitIndices %>% 
  mutate(scale = -sign(slope)*SFT.R.sq) %>%
  select(Power, mean.k., scale) %>% pivot_longer(-Power) %>%
  mutate(name = ifelse(name %in% 'scale', title1, title2)) %>%
  ggplot(aes(y = Power, x = value)) +
  facet_grid(~name, scales = 'free_x', switch = "x") +
  geom_text(aes(label = Power), size = 3.5, family = "GillSans") +
  geom_abline(slope = 0, intercept = softPower, linetype="dashed", alpha=0.5) +
  # geom_vline(xintercept = min(meanK), linetype="dashed", alpha=0.5) +
  labs(y = 'Soft Threshold (power)', x = '', 
    caption = caption) +
  # scale_x_continuous(position = "top") +
  theme_light(base_family = "GillSans",base_size = 16) -> psave

psave

allowWGCNAThreads()

wd <- paste0(dir,"/",Sys.Date(),"_wgcna")

system(paste0('mkdir ', wd))

setwd(wd)

# The variable datExpr now contains the expression data ready for network analysis.

getwd()

bwnet <- blockwiseModules(datExpr, 
  maxBlockSize = 2000,
  power = 29, # softPower, 
  TOMType = "unsigned", 
  networkType = "unsigned",
  minModuleSize = 50,
  corType = "bicor",
  reassignThreshold = 0, 
  mergeCutHeight = 0.3,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  saveTOMFileBase = "TOM-blockwise",
  verbose = 3)

saveRDS(bwnet, "bwnet.rds")
