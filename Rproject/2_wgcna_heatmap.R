
# label Conotoxins

library(WGCNA)
library(flashClust)
library(tidyverse)


rm(list = ls());

if(!is.null(dev.list())) dev.off()

pub_dir <- "/Users/cigom/Documents/GitHub/conopeptides/PUBLICATION_DIR/"

dir <- "/Users/cigom/Documents/GitHub/conopeptides/06.Quantification/MATRIX_RSEM_dir/"


DB <- read_rds(paste0(pub_dir, "/structured_db.rds"))


.datTraits <- list.files(dir, pattern = "Manifest", full.names = T) 

.datTraits <- read_tsv(.datTraits) %>% 
  # select_at(vars(starts_with("Group"))) %>%
  # rename("LIBRARY_ID" = "Group2", "Population" = "Group3", "Condition" = "Group4") %>%
  mutate(datTraits = paste0(Diatery, "-", Time))

datTraits <- table(.datTraits$LIBRARY_ID, .datTraits$datTraits)

# 
# .datTraits <- .datTraits %>% 
#   mutate(values_from = 1) %>%
#   pivot_wider(names_from = LIBRARY_ID, values_from = values_from, values_fill = 0) 
# 
# datTraits <- .datTraits %>% 
#   mutate_if(is.character, as.factor) %>% 
#   mutate_if(is.factor, as.numeric) %>%
#   data.frame(row.names = datTraits$LIBRARY_ID)

f <- list.files(path = dir, pattern = "wgcna.input.rds", full.names = T)

dim(datExpr <- readRDS(f))


bwnnet_dir <- paste0(dir, "/2025-03-25_wgcna/")

bwnnet_f <- list.files(path = bwnnet_dir, pattern = "bwnet.rds", full.names = T)

bwnet <- readRDS(bwnnet_f)

bwmodules <- labels2colors(bwnet$colors)

names(bwmodules) <- names(bwnet$colors)

barplot(table(bwmodules))

length(table(bwmodules))

# reads <- colSums(datExpr)
sum(keep <- colnames(datExpr) %in% names(bwmodules))
reads <- colSums(datExpr[,keep])
Total <- sum(reads)

data.frame(reads, bwmodules) %>% 
  as_tibble() %>% 
  group_by(bwmodules) %>% 
  summarise(n = n(), reads = sum(reads)) %>%
  dplyr::rename('WGCNA' = 'bwmodules') -> stats


as_tibble(bwmodules, rownames = "gene_id") %>% 
  dplyr::rename("WGCNA" = "value") %>% 
  write_tsv(file = paste0(pub_dir, "/WGCNA.tsv"))

str(unique(DB$gene_id)) # 66, 528 unique genes

sum(sort(unique(DB$gene_id)) %in% sort(names(bwmodules)))

WGCNADB <- read_tsv(paste0(pub_dir, "/WGCNA.tsv")) %>% right_join(DB)

# split EGGnog by modules

WGCNADB %>%
  filter(COG_category != "-") %>%
  select(protein_id, WGCNA, COG_category) %>%
  mutate(COG_category = strsplit(COG_category, "")) %>%
  unnest(COG_category) %>%
  right_join(NOG.col, by = c("COG_category" = "EggNM.COG_category")) %>%
  mutate(COG_category = paste0(COG_category, ", ", COG_name)) %>%
  dplyr::count(COG_category, WGCNA, sort = T) %>%
  group_by(COG_category) %>% mutate(frac = n/sum(n)) %>%
  arrange(desc(frac)) %>%
  mutate(COG_category = factor(COG_category, levels = unique(COG_category))) %>%
  # mutate(facet = "A) Transcriptome") %>%
  drop_na(WGCNA) %>%
  ggplot(aes(y = WGCNA, x = COG_category, fill = frac)) +
  # labs(y = "Nested Orthologous Gene Group (NOGS)", x = "Enrichment ratio (N transcripts/Total transcripts)") +
  geom_tile() +
  # scale_x_continuous(labels = scales::percent_format()) +
  ggh4x::scale_y_dendrogram(hclust = hclust) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  # xlim(0,0.4) +
  theme(legend.position = "top", 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    # panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = 'grey95', color = 'white')) 
     geom_text(data = data_text, 
    aes(label = n), size = 3.5,
    hjust = -0.1, vjust = 0, 
    family = "GillSans", position = position_dodge(width = 1)) 

# find toxins by modules
# as datExp was filtered to low, some peptides was dropped during the Wgcna analysis, therefore, focus the analysis to conopeptides?

conopep_module <- WGCNADB %>%
  drop_na(tab) %>%
  distinct(tab, Signalp_class, Superfamily, protein_id, WGCNA) %>%
  filter(Signalp_class == "SP") %>%
  # filter(is.na(WGCNA))
  dplyr::count(WGCNA, sort = T) %>%
  drop_na(WGCNA) %>%
  mutate(n = paste0("sf: ",n,""))  %>% dplyr::rename("n_conopeptides" = "n")

Thioredoxin_module <- WGCNADB %>%
  filter(if_any(where(is.character), ~ grepl(pattern = 'Thioredoxin', x = .x, ignore.case = T))) %>%
  distinct(COG_category, Preferred_name, protein_id, PFAMs, WGCNA) %>%
  dplyr::count(WGCNA, sort = T) %>%
  drop_na(WGCNA) %>%
  mutate(n = paste0("pdi: ",n,"")) %>% dplyr::rename("n_pdi_domain" = "n")


# Recalculate MEs with color labeLIBRARY_ID# Recalculate MEs with color labels

MEs0 = moduleEigengenes(datExpr, bwmodules)$eigengenes

MEs = orderMEs(MEs0)

names(MEs) <- str_replace_all(names(MEs), '^ME', '')

identical(sort(rownames(datTraits)), sort(rownames(MEs)))

# datTraits <- datTraits[match(rownames(datTraits), rownames(MEs)),]

# identical(rownames(datTraits), rownames(MEs))

moduleTraitCor = cor(MEs, datTraits, use= "p")

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(datTraits))

moduleTraitCor %>% as_tibble(rownames = 'WGCNA') %>% 
  pivot_longer(-WGCNA, values_to = 'moduleTraitCor') -> df1

moduleTraitPvalue %>% as_tibble(rownames = 'WGCNA') %>% 
  pivot_longer(-WGCNA, values_to = 'corPvalueStudent') %>%
  right_join(df1) -> df1

hclust <- hclust(dist(moduleTraitCor), "complete")

plot(hclust)
# up_df %>% distinct(transcript) %>% pull() -> updegs
# down_df %>% distinct(transcript) %>% pull() -> dwndegs

bwmodules %>% 
  as_tibble(., rownames = 'transcript') %>%
  # mutate(degs = ifelse(transcript %in% updegs, 'up', 
  # ifelse(transcript %in% dwndegs, 'down', ''))) %>%
  dplyr::rename('WGCNA' = 'value') -> bwModuleCol

bwModuleCol %>% group_by(WGCNA) %>% dplyr::count(sort = T) %>% left_join(stats)

bwModuleCol %>% group_by(WGCNA) %>% dplyr::count(sort = T) -> bwModuleDF

bwModuleDF %>% mutate(WGCNA = factor(WGCNA, levels = hclust$labels[hclust$order])) -> bwModuleDF

bwModuleDF %>% group_by(WGCNA) %>% mutate(pct = n / sum(n)) -> bwModuleDF

bwModuleDF

bwModuleDF <- bwModuleDF %>% left_join(Thioredoxin_module) %>% left_join(conopep_module)


df1 %>%
  mutate(star = ifelse(corPvalueStudent <.001, "***", 
    ifelse(corPvalueStudent <.01, "**",
      ifelse(corPvalueStudent <.05, "*", "")))) -> df1

df1 <- df1 %>% separate(name, into = c("Diatery", "Time"), sep = "-", remove = F)

lo = floor(min(df1$moduleTraitCor))
up = ceiling(max(df1$moduleTraitCor))
mid = (lo + up)/2

time_levs <- c("Ctrl", "2", "4", "6")
recode_time <- structure(c("Control", "2 months", "4 months", "6 months"), names = time_levs)

Diatery_levs <- c("Ctrl","Cam", "Lit", "Pol", "Mix")
recode_Diatery <- structure(c("Control","Shrimp", "Mollusk", "Polychaete", "Mixed"), names = Diatery_levs)


y_labels <- bwModuleDF %>% 
  mutate(label = paste0(WGCNA, " (", n,")")) %>% 
  mutate(label = ifelse(!is.na(n_pdi_domain), paste0(label, " ", n_pdi_domain), label)) %>%
  mutate(label = ifelse(!is.na(n_conopeptides), paste0(label, " ", n_conopeptides), label)) %>%
  pull(label, name = WGCNA)

# y_labels <- bwModuleDF %>% mutate(label = paste0(WGCNA, " (", n,")")) %>% pull(label, name = WGCNA)

# library(ggh4x)

df1 %>%
  mutate(Time = recode_factor(Time, !!!recode_time, .ordered = T)) %>%
  mutate(Diatery = recode_factor(Diatery, !!!recode_Diatery, .ordered = T)) %>%
  mutate(moduleTraitCor = round(moduleTraitCor, 2)) %>%
  mutate(star = ifelse(star != '', paste0(moduleTraitCor, '(', star,')'), '')) %>%
  ggplot(aes(y = WGCNA, x = Time, fill = moduleTraitCor)) +
  facet_grid(~ Diatery, scales = "free_x", space = "free_x") +
  # geom_tile(color = 'black', size = 0.5, width = 0.7) + 
  geom_raster() +
  geom_text(aes(label = star),  vjust = 0.5, hjust = 0.5, size= 1.5, family =  "GillSans") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
    na.value = "white", midpoint = mid, limit = c(lo, up),
    name = NULL) +
  ggh4x::scale_y_dendrogram(hclust = hclust, labels = NULL) +
  labs(x = '', y = 'Module') +
  theme_classic(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top",
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 0, size = 7),
    axis.line.y = element_line(color = 'white'),
    axis.text.y = element_text(hjust = 1),
    axis.ticks.length = unit(5, "pt")) -> p1 

p1 <- p1 +   guides(
  fill = guide_colorbar(barwidth = unit(1.5, "in"),
    barheight = unit(0.1, "in"), label.position = "bottom",
    alignd = 0.5,
    title = "Trait correlation",
    title.position  = "top",
    title.theme = element_text(size = 10, family = "GillSans", hjust = 0),
    ticks.colour = "black", ticks.linewidth = 0.35,
    frame.colour = "black", frame.linewidth = 0.35,
    label.theme = element_text(size = 10, family = "GillSans")),
  
  y.sec = ggh4x::guide_axis_manual(labels = y_labels, label_size = 8, label_family = "GillSans")
  
) 

p1 <- p1 + theme(axis.ticks.x = element_blank(), 
  axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
  # axis.text.x = element_blank(), 
  axis.line.x = element_blank())
# 
p1 <- p1 + theme(panel.spacing.x = unit(0, "mm"))

# p1

ggsave(p1, filename = 'ModuleTraitRelationship.png', 
  path = pub_dir, width = 5.5, height = 10, dpi = 1000, device = png)



# Intramodular connectivity ----
# We would find modules containing genes w/ high positive / negative correlations in spite of a variable (for example morphotype awa site) previosly correlated with plotEigengeneNetworks

# Calculate the correlations between modules

geneModuleMembership <- as.data.frame(WGCNA::cor(datExpr, MEs, use = "p"))

# pca of modules

PCA <- prcomp(t(geneModuleMembership), scale. = FALSE)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
# sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dtvis <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])


dtvis %>%
  mutate(id = rownames(.)) %>%
  ggplot(., aes(PC1, PC2)) +
  # geom_point(size = 5, alpha = 0.9) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  geom_text(aes(label = id), alpha = 0.9) +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) 
# theme_bw(base_family = "GillSans", base_size = 16) +
# theme(plot.title = element_text(hjust = 0.5), legend.position = 'top')


# What are the p-values for each correlation?
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr)))

# What's the correlation for the trait?
geneTraitSignificance <- as.data.frame(cor(datExpr,datTraits, use = "p"))

# What are the p-values for each correlation?
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples = nrow(datExpr)))
names(geneTraitSignificance) <- paste("GS.", colnames(datTraits), sep = "")
names(GSPvalue) <- paste("p.GS.", colnames(datTraits), sep = "")


identical(rownames(geneTraitSignificance), rownames(GSPvalue))

bind <- cbind(geneTraitSignificance, GSPvalue) %>% as_tibble(rownames = "Transcript")

as_tibble(bwmodules, rownames = "Transcript") %>% 
  dplyr::rename("WGCNA" = "value") %>% 
  left_join(bind) %>% 
  write_tsv(file = paste0(dir, "GeneTraitSignificance.tsv"))


GSpval <- GSPvalue %>%
  tibble::rownames_to_column(var = "ids")

gMM_df <- geneModuleMembership %>%
  tibble::rownames_to_column(var = "ids") %>%
  gather(key = "module", value = "moduleMemberCorr", -ids) 

gMM_df %>% 
  mutate(module = str_replace_all(module, '^ME', '')) %>%
  # mutate(module = factor(module, levels = yLabels)) %>% 
  as_tibble() -> gMM_df

# Prepare gene significance df
GS_bacprod_df <- geneTraitSignificance %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "ids")


# Put everything together 
allData_df <- gMM_df %>%
  left_join(GS_bacprod_df, by = "ids") %>%
  left_join(GSpval, by = "ids") %>%
  as_tibble()

# Filter by significance genes 

# allData_df %>% filter(ids %in% names(bwmodules_sig)) -> allData_df

bwmodules %>% as_tibble(rownames = 'ids') %>% rename('module' = 'value') %>%
  left_join(allData_df) -> allData_df

saveRDS(allData_df, file = paste0(dir,'allData_df.rds'))

# Sanity check
allData_df %>% group_by(module) %>% count() 

# Longer
allData_df %>% select_at(vars(starts_with('GS'))) %>% names() -> colName

allData_df %>%
  select_at(vars(starts_with('GS'), 'module')) %>%
  # mutate(n = 1:nrow(.)) %>%
  pivot_longer(cols = all_of(colName), names_to = 'group', values_to = 'GS') %>%
  mutate(group = str_replace_all(group, '^GS.', '')) %>% # pairs = paste0(module,'-', group)
  mutate(GS = abs(GS)) -> gs_df

allData_df %>% select_at(vars(starts_with('p.GS'))) %>% names() -> colName

allData_df %>%
  select_at(vars(starts_with('p.GS'), 'moduleMemberCorr')) %>%
  # mutate(n = 1:nrow(.)) %>%
  pivot_longer(cols = colName, names_to = 'group', values_to = 'p.GS') %>%
  mutate(group = str_replace_all(group, '^p.GS.', '')) %>% 
  select(p.GS, moduleMemberCorr) %>%
  cbind(gs_df, .) -> df_viz

plotscatterSig(df_viz, M = psME)

colVal <- c("#54278f", "grey30")

df_viz %>%
  mutate(color = ifelse(p.GS < 0.05, '*', 'ns')) %>%
  filter(module %in% "blue") %>%
  # filter(group != 'HC') %>%
  ggplot(aes(x = moduleMemberCorr, y = abs(GS), color = color, alpha = -log10(p.GS))) +  
  # geom_vline(color = "black", xintercept = 0.8, linetype="dashed", alpha=0.5) +
  # geom_hline(color = "black", yintercept = 0.8, linetype="dashed", alpha=0.5) +
  # geom_hline(color = "black", yintercept = -0.8, linetype="dashed", alpha=0.5) +
  geom_point(size = 1) +
  labs(x = 'Module membership', y = 'Gene significance') + # 'Gene Correlation'
  scale_color_manual(name = '', values = colVal) + # breaks = c(0, 0.25 ,0.5, 0.75, 1), limits= c(0, 1)
  scale_alpha_continuous(name = expression(-Log[10] ~ "P")) +
  facet_grid(module~group) +
  # ggh4x::facet_nested_wrap(module+group ~., nest_line = T) +
  theme_bw(base_family = "GillSans", base_size = 16) +
  # geom_smooth(method = "lm", linetype="dashed", size = 0.5, alpha=0.5, color = 'white',
  #   se = TRUE, na.rm = TRUE) +
  theme(legend.position = 'top') -> psave

allData_df %>%
  ggplot(aes(moduleMemberCorr)) +
  geom_histogram() +
  facet_wrap(~ module)
