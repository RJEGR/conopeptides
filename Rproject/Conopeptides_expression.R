


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

dir <- "C://Users//cinai/OneDrive/Documentos/PUBLICATION_DIR/03.Coverage/"


f <- list.files(dir, pattern = "counts.txt", full.names = T) 

Manifest <- list.files(dir, pattern = "Manifest.tsv", full.names = T)

Manifest <- read_tsv(Manifest) %>% distinct()

cols <- read_tsv(f[1], skip = 1) %>% select(contains(".sorted.bam")) %>% names()

df <- read_tsv(f[1], skip = 1) 

# BLAST
f <- list.files(dir, pattern = "outfmt6", full.names = T) 

read_outfmt6 <- function(f) {
  
  # seqid = transcript_id
  outfmt6.names <- c("transcript_id", "subject", "identity", "coverage", "mismatches", "gaps", "seq_start", "seq_end", "sub_start", "sub_end", "e", "score")
  
  
  df <- read_tsv(f, col_names = F) 
  
  colnames(df) <- outfmt6.names
  
  df <- df %>% mutate(db = basename(f))
  
  return(df)
  
  
}

blast <- do.call(rbind, lapply(f[1], read_outfmt6)) %>%
  mutate(subject =  sapply(strsplit(subject, "[|]"), `[`, 1))

dir <- "~/Documents/GitHub/conopeptides/"

conoserverDB <- read_tsv(file.path(dir, "conoserver_protein.tsv")) %>%
  rename("subject" = "identifier") 

data <- df %>% 
  rename("transcript_id" = "Chr") %>%
  right_join(distinct(blast, transcript_id, subject)) %>%
  select(all_of(c("subject",cols))) %>%  
  left_join(conoserverDB) %>%
  mutate(subject = ifelse(is.na(`gene superfamily`), name, `gene superfamily`)) %>%
  mutate_if(is.numeric, replace_na, replace = 0) %>%
  group_by(subject) %>%
  summarise_at(vars(cols), sum)


rowNames <- data$subject

data <- data %>% select(all_of(cols)) %>% as("matrix") 
data <- DESeq2::varianceStabilizingTransformation(round(data)+1)
rownames(data) <- rowNames

colNames <- colnames(data)

colNames <- sapply(strsplit(basename(colNames), "_L"), `[`, 1)

colnames(data) <- gsub(".sorted.bam", "", colNames)

z_scores <- function(x) {(x-mean(x))/sd(x)}

datExpr <- t(apply(data, 1, z_scores))


h <- heatmap(datExpr, col = cm.colors(12), keep.dendro = T)

hc_samples <- as.hclust(h$Colv)
plot(hc_samples)
hc_order <- hc_samples$labels[h$colInd]

hc_genes <- as.hclust(h$Rowv)
plot(hc_genes)
order_genes <- hc_genes$labels[h$rowInd]


Dataviz <- datExpr %>% 
  as_tibble(rownames = 'subject') %>%
  pivot_longer(cols = all_of(colnames(datExpr)), values_to = "fill", names_to = "LIBRARY_ID") %>%
  # right_join(distinct(blast, transcript_id, subject)) %>%
  left_join(Manifest) %>%
  left_join(conoserverDB) 
  
Dataviz <- Dataviz %>% drop_na(fill)

lo = floor(min(Dataviz$fill))
up = ceiling(max(Dataviz$fill))
mid = (lo + up)/2

gene_names <- Dataviz %>% 
  mutate(`gene superfamily` = ifelse(is.na(`gene superfamily`), `protein type`, `gene superfamily`)) %>%
  mutate(`cysteine framework` = ifelse(is.na(`cysteine framework`), `protein type`, `cysteine framework`)) %>%
  distinct(subject, `gene superfamily`, `cysteine framework`) %>%
  mutate(relabel = paste0(subject, " (",`gene superfamily`," | ", `cysteine framework`,")")) %>% 
  pull(relabel, name = subject)

gene_names <- gene_names[match(order_genes, names(gene_names))]

identical(names(gene_names), order_genes) # TRUE

Dataviz %>%
  mutate(Feeding = factor(Feeding, levels = c("Ctrl", "Cam", "Lit", "Mix", "Pol"))) %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = hc_order)) %>%
  mutate(subject = factor(subject, levels = order_genes)) %>%
  ggplot(aes(x = LIBRARY_ID, y = subject, fill = fill)) + 
  # facet_grid(gene_wrap ~ CONTRASTE_D, scales = "free", space = "free") + 
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
    na.value = "white", midpoint = mid, limit = c(lo, up), 
    breaks = seq(lo, up, by = 3),
    name = NULL) +
  labs(x = "", y = "") +
  scale_y_discrete(position = 'left', labels = gene_names) +
  scale_x_discrete(position = "bottom") +
  # ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top') +
  # ggh4x::scale_y_dendrogram(hclust = hc_genes, position = "left", labels = NULL) +
  guides(
    fill = guide_colorbar(barwidth = unit(1.5, "in"),
      barheight = unit(0.1, "in"), label.position = "bottom",
      alignd = 0.5,
      title = "Row Z-score",
      title.position  = "top",
      title.theme = element_text(size = 10, family = "GillSans", hjust = 1),
      ticks.colour = "black", ticks.linewidth = 0.35,
      frame.colour = "black", frame.linewidth = 0.35,
      label.theme = element_text(size = 12, family = "GillSans"))
    # y.sec = ggh4x::guide_axis_manual(labels = gene_names, label_size = 6)
  ) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top",
    axis.line.y = element_line(color = 'white'),
    axis.text.y = element_text(hjust = 1, size = 3),
    # axis.ticks.y = element_blank(),
    axis.ticks.length = unit(2.5, "pt")) -> p

p <- p +  theme(
  # axis.ticks.x = element_blank(), 
  axis.text.x = element_text(angle = -90, hjust = 0, size = 10),
  # axis.text.x = element_blank(), 
  axis.line.x = element_blank()
)

p <- p + theme(panel.spacing.x = unit(0, "mm"))

# `gene superfamily` `cysteine framework`

p <- p + ggh4x::facet_nested( ~ Feeding+Time, scales = "free", space = "free") + 
  theme(strip.background = element_rect(fill = 'grey89', color = 'white'),
  strip.text.x.top = element_text(color = "black",hjust = 0, size = 8),
    strip.text.y.right = element_text(color = "black",hjust = 0, size = 8))

ggsave(p, filename = 'known_protein_expression_.png', path = dir, 
  width = 5, height = 10, device = png, dpi = 300)


p <- Dataviz %>% 
  filter(fill > 0) %>%
  mutate(`gene superfamily` = ifelse(is.na(`gene superfamily`), `protein type`, `gene superfamily`)) %>%
  mutate(`cysteine framework` = ifelse(is.na(`cysteine framework`), `protein type`, `cysteine framework`)) %>%
  mutate(Feeding = factor(Feeding, levels = c("Ctrl", "Cam", "Lit", "Mix", "Pol"))) %>%
  mutate(Time = factor(Time, levels = c("Ctrl", "2", "4", "6"))) %>%
  count(Time, Feeding, `protein type`, `cysteine framework`,`gene superfamily`, sort = T) %>%
  mutate(relabel = paste0("(",`gene superfamily`," | ", `cysteine framework`,")")) %>%
  group_by(relabel, Time) %>%
  # arrange(desc(n)) %>%
  mutate(freq = n / sum(n)) %>%
  mutate(relabel = factor(relabel, levels=unique(relabel))) %>%
  ggplot(aes(y = relabel, x = freq, fill = Feeding)) + 
  geom_col() +
  scale_x_continuous(breaks = seq(0,1, by = 0.5)) +
  facet_grid(~ Time) +
  ggthemes::scale_fill_calc() + theme_bw(base_size = 12, base_family = "GillSans") +
  theme(strip.background.y = element_rect(fill = 'grey89', color = 'white'),
    strip.text.y = element_text(angle = 0, size = 12, hjust = 0)) +
  theme(legend.position = "top")
  
ggsave(p, filename = 'known_protein_expression_bar.png', path = dir, 
  width = 7, height = 10, device = png, dpi = 300)


# Diversity plot


prevelancedf = apply(X = data,
  MARGIN = 2,
  FUN = function(x){sum(x > 10)})

hist(prevelancedf)

df <- data.frame(Diversity = prevelancedf, 
  TotalAbundance = colSums(data))

df %>% 
  as_tibble(rownames = 'LIBRARY_ID') %>%
  left_join(Manifest) %>%
  mutate(Feeding = factor(Feeding, levels = c("Ctrl", "Cam", "Lit", "Mix", "Pol"))) %>%
  mutate(Time = factor(Time, levels = c("Ctrl", "2", "4", "6"))) %>%
  # left_join(conoserverDB) %>%
  ggplot(aes(Diversity, LIBRARY_ID, fill = Feeding)) + geom_col() +
  facet_grid(Time ~., scales = "free_y", space = "free_y") +
  ggthemes::scale_fill_calc()


df %>% 
  as_tibble(rownames = 'LIBRARY_ID') %>%
  left_join(Manifest) %>%
  mutate(Feeding = factor(Feeding, levels = c("Ctrl", "Cam", "Lit", "Mix", "Pol"))) %>%
  mutate(Time = factor(Time, levels = c("Ctrl", "2", "4", "6"))) %>%
  # left_join(conoserverDB) %>%
  ggplot(aes(y = Diversity, x = Feeding, fill = Feeding)) + 
  geom_boxplot() +   ggthemes::scale_fill_calc()
