# analyze False or redundant conotoxin trasncript
# Based on contig impact score, SP signal, and 
# Redundant transcript encoding same CDS probably correspond to:
# 1) allelic variation encoding same peptide (contig_impact_score > 0) OR
# 2) Missassembled transcripts (contig_impact_score < 0)
# Missassembled (ie.e Incomplete or chimeric transcripts) not include SP signal (Signalp_class == SP)
# Those display a smooth decay in contig_impact_Score 
# Conflicts in classification correlate to some degee with low conig impact_Score
# 

# Previz filtering (facet of filters, y= sf and x=n transcripts)

df1 <- DB %>%
  filter(nchar(pep_seq) < 200) %>%
  dplyr::count(Superfamily, prediction_tool, sort = T) %>%
  drop_na(Superfamily) %>%
  mutate(facet = "A) Total: 8,861")

sf_levs <- df1 %>% distinct(Superfamily) %>% pull()

df2 <- DB %>%
  filter(nchar(pep_seq) < 200 & Signalp_class == "SP") %>%
  dplyr::count(Superfamily, prediction_tool, sort = T) %>%
  mutate(facet = "B) Signalp_class (SP): 2,461")

df3 <- DB %>%
  filter(nchar(pep_seq) < 200 & Signalp_class == "SP" & contig_impact_score > 0) %>%
  dplyr::count(Superfamily, prediction_tool, sort = T) %>%
  mutate(facet = "C) SP & contig_impact_score (CQ): 2,010")


df4 <- DB %>%
  filter(nchar(pep_seq) < 200 & Signalp_class == "SP" & contig_impact_score > 0 & prediction_tool == "BOTH") %>%
  dplyr::count(Superfamily, prediction_tool, sort = T) %>%
  mutate(facet = "D) SP & CQ & prediction_tool: 625")

rbind(df1, df2, df3, df4) %>%
  drop_na(Superfamily) %>%
  mutate(Superfamily = factor(Superfamily, levels = rev(sf_levs))) %>%
  ggplot(aes(x = n, y = Superfamily, fill = prediction_tool)) +
  facet_grid(~ facet, scales = "free_x") +
  geom_col() +
  # geom_text(aes(label= n), hjust= 1, vjust = 0.5, size = 3, family = "GillSans", color = "black") +
  theme_bw(base_family = "GillSans", base_size = 12) +
  labs(x = "", y = "Superfamily") +
  theme(
    legend.position = "top",
    # panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    # axis.text.y.right = element_text(angle = 0, hjust = 1, vjust = 0, size = 2.5),
    axis.text.y = element_text(angle = 0, size = 7),
    axis.text.x = element_text(angle = 0),
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 0, size = 10)) -> p 
p

ggsave(p, filename = 'Superfamilies_by_filters.png', 
  path = pub_dir, width = 12, height = 10, dpi = 500, device = png)

DataViz <- DB %>%
  # filter(Signalp_class == "SP" & contig_impact_score > 0 & prediction_tool == "BOTH") %>%
  mutate(uniprotkb_toxprot = sapply(strsplit(uniprotkb_toxprot, "[|]"), `[`, 2)) %>%
  mutate(uniprotkb_toxprot = ifelse(is.na(uniprotkb_toxprot), "uID", uniprotkb_toxprot)) %>%
  mutate(conoserver_protein = ifelse(is.na(conoserver_protein), "uID", conoserver_protein)) %>%
  mutate(conoserver_protein = ifelse(is.na(conoserver_protein), "uID", conoserver_protein)) %>%
  mutate(hmm_pred_conodictor = stringr::str_to_sentence(hmm_pred_conodictor)) %>%
  mutate(hmm_pred_conodictor = ifelse(is.na(hmm_pred_conodictor), "uID", hmm_pred_conodictor)) %>%
  mutate(Superfamily = ifelse(is.na(Superfamily), "uID", Superfamily)) %>%
  mutate(prediction_tool = ifelse(tab %in% "pHMM", paste0(prediction_tool,"_",tab), prediction_tool)) %>%
  mutate(sf = Superfamily) %>%
  select(protein_id, pep_seq, sf, Signalp_class, prediction_tool, hmm_pred_conodictor, Superfamily, uniprotkb_toxprot, conoserver_protein)

Redundancydf <- DataViz %>% count(pep_seq, sf, sort = T) %>% mutate(pep_len = nchar(pep_seq)-1)

Redundancydf

DB %>% mutate(dna_seq = nchar(dna_seq)) %>%
  mutate(x = contig_impact_score, x = sign(x) * log(1+abs(x)) ) %>% 
  ggplot(aes(y = dna_seq, x = x)) + geom_point(shape = 1, aes(alpha = effective_length_frac))

DB %>% 
  mutate(x = contig_impact_score, x = sign(x) * log(1+abs(x)) ) %>% 
  ggplot(aes(y = x, x = Superfamily)) + geom_jitter(aes(color = prediction_tool, size = effective_length_frac), shape = 1)

DB %>% 
  ggplot(aes(y = effective_length_frac, x = Signalp_class)) + geom_jitter(aes(color = prediction_tool))
