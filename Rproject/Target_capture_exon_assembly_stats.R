#  

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)


read_outfmt6 <- function(f) {
  
  # 1	qseqid	Query sequence ID
  # 2	sseqid	Subject sequence ID
  # 3	pident	Percentage of identical matches
  # 4	length	Alignment length
  # 5	mismatch	Number of mismatches
  # 6	gapopen	Number of gap openings
  # 7	qstart	Start of alignment in query sequence
  # 8	qend	End of alignment in query sequence
  # 9	sstart	Start of alignment in subject sequence
  # 10	send	End of alignment in subject sequence
  # 11	evalue	Expect value (E-value) for the alignment
  # 12	bitscore	Bit score of the alignment
  
  # seqid = transcript_id
  outfmt6.names <- c("transcript_id", "subject", "identity", "aln_length", "mismatches", "gaps", "seq_start", "seq_end", "sub_start", "sub_end", "e", "score")
  
  
  df <- read_tsv(f, col_names = F) 
  
  colnames(df) <- outfmt6.names
  
  df <- df %>% mutate(db = basename(f))
  
  return(df)
  
  
}

library(tidyverse)

dir <- "/Users/cigom/Documents/GitHub/conopeptides/07.Reference/"

dir <- "/Users/cigom/Documents/GitHub/conopeptides/03.Coverage/"

dir_out <- dir

f <- list.files(dir, pattern = "blastx.outfmt6", full.names = T)

dbname <- gsub("_vs_conoserver_protein.diamond.blastx.outfmt6","",basename(f))

DB1 <- do.call(rbind, lapply(f, read_outfmt6)) %>% mutate(db = dbname) %>%
  dplyr::rename("identifier" = "subject") %>%
  mutate(identifier =  sapply(strsplit(identifier, "[|]"), `[`, 1)) %>%
  # select(transcript_id, identifier, identity, e) %>%
  dplyr::rename("dna_identity" = "identity", "dna_e" = "e") 

# DB1 %>% count(identifier, sort = T)

# Load conoserverDB

dir <- "~/Documents/GitHub/conopeptides/"

conoserverDB <- read_tsv(file.path(dir, "conoserver_protein.tsv")) 

names(conoserverDB) <- paste0(gsub(" ","_", names(conoserverDB)), "_conoserver")

names(conoserverDB)[1] <- "identifier"

# Bind conoserver + conodictor|conosorter

DB1 <- DB1 %>% 
  # distinct(transcript_id, identifier) %>% 
  left_join(conoserverDB)

# A higher Bit Score indicates a more significant alignment between the query and subject sequences.

DB1 %>% 
  ggplot(aes(aln_length, score)) +
  # ggplot(aes(aln_length, mismatches)) + 
  geom_point(aes(color = -log10(dna_e))) + 
  facet_grid(~ protein_type_conoserver) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()) 


# Contrast to our assembly

dir <- "~/Documents/GitHub/conopeptides/04.Merge/"


f <- list.files(dir, pattern = "outfmt6", full.names = T) 

f <- f[grepl("Merged_clusters_vs", f)]

DB2 <- do.call(rbind, lapply(f, read_outfmt6))

DB2 <- DB2 %>% 
  dplyr::rename("identifier" = "subject") %>%
  mutate(identifier =  sapply(strsplit(identifier, "[|]"), `[`, 1)) %>%
  dplyr::rename("dna_identity" = "identity", "dna_e" = "e") %>%
  mutate(db = gsub("vs_conoserver_protein.fa.diamond.blastx.outfmt6","", db)) %>%
  left_join(conoserverDB)

rbind(DB1, DB2) %>%
  ggplot(aes(dna_identity, score)) +
  # ggplot(aes(aln_length, mismatches)) + 
  geom_point(aes(color = -log10(dna_e))) + 
  facet_grid(db ~ protein_type_conoserver) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank())


rbind(DB1, DB2) %>%
  # filter(protein_type_conoserver %in% "Precursor") %>%
  count(db, protein_type_conoserver) 
  

rbind(DB1, DB2) %>%
  ggplot(aes(score)) +
  geom_histogram(aes(fill = db)) + 
  # geom_point(aes(color = -log10(dna_e))) + 
  # facet_grid(db ~ protein_type_conoserver) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank())


# Coverage

# blast distribution hits plot ----


read_outfmt6 <- function(f) {
  
  # 1	qseqid	Query sequence ID
  # 2	sseqid	Subject sequence ID
  # 3	pident	Percentage of identical matches
  # 4	length	Alignment length
  # 5	mismatch	Number of mismatches
  # 6	gapopen	Number of gap openings
  # 7	qstart	Start of alignment in query sequence
  # 8	qend	End of alignment in query sequence
  # 9	sstart	Start of alignment in subject sequence
  # 10	send	End of alignment in subject sequence
  # 11	evalue	Expect value (E-value) for the alignment
  # 12	bitscore	Bit score of the alignment
  
  # seqid = transcript_id
  outfmt6.names <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "seq_end", "qend", "send", "evalue", "bitscore")
  
  
  df <- read_tsv(f, col_names = F) 
  
  colnames(df) <- outfmt6.names
  
  df <- df %>% mutate(db = basename(f))
  
  return(df)
  
  
}

ggblast <- function(blast_out, num_hit){
  require(zoo)
  require(RColorBrewer)
  require(dplyr)
  # usage: ggblast(blastout=bnout, 　     # blast結果data.frame
  #                num_query=c(1,3,5),    # uniqueなqueryの番号、もしくはid
  #                num_hit=10             # top10
  #                 ) 
  pacman::p_load(RColorBrewer, zoo, dplyr, ggplot2)
  
  # Etiqueta de leyenda
  vlab <- c(paste0(zoo::rollapply(c(0,seq(70,100,5)), width=2, by=1, 
    function(x){paste(x, collapse = ":")}), "(-)"),
    zoo::rollapply(c(0,seq(70,100,5)), width=2, by=1, 
      function(x){paste(x, collapse = ":")}))
  # identity - Codigo de color
  vcol <- c(brewer.pal(7, "Reds"), brewer.pal(7, "Blues"))
  if(missing(num_hit)){ num_hit <- sum(blast_out$qseqid %in% x) }
  # queryごとにplot -----
  gg_blast <-
    lapply(unique(blast_out$qseqid),
      function(x){
        # queryごとのデータフレーム
        qdat <-blast_out %>%
          # identityを離散値に変換
          mutate(pident=as.integer(cut(pident, breaks=c(0,seq(70,100,5)), 1:7))) %>%
          # hit position
          mutate(pident=ifelse((qstart-qend) > 0, pident, pident+7)) %>%
          # query selection
          filter(qseqid==x) %>%
          # top hitの数を指定
          filter(row_number() %in% 
              1:ifelse(is.null(num_hit), sum(blast_out$qseqid %in% x), num_hit)) %>%
          # identity, hitをfactorにして水準を逆に
          mutate(pident=factor(pident, levels=rev(levels(factor(pident)))),
            sseqid=factor(sseqid, levels=rev(unique(sseqid))))
        
        # queryごとのカラーコードとレジェンドラベル
        col_qdat <- vcol[as.integer(levels(factor(qdat$pident)))]
        lab_qdat <- vlab[as.integer(levels(factor(qdat$pident)))]
        
        # 
        ggplot(qdat, aes(x=qstart,xend=qend, y=sseqid, yend=sseqid, colour=pident)) +
          theme_bw() + labs(title=x) +
          geom_segment(size=3) +
          scale_color_manual(values = col_qdat, name='Key of\nAligment\nScore', labels = lab_qdat) +
          labs(title = paste0("Distribution of top ", num_hit, " blast hits")) + 
          labs(x = 'Scafold Position (Genome)', y = 'Consensus hairpin sequence') +
          theme_classic(base_size = 12) 
      }
    )
  return(gg_blast)
}


blast_out <- do.call(rbind, lapply(f, read_outfmt6)) %>% mutate(db = dbname)

ggblast(blast_out = blast_out, num_hit = 1)


blast_out %>%
  filter(grepl("P04038", sseqid)) %>%
  mutate(y = row_number()) %>%
  ggplot(aes(x=qstart, xend=qend, y=y, yend=y, colour=bitscore)) +
  theme_bw() + 
  # labs(title=x) +
  geom_segment(size=3)


# DB1 %>% count(protein_type_conoserver)
# DB1 %>% filter(protein_type_conoserver %in% "Precursor")
