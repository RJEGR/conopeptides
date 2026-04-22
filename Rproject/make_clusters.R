

library(tidyverse)

rm(list = ls());

if(!is.null(dev.list())) dev.off()

dir <- "~/Documents/GitHub/conopeptides/04.Merge/"

psl_file <- list.files(path = dir, pattern = "Trinity_vs_transcripts.blat.psl", full.names = T)

##### inputs
tc_psl_file="trinity.chicken.psl"
ch_psl_file="chicken.human.psl"
th_psl_file="trinity.human.psl"

##### output
out_file=paste0(gsub(".psl", "", basename(psl_file)), ".clusters.txt")


out_file <- file.path(dir, out_file)

####################################################
## Processing the pairwise alignment between assemblies.
## Overlapping sequences will be merged based on the blat output
####################################################

tc_psl=unique(read.delim(psl_file,skip=5,stringsAsFactors=F,header=F)[,c(10,14)])

#each gene starts as a cluster
gene_clusters<-1:length(unique(tc_psl$V14))
names(gene_clusters)<-unique(tc_psl$V14)

#work out when a contig matches two or more genes
c_split=split(tc_psl$V14,tc_psl$V10)
#remove contigs that match multiple genes (these are possible chimeras)
c_split=c_split[sapply(c_split,length)==1]

contig_clusters=gene_clusters[unlist(c_split)]
gene_clusters=gene_clusters[unique(names(contig_clusters))]
names(contig_clusters)<-names(unlist(c_split))

clusters=c(gene_clusters,contig_clusters)
clusters=sort(clusters)

write.table(clusters,out_file,quote=F,sep="\t",col.names=F)
