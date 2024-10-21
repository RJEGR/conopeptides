#!/bin/bash

#SBATCH -p cicese
#SBATCH --job-name=corset1
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 06-00:00:00

# 
 

## blat pairwise alignments between transcripts from transcriptome-transcriptome 
# From literature: Transcripts are aligned against one another using blat and the nodes of shared bases are merged. Lace then simplifies the graph by compacting unforked edges. The graph is topologically sorted and the resulting superTranscript annotated with transcripts and blocks (Davidson et al., 2017)
# https://github.com/Oshlack/superTranscript_paper_code/blob/master/chicken_superTranscriptome/get_clusters.sh

refone=$1
reftwo=$2


output=${refone%.*}_vs_${reftwo%.*}.blat.psl

blat $refone $reftwo -minScore=200 -minIdentity=98 $output

# add step 2, make_cluster.R and then corset and lace

exit


corset 


# Here we propose an alternative representation for each gene, which we refer to as a superTranscript. SuperTranscripts contain the sequence of all exons of a gene without redundancy (Figure 1A). They can be constructed from any set of transcripts including de novo assemblies and we have developed a python program to build them called Lace (available from https://github.com/Oshlack/Lace/wiki). Lace works by building a splice graph[13] for each gene, then topologically sorting the graph using Kahn’s algorithm[14] (Figure 1B). Building superTranscripts is a simple post-assembly step that promises to unlock numerous analytical approaches for non-model organisms 
# https://www.biorxiv.org/content/10.1101/077750v3

exit
