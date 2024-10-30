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

exit
