#!/bin/bash

#SBATCH -p cicese
#SBATCH --job-name=corset
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 06-00:00:00

# 

EXPORT=/LUSTRE/apps/bioinformatica/corset-1.09-linux64/
export PATH=$PATH:$EXPORT


# corset -f true -r true *.bam
# corset  -f true -i corset *.corset-reads


for f in $(ls *.sorted.bam.corset-reads)
do
basename=${f##*/}
bs="${basename%.sorted.bam.corset-reads}"

if [ ! -f "${bs}.sorted.bam.corset-reads" ]; then

    call="corset -f true -r true ${bs}.sorted.bam"

    echo "Running corst in '${bs}' sample."

    echo $call

    eval $call

else
    echo "'${bs}.sorted.bam.corset-reads' already exists."
    #echo "Continue with next sample."
fi

    echo "Continue with next step for corset."

done

corset  -f true -i corset *.corset-reads

exit

# https://github.com/Adamtaranto/Corset-tools
# error in Target cluster not in Map file: TRINITY_***
# see NotFound_clusters.txt.log

EXPORT=~/Documents/Tools/Corset-tools
# EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/Corset-tools-master/
export PATH=$PATH:$EXPORT


# module load conda-2024

fetchClusterSeqs.py -i ../Merged_clusters.fasta -t clusters.txt -o Merged_clusters_fetched.fasta -c clusters.txt


# Here we propose an alternative representation for each gene, which we refer to as a superTranscript. SuperTranscripts contain the sequence of all exons of a gene without redundancy (Figure 1A). They can be constructed from any set of transcripts including de novo assemblies and we have developed a python program to build them called Lace (available from https://github.com/Oshlack/Lace/wiki). Lace works by building a splice graph[13] for each gene, then topologically sorting the graph using Kahn’s algorithm[14] (Figure 1B). Building superTranscripts is a simple post-assembly step that promises to unlock numerous analytical approaches for non-model organisms 
# https://www.biorxiv.org/content/10.1101/077750v3
