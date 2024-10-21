# algunas mejoras en desempeno para unir ensambles fueron reportadas para transfuse, seguido de evidentialgene (https://pmc.ncbi.nlm.nih.gov/articles/PMC7188282)

#!/bin/sh
## Directivas
#SBATCH --job-name=MMseqs
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20


EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/mmseqs/bin
export PATH=$PATH:$EXPORT

thread_count=$SLURM_NPROCS

cat $f1 $f2 > Merged.fasta

# s1)
mmseqs createdb Merged.fasta DB 

# s2
# How to redundancy filter sequences with identical length and 100% length overlap?
# To redundancy filter sequences of identical length and 100% overlap mmseqs clusthash can be used. It reduces each sequence to a five-letter alphabet, computes a 64 bit CRC32 hash value for the full-length sequences, and places sequences with identical hash code that satisfy the sequence identity threshold into the same cluster.

#mmseqs clusthash sequenceDB resultDB --min-seq-id 0.9
#mmseqs clust sequenceDB resultDB clusterDB

# How to add sequence identities and other alignment information to a clustering result
# We can add sequence identities and other alignment information to the clustering result outDB by running an additional align step

mmseqs cluster DB Merged_clusters tmp --threads $thread_count --min-seq-id 0.98

# s3 create tsv format
mmseqs createtsv DB DB Merged_clusters Merged_clusters.tsv

exit
