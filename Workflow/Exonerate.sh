#!/bin/bash
#SBATCH --job-name=exonerate
#SBATCH -N 2
#SBATCH --mem=80GB
#SBATCH -t 6-00:00:00
#SBATCH --ntasks-per-node=20
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=rgomez@uabc.edu.mx
#SBATCH --error=slurm-%j.err

#EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/exonerate-2.2.0-x86_64/bin
EXPORT=/LUSTRE/apps/bioinformatica/exonerate-2.2.0-x86_64/bin
export PATH=$PATH:$EXPORT

QUERY=$1


# Step 1, prep exonerate sequence database file
# --softmask option allows you to specify whether the input FASTA sequences should be treated as softmasked. 
# Softmasking typically means that lower-case letters represent masked regions, which are ignored during alignment.

REFDIR=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/californicus/07.Reference/EXONERATE/TARGETMODEL/Lottia
REFERENCE=${REFDIR}/Lottia_gigantea.Lotgi1.dna_sm.toplevel.fa

ESDREF=${REFERENCE%.*}.esd

fasta2esd --fasta $REFERENCE --softmask yes --alphabet DNA --output $ESDREF

# Step 2,

# -ryo to report all the sections of query sequences which feature in alignments in fasta format

exonerate --model est2genome $QUERY $REFERENCE --softmasktarget yes --showvulgar FALSE --subopt FALSE --showalignment FALSE --ryo ">%qi %qd\n%qas\n"

# exonerate --model est2genome $QUERY $REFERENCE --softmasktarget yes --showalignment yes > ${QUERY.*}.est2genome.fa

exit