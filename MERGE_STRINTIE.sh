#!/bin/bash
#SBATCH --job-name=stringtie
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20

EXPORT=/LUSTRE/apps/bioinformatica/stringtie/
export PATH=$PATH:$EXPORT


CPU=$SLURM_NPROCS

mkdir -p S2_STRINGTIE_MERGE

for i in $(ls *.sorted.bam)
do
withpath="${i}"
filename=${withpath##*/}
bs="${filename%*.sorted.bam}"
stringtie --rf -p $CPU -o S2_STRINGTIE_MERGE/${bs}_transcripts.gtf $i
done

exit
```

Then merge
```bash
cd S2_STRINGTIE_MERGE

ls -1 *_transcripts.gtf > stringtie_gtf_list.txt

stringtie --rf --merge -p 24 -o transcripts.gtf stringtie_gtf_list.txt

# 3) generate FASTA
# Using concat from both genomes, red and green
RNA_REF_FASTA=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Haliotis/Hybryd/RF_Ref.fna


gffread -w transcripts.fa -g $RNA_REF_FASTA transcripts.gtf

grep "^>" -c transcripts.fa # 76, 658 scaffolds from ctrl samples