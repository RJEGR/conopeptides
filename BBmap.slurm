#!/bin/bash
#SBATCH --job-name=bbmap
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=rgomez@uabc.edu.mx
#SBATCH -t 6-00:00:00


EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/bbmap
export PATH=$PATH:$EXPORT


EXPORT=/LUSTRE/apps/bioinformatica/RSEM/bin/samtools-1.3/
export PATH=$PATH:$EXPORT

Manifest=$1 # samples.txt
ref=$2

REF_PREFIX=${ref%.*}

# ...to build an index in memory and map paired reads to it in a single command

WD=S1_BBMAP_SAM_BAM_FILES_${REF_PREFIX}_DIR

mkdir -p $WD

while read -r line; do
    in1=$(echo "$line" | awk '{print $3}')
    # in2=$(echo "$line" | awk '{print $4}')
    bs=${in1##*/}
    call="bbmap.sh in1="$in1" out="${bs%_1.fastq}.sam" mappedonly=t ref=$ref  threads=$SLURM_NPROCS"
    eval $call

    call="samtools sort -@ $SLURM_NPROCS -o ${bs}.sorted.bam ${bs}.sam"
    eval $call

    mv ${bs}.sorted.bam $WD

done < "$Manifest"

rm *.sam

exit