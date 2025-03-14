#!/bin/bash

#SBATCH -p cicese
#SBATCH --job-name=STAR
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 06-00:00:00

# increase the max number of open files with Linux command ulimit before running STAR
ulimit -s unlimited

# 
# Continue with Differential Transcript Usage on a non model organism
# Map reads to superTranscriptome using STAR


EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/STAR-2.7.11b/bin/Linux_x86_64_static
export PATH=$PATH:$EXPORT


limitRAM=100000000000 # $SLURM_MEM_PER_NODE

INPUT=SuperDuper.fasta
GFF=SuperDuper.gff

REF_PREFIX=${INPUT%.*}

HASH=`md5sum $INPUT | awk '{print $1}'`

genomeDir=INDEX_${HASH}

mkdir -p $genomeDir

if [ ! -f "$genomeDir/genomeGenerate.chkpt" ]; then
    STAR --runMode genomeGenerate --runThreadN $SLURM_NPROCS --limitGenomeGenerateRAM $limitRAM --genomeDir $genomeDir --genomeFastaFiles $INPUT --sjdbGTFfile $GFF --sjdbGTFtagExonParentTranscript gene_id

    touch $genomeDir/genomeGenerate.chkpt
else
    echo "index already exists."
    echo "Continue with step two"
fi



WD=S1_STAR_FILES_${REF_PREFIX}_DIR

mkdir -p $WD
mkdir -p CHKPNT_DIR

for i in $(ls *_R1.fq)
do
bs="${i%*_R1.fq}"

left_file=${bs}_R1.fq
right_file=${bs}_R2.fq


outFileNamePrefix=$WD/${bs}


if [ ! -f "CHKPNT_DIR/${bs}_STAR.chkpt" ]; then

    call="STAR --runMode alignReads --twopassMode Basic --runThreadN $SLURM_NPROCS --outSAMtype BAM SortedByCoordinate --genomeDir $genomeDir --readFilesIn $left_file $right_file --outFileNamePrefix $outFileNamePrefix"

    eval $call

    touch  CHKPNT_DIR/${bs}_STAR.chkpt

else
    echo "STAR for $bs already exists"

fi

done