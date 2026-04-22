#!/bin/bash

#SBATCH -p cicese
#SBATCH --job-name=Conotrace
#SBATCH -N 2
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 06-00:00:00


EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/FastQ-Screen-0.16.0/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/bowtie2/bin/
export PATH=$PATH:$EXPORT

NPROCS=$SLURM_NPROCS

mkdir -p CONOTRACE_DIR
mkdir -p CHKPNT_DIR


for f in $(ls *_merged.fq.gz)
do
basename=${f##*/}
bs="${basename%_merged.fq.gz}"


if [ ! -f "CHKPNT_DIR/${bs}_trace.chkpt" ]; then

    call="fastq_screen --threads $NPROCS --conf fastq_screen.conf --aligner bowtie2 --subset 0 --filter 33 --tag --force --outdir CONOTRACE_DIR $f"

    echo "Running\n."

    echo $call

    eval $call

else
    echo "'${bs}_trace.chkpt' already exists."
fi

    touch  CHKPNT_DIR/${bs}_trace.chkpt


done

exit

# --aligner bowtie --subset 100000”. The “one hit/one genome” and “multiple hits/one genome” reads were used for the species analysis
# _merged.fq.gz
# ABC (genome order)
#filter --3

# Prepare indexes 
cd /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/DataBases/Conopeptides/
mkdir -p INDEX
REFERENCE=all_definedseqs_v6.fa
bowtie2-build --threads 24 $REFERENCE INDEX/${REFERENCE%.*}