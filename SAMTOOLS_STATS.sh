#!/bin/bash
#SBATCH --job-name=stats
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

thread_count=$SLURM_NPROCS


EXPORT=/LUSTRE/apps/bioinformatica/RSEM/bin/samtools-1.3/samtools
export PATH=$PATH:$EXPORT

for f in $(ls *.bam)
do
bs="${f%*.bam}"

input_bam=${bs}.bam

bam_file=${bs}.tmp.bam

touch_file=${bs}.chkpt



if [ ! -f "$touch_file" ]; then

    samtools sort -@ $thread_count -n -o $input_bam > $bam_file

    call="samtools flagstat $bam_file > ${bs}.flagstats.txt"

    eval $call
    
    call="samtools stats $bam_file > ${bs}.stats.txt"

    eval $call

    call="samtools depth $bam_file > ${bs}.depth.txt" #  will contain three columns: chromosome, base position, and read depth at each position.

    eval $call


    touch $touch_file

    
else
    echo "File $touch_file already exists"

fi

done