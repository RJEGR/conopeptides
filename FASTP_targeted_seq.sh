#!/bin/bash
#SBATCH --job-name=FASTP
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

FASTQC=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software
export PATH=$PATH:$FASTQC

NPROCS=$SLURM_NPROCS


mkdir -p MULTIQC_VIZ_DIR

mkdir -p FASTP_OUT_DIR

mkdir -p CHKPNT_DIR

for f in $(ls *.fastq)
do
basename=${f##*/}
bs="${basename%_1.fastq}"
infile="${f%_1.fastq}"

left_file=${infile}_1.fastq
right_file=${infile}_2.fastq

# this configure able the Data Filtration and Initial Assembly step from Phuong, M et al 2019

if [ ! -f "CHKPNT_DIR/${bs}_fastp.chkpt" ]; then

    call="fastp --thread $NPROCS --detect_adapter_for_pe \
    --json MULTIQC_VIZ_DIR/${bs}_fastp.json \
    --html MULTIQC_VIZ_DIR/${bs}_fastp.html \
    -i $left_file -I $right_file \
    --length_required 36
    --merge --merged_out FASTP_OUT_DIR/${bs}_merged.fq.gz \
    --dedup --dup_calc_accuracy 3 \
    --low_complexity_filter --complexity_threshold 60 \
    -o FASTP_OUT_DIR/${bs}_R1.fq.gz  -O FASTP_OUT_DIR/${bs}_R2.fq.gz"

    echo "Running fastp in '${bs}' sample."

    echo $call

    eval $call

    unlink $left_file
    unlink $right_file

    touch  CHKPNT_DIR/${bs}_fastp.chkpt

else
    echo "'CHKPNT_DIR/${bs}_fastp.chkpt' already exists."
    echo "Continue with next sample."
fi


done

WDM=/LUSTRE/apps/Anaconda/2023/miniconda3/bin/

$WDM/multiqc MULTIQC_VIZ_DIR/*_fastp.json  -o MULTIQC_VIZ_DIR --config multiqc_info.conf


exit