#!/bin/bash
#SBATCH --job-name=Conotrace
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

# All-in-one workflow to trimm, filter lowQ and conotoxin reads from Paired-fastq files

FASTQC=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software
export PATH=$PATH:$FASTQC

EXPORT=/LUSTRE/apps/bioinformatica/diamond_v2.1.8/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/ncbi-blast-2.14.0+/bin/
export PATH=$PATH:$EXPORT

NPROCS=$SLURM_NPROCS


mkdir -p MULTIQC_VIZ_DIR

mkdir -p FASTP_OUT_DIR

mkdir -p CHKPNT_DIR

for f in $(ls *_1.fq.gz)
do
basename=${f##*/}
bs="${basename%_1.fq.gz}"
infile="${f%_1.fq.gz}"

left_file=${infile}_1.fq.gz
right_file=${infile}_2.fq.gz

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

# Filtering step using Diamond agains conoserver

QUERY=FASTP_OUT_DIR/${bs}_merged.fq.gz
DBDIR=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/DataBases/Conopeptides
REF=$DBDIR/conoserver_protein.fa
BS=${REF##*/}

outfmt=${QUERY%.*}_vs_${BS%.*}.diamond.blastx.outfmt6
alfmt=${QUERY%.*}_vs_${BS%.*}.diamond.blastx.fastq

if [ ! -f "CHKPNT_DIR/${bs}_diamond.chkpt" ]; then

    #call="diamond blastx -d ${REF%.*} -q $QUERY -p $NPROCS --min-orf 1 -k 1 -e 1e-3 --outfmt 6 qseqid evalue pident -o $outfmt --alfmt fastq --al $alfmt"
    
    # If nucleotide 

    call="blastn -db ${REF%.*} -query $QUERY -num_threads $NPROCS -e 1e-3 --outfmt 6 -o $outfmt"
   
    echo $call

    eval $call

    touch  CHKPNT_DIR/${bs}_diamond.chkpt

else
    echo "Creating ${bs}_diamond.chkpt"
    echo "Continue with next sample."
fi
    
done

WDM=/LUSTRE/apps/Anaconda/2023/miniconda3/bin/

$WDM/multiqc MULTIQC_VIZ_DIR/*_fastp.json  -o MULTIQC_VIZ_DIR --config multiqc_info.conf --force


exit

# improve w/ https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options
# https://gensoft.pasteur.fr/docs/diamond/0.9.19/diamond_manual.pdf

--min-orf                this is the parameter to available short conopeptides in the search
--min-orf   Ignore translated sequences that do not contain an open reading frame of at least this length. 
 By default this feature is disabled for sequences of length below 30, 
 set to 20 for sequences of length below 100, and set to 40 otherwise. Setting this option to 1 will disable this feature.
--seed-cut               cutoff for seed complexity
