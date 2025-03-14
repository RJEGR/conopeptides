#!/bin/bash
#SBATCH --job-name=Coverage
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=rgomez@uabc.edu.mx


EXPORT=/LUSTRE/apps/bioinformatica/hisat2/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/RSEM/bin/samtools-1.3
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/stringtie/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/apps/gffread
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/TransDecoder-v5.7.0/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/subread-2.0.6-Linux-x86_64/bin
export PATH=$PATH:$EXPORT

# Vars

REFERENCE=$1 # Assembly file

REF_PREFIX=`basename ${REFERENCE%.f*}`

CPU=$SLURM_NPROCS


# 1) Create index if does not exist

mkdir -p INDEX

if [ ! -f "INDEX/${REF_PREFIX}.1.ht2" ]; then
hisat2-build -p $CPU $REFERENCE INDEX/$REF_PREFIX

else
    echo "index for '$REFERENCE' already exists."
fi

# 2) Align reads back to the reference

WD=`pwd`
WD=S1_HISAT2_SAM_BAM_FILES_${REF_PREFIX}_DIR

mkdir -p $WD

mkdir -p $WD/STATS

for i in $(ls *.fq.gz)
do
bs="${i%*.fq.gz}"

bam_file=$WD/${bs}.sorted.bam

if [ ! -f "$WD/${bs}.sorted.bam.chkpt" ]; then

merged_reads=${bs}.fq.gz

    hisat2  --phred33 --dta -p $CPU \
        -x INDEX/$REF_PREFIX -U $merged_reads \
        --rg-id=${bs} --rg SM:${bs} -S $WD/${bs}.sam \
        --summary-file $WD/${bs}.summary.txt
        #--met-file $WD/${bs}.hisat.met.txt

    samtools sort -@ $CPU -o $bam_file $WD/${bs}.sam

    touch $WD/${bs}.sorted.bam.chkpt

    call="samtools flagstat $bam_file > $WD/STATS/${bs}.flagstats.txt"

    eval $call
    
    call="samtools stats $bam_file > $WD/STATS/${bs}.stats.txt"

    eval $call

    call="samtools depth $bam_file > $WD/STATS/${bs}.depth.txt"

    eval $call

else 
    echo "'$WD/${bs}.sorted.bam' already exists."
    echo "Continue with sorting gtf."
fi

WD2=`pwd`
WD2=S2_STRINGTIE_DENOVO_MODE_${REF_PREFIX}_DIR


mkdir -p $WD2

if [ ! -f "$WD2/${bs}.gtf.chkpt" ]; then

    stringtie --rf -p $CPU -l $bs -o $WD2/${bs}.gtf $WD/${bs}.sorted.bam

    touch $WD2/${bs}.gtf.chkpt

else
 echo "$WD2/${bs}.gtf.chkpt already exists"
fi

done

WDM=/LUSTRE/apps/Anaconda/2023/miniconda3/bin/

MULTIQCDIR=MULTIQC_VIZ_${REF_PREFIX}_DIR

mkdir -p $MULTIQCDIR

$WDM/multiqc $WD/*.summary.txt -o $MULTIQCDIR

# Merging

echo "Continue with next sample."

ls -d -1 $WD2/*.gtf > ${REF_PREFIX}_stringtie_gtf_list.txt

if [ ! -f "${REF_PREFIX}_transcripts.gtf" ]; then
 stringtie --rf --merge -p $CPU -o ${REF_PREFIX}_transcripts.gtf ${REF_PREFIX}_stringtie_gtf_list.txt

else
    echo "Merging step already exists. Continue w/ Feature Count"
fi

# 2) Count 
# ls -d -1 $WD/*.sorted.bam > ${REF_PREFIX}_sorted_bam_list.txt

if [ ! -f "${REF_PREFIX}_counts.txt" ]; then
 featureCounts -T $CPU -a ${REF_PREFIX}_transcripts.gtf -o $WD2/${REF_PREFIX}_counts.txt $WD/${bs}.sorted.bam

else
    echo "${REF_PREFIX}_counts.txt file already exists"
fi


exit
