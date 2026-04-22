#!/bin/bash
#SBATCH --job-name=Bwt2rnaseq
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00


EXPORT=/LUSTRE/apps/bioinformatica/bowtie2/bin/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/RSEM/bin/samtools-1.3/
export PATH=$PATH:$EXPORT

# Parameters for alignment

# @ Adjust parameters like allowed mismatches (--N) and seed length (--L) to optimize sensitivity and specificity
#  By default, Bowtie2 performs end-to-end alignment. For RNA-seq, you may want to use local alignment to allow soft-clipping of reads:
# Using --local to soft clipped the reads. The primary goal of soft-clipping is to improve alignment accuracy by ignoring low-quality or mismatched portions of the read that do not correspond to any part of the reference genome. This is particularly useful in RNA-seq data, where reads may contain non-coding regions or adapter sequences that do not align with the reference

# aligner_param is equivalent to --sensitive-local flag

aligner_params="--no-mixed --no-discordant -N 1 -L 20 -D 15 -R 2"
read_type="-q"
max_ins_size=500 # If paired-end include -X $max_ins_size


# 1) INDEX THE REFERENCE

REFERENCE=$1 # reference.fasta


REF_PREFIX=${REFERENCE%.*}

thread_count=$SLURM_NPROCS

mkdir -p INDEX

# Build Bowtie2 index if not already present
if [ ! -f "INDEX/$REF_PREFIX.1.bt2" ]; then
    bowtie2-build --threads $thread_count $REFERENCE INDEX/$REF_PREFIX
else
    echo "index '$REF_PREFIX' already exists."
fi

echo "Continue with step two"

WD=S1_BOWTIE2_BAM_FILES_${REF_PREFIX}_DIR

mkdir -p $WD
mkdir -p $WD/STATS

for i in $(ls *_R1.fq.gz)
do
bs="${i%*_R1.fq.gz}"


# Input files
# reads_file=${bs}.fq.gz

left_file=${bs}_R1.fq.gz
right_file=${bs}_R2.fq.gz

# Output files

bam_file=$WD/${bs}.sorted.bam

met_file=$WD/STATS/${bs}.met.txt
err_file=$WD/STATS/${bs}.stderr 
# # Test if the alignment was previously done!

if [ ! -f "$bam_file" ]; then
    
    echo "Aligning reads back to reference"

    # If merged or single end
    call="bowtie2 --local --met-file $met_file $aligner_params $read_type -x INDEX/$REF_PREFIX -U $reads_file -p $thread_count 2> $err_file | samtools view -@ $thread_count -F 4 -S -b | samtools sort -@ $thread_count -o $bam_file"

    # If paired-end
    #call="bowtie2 --local --met-file $met_file $aligner_params $read_type -x INDEX/$REF_PREFIX -1 $left_file -2 $right_file -p $thread_count 2> $err_file| samtools view -@ $thread_count -F 4 -S -b | samtools sort -n -@ $thread_count -o $bam_file"


    echo $call

    eval $call
    
    echo "Alignment $bam_file was done"

    echo "Summarizing stats for  $bam_file"

    # samtools_0.1.18 flagstat $bam_file > $WD/STATS/${bs}.flagstats.txt

    call="samtools flagstat $bam_file > $WD/STATS/${bs}.flagstats.txt"

    eval $call
    
    call="samtools stats $bam_file > $WD/STATS/${bs}.stats.txt"

    eval $call

    call="samtools depth $bam_file > $WD/STATS/${bs}.depth.txt"

    eval $call

else
    echo "File $bam_file already exists"

fi

#unlink $left_file
#unlink $right_file

done


samtools merge -@ $thread_count $WD/MERGED.bam $WD/*.sorted.bam

WDM=/LUSTRE/apps/Anaconda/2023/miniconda3/bin/

MULTIQCDIR=MULTIQC_VIZ_${REF_PREFIX}_DIR

mkdir -p $MULTIQCDIR

$WDM/multiqc $WD/STATS/* -o $MULTIQCDIR

# samtools merge -@ 12 *.sorted.bam | samtools sort -@ 12 -n -o MERGED.bam 

exit

# Axtell (Shorstacks )
#   - c_file : file path to a condensed FASTA file. Might be .gz compressed.

# bt_args = (['bowtie', input_option, '-p', str(args.threads), '-v', '1', '-k', '20', '-S', '--best', '--strata', '-x', args.genomefile, c_file]) 
# "score" set to zero, reflecting unknown alignments at this point mir_samfile = args.outdir + '/user_mir.sam' un_file = args.outdir + '/known_miRNAs_unaligned.fasta' 
# bt_args = ['bowtie', '-f', '-a', '-S', '-v', '0', '-p', str(args.threads), '--mapq', '0', '--un', un_file, '-x', args.genomefile, sane_known_miRNAs, mir_samfile]


thread_count=12
i=cam2_E_CKDL230016188-1A_H75JVDSX7_L3_merged.fq.gz
bowtie2 -q -p $thread_count -v 1 -k 20 -S --best --strata -x INDEX/$REF_PREFIX $i

# Calculate coverage
# Count number of records (of alignments ) by:
samtools view -c -F 256 cam2_E_CKDL230016188-1A_H75JVDSX7_L3.sorted.bam