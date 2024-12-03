#!/bin/bash
#SBATCH --job-name=Bowtie2ref
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


for i in $(ls *.fq.gz)
do
bs="${i%*.fq.gz}"

reads_file=${bs}.fq.gz
#If paired-end use -1 $left_file -2 $right_file
#left_file=${bs}_R1.fq.gz
#right_file=${bs}_R2.fq.gz

bam_file=${bs}.sorted.bam
met_file=${bs}.met.txt

# # Test if the alignment was previously done!

if [ ! -f "$WD/$bam_file" ]; then
    
    echo "Aligning reads back to reference"

    call="bowtie2 --local --met-file $met_file $aligner_params $read_type \
    -x INDEX/$REF_PREFIX -U $reads_file -p $thread_count 2> ${bs}.stderr | \
    samtools view -@ $thread_count -F 4 -S -b | samtools sort -@ $thread_count -n -o $WD/$bam_file"

    echo $call

    eval $call


else
    echo "File $bam_file already exists in $WD directory"

fi

# unlink $reads_file

done

mkdir -p stats
mv *.met.txt stats
mv *.stderr stats
mv stats $WD

exit