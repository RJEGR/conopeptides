#!/bin/bash
#SBATCH --job-name=detonate
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00


EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/detonate-master/rsem-eval
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/bowtie2/bin/
export PATH=$PATH:$EXPORT


#  RSEM-EVAL is a reference-free evaluation method based on a novel probabilistic model that depends only on an assembly and the RNA-Seq reads used for its construction.

assembly_fasta_f=$1 # fasta file input.fasta is a multi-FASTA file contains all transcript sequences used to learn the true transcript length distribution's parameters.
upstream_read_file=$2
downstream_read_file=$3

sample_name=${assembly_fasta_f%.*}

#param_dir=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/detonate-master/rsem-eval/true_transcript_length_distribution
#parameter_file=$param_dir/conoserver_jan_2025.txt

parameter_file=${assembly_fasta_f%.*}_length_distribution.txt

rsem-eval-estimate-transcript-length-distribution $assembly_fasta_f $parameter_file

rsem-eval-calculate-score --bowtie2 -p $SLURM_NPROCS --transcript-length-parameters $parameter_file --paired-end $upstream_read_file $downstream_read_file $assembly_fasta_f  $sample_name 150


exit

# L represents the average fragment length. It should be a positive integer (real value will be rounded to the nearest integer).
#rsem-eval-calculate-score [options] upstream_read_file(s) assembly_fasta_file sample_name L
#     rsem-eval-calculate-score [options] --paired-end upstream_read_file(s) downstream_read_file(s) assembly_fasta_file sample_name L
#     rsem-eval-calculate-score [options] --sam/--bam [--paired-end] input assembly_fasta_file sample_name L


bamfile=cam2_E_CKDL230016188-1A_H75JVDSX7_L3.sorted.bam

rsem-eval-calculate-score --bam $bamfile $assembly_fasta_f $sample_name 150