#!/bin/bash
#SBATCH --job-name=evals
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00


# 1) using blastn, blat and detonate, evaluate the assembly generated from the randomized paired-end reads


EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/detonate-master/ref-eval
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/RSEM/bin/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/RSEM/bin/samtools-1.3/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/bowtie2/bin/
export PATH=$PATH:$EXPORT

# 1.1) Detonate:
# From now on, we will call the estimated "true" assembly or the
# collection of full-length reference sequences (whichever you choose
# to use) the reference. Let's assume that the assembly of interest is
# in A.fa, and the reference is in B.fa.

contigs_a=$1 # Assembly of interest A.fa
true_ref_b=$2 # Reference B.fa
Manifest=$3

# Concat reads prior to map back to the assemblt and reference

fqreads=${Manifest%.txt}_sets_1.fastq

call=`awk '{print $3}' $Manifest | tr "\n" " "`

#cat $call > $fqreads

cat $call | awk 'NR%4==1 {print "@new_id_" ++i} NR%4!=1 {print}' > $fqreads

blatout1=${contigs_a}_to_${true_ref_b}.psl
blatout2=${true_ref_b}_to_${contigs_a}.psl

blat -minIdentity=80 $true_ref_b $contigs_a $blatout1
blat -minIdentity=80 $contigs_a $true_ref_b $blatout2

A_ref=${contigs_a%.*}_prepref
B_ref=${true_ref_b%.*}_prepref

A_expr=${contigs_a%.*}
B_expr=${true_ref_b%.*}
 
rsem-prepare-reference $contigs_a $A_ref
rsem-prepare-reference $true_ref_b $B_ref

rsem-calculate-expression --bowtie2 -p $SLURM_NPROCS --no-bam-output $fqreads $A_ref $A_expr

rsem-calculate-expression --bowtie2 -p $SLURM_NPROCS --no-bam-output $fqreads $B_ref $B_expr

ref-eval --scores=nucl,pair,contig,kmer,kc \
              --weighted=both \
              --A-seqs $contigs_a \
              --B-seqs $true_ref_b \
              --A-expr ${A_expr}.isofr.results \
              --B-expr ${B_expr}.isofr.results \
              --A-to-B $blatout1 \
              --B-to-A $blatout2 \
              --num-reads 5000000 \
              --readlen 76 \
              --kmerlen 76 \
              | tee ${contigs_a}_and_${true_ref_b}.scores.txt


 # 1.2) Tranrate:

TRANSRATE=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/bin
export PATH=$PATH:$TRANSRATE

export PATH=/LUSTRE/apps/bioinformatica/.local/bin:$PATH
export PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/bin:$PATH
export LD_LIBRARY_PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/lib:$LD_LIBRARY_PATH


output=${contigs_a%.*}_to_${true_ref_b%.*}_transrate_dir

transrate --assembly $contigs_a --reference $true_ref_b --threads $SLURM_NPROCS --output $output


exit

# 1.3) BLASTN
EXPORT=/LUSTRE/apps/bioinformatica/ncbi-blast-2.14.0+/bin/
export PATH=$PATH:$EXPORT


DBDIR=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/californicus/true_data_dir/INDEX
REF=$DBDIR/conoserver_nucleic_californicus_length_100.fa
BS=${REF##*/}

QUERY=$fqreads

outfmt=${QUERY%.*}_vs_${BS%.*}.blastn.outfmt6
alfmt=${QUERY%.*}_vs_${BS%.*}.blastn.fastq
chkpt=${QUERY%.*}_vs_${BS%.*}.blastn.chkpt

if [ ! -f "$chkpt" ]; then

    #call="diamond blastx -d ${REF%.*} -q $QUERY -p $NPROCS --min-orf 1 -k 1 -e 1e-3 --outfmt 6 qseqid evalue pident -o $outfmt --alfmt fastq --al $alfmt"
    
    # If nucleotide 

    call="blastn -db ${REF%.*} -query $QUERY -num_threads $CPU -evalue 1e-3 -outfmt 6 -out $outfmt"
   
    echo $call

    eval $call

    touch  $chkpt

else
    echo "Continue with next sample."
fi
    
done