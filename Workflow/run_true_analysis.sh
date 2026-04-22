#!/bin/bash
#SBATCH --job-name=TrueSet
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

# sbatch run_true_analysis.sh randomize_samples.txt

#1)  Test Assembly
# 1.1) spades

EXPORT=/LUSTRE/apps/bioinformatica/SPAdes-3.15.5-Linux/bin/
export PATH=$PATH:$EXPORT

# for i in $(seq 1 5);do call="python randomize_paired.py conoserver_nucleic_californicus_length_100.fa randomize_set${i} 1000000 100 --phred_score 40"; eval $call; done

# for file in $(ls *_1.fastq); do basename=${file##*/}; base="${basename%*_1.fastq}"; bs="${base}"; fwr=${base}_1.fastq; rev=${base}_2.fastq; echo "$bs" "$bs" `printf "$PWD/$fwr"` `printf "$PWD/$rev"`; done > randomize_samples.txt
# for file in $(ls *_1.fastq); do basename=${file##*/}; base="${basename%*_1.fastq}"; bs="${base}"; fwr=${base}_1.fastq; rev=${base}_2.fastq; echo "$bs" "$bs" `printf "$PWD/$fwr"`; done > randomize_samples.txt

CPU=$SLURM_NPROCS
MEM=$SLURM_MEM_PER_NODE

Manifest=$1 # samples.txt

left_reads=`awk '{print "--pe1-1",$3}' $Manifest`
right_read=`awk '{print "--pe1-2", $4}' $Manifest`

single_reads=`awk '{print "--s" NR, $3}' $Manifest`

OUTDIR=${Manifest%.txt}_spades_dir

mkdir -p $OUTDIR


if [ ! -f "1_spades.chkp" ]; then

    #call="rnaspades.py $left_reads $right_read --pe1-fr -t $CPU -m $MEM -o $OUTDIR"

    call="rnaspades.py $single_reads -t $CPU -m $MEM -o $OUTDIR"

    eval $call

    ln -s ${OUTDIR}/transcripts.fasta  Spades.fasta

    touch 1_spades.chkp

else
    echo "..."
fi


# 1.2) Trinity

module load trinityrnaseq-v2.15.1

OUTDIR=${Manifest%.txt}_trinity_dir

mkdir -p $OUTDIR

if [ ! -f "1_trinity.chkp" ]; then

    # call="Trinity --seqType fq --max_memory 100G --samples_file $Manifest --no_normalize_reads --CPU $CPU --output $OUTDIR"
    
    call="Trinity --seqType fq --max_memory 100G --samples_file $Manifest --no_normalize_reads --CPU $CPU --output $OUTDIR"
    
    eval $call

    ln -s ${OUTDIR}.Trinity.fasta Trinity.fasta

    touch 1_trinity.chkp

else
    echo "..."
fi


# From raw assembly match the false/true number of nucleide sequences

# Concat assemblies and map reads back (individually) (run supertranscript.sh script)

cat Spades.fasta Trinity.fasta > all_assemblies.fasta

sbatch SuperTranscript.sh all_assemblies.fasta 

# Insilico prediction of orfs (also run for conoserver_nucleic.fa.gz) using transdecoderPredict.sh


# ConoSorter

# Count the number of conotoxin used to create the simulated data:

EXPORT=/LUSTRE/apps/bioinformatica/ncbi-blast-2.14.0+/bin/
export PATH=$PATH:$EXPORT


QUERY=${Manifest%.txt}_sets_1.fastq.tmp

call=`awk '{print $3}' $Manifest | tr "\n" " "`

cat $call > $QUERY

# Convert FASTQ to FASTA using awk
awk 'NR%4==1 {print ">" substr($0, 2)} NR%4==2 {print}' $QUERY > ${QUERY%.fastq.tmp}.fasta

QUERY=${QUERY%.fastq}.fasta

rm *fastq.tmp

DBDIR=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/californicus/true_data_dir/INDEX
REF=$DBDIR/conoserver_nucleic_californicus_length_100.fa
BS=${REF##*/}

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

exit