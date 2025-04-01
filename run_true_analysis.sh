#!/bin/bash
#SBATCH --job-name=TrueSet
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

#1)  Test Assembly
# 1.1) spades

EXPORT=/LUSTRE/apps/bioinformatica/SPAdes-3.15.5-Linux/bin/
export PATH=$PATH:$EXPORT

# for file in $(ls *_1.fasta); do basename=${file##*/}; base="${basename%*_1.fasta}"; bs="${base}"; fwr=${base}_1.fasta; rev=${base}_2.fasta; echo "$bs" "$bs" `printf "$PWD/$fwr"` `printf "$PWD/$rev"`; done > samples.txt

CPU=$SLURM_NPROCS
MEM=$SLURM_MEM_PER_NODE

Manifest=$1 # samples.txt

left_reads=`awk '{print "--pe1-1",$3}' $Manifest`
right_read=`awk '{print "--pe1-2", $4}' $Manifest`

OUTDIR=${Manifest%.txt}_spades_dir

mkdir -p $OUTDIR


if [ ! -f "1_spades.chkp" ]; then

    call="rnaspades.py $left_reads $right_read --pe1-fr -t $CPU -m $MEM -o $DIR/$OUTDIR"

    eval $call

    touch 1_spades.chkp

else
    echo "..."
fi


# 1.2) Trinity

module load trinityrnaseq-v2.15.1

OUTDIR=${Manifest%.txt}_trinity_dir

mkdir -p $OUTDIR

if [ ! -f "1_trinity.chkp" ]; then

    call="Trinity --seqType fa --max_memory 100G --samples_file $Manifest --no_normalize_reads --CPU $CPU --output $OUTDIR"

    eval $call

    touch 1_trinity.chkp

else
    echo "..."
fi




# From raw assembly match the false/true number of nucleide sequences

# Concat assemblies and map reads back (individually) (run supertranscript.sh script)



# Insilico prediction of orfs (also run for conoserver_nucleic.fa.gz) using transdecoderPredict.sh


exit