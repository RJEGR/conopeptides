#!/bin/bash
#SBATCH -p cicese
#SBATCH -t 6-00:00:00
#SBATCH --job-name=splitSpades
#SBATCH -N 2
#SBATCH --mem=80GB
#SBATCH --ntasks-per-node=20
#SBATCH --error=slurm-%j.err

# after run split_samples.pl, then

EXPORT=/LUSTRE/apps/bioinformatica/SPAdes-3.15.5-Linux/bin/
export PATH=$PATH:$EXPORT

CPU=$SLURM_NPROCS
MEM=$SLURM_MEM_PER_NODE
DIR=$SLURM_SUBMIT_DIR

mkdir -p SPADES_CHKPNT_DIR

for f in $(ls samples_part_*.txt)
do

chkpt_file=${f%.txt}.chkpt

if [ ! -f "SPADES_CHKPNT_DIR/$chkpt_file" ]; then

Manifest=$f

left_reads=`awk '{print "--pe1-1",$3}' $Manifest`
right_read=`awk '{print "--pe1-2", $4}' $Manifest`

OUTDIR=${f%.txt}_spades_dir

mkdir -p $OUTDIR

   call="rnaspades.py $left_reads $right_read --pe1-fr -t $CPU -m $MEM -o $DIR/$OUTDIR"

   echo "Running spades for list of samples: '${f%.txt}'."

   echo $call

   #eval $call

   touch  SPADES_CHKPNT_DIR/$chkpt_file

else
    echo "'SPADES_CHKPNT_DIR/$chkpt_file' already exists."
    echo "Continue with next list of samples."
fi

done

exit