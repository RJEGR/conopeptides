#!/bin/bash
#SBATCH --job-name=Conotrace
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=rgomez@uabc.edu.mx

# Continue w/ /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/californicus/01.Preprocess/Conotrace

EXPORT=/LUSTRE/apps/bioinformatica/diamond_v2.1.8/
export PATH=$PATH:$EXPORT

NPROCS=$SLURM_NPROCS

DBDIR=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/DataBases/Conopeptides

REF=$DBDIR/conoserver_protein.fa

BS=${REF##*/}

QUERY=$1

output1=${QUERY%.*}_vs_${BS%.*}.diamond.blastx.outfmt6
output2=${QUERY%.*}_vs_${BS%.*}.diamond.blastx.fastq

START_TIME=$SLURM_JOB_START_TIME


if [ ! -f "${REF%.*}.dmnd" ]; then

    diamond makedb --in $REF --db ${REF%.*}

else
     echo "${REF%.*}.dmnd already exists."
fi

# diamond blastx -d ${REF%.*} -q $QUERY -p $NPROCS -k 1 -e 1e-5 -o $output --outfmt 6 

diamond blastx -d ${REF%.*} -q $QUERY -p $NPROCS -k 1 -e 1e-5 -o $output1 --outfmt 6 qseqid evalue  --alfmt fastq --al $output2

END_TIME=$SLURM_JOB_END_TIME

echo "Start job: $START_TIME"
echo "End job: $END_TIME"

exit

exit