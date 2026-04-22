#!/bin/bash
#SBATCH --job-name=Diamond
#SBATCH -N 1
#SBATCH --mem=50GB
#SBATCH --ntasks-per-node=20
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=rgomez@uabc.edu.mx

EXPORT=/LUSTRE/apps/bioinformatica/diamond_v2.1.8/
export PATH=$PATH:$EXPORT

NPROCS=$SLURM_NPROCS

DBDIR=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/DataBases/Conopeptides

REF=$DBDIR/conoserver_protein.fa

BS=${REF##*/}

QUERY=$1

output=${QUERY%.*}_vs_${BS%.*}.diamond.blastx.outfmt6

START_TIME=$SLURM_JOB_START_TIME


if [ ! -f "${REF%.*}.dmnd" ]; then

    diamond makedb --in $REF --db ${REF%.*}

else
     echo "${REF%.*}.dmnd already exists."
fi

diamond blastx -d ${REF%.*} -q $QUERY -p $NPROCS -k 1 -e 1e-5 -o $output --outfmt 6

END_TIME=$SLURM_JOB_END_TIME

echo "Start job: $START_TIME"
echo "End job: $END_TIME"

exit

seqkit grep -f IDs.txt read_1.fq.gz -o dir/read_1.fq.gz
