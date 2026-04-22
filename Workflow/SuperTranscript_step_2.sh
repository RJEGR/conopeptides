#!/bin/bash
#SBATCH --job-name=supertr
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=rgomez@uabc.edu.mx
#SBATCH -t 6-00:00:00

# After running BBMap.slurm, the script will run Corset and Lace

####
# Run Corset
####

REFERENCE=$1 # Reference 
REF_PREFIX=${REFERENCE%.*}

WDD=S2_CORSET_FILES_${REF_PREFIX}_DIR

mkdir -p $WDD

EXPORT=/LUSTRE/apps/bioinformatica/corset-1.09-linux64/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/RSEM/bin/samtools-1.3/
export PATH=$PATH:$EXPORT


for f in $(ls *.bam)
do
basename=${f##*/}
bs="${basename%.bam}"

if [ ! -f "${bs}.corset-reads.ok" ]; then

    #call="samtools sort -@ $SLURM_NPROCS -o ${bs}.sorted.bam ${bs}.bam"

    #eval $call

    call="corset -f true -r true ${bs}.sorted.bam"

    echo "Running corset in '${bs}' sample."

    echo $call

    eval $call

    touch ${bs}.corset-reads.ok


else
    echo "'${bs}.corset-reads' already exists."
    echo "Continue with next sample."
fi

    echo "Continue with next step."

done

corset  -f true -i corset *.corset-reads

mv *corset-reads $WDD
mv *txt $WDD

# Run Lace


module load conda-2024-2
source activate lace
export PATH=/LUSTRE/apps/bioinformatica/Lace-1.14.1/Lace:$PATH

OUTDIR=S3_LACE_${REFERENCE%.*}_DIR

mkdir -p $OUTDIR

Lace $REFERENCE $WDD/clusters.txt -o $OUTDIR

exit