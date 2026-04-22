#!/bin/sh
## Directivas
#SBATCH --job-name=Trinotate
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20

FILE=$1 #rnaspades.fasta
# prefix=`basename ${FILE%.f*}`

EXPORT=/LUSTRE/apps/bioinformatica/TransDecoder-v5.7.0/
export PATH=$PATH:$EXPORT

module load conda-2024

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/eggnog-mapper-master
export PATH=$PATH:$EXPORT

emapper.py -i rnaspades.fasta.transdecoder.pep --cpu 20 --data_dir /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Trinotate/TRINOTATE_DB//EGGNOG_DATA_DIR -o eggnog_mapper --override

