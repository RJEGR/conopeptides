# DO NOT USE, ERRORS WITH CONOSORTER SCRIPTS
#!/bin/bash
#SBATCH -p cicese
#SBATCH --job-name=ConoSort
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 06-00:00:00

# The program requires only the assembled transcriptome or raw reads file, 
# in either DNA or amino acid format. ConoDictor 2 automatically recognizes the alphabet used.
NPROCS=$SLURM_NPROCS


# Run conosorter

FASTA=$1

export PATH=/LUSTRE/apps/bioinformatica/hmmer-3.3.2/bin:$PATH
export PATH=/LUSTRE/apps/bioinformatica/ConoSorter_v1.1:$PATH

ConoSorter -d $FASTA

exit