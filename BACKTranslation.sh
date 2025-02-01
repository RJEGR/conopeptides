#!/bin/bash
#SBATCH --job-name=backtranseq
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=rgomez@uabc.edu.mx

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/EMBOSS-6.6.0/emboss
export PATH=$PATH:$EXPORT

WD=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/DataBases/Conopeptides

cd $WD

#backtranambig -sequence conoserver_protein.curated.fa -outfile conoserver_protein.curated_backtranambig.fa

backtranseq -sequence conoserver_protein.curated.fa -outfile conoserver_protein.curated_backtranseq.fa -sprotein1 -auto -stdout

exit

 -warning  -error   -fatal  -die