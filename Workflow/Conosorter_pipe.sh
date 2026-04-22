#!/bin/bash
#SBATCH -p cicese
#SBATCH --job-name=CnSrt
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 06-00:00:00

export PATH=/LUSTRE/apps/bioinformatica/hmmer-3.3.2/bin:$PATH

NPROCS=$SLURM_NPROCS

file=$1

B=`basename $file .fasta`


ConoSorterDir=/LUSTRE/apps/bioinformatica/ConoSorter_v1.1

echo -e "Input file contains DNA sequences:\n\nStep 1/3: Formatting & Translating..."
${ConoSorterDir}/transl $file
## REGEX ##
echo "Step 2/3: Matching regular expressions..."
${ConoSorterDir}/regex_d translated.temp

## pHMM #######
echo "Step 3/3: Matching Markov models..."
if [ -f HMMER_input_superfam.fasta ]
   then
        hmmscan --cpu $NPROCS --tblout HMMER_output_superfam.temp --notextw ${ConoSorterDir}/models/HMM_superfam_db HMMER_input_superfam.fasta &> /dev/null
        ${ConoSorterDir}/pHMM HMMER_output_superfam.temp | sort -t'|' -k 1 > sorted_pHMM_superfam.temp
        ${ConoSorterDir}/pHMM HMMER_output_superfam.temp | sort -t'|' -k 1 > sorted_pHMM_superfam.temp
        else
        > sorted_pHMM_superfam.temp
        fi

if [ -f HMMER_input_class.fasta ]
    then
    #hmmscan --tblout HMMER_output_class.temp --notextw ./models/HMM_class_db HMMER_input_class.fasta &> /dev/null
    hmmscan --cpu $NPROCS  --tblout HMMER_output_class.temp --notextw ${ConoSorterDir}/models/HMM_class_db HMMER_input_class.fasta &> /dev/null
    ${ConoSorterDir}/pHMM HMMER_output_class.temp | sort -t'|' -k 1 > sorted_pHMM_class.temp
    else
    > sorted_pHMM_class.temp
    fi

        join -t'|' -a1 -a2 -1 1 -2 1 sorted_pHMM_superfam.temp sorted_pHMM_class.temp | ${ConoSorterDir}/layout


        ## Rename output files ####

        mv Regex.tab ${B}_Regex.tab
        mv pHMM.tab ${B}_pHMM.tab
