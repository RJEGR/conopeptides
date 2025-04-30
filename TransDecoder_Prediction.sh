#!/bin/sh
## Directivas
#SBATCH -t 6-00:00:00
#SBATCH --error=slurm-%j.err
#SBATCH --job-name=Decoder
#SBATCH -N 2
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20

# Trinotate

EXPORT=/LUSTRE/apps/bioinformatica/TransDecoder-v5.7.0/
export PATH=$PATH:$EXPORT


EXPORT=/LUSTRE/apps/bioinformatica/diamond_v2.1.8/
export PATH=$PATH:$EXPORT


EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/eggnog-mapper-master
export PATH=$PATH:$EXPORT

EXPORT=/home/rgomez/.local/bin #signalP
export PATH=$PATH:$EXPORT

FILE=$1 #transcriptome.fasta

bs=`basename ${FILE%.f*}`

output_dir=${bs}.transdecoder_dir

if [ ! -f "${bs}.TransDecoder.LongOrfs.chkp" ]; then
    call="TransDecoder.LongOrfs -t $FILE --output_dir $output_dir --complete_orfs_only -m 50"
    
    eval $call

    touch ${bs}.TransDecoder.LongOrfs.chkp

else
    echo "TransDecoder.LongOrfs for '$FILE' already exists."
fi


ln -s ${output_dir}/longest_orfs.pep ${bs}_longest_orfs.pep


# if not homology evidence is provided:


if [ ! -f "${bs}.TransDecoder.Predict.chkp" ]; then

    call="TransDecoder.Predict -t $FILE --cpu $SLURM_NPROCS --output_dir $output_dir"
    
    #echo $call
    
    # eval $call # omit 

    touch ${bs}.TransDecoder.Predict.chkp

else
    echo "TransDecoder.Predict for '$FILE' already exists."
fi

# Generate homology evidence for all _longest_orfs.pep :

# Evidence of known conopeptides
DBDIR=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/DataBases/Conopeptides

REF1=$DBDIR/conoserver_protein_curated.fa

REF2=$DBDIR/uniprotkb_taxonomy_id_33208_AND_cc_tiss_2025_02_19.fasta


BS1=${REF1##*/}
BS2=${REF2##*/}

QUERY=${bs}_longest_orfs.pep

output1=${QUERY%.*}_vs_${BS1%.*}.diamond.blastp.outfmt6
output2=${QUERY%.*}_vs_${BS2%.*}.diamond.blastp.outfmt6

chkp1=${QUERY%.*}_vs_${BS1%.*}.diamond.blastp.ok
chkp2=${QUERY%.*}_vs_${BS2%.*}.diamond.blastp.ok

if [ ! -f "${REF1%.*}.dmnd" ]; then
    diamond makedb --in $REF1 --db ${REF1%.*}
else
     echo "${REF1%.*}.dmnd already exists."
fi


if [ ! -f "${REF2%.*}.dmnd" ]; then
    diamond makedb --in $REF2 --db ${REF2%.*}
else
     echo "${REF2%.*}.dmnd already exists."
fi

# Run blastp individually
# 1)
if [ ! -f "$chkp1" ]; then

    diamond blastp --sensitive -d ${REF1%.*} -q $QUERY -p $SLURM_NPROCS -k 3 -e 1e-3 -o $output1 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp

    touch $chkp1

else
    echo "$output1 already exists."
fi

# 2)

if [ ! -f "$chkp2" ]; then

    diamond blastp --sensitive -d ${REF2%.*} -q $QUERY -p $SLURM_NPROCS -k 3 -e 1e-3 -o $output2 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp

    touch $chkp2

else
    echo "$output2 already exists."
fi

# Evidence of other proteins for all _longest_orfs.pep :

module load conda-2024

# 3)

eggnog_mapper=${bs}_eggnog_mapper

outfmt_evidence=${bs}_blastp_evidence.hits


if [ ! -f "${bs}.eggnog_mapper.chkp" ]; then

    emapper.py -i $QUERY --itype proteins --cpu $SLURM_NPROCS --data_dir /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Trinotate/TRINOTATE_DB//EGGNOG_DATA_DIR -o $eggnog_mapper --override

    touch "${bs}.eggnog_mapper.chkp"

else
    echo "Emapper already exists."
fi

cat ${eggnog_mapper}.emapper.hits $output1 $output2 > $outfmt_evidence


if [ ! -f "${bs}.TransDecoder.homology.Predict.chkp" ]; then

    call="TransDecoder.Predict -t $FILE --cpu $SLURM_NPROCS --retain_blastp_hits $outfmt_evidence --output_dir $output_dir"
    
    #echo $call
    
    eval $call

    touch ${bs}.TransDecoder.homology.Predict.chkp

else
    echo "TransDecoder.Predict for '$FILE' using homology evidence already exists."
fi

# 4) Predict signalP ()
output_dir="S3_signalp6_${FILE%.*}_output_dir"

mkdir -p $output_dir

signalp6 --fastafile ${FILE}.transdecoder.pep --organism other --output_dir $output_dir --format txt --mode fast --torch_num_threads 20 --write_procs 8

# 5) Sort pHMM and RegExp conosorter

# 6) rerun emapper


exit

# Load module for mapper.py module
module load conda-2024
source activate base
conda activate trinotate2

#* to run ALL supported computes
TRINOTATE_DATA_DIR=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Trinotate/TRINOTATE_DB/

Trinotate --db sqlite.db --run ALL --CPU $SLURM_NPROCS --transcript_fasta $FILE --transdecoder_pep ${FILE}.transdecoder.pep --trinotate_data_dir  $TRINOTATE_DATA_DIR --use_diamond

exit