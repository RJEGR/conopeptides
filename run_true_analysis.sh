#1)  Test Assembly
# 1.1) spades

EXPORT=/LUSTRE/apps/bioinformatica/SPAdes-3.15.5-Linux/bin/
export PATH=$PATH:$EXPORT


Manifest=$f

left_reads=`awk '{print "--pe1-1",$3}' $Manifest`
right_read=`awk '{print "--pe1-2", $4}' $Manifest`

OUTDIR=${f%.txt}_spades_dir

mkdir -p $OUTDIR

rnaspades.py $left_reads $right_read --pe1-fr -t $CPU -m $MEM -o $DIR/$OUTDIR


# 1.2) Trinity


module load trinityrnaseq-v2.15.1


# From raw assembly match the false/true number of nucleide sequences

# Concat assemblies and map reads back (individually) (run supertranscript.sh script)



# Insilico prediction of orfs (also run for conoserver_nucleic.fa.gz) using transdecoderPredict.sh


