
#!/bin/bash
#SBATCH --job-name=spades
#SBATCH -N 2
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20
#SBATCH --error=slurm-%j.err

EXPORT=/LUSTRE/apps/bioinformatica/SPAdes-3.15.5-Linux/bin/
export PATH=$PATH:$EXPORT

CPU=$SLURM_NPROCS
MEM=$SLURM_MEM_PER_NODE
DIR=$SLURM_SUBMIT_DIR

Manifest=$1 # samples.txt

left_reads=`awk '{print "--pe1-1",$3}' $Manifest`
right_read=`awk '{print "--pe1-2", $4}' $Manifest`

rnaspades.py $left_reads $right_read \
  -k 21,25,49,73 \
  --pe1-fr -t $CPU -m $MEM \
  --checkpoints last -o $DIR/Rnaspades_out

exit

# samples txt

for file in $(ls *_R1.fq.gz | grep fq)
do
basename=${file##*/}
base="${basename%*_R1.fq.gz}"
bs="${base%*_E_*}"
fwr=${base}_R1.fq.gz
rev=${base}_R2.fq.gz
echo "$bs" "$bs" `printf "$PWD/$fwr"` `printf "$PWD/$rev"`
done > samples.txt

# single-read

for file in $(ls *_merged.fq.gz| grep fq)
do
basename=${file##*/}
bs="${basename%*_merged.fq.gz}"
single_read=${bs}_merged.fq.gz
rev=${base}_R2.fq.gz
echo "$bs" "$bs" `printf "$PWD/$single_read"` 
done > samples_merged.txt

# OR
cd /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/californicus/07.Reference/FASTP_OUT_DIR

cat *_merged.fq.gz > MERGED.fq.gz

#!/bin/bash
#SBATCH --job-name=spades
#SBATCH -N 2
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20
#SBATCH --error=slurm-%j.err

EXPORT=/LUSTRE/apps/bioinformatica/SPAdes-3.15.5-Linux/bin/
export PATH=$PATH:$EXPORT

CPU=$SLURM_NPROCS
MEM=$SLURM_MEM_PER_NODE
DIR=$SLURM_SUBMIT_DIR


rnaspades.py --s1 MERGED.fq.gz -t $CPU -m $MEM -o $DIR/MERGED_Rnaspades_out

exit