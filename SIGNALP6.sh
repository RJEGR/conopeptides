#!/bin/sh
## Directivas
#SBATCH --job-name=SIGNALP
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24

module load conda-2024

EXPORT=/home/rgomez/.local/bin #signalP
export PATH=$PATH:$EXPORT

NPROCS=$SLURM_NPROCS

# after installation, edit:
# SIGNALP_DIR=$(python3 -c "import signalp; import os; print(os.path.dirname(signalp.__file__))" )
# cp -r signalp-6-package/models/* $SIGNALP_DIR/model_weights/
# vi line 144 to torch.set_num_threads(int(args.torch_num_threads)) from signalp/predict.py 

# signalp6  --torch_num_threads $NPROCS --write_procs 8
# example

INPUT=$1

output_dir=${PWD}/signalp6_${INPUT%.*}_output_dir

mkdir -p $output_dir

signalp6 --fastafile $INPUT --organism other --output_dir $output_dir --format txt --mode fast --torch_num_threads 20 --write_procs 8

exit

# Run prop 1.0 (use signalp 3.0)

$SIGNALP -t euk -m nn -trunc 100 $f | grep '^# Most'`

# signalp -od None -org euk -m fast -bs 100 -ff $f | grep '^# Most'

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/prop-1.0c
export PATH=$PATH:$EXPORT

prop -s four.fsa

 4. Test the package:

      ./prop test/EDA_HUMAN.fsa         # one sequence
      ./prop test/FGF2_HUMAN.fsa
      ./prop -g -s test/GDNF_HUMAN.fsa  # one sequence, graphics, signalp
      ./prop -g -s test/P53_HUMAN.fsa
      ./prop -g -s test/four.fsa        # many sequences, graphics, signalp