#! /bin/bash -l
 
#SBATCH --partition=scu-gpu
#SBATCH --job-name=pymethrw
#SBATCH --time=5:40:00
#SBATCH --mem=7G   # memory requested, units available: K,M,G,T
#SBATCH --cpus-per-task=1
#SBATCH --output pymethrw-%j.out
#SBATCH --error pymethrw-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

source ~/.bashrc
mamba activate bc_benchmark
echo "conda activated?"
## PLEASE MAKE SURE THE UNCOMMENT THE SECTION STARTING WITH
## ## RUN HARMONY/PERCENTILE_NORM FOR RW
## AND END WITH
## ## EVALUATE METHODS ON REAL-WORLD DATASET
## IN evaluate.py before running the following
python3 evaluate.py -o 4
exit
