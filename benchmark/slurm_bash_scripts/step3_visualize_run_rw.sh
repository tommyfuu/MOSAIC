#! /bin/bash -l
 
#SBATCH --partition=scu-cpu
#SBATCH --job-name=ev_rw
#SBATCH --time=10:00:00
#SBATCH --mem=5G   # memory requested, units available: K,M,G,T
#SBATCH --cpus-per-task=2
#SBATCH --output ev_rw-%j.out
#SBATCH --error ev_rw-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

source ~/.bashrc
mamba activate bc_benchmark
echo "conda activated?"
python3 evaluate.py -o 5
exit


