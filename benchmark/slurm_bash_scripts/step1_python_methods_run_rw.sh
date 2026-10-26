#! /bin/bash -l
 
#SBATCH --partition=scu-gpu
#SBATCH --job-name=pymethrw
#SBATCH --time=5:40:00
#SBATCH --mem=7G   # memory requested, units available: K,M,G,T
#SBATCH --cpus-per-task=1
#SBATCH --output ../logs/slurm/pymethrw-%j.out
#SBATCH --error ../logs/slurm/pymethrw-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

source ~/.bashrc
mamba activate bc_benchmark
echo "conda activated?"
# Before running, uncomment the block in archive/evaluate_legacy_pre052626.py from
# "## RUN HARMONY/PERCENTILE_NORM FOR RW" to "## EVALUATE METHODS ON REAL-WORLD DATASET".
# Legacy entrypoint; current runs use benchmark/evaluate.py -o 4.
python3 archive/evaluate_legacy_pre052626.py -o 4
exit
