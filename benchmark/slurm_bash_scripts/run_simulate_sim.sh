#! /bin/bash -l
 
#SBATCH --partition=scu-cpu
#SBATCH --array=1-1000
#SBATCH --job-name=midasgen
#SBATCH --time=72:00:00
#SBATCH --mem=20G   # memory requested, units available: K,M,G,T
#SBATCH --cpus-per-task=10
#SBATCH --output midasgen-%j.out
#SBATCH --error midasgen-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

ITER_ARRAY=( $(seq 1000 ) )
i=$SLURM_ARRAY_TASK_ID

source ~/.bashrc
mamba activate bc_benchmark
echo "conda activated?"
Rscript generate_data_MIDAS.R ${ITER_ARRAY[$i-1]}
exit