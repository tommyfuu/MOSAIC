#! /bin/bash -l
 
#SBATCH --partition=scu-cpu
#SBATCH --array=344,345,347,361
#SBATCH --job-name=ogsim
#SBATCH --time=5:30:00
#SBATCH --mem=3G   # memory requested, units available: K,M,G,T
#SBATCH --cpus-per-task=1
#SBATCH --output ../logs/slurm/ogsim-%j.out
#SBATCH --error ../logs/slurm/ogsim-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

ITER_ARRAY=( $(seq 1000 ) )
i=$SLURM_ARRAY_TASK_ID

source ~/.bashrc
mamba activate mic_data_integration
echo "conda activated?"
# Rscript methods_benchmarking_sim.R ${ITER_ARRAY[$i-1]}
Rscript methods_benchmarking_sim_052626.R ${ITER_ARRAY[$i-1]}
exit