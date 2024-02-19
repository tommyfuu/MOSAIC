#! /bin/bash -l
 
#SBATCH --partition=scu-cpu
#SBATCH --array=951-1000
#SBATCH --job-name=rtcountyes
#SBATCH --time=6:40:00
#SBATCH --mem=20G   # memory requested, units available: K,M,G,T
#SBATCH --cpus-per-task=2
#SBATCH --output rtcount_yes-%j.out
#SBATCH --error rtcount_yes-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

ITER_ARRAY=( $(seq 1000 ) )
i=$SLURM_ARRAY_TASK_ID

source ~/.bashrc
mamba activate bc_benchmark
echo "conda activated?"
Rscript methods_benchmarking_sim.R ${ITER_ARRAY[$i-1]}
exit