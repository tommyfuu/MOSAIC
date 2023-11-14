#! /bin/bash -l
 
#SBATCH --partition=scu-cpu
#SBATCH --array=111
#SBATCH --job-name=mtrelyes
#SBATCH --time=12:00:00
#SBATCH --mem=3G   # memory requested, units available: K,M,G,T
#SBATCH --cpus-per-task=1
#SBATCH --output mtrelyes-%j.out
#SBATCH --error mtrelyes-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

ITER_ARRAY=( $(seq 1000 ) )
i=$SLURM_ARRAY_TASK_ID

source ~/.bashrc
mamba activate bc_benchmark
echo "conda activated?"
Rscript methods_benchmarking_sim_yes.R ${ITER_ARRAY[$i-1]}
exit