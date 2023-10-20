#! /bin/bash -l
 
#SBATCH --partition=scu-cpu
#SBATCH --array=1-5
#SBATCH --job-name=methods_benchmark_sim
#SBATCH --time=00:15:00
#SBATCH --mem=3G   # memory requested, units available: K,M,G,T
#SBATCH --cpus-per-task=10
#SBATCH --output methods_benchmark_sim-%j.out
#SBATCH --error methods_benchmark_sim-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

ITER_ARRAY=( $(seq 100 ) )
i=$SLURM_ARRAY_TASK_ID

source ~/.bashrc
mamba activate bc_benchmark
echo "conda activated?"
Rscript methods_benchmarking_sim.R ${ITER_ARRAY[$i-1]}
exit