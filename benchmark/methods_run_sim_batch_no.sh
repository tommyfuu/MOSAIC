#! /bin/bash -l
 
#SBATCH --partition=scu-cpu
#SBATCH --array=427
#SBATCH --job-name=rtcountno
#SBATCH --time=2:30:00
#SBATCH --mem=20G   # memory requested, units available: K,M,G,T
#SBATCH --cpus-per-task=2
#SBATCH --output rtcount_no-%j.out
#SBATCH --error rtcount_no-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

ITER_ARRAY=( $(seq 1000 ) )
i=$SLURM_ARRAY_TASK_ID

source ~/.bashrc
mamba activate bc_benchmark
echo "conda activated?"
Rscript methods_benchmarking_sim_no.R ${ITER_ARRAY[$i-1]}
exit