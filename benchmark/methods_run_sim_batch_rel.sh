#! /bin/bash -l
 
#SBATCH --partition=scu-gpu
#SBATCH --array=771-1000
#SBATCH --job-name=rtrelno
#SBATCH --time=8:00:00
#SBATCH --mem=20G   # memory requested, units available: K,M,G,T
#SBATCH --cpus-per-task=2
#SBATCH --output rtrel_no-%j.out
#SBATCH --error rtrel_no-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

ITER_ARRAY=( $(seq 1000 ) )
i=$SLURM_ARRAY_TASK_ID

source ~/.bashrc
mamba activate bc_benchmark
echo "conda activated?"
Rscript methods_benchmarking_sim_rel.R ${ITER_ARRAY[$i-1]}
exit