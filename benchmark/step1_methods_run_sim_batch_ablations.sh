#! /bin/bash -l

#SBATCH --partition=scu-cpu
#SBATCH --array=344,345,347,361
#SBATCH --job-name=ablations
#SBATCH --time=10:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=10
#SBATCH --output ../logs/slurm/ablations-%j.out
#SBATCH --error ../logs/slurm/ablations-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

ITER_ARRAY=( $(seq 1000 ) )
i=$SLURM_ARRAY_TASK_ID

source ~/.bashrc
mamba activate mic_data_integration
echo "conda activated?"
Rscript methods_benchmarking_sim_052626.R ${ITER_ARRAY[$i-1]}
exit
