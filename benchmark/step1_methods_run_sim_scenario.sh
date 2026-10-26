#! /bin/bash -l
#SBATCH --partition=scu-cpu,scu-gpu
#SBATCH --qos=normal
#SBATCH --array=1-100
#SBATCH --job-name=scenR
#SBATCH --time=6:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=10
#SBATCH --output ../logs/slurm/scenR-%j.out
#SBATCH --error ../logs/slurm/scenR-%j.err
# R batch-correction methods for a named ablation scenario. Usage: sbatch step1_methods_run_sim_scenario.sh <scenario> [count|relab]
source ~/.bashrc
mamba activate mic_data_integration
Rscript methods_benchmarking_sim_scenario.R $SLURM_ARRAY_TASK_ID $1 ${2:-count}
exit
