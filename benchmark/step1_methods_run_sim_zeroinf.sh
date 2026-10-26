#! /bin/bash -l
#SBATCH --partition=scu-cpu,scu-gpu
#SBATCH --qos=normal
#SBATCH --array=1-50
#SBATCH --job-name=zeroinfR
#SBATCH --time=6:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=10
#SBATCH --output ../logs/slurm/zeroinfR-%j.out
#SBATCH --error ../logs/slurm/zeroinfR-%j.err
# R batch-correction methods on the zero-inflation ablation data. Usage: sbatch step1_methods_run_sim_zeroinf.sh <75|80|85|90>
source ~/.bashrc
mamba activate mic_data_integration
Rscript methods_benchmarking_sim_zeroinf.R $SLURM_ARRAY_TASK_ID ${1:-all}
exit
