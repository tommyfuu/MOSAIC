#! /bin/bash -l
#SBATCH --partition=scu-cpu,scu-gpu
#SBATCH --qos=normal
#SBATCH --array=1-100
#SBATCH --job-name=imbgen
#SBATCH --time=4:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --output ../logs/slurm/imbgen-%j.out
#SBATCH --error ../logs/slurm/imbgen-%j.err
# Imbalance ablation data (batchfrac20/10, condfrac20/10), one iteration per task.
source ~/.bashrc
mamba activate mic_data_integration
Rscript generate_data_MIDASim_imbalance.R $SLURM_ARRAY_TASK_ID
exit
