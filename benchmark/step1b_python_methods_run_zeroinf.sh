#! /bin/bash -l
#SBATCH --partition=scu-cpu,scu-gpu
#SBATCH --qos=normal
#SBATCH --array=1-50
#SBATCH --job-name=zeroinfpy
#SBATCH --time=3:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --output ../logs/slurm/zeroinfpy-%j.out
#SBATCH --error ../logs/slurm/zeroinfpy-%j.err
# Python batch-correction methods (harmony, percentile_norm, debiasm) on the zero-inflation ablation data.
source ~/.bashrc
mamba activate mic_data_integration
python3 evaluate.py -o 13 -i $SLURM_ARRAY_TASK_ID
exit
