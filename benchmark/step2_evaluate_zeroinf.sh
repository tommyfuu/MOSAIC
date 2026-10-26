#! /bin/bash -l
#SBATCH --partition=scu-cpu,scu-gpu
#SBATCH --qos=normal
#SBATCH --array=1-50
#SBATCH --job-name=zeroinfE
#SBATCH --time=6:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --output ../logs/slurm/zeroinfE-%j.out
#SBATCH --error ../logs/slurm/zeroinfE-%j.err
# Evaluate the zero-inflation ablation (regular and nocovar RF). Usage: sbatch step2_evaluate_zeroinf.sh [75|80|85|90]
source ~/.bashrc
mamba activate mic_data_integration
python3 evaluate.py -o 14 -i $SLURM_ARRAY_TASK_ID ${1:+--zeroinf-level $1}
python3 evaluate.py -o 14 -i $SLURM_ARRAY_TASK_ID --rerun-rf ${1:+--zeroinf-level $1}
exit
