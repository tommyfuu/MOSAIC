#! /bin/bash -l
#SBATCH --partition=scu-cpu
#SBATCH --qos=cpu-limited
#SBATCH --array=1-100
#SBATCH --job-name=scenE
#SBATCH --time=6:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --output ../logs/slurm/scenE-%j.out
#SBATCH --error ../logs/slurm/scenE-%j.err
# Python methods (skipped if present) then evaluation (regular and nocovar RF) for a named ablation scenario.
# Usage: sbatch step2_evaluate_scenario.sh <scenario> [count|relab]
source ~/.bashrc
mamba activate mic_data_integration
python3 evaluate.py -o 13 -i $SLURM_ARRAY_TASK_ID --scenario $1 -d ${2:-count}
python3 evaluate.py -o 14 -i $SLURM_ARRAY_TASK_ID --scenario $1 -d ${2:-count}
python3 evaluate.py -o 14 -i $SLURM_ARRAY_TASK_ID --scenario $1 -d ${2:-count} --rerun-rf
exit
