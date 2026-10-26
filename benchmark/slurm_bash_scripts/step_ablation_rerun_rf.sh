#!/bin/bash
#SBATCH --partition=preempt_cpu
#SBATCH --job-name=ablation_rerun_rf
#SBATCH --array=1-1000
#SBATCH --output=../logs/slurm/rf_%a.out
#SBATCH --error=../logs/slurm/rf_%a.err
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2

# --refresh-rf together with --rerun-rf refreshes both the regular and nocovar RF summaries.
source ~/.bashrc
conda activate mic_data_integration


# plot_PCOA_multiple() (called on every 50th "keep" iteration) looks for
# PERMANOVA_supporting.R relative to the working directory, so run from benchmark/.
cd "${SLURM_SUBMIT_DIR:-.}"
python3 -u evaluate.py \
    -o 7 -i ${SLURM_ARRAY_TASK_ID} --refresh-rf --rerun-rf
