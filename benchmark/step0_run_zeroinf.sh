#! /bin/bash -l
#SBATCH --partition=scu-cpu,scu-gpu
#SBATCH --qos=normal
#SBATCH --array=1-50
#SBATCH --job-name=zeroinfgen
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output ../logs/slurm/zeroinfgen-%j.out
#SBATCH --error ../logs/slurm/zeroinfgen-%j.err
# Zero-inflation ablation data: raise the zero fraction of the main-grid count data (no relation, OR=1, cond 0.5, batch 0.5) to 0.75-0.90.
source ~/.bashrc
mamba activate mic_data_integration
python3 generate_zeroinf_ablation.py -i $SLURM_ARRAY_TASK_ID
exit
