#! /bin/bash -l

#SBATCH --partition=scu-gpu
#SBATCH --array=51-1000
#SBATCH --job-name=condsweepE-ryfull
#SBATCH --time=3:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --output ../logs/slurm/condsweepE-ry-%j.out
#SBATCH --error ../logs/slurm/condsweepE-ry-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

# Evaluate condsweep results, relab/yes-relation only (regular and nocovar RF).
i=$SLURM_ARRAY_TASK_ID

source ~/.bashrc
mamba activate mic_data_integration
echo "conda activated?"
python3 evaluate.py -o 7 -i $i -d relab --condsweep-relation yes_relation
python3 evaluate.py -o 7 -i $i -d relab --condsweep-relation yes_relation --rerun-rf
exit
