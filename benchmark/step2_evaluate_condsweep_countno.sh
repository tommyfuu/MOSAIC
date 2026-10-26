#! /bin/bash -l

#SBATCH --partition=scu-cpu
#SBATCH --array=1-50
#SBATCH --job-name=condsweepE-cn
#SBATCH --time=3:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --output ../logs/slurm/condsweepE-cn-%j.out
#SBATCH --error ../logs/slurm/condsweepE-cn-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

# Evaluate condsweep results, count/no-relation only (regular and nocovar RF).
i=$SLURM_ARRAY_TASK_ID

source ~/.bashrc
mamba activate mic_data_integration
echo "conda activated?"
python3 evaluate.py -o 7 -i $i -d count --condsweep-relation no_relation
python3 evaluate.py -o 7 -i $i -d count --condsweep-relation no_relation --rerun-rf
exit
