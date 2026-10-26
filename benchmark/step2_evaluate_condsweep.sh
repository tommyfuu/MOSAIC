#! /bin/bash -l

#SBATCH --partition=scu-cpu
#SBATCH --array=1-50
#SBATCH --job-name=condsweepeval
#SBATCH --time=3:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --output ../logs/slurm/condsweepeval-%j.out
#SBATCH --error ../logs/slurm/condsweepeval-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

# Evaluate all methods on condsweep data (regular and nocovar RF). Run after step1 and step1b finish.
i=$SLURM_ARRAY_TASK_ID

source ~/.bashrc
mamba activate mic_data_integration
echo "conda activated?"
python3 evaluate.py -o 7 -i $i -d count
python3 evaluate.py -o 7 -i $i -d count --rerun-rf
python3 evaluate.py -o 7 -i $i -d relab
python3 evaluate.py -o 7 -i $i -d relab --rerun-rf
exit
