#! /bin/bash -l

#SBATCH --partition=scu-cpu
#SBATCH --job-name=condsweepViz-cn
#SBATCH --time=1:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --output ../logs/slurm/condsweepViz-cn-%j.out
#SBATCH --error ../logs/slurm/condsweepViz-cn-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

# Condsweep line plots for count/no-relation, 1000 iterations (regular and nocovar).
source ~/.bashrc
mamba activate mic_data_integration
echo "conda activated?"
python3 evaluate.py -o 9 -d count -r no -i 1000
python3 evaluate.py -o 9 -d count -r no -i 1000 --rerun-rf
exit
