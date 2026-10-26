#! /bin/bash -l

#SBATCH --partition=scu-cpu
#SBATCH --array=1-1
#SBATCH --job-name=plotrnnc
#SBATCH --time=6:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --output ../logs/slurm/plotrnnc-%j.out
#SBATCH --error ../logs/slurm/plotrnnc-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

# Main-grid line plots: relab, no-relation, nocovar RF, 1000 iterations.
source ~/.bashrc
mamba activate mic_data_integration
echo "conda activated?"
python3 evaluate.py -o 5 -d relab -r no -i 1000 --rerun-rf
exit
