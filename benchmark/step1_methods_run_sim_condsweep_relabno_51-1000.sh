#! /bin/bash -l

#SBATCH --partition=scu-gpu
#SBATCH --array=51-1000
#SBATCH --job-name=condsweepR-rnfull
#SBATCH --time=4:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=10
#SBATCH --output ../logs/slurm/condsweepR-rn-%j.out
#SBATCH --error ../logs/slurm/condsweepR-rn-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

# R batch-correction methods on condsweep data, relab/no-relation only.
source ~/.bashrc
mamba activate mic_data_integration
echo "conda activated?"
Rscript methods_benchmarking_sim_condsweep.R $SLURM_ARRAY_TASK_ID relab_no
exit
