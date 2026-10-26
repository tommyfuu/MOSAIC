#! /bin/bash -l

#SBATCH --partition=scu-cpu
#SBATCH --array=1-50
#SBATCH --job-name=condsweepR
#SBATCH --time=10:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=10
#SBATCH --output ../logs/slurm/condsweepR-%j.out
#SBATCH --error ../logs/slurm/condsweepR-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

# R batch-correction methods on condsweep data, all 4 datatype/relation combos, covar and nocovar passes.
source ~/.bashrc
mamba activate mic_data_integration
echo "conda activated?"
Rscript methods_benchmarking_sim_condsweep.R $SLURM_ARRAY_TASK_ID
exit
