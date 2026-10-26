#! /bin/bash -l

#SBATCH --partition=scu-cpu
#SBATCH --array=1-50
#SBATCH --job-name=condsweepgen
#SBATCH --time=4:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --output ../logs/slurm/condsweepgen-%j.out
#SBATCH --error ../logs/slurm/condsweepgen-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

# Generate condsweep data (OR=1.25, cond_effect in {0.005,0.01,0.02,0.05,0.1}, batch=0.5, no/yes relation), iterations 1-50.
source ~/.bashrc
mamba activate mic_data_integration
echo "conda activated?"
Rscript generate_data_MIDASim_condsweep.R $SLURM_ARRAY_TASK_ID
exit
