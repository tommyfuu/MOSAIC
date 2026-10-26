#! /bin/bash -l

#SBATCH --partition=scu-cpu
#SBATCH --array=1-50
#SBATCH --job-name=condsweeppy
#SBATCH --time=2:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --output ../logs/slurm/condsweeppy-%j.out
#SBATCH --error ../logs/slurm/condsweeppy-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

# Python batch-correction methods (harmony, percentile_norm, debiasm) on condsweep data, iterations 1-50. Can run alongside the R methods.
i=$SLURM_ARRAY_TASK_ID

source ~/.bashrc
mamba activate mic_data_integration
echo "conda activated?"
python3 evaluate.py -o 6 -i $i -d count
python3 evaluate.py -o 6 -i $i -d relab
exit
