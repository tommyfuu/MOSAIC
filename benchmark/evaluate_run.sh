#! /bin/bash -l
 
#SBATCH --partition=scu-cpu
#SBATCH --job-name=evaluate
#SBATCH --time=72:00:00
#SBATCH --mem=20G   # memory requested, units available: K,M,G,T
#SBATCH --cpus-per-task=10
#SBATCH --output evaluate-%j.out
#SBATCH --error evaluate-%j.err


source ~/.bashrc
mamba activate bc_benchmark
echo "conda activated?"
# python /home/chf4012/camp_short-read-assembly/workflow/short-read-assembly.py -slurm --cores 24 -d /athena/masonlab/scratch/users/chf4012/fairbanks/fairbanks_2018_metagenome/short_read_qc/3_error_removal -s /home/chf4012/tom_fairbanks/2018_assembly_samples_filtered.csv --unlock
python3 evaluate.py
exit