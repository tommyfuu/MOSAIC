#! /bin/bash -l
 
#SBATCH --partition=scu-cpu
#SBATCH --array=1-100
#SBATCH --job-name=yes_evaluate_sim
#SBATCH --time=24:00:00
#SBATCH --mem=40G   # memory requested, units available: K,M,G,T
#SBATCH --cpus-per-task=10
#SBATCH --output evaluate-%j.out
#SBATCH --error evaluate-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

ITER_ARRAY=( $(seq 100 ) )
i=$SLURM_ARRAY_TASK_ID

source ~/.bashrc
mamba activate bc_benchmark
echo "conda activated?"
# python /home/chf4012/camp_short-read-assembly/workflow/short-read-assembly.py -slurm --cores 24 -d /athena/masonlab/scratch/users/chf4012/fairbanks/fairbanks_2018_metagenome/short_read_qc/3_error_removal -s /home/chf4012/tom_fairbanks/2018_assembly_samples_filtered.csv --unlock
# python3 evaluate.py -o 1 -i ${ITER_ARRAY[$i-1]} -r yes -d count
python3 evaluate.py -o 2 -i ${ITER_ARRAY[$i-1]} -r yes -d count
# python3 evaluate.py -o 1 -i ${ITER_ARRAY[$i-1]} -r yes -d relab
python3 evaluate.py -o 2 -i ${ITER_ARRAY[$i-1]} -r yes -d relab
exit