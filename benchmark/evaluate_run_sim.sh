#! /bin/bash -l
 
#SBATCH --partition=scu-cpu
#SBATCH --array=1-100
#SBATCH --job-name=evaluate_sim
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
# python3 evaluate.py -o 1 -i ${ITER_ARRAY[$i-1]} -r no -d count
python3 evaluate.py -o 2 -i ${ITER_ARRAY[$i-1]} -r no -d count
# python3 evaluate.py -o 1 -i ${ITER_ARRAY[$i-1]} -r no -d relab
python3 evaluate.py -o 2 -i ${ITER_ARRAY[$i-1]} -r no -d relab
exit