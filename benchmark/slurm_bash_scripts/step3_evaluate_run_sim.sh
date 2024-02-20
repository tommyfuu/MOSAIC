#! /bin/bash -l
 
#SBATCH --partition=scu-cpu
#SBATCH --array=295,874
#SBATCH --job-name=fixryes
#SBATCH --time=10:00:00
#SBATCH --mem=5G   # memory requested, units available: K,M,G,T
#SBATCH --cpus-per-task=2
#SBATCH --output fixryes-%j.out
#SBATCH --error fixryes-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

ITER_ARRAY=( $(seq 1000 ) )
i=$SLURM_ARRAY_TASK_ID

source ~/.bashrc
mamba activate bc_benchmark
echo "conda activated?"
# python3 evaluate.py -o 2 -i ${ITER_ARRAY[$i-1]} -r yes -d count -a cond_1
python3 evaluate.py -o 2 -i ${ITER_ARRAY[$i-1]} -r yes -d relab -a cond_1
exit