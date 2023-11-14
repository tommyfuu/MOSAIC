#! /bin/bash -l
 
#SBATCH --partition=scu-cpu
#SBATCH --array=1-1000
#SBATCH --job-name=sim2
#SBATCH --time=0:02:00
#SBATCH --mem=1G   # memory requested, units available: K,M,G,T
#SBATCH --cpus-per-task=1
#SBATCH --output evaluate-%j.out
#SBATCH --error evaluate-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

ITER_ARRAY=( $(seq 1000 ) )
i=$SLURM_ARRAY_TASK_ID

source ~/.bashrc
mamba activate bc_benchmark
echo "conda activated?"
# python3 evaluate.py -o 1 -i ${ITER_ARRAY[$i-1]} -r no -d count
python3 evaluate.py -o 1 -i ${ITER_ARRAY[$i-1]} -r no -d relab
# python3 evaluate.py -o 1 -i ${ITER_ARRAY[$i-1]} -r yes -d count
python3 evaluate.py -o 1 -i ${ITER_ARRAY[$i-1]} -r yes -d relab
exit