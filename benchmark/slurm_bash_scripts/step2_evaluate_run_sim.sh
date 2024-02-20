#! /bin/bash -l
 
#SBATCH --partition=scu-cpu
#SBATCH --array=206,753,861,963
#SBATCH --job-name=fixrno
#SBATCH --time=10:00:00
#SBATCH --mem=2G   # memory requested, units available: K,M,G,T
#SBATCH --cpus-per-task=2
#SBATCH --output fixrno-%j.out
#SBATCH --error fixrno-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

ITER_ARRAY=( $(seq 1000 ) )
i=$SLURM_ARRAY_TASK_ID

source ~/.bashrc
mamba activate bc_benchmark
echo "conda activated?"
# python3 evaluate.py -o 2 -i ${ITER_ARRAY[$i-1]} -r no -d count -a cond_1
python3 evaluate.py -o 2 -i ${ITER_ARRAY[$i-1]} -r no -d relab -a cond_1
exit