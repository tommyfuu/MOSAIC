#! /bin/bash -l
 
#SBATCH --partition=scu-cpu
#SBATCH --array=1-1000
#SBATCH --job-name=fixrno
#SBATCH --time=10:00:00
#SBATCH --mem=2G   # memory requested, units available: K,M,G,T
#SBATCH --cpus-per-task=2
#SBATCH --output ../logs/slurm/fixrno-%j.out
#SBATCH --error ../logs/slurm/fixrno-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

ITER_ARRAY=( $(seq 1000 ) )
i=$SLURM_ARRAY_TASK_ID

source ~/.bashrc
mamba activate bc_benchmark
echo "conda activated?"
python3 archive/evaluate_legacy_pre052626.py -o 2 -i ${ITER_ARRAY[$i-1]} -r no -d count -a cond_1
python3 archive/evaluate_legacy_pre052626.py -o 2 -i ${ITER_ARRAY[$i-1]} -r no -d relab -a cond_1
exit