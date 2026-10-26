#! /bin/bash -l
 
#SBATCH --partition=scu-cpu
#SBATCH --array=101-1000
#SBATCH --job-name=evalyesn
#SBATCH --time=8:30:00
#SBATCH --mem=2G   # memory requested, units available: K,M,G,T
#SBATCH --cpus-per-task=1
#SBATCH --output ../logs/slurm/evalyesn-%j.out
#SBATCH --error ../logs/slurm/evalyesn-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

ITER_ARRAY=( $(seq 1000 ) )
i=$SLURM_ARRAY_TASK_ID

source ~/.bashrc
mamba activate mic_data_integration
echo "conda activated?"
# python3 evaluate.py -o 2 -i ${ITER_ARRAY[$i-1]} -r no -d count  -a cond_1
# python3 evaluate.py -o 2 -i ${ITER_ARRAY[$i-1]} -r no -d relab -a cond_1
python3 evaluate.py -o 2 -i ${ITER_ARRAY[$i-1]} -r yes -d count -a cond_1
python3 evaluate.py -o 2 -i ${ITER_ARRAY[$i-1]} -r yes -d relab -a cond_1
exit

