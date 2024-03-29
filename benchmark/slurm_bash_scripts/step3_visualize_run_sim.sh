#! /bin/bash -l
 
#SBATCH --partition=scu-cpu
#SBATCH --array=1-1000
#SBATCH --job-name=ev_sim
#SBATCH --time=10:00:00
#SBATCH --mem=5G   # memory requested, units available: K,M,G,T
#SBATCH --cpus-per-task=2
#SBATCH --output ev_sim-%j.out
#SBATCH --error ev_sim-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

ITER_ARRAY=( $(seq 1000 ) )
i=$SLURM_ARRAY_TASK_ID

source ~/.bashrc
mamba activate mic_data_integration
echo "conda activated?"
python3 evaluate.py -o 3 -i ${ITER_ARRAY[$i-1]} -r no -d count -a cond_1
python3 evaluate.py -o 3 -i ${ITER_ARRAY[$i-1]} -r no -d relab -a cond_1
python3 evaluate.py -o 3 -i ${ITER_ARRAY[$i-1]} -r yes -d count -a cond_1
python3 evaluate.py -o 3 -i ${ITER_ARRAY[$i-1]} -r yes -d relab -a cond_1
exit