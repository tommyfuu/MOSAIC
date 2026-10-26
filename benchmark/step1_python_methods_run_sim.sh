#! /bin/bash -l
 
#SBATCH --partition=scu-cpu
#SBATCH --array=1-1000
#SBATCH --job-name=pymethsim
#SBATCH --time=7:00:00
#SBATCH --mem=7G   # memory requested, units available: K,M,G,T
#SBATCH --cpus-per-task=1
#SBATCH --output ../logs/slurm/pymethsim-%j.out
#SBATCH --error ../logs/slurm/pymethsim-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

ITER_ARRAY=( $(seq 1000 ) )
i=$SLURM_ARRAY_TASK_ID

source ~/.bashrc
mamba activate mic_data_integration
echo "conda activated?"
python3 evaluate.py -o 1 -i ${ITER_ARRAY[$i-1]} -r yes -d count
python3 evaluate.py -o 1 -i ${ITER_ARRAY[$i-1]} -d count
python3 evaluate.py -o 1 -i ${ITER_ARRAY[$i-1]} -d relab
python3 evaluate.py -o 1 -i ${ITER_ARRAY[$i-1]} -r yes -d relab
exit



# python3 archive/evaluate_legacy_pre052626.py -o 1 -i 101 -r yes -d count
# python3 archive/evaluate_legacy_pre052626.py -o 2 -i 101 -r yes -d count