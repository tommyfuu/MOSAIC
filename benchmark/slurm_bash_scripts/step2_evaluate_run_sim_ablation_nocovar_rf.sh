#! /bin/bash -l

#SBATCH --partition=scu-cpu
#SBATCH --array=1-1000
#SBATCH --job-name=evalabl
#SBATCH --time=10:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=2
#SBATCH --output ../logs/slurm/evalabl-%j.out
#SBATCH --error ../logs/slurm/evalabl-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

ITER_ARRAY=( $(seq 1000 ) )
i=$SLURM_ARRAY_TASK_ID

source ~/.bashrc
mamba activate mic_data_integration
echo "conda activated?"
cd "${SLURM_SUBMIT_DIR:-.}"
python3 evaluate.py -o 7 -i ${ITER_ARRAY[$i-1]}
python3 evaluate.py -o 7 -i ${ITER_ARRAY[$i-1]} --rerun-rf
exit
