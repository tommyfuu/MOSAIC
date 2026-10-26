#! /bin/bash -l

#SBATCH --partition=scu-gpu
#SBATCH --array=1-1000
#SBATCH --job-name=pyabl
#SBATCH --time=0:40:00
#SBATCH --mem=7G
#SBATCH --cpus-per-task=1
#SBATCH --output ../logs/slurm/pyabl-%j.out
#SBATCH --error ../logs/slurm/pyabl-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

ITER_ARRAY=( $(seq 1000 ) )
i=$SLURM_ARRAY_TASK_ID

source ~/.bashrc
mamba activate mic_data_integration
echo "conda activated?"
cd "${SLURM_SUBMIT_DIR:-.}"
python3 evaluate.py -o 6 -i ${ITER_ARRAY[$i-1]} -d count
exit
