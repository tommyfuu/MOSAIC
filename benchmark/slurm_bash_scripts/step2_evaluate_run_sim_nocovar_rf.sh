#! /bin/bash -l

#SBATCH --partition=scu-gpu
#SBATCH --array=816-821,896
#SBATCH --job-name=evalnc
#SBATCH --time=5:30:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=2
#SBATCH --output ../logs/slurm/evalnc-%j.out
#SBATCH --error ../logs/slurm/evalnc-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

ITER_ARRAY=( $(seq 1000 ) )
i=$SLURM_ARRAY_TASK_ID

source ~/.bashrc
mamba activate mic_data_integration
echo "conda activated?"
cd "${SLURM_SUBMIT_DIR:-.}"
python3 evaluate.py -o 2 -i ${ITER_ARRAY[$i-1]} -d count --rerun-rf
python3 evaluate.py -o 2 -i ${ITER_ARRAY[$i-1]} -d count -r yes --rerun-rf
python3 evaluate.py -o 2 -i ${ITER_ARRAY[$i-1]} -d relab --rerun-rf
python3 evaluate.py -o 2 -i ${ITER_ARRAY[$i-1]} -d relab -r yes --rerun-rf
exit
