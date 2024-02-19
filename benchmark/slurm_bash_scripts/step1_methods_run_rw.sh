#! /bin/bash -l
 
#SBATCH --partition=scu-cpu
#SBATCH --job-name=rwmeth
#SBATCH --time=9:00:00
#SBATCH --mem=20G   # memory requested, units available: K,M,G,T
#SBATCH --cpus-per-task=10
#SBATCH --output rwmeth-%j.out
#SBATCH --error rwmeth-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
# ITER_ARRAY=( $(seq 5 ) )
source ~/.bashrc
mamba activate bc_benchmark
echo "conda activated?"
Rscript methods_benchmarking.R 3
exit