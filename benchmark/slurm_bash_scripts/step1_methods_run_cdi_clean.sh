#! /bin/bash -l

#SBATCH --partition=scu-cpu
#SBATCH --job-name=cdicleanmeth
#SBATCH --time=9:00:00
#SBATCH --mem=20G   # memory requested, units available: K,M,G,T
#SBATCH --cpus-per-task=10
#SBATCH --output ../logs/slurm/cdicleanmeth-%j.out
#SBATCH --error ../logs/slurm/cdicleanmeth-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu

source ~/.bashrc
mamba activate mic_data_integration   # has the ruv package needed for ruv3
echo "conda activated?"

cd "${SLURM_SUBMIT_DIR:-.}/.."   # submit from benchmark/; config paths (src: ./data, ...) are relative to the repo root

# methods_benchmarking_fixed_052626.R has no --config flag: it always reads the
# Snakefile's first `configfile:` line. Point it at config_cdi_clean.yml, run, restore.
sed -i '1s#.*#configfile: "./config_cdi_clean.yml"#' Snakefile
Rscript benchmark/methods_benchmarking_fixed_052626.R
sed -i '1s#.*#configfile: "./config_simulate_052626.yml"#' Snakefile

exit
