#! /bin/bash -l

#SBATCH --partition=scu-cpu
#SBATCH --job-name=pubtator-search
#SBATCH --time=1:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output ../logs/slurm/pubtator_search_job.out
#SBATCH --error ../logs/slurm/pubtator_search_job.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

# Sequential on purpose: PubTator3 allows ~3 requests/sec, so a job array would not be faster.

source ~/.bashrc
mamba activate bc_benchmark
cd "${SLURM_SUBMIT_DIR:-.}"   # submit from pubtator/

REPO="$(cd .. && pwd)"

python3 pubtator_api_search.py \
  --input "$REPO/trial_rw/outputs/autism_2_microbiomeHD_taxa_union_wilcoxon.csv" \
  --taxa-col cleaned_taxa --condition "ASD" \
  --output "$REPO/trial_rw/outputs/autism_2_microbiomeHD_pubtator_counts.csv"

python3 pubtator_api_search.py \
  --input "$REPO/trial_rw/outputs/cdi_3_microbiomeHD_clean_taxa_union_wilcoxon.csv" \
  --taxa-col cleaned_taxa --condition "CDI" \
  --output "$REPO/trial_rw/outputs/cdi_3_microbiomeHD_clean_pubtator_counts.csv"

python3 pubtator_api_search.py \
  --input "$REPO/outputs/ibd_3_CMD_taxa_union_wilcoxon.csv" \
  --taxa-col cleaned_taxa --condition "IBD" \
  --output "$REPO/outputs/ibd_3_CMD_pubtator_counts.csv"

python3 pubtator_api_search.py \
  --input "$REPO/outputs/crc_11_CMD_taxa_union_wilcoxon.csv" \
  --taxa-col cleaned_taxa --condition "CRC" \
  --output "$REPO/outputs/crc_11_CMD_pubtator_counts.csv"

echo "ALL DATASETS DONE"
