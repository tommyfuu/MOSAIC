#! /bin/bash -l
#SBATCH --partition=scu-cpu
#SBATCH --job-name=pubtator-keyterms
#SBATCH --time=6:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output ../logs/slurm/pubtator_keyterms-%j.out
#SBATCH --error ../logs/slurm/pubtator_keyterms-%j.err

# Searches each dataset's significant taxa against its clinical key terms, read from the pubtator_key_terms
# section of the dataset's config YAML (key_terms.py; at least five terms). Sequential on purpose: PubTator3
# allows ~3 requests/sec. Submit from pubtator/, optionally with dataset names as arguments
# (default: the four real-world datasets). Then run pubtator_relevance_tables.py.
source ~/.bashrc
mamba activate bc_benchmark
cd "${SLURM_SUBMIT_DIR:-.}"
REPO="$(cd .. && pwd)"
IN="$REPO/trial_rw/outputs"
OUT="$IN/pubtator_keyterms"
mkdir -p "$OUT"

search() {  # dataset, term
  slug=$(echo "$2" | tr ' ' '_' | tr -cd '[:alnum:]_-')
  python3 pubtator_api_search.py --input "$IN/$1_taxa_union_wilcoxon.csv" --taxa-col cleaned_taxa \
    --condition "$2" --output "$OUT/$1__${slug}.csv"
}

if [ $# -gt 0 ]; then DATASETS=("$@")
else DATASETS=(autism_2_microbiomeHD cdi_3_microbiomeHD_clean ibd_3_CMD crc_11_CMD); fi
# key-term history: IBD terms revised in response to review (previously "Irritable Bowel Syndrome", "FODMAP diet",
# "inflammation"); CRC "diarrhea" replaced by "tumorigenesis" (Oct 2026; non-specific and already a CDI term)
for ds in "${DATASETS[@]}"; do
  terms=$(python3 key_terms.py "$ds") || exit 1
  while IFS= read -r t; do search "$ds" "$t"; done <<< "$terms"
done
# taxon-only publication counts ("popularity"), for the sensitivity analysis (pubtator_sensitivity.py)
for ds in "${DATASETS[@]}"; do
  python3 pubtator_api_search.py --input "$IN/${ds}_taxa_union_wilcoxon.csv" --taxa-col cleaned_taxa \
    --condition "" --output "$OUT/${ds}_taxon_total.csv"
done
echo "ALL KEY-TERM SEARCHES DONE"
