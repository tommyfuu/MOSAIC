#! /bin/bash -l
#SBATCH --partition=scu-cpu
#SBATCH --job-name=pubtator-keyterms
#SBATCH --time=6:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output ../logs/slurm/pubtator_keyterms-%j.out
#SBATCH --error ../logs/slurm/pubtator_keyterms-%j.err

# Searches each dataset's significant taxa against the five clinical key terms per disease
# used in the MOSAIC preprint (Methods). Sequential on purpose: PubTator3 allows ~3 requests/sec.
# Submit from pubtator/. Then run pubtator_relevance_tables.py.
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

for t in "autism" "autism spectrum disorders" "gut-brain axis" "neurological disorders" "nervous systems"; do
  search autism_2_microbiomeHD "$t"; done
for t in "Clostridioides difficile" "diarrhea" "antibiotics" "toxins" "infection"; do
  search cdi_3_microbiomeHD_clean "$t"; done
# IBD terms revised in response to review (previously "Irritable Bowel Syndrome", "FODMAP diet", "inflammation")
for t in "inflammatory bowel disease" "Crohn's disease" "ulcerative colitis" "abdominal pain" "intestinal inflammation"; do
  search ibd_3_CMD "$t"; done
# CRC: "diarrhea" replaced by "tumorigenesis" (Oct 2026; diarrhea is non-specific and already a CDI term)
for t in "colon cancer" "colorectal cancer" "tumorigenesis" "rectal bleeding" "MSI subtype"; do
  search crc_11_CMD "$t"; done
# taxon-only publication counts ("popularity"), for the sensitivity analysis (pubtator_sensitivity.py)
for ds in autism_2_microbiomeHD cdi_3_microbiomeHD_clean ibd_3_CMD crc_11_CMD; do
  python3 pubtator_api_search.py --input "$IN/${ds}_taxa_union_wilcoxon.csv" --taxa-col cleaned_taxa \
    --condition "" --output "$OUT/${ds}_taxon_total.csv"
done
echo "ALL KEY-TERM SEARCHES DONE"
