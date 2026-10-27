#!/bin/bash
# Re-evaluate DEBIAS-M random forest metrics after the within-fold input fix (evaluate.py --refresh-debiasm-rf).
# One task = one line of the task file: <datatype> <eval_dir> <data_dir> <output_dir> <or_l> <cond_l> <batch_l> <iter>.
# The original (pre-fix) DEBIAS-M RF files and the per-replicate global stats are tarred to BACKUP before overwriting.
# Usage: sbatch --array=1-N%M rerun_debiasm_rf_fix.sh <task_file> [line_offset]
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --exclude=c12,c41
#SBATCH -o /athena/linglab/scratch/chf4012/simulation_outputs/backup_debiasm_rf_bug_20261026/logs/%A_%a.out
#SBATCH -e /athena/linglab/scratch/chf4012/simulation_outputs/backup_debiasm_rf_bug_20261026/logs/%A_%a.out

export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ENV=/athena/linglab/scratch/chf4012/miniconda3/envs/mic_data_integration
export PATH=$ENV/bin:$PATH
PY=$ENV/bin/python
BENCH=/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark
SIM=/athena/linglab/scratch/chf4012/simulation_outputs
BACKUP=$SIM/backup_debiasm_rf_bug_20261026
TASKS=$1
LINE=$(( ${SLURM_ARRAY_TASK_ID:-1} + ${2:-0} ))
read -r DT EVAL DATA OUT ORL CONDL BATCHL ITER <<< "$(sed -n "${LINE}p" "$TASKS")"
[ -z "$ITER" ] && { echo "no task on line $LINE"; exit 1; }

DONE=$BACKUP/done/$(basename "$EVAL")
mkdir -p "$DONE" "$BACKUP/$(basename "$EVAL")"
TAG=or${ORL}_iter_$ITER
[ -e "$DONE/$TAG" ] && { echo "already done: $EVAL iter $ITER"; exit 0; }

# back up once; a rerun after a failure keeps the first (original) backup
TAR=$BACKUP/$(basename "$EVAL")/$TAG.tar.gz
if [ ! -s "$TAR" ]; then
    cd "$EVAL" || exit 1
    ls -d out_*_iter_"$ITER" > /dev/null 2>&1 || { echo "no outputs for iter $ITER"; exit 1; }
    find $(for o in ${ORL//,/ }; do ls -d out_${o}_*_iter_"$ITER"; done) -maxdepth 2 \( -path '*/debiasm/debiasm_summary*.csv' -o -path '*/debiasm/debiasm_cond_rf_*' \
        -o -name 'global_benchmarking_stats*_batch_0.csv' \) -print0 | tar --null -czf "$TAR.part" -T - && mv "$TAR.part" "$TAR" || exit 1
fi

cd "$BENCH" || exit 1
$PY evaluate.py -o 17 -d "$DT" -i "$ITER" --data-dir "$DATA" --output-dir "$OUT" --eval-dir "$EVAL" \
    --or-l "$ORL" --cond-l "$CONDL" --batch-l "$BATCHL" --rerun-rf --refresh-debiasm-rf \
    && touch "$DONE/$TAG"
