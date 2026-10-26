#!/bin/bash
# Regression tests for the simulation workflow (Snakefile_simulation): the shared code must reproduce the scripts
# used for the manuscript. Run from the repository root inside the mic_data_integration environment:
#   bash benchmark/tests/test_simulation_workflow.sh [test_dir]          # quick (about 3 min)
#   FULL=1 bash benchmark/tests/test_simulation_workflow.sh [test_dir]   # also integration + evaluation (about 1 h)
# Checks
#   1. simulate_scenario.R == generate_data_MIDASim_imbalance.R (batchfrac20, replicate 1): byte-identical tables
#   2. apply_zero_inflation.py: reaches the target zero fraction, never empties a sample, is idempotent
#   3. (FULL) run_sim_methods.R == methods_benchmarking_sim_scenario.R (all R strategies byte-identical), and
#      evaluate.py -o 16/17 == -o 13/14 for every reported metric (DEBIAS-M and permutation p-values are unseeded)
#   4. Snakefile_simulation dry run resolves every scenario in the config
set -uo pipefail
cd "$(dirname "$0")/../.."
REPO=$(pwd)
T=${1:-$(mktemp -d -p "${SIM_OUTPUT_ROOT:-$REPO}" mosaic_test_XXXX)}
mkdir -p "$T"
fail=0
ok()  { echo "PASS  $1"; }
bad() { echo "FAIL  $1"; fail=1; }

echo "test folder: $T"
# 1. simulation
( cd benchmark && SIM_OUTPUT_ROOT="$T/orig" Rscript generate_data_MIDASim_imbalance.R 1 > "$T/orig_sim.log" 2>&1 )
( cd benchmark && Rscript simulate_scenario.R --out="$T/new/data" --iter=1 --p_batch=0.2 > "$T/new_sim.log" 2>&1 )
O="$T/orig/scenarios_simulation_data_052626/count_batchfrac20"
for f in count_1_0.5_0.5_iter_1.csv relab_1_0.5_0.5_iter_1.csv meta_1_0.5_0.5_iter_1.csv id_cond_1_0.5_0.5_iter_1.txt id_batch_1_0.5_0.5_iter_1.txt; do
    cmp -s "$O/ibd_150_$f" "$T/new/data/ibd_150_$f" && ok "simulation identical: $f" || bad "simulation differs: $f"
done

# 2. zero inflation
mkdir -p "$T/zero"; cp "$T/new/data/ibd_150_count_1_0.5_0.5_iter_1.csv" "$T/zero/"
python3 benchmark/apply_zero_inflation.py --data-dir "$T/zero" --iter 1 --target 0.9 > /dev/null
cp "$T/zero/ibd_150_count_1_0.5_0.5_iter_1.csv" "$T/zero/first.csv"
python3 benchmark/apply_zero_inflation.py --data-dir "$T/zero" --iter 1 --target 0.9 > /dev/null
cmp -s "$T/zero/first.csv" "$T/zero/ibd_150_count_1_0.5_0.5_iter_1.csv" && ok "zero inflation idempotent" || bad "zero inflation not idempotent"
python3 - "$T/zero/ibd_150_count_1_0.5_0.5_iter_1.csv" <<'EOF' && ok "zero fraction 0.90, no empty sample" || bad "zero inflation check"
import sys, pandas as pd
x = pd.read_csv(sys.argv[1]).values
sys.exit(0 if abs((x == 0).mean() - 0.9) < 1e-4 and (x.sum(axis=1) > 0).all() else 1)
EOF

# 3. integration and evaluation (slow)
if [ "${FULL:-0}" = 1 ]; then
    M=combat_seq,limma,MMUPHin,ConQuR,ConQuR_libsize,ruv3,bmc,bat,rbe,svd
    ( cd benchmark && SIM_OUTPUT_ROOT="$T/orig" Rscript methods_benchmarking_sim_scenario.R 1 batchfrac20 count > "$T/orig_methods.log" 2>&1 )
    ( cd benchmark && Rscript run_sim_methods.R --data="$T/new/data" --out="$T/new/output_count" --iter=1 --datatype=count --methods=$M > "$T/new_methods.log" 2>&1 )
    for v in "" _nocovar; do
        A="$O/../output_count_batchfrac20$v/out_1_0.5_0.5_iter_1"; B="$T/new/output_count$v/out_1_0.5_0.5_iter_1"
        for f in $(ls "$A" | grep -v runtime); do cmp -s "$A/$f" "$B/$f" && ok "R integration identical: output_count$v/$f" || bad "R integration differs: output_count$v/$f"; done
    done
    ( cd benchmark && SIM_OUTPUT_ROOT="$T/orig" python3 evaluate.py -o 13 -i 1 -d count --scenario batchfrac20 > "$T/orig_py.log" 2>&1
      SIM_OUTPUT_ROOT="$T/orig" python3 evaluate.py -o 14 -i 1 -d count --scenario batchfrac20 > "$T/orig_eval.log" 2>&1
      SIM_OUTPUT_ROOT="$T/orig" python3 evaluate.py -o 14 -i 1 -d count --scenario batchfrac20 --rerun-rf >> "$T/orig_eval.log" 2>&1 )
    G="--config config_simulation_scenarios.yml -i 1 -d count --data-dir $T/new/data --output-dir $T/new/output_count --eval-dir $T/new/eval_count"
    ( cd benchmark && python3 evaluate.py -o 16 $G > "$T/new_py.log" 2>&1 && python3 evaluate.py -o 17 $G > "$T/new_eval.log" 2>&1 \
      && python3 evaluate.py -o 17 $G --rerun-rf >> "$T/new_eval.log" 2>&1 )
    for s in "" _nocovar; do
        python3 - "$O/../eval_count_batchfrac20/out_1_0.5_0.5_iter_1/global_benchmarking_stats${s}_batch_0.csv" \
                  "$T/new/eval_count/out_1_0.5_0.5_iter_1/global_benchmarking_stats${s}_batch_0.csv" <<'EOF' \
            && ok "evaluation identical (reported metrics): global_benchmarking_stats${s}" || bad "evaluation differs: global_benchmarking_stats${s}"
import sys, numpy as np, pandas as pd
# runtime, PERMANOVA permutation p-values (unseeded permutations) and DEBIAS-M (unseeded training) are not compared
skip = ['runtime', 'aitch_pval_batch', 'aitch_pval_biovar', 'bray_pval_batch', 'bray_pval_biovar']
a, b = (pd.read_csv(f, index_col=0).drop(index=skip, errors='ignore') for f in sys.argv[1:3])
a = a.drop(columns=[c for c in ['debiasm'] if c in a])
b = b[a.columns].loc[a.index]
diff = (a.apply(pd.to_numeric, errors='coerce') - b.apply(pd.to_numeric, errors='coerce')).abs().max().max()
print('max abs difference', diff); sys.exit(0 if diff < 1e-9 else 1)
EOF
    done
fi

# 4. workflow dry run
snakemake -s Snakefile_simulation -n -q --config output_root="$T/smk" iterations=2 \
    'run_scenarios=[main_grid_no_libsize_confounding,main_grid_libsize_confounded,condsweep,baseline,n150,halftaxa,zeroinf75,zeroinf80,zeroinf85,zeroinf90,batchfrac40,batchfrac20,batchfrac10,condfrac40,condfrac20,condfrac10]' \
    'datatypes=[count,relab]' > "$T/dryrun.log" 2>&1 && ok "workflow dry run (all scenarios, count + relab)" || bad "workflow dry run (see $T/dryrun.log)"

[ $fail = 0 ] && echo "ALL TESTS PASSED" || echo "SOME TESTS FAILED"
exit $fail
