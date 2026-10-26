#!/usr/bin/env python3
"""
Check which iterations have incomplete evaluation results.
Usage: python3 check_eval_completeness_052626.py <eval_dir> <expected_methods_csv> [max_iter]

Example:
  python3 check_eval_completeness_052626.py \
    /athena/linglab/scratch/chf4012/simulation_outputs/082924check/simulation_data_eval_new_count_norelation_102023 \
    "nobc,combat_seq,limma,MMUPHin,ConQuR,ConQuR_libsize,ruv3,bmc,bat,rbe,svd,harmony,percentile_norm,debiasm" \
    100
"""
import os, sys, csv
from collections import defaultdict

eval_dir   = sys.argv[1]
expected   = set(sys.argv[2].split(',')) if len(sys.argv) > 2 else set()
max_iter   = int(sys.argv[3]) if len(sys.argv) > 3 else None

# # for counts
# expected = "nobc,combat_seq,limma,MMUPHin,ConQuR,ConQuR_libsize,ruv3,bmc,bat,rbe,svd,harmony,percentile_norm,debiasm".split(',') if not expected else expected
# for relab
# expected = "nobc,combat,limma,MMUPHin,ConQuR_rel,ruv3,bmc,bat,rbe,svd,harmony,percentile_norm,debiasm".split(',') if not expected else expected
# (mic_data_integration) [chf4012@g40 benchmark]$ python3 check_eval_completeness_052626.py /athena/linglab/scratch/chf4012/simulation_outputs/082924check/simulation_data_eval_new_count_yesrelation_102023 "nobc,combat_seq,limma,MMUPHin,ConQuR,ConQuR_libsize,ruv3,bmc,bat,rbe,svd,harmony,percentile_norm,debiasm" 100
# python3 check_eval_completeness_052626.py /athena/linglab/scratch/chf4012/simulation_outputs/082924check/simulation_data_eval_new_count_norelation_102023 "nobc,combat_seq,limma,MMUPHin,ConQuR,ConQuR_libsize,ruv3,bmc,bat,rbe,svd,harmony,percentile_norm,debiasm" 100
# ITER     N_COMBOS_MISSING   MISSING_METHODS
# ------------------------------------------------------------
# ?        4                  ConQuR, ConQuR_libsize, MMUPHin, bat, bmc, combat_seq, debiasm, harmony, limma, nobc, percentile_norm, rbe, ruv3, svd

# Total incomplete iterations: 1
# (mic_data_integration) [chf4012@g40 benchmark]$ python3 check_eval_completeness_052626.py /athena/linglab/scratch/chf4012/simulation_outputs/082924check/simulation_data_eval_new_relab_yesrelation_102023 "nobc,combat,limma,MMUPHin,ConQuR_rel,ruv3,bmc,bat,rbe,svd,harmony,percentile_norm,debiasm" 100
# python3 check_eval_completeness_052626.py /athena/linglab/scratch/chf4012/simulation_outputs/082924check/simulation_data_eval_new_relab_norelation_102023 "nobc,combat,limma,MMUPHin,ConQuR_rel,ruv3,bmc,bat,rbe,svd,harmony,percentile_norm,debiasm" 100
# iter -> {combo -> set of missing methods}
missing_by_iter = defaultdict(dict)

for combo_dir in sorted(os.listdir(eval_dir)):
    # filter by max_iter if specified
    if max_iter is not None and '_iter_' in combo_dir:
        iter_num_check = combo_dir.split('_iter_')[-1]
        if iter_num_check.isdigit() and int(iter_num_check) > max_iter:
            continue
    csv_path = os.path.join(eval_dir, combo_dir, 'global_benchmarking_stats_batch_0.csv')
    if not os.path.isfile(csv_path):
        # whole combo missing — extract iter from dir name
        iter_num = combo_dir.split('_iter_')[-1] if '_iter_' in combo_dir else '?'
        missing_by_iter[iter_num][combo_dir] = expected or {'ENTIRE CSV MISSING'}
        continue
    with open(csv_path) as f:
        cols = set(next(csv.reader(f)))
    cols.discard('')
    methods_to_check = expected if expected else cols  # if no expected given, just report what's there
    missing = expected - cols if expected else set()
    if missing:
        iter_num = combo_dir.split('_iter_')[-1] if '_iter_' in combo_dir else '?'
        missing_by_iter[iter_num][combo_dir] = missing

if not missing_by_iter:
    print("All iterations complete.")
else:
    # summarise: for each iter, which methods are missing across any combo
    print(f"{'ITER':<8} {'N_COMBOS_MISSING':<18} {'MISSING_METHODS'}")
    print('-'*60)
    for iter_num in sorted(missing_by_iter, key=lambda x: int(x) if x.isdigit() else 0):
        combos = missing_by_iter[iter_num]
        all_missing = set().union(*combos.values())
        print(f"{iter_num:<8} {len(combos):<18} {', '.join(sorted(all_missing))}")
    print(f"\nTotal incomplete iterations: {len(missing_by_iter)}")
