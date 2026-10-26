"""Summary table of all ablation scenarios: mean of each metric per method over iterations (count data).

Rows: original main-grid setting (no relation, OR=1, cond 0.5, batch 0.5) and each ablation of it (n150, halftaxa,
batchfrac40, condfrac40, batchfrac10, condfrac10, zero inflation 0.75-0.90), plus the condsweep (OR=1.25, batch 0.5, no relation).
Uses global_benchmarking_stats_nocovar_batch_0.csv (RF evaluated on data corrected without covariates).
Writes ablation_summary_long.csv, ablation_summary_<metric>.csv and ablation_summary.xlsx (if openpyxl is available)
to <SIM_OUTPUT_ROOT>/scenarios_simulation_data_052626/ablation_summary/.

Usage: python3 ablation_summary_table.py [-i 1000]
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import SIM_OUTPUT_ROOT

ap = argparse.ArgumentParser()
ap.add_argument('-i', '--iterations', type=int, default=1000)
N = ap.parse_args().iterations

B = os.path.join(SIM_OUTPUT_ROOT, 'scenarios_simulation_data_052626')
ORIG = os.path.join(SIM_OUTPUT_ROOT, '082924check', 'simulation_data_eval_new_count_norelation_102023')
SETUPS = [('original (n=450, 301 taxa, 50/50, zeros 0.70)', ORIG, '1_0.5_0.5'),
          ('n150', f'{B}/eval_count_n150_no_relation', '1_0.5_0.5'),
          ('halftaxa', f'{B}/eval_count_halftaxa', '1_0.5_0.5'),
          ('batchfrac40', f'{B}/eval_count_batchfrac40', '1_0.5_0.5'),
          ('condfrac40', f'{B}/eval_count_condfrac40', '1_0.5_0.5'),
          ('batchfrac10', f'{B}/eval_count_batchfrac10', '1_0.5_0.5'),
          ('condfrac10', f'{B}/eval_count_condfrac10', '1_0.5_0.5')] + \
         [(f'zero inflation {l / 100:.2f}', f'{B}/eval_count_zeroinf{l}', '1_0.5_0.5') for l in [75, 80, 85, 90]] + \
         [(f'condsweep cond={c} (OR=1.25)', f'{B}/eval_count_condsweep_no_relation', f'1.25_{c}_0.5') for c in [0.005, 0.01, 0.02, 0.05, 0.1]]
METRICS = {'batch_aitch_r2': 'batch R2 (Aitchison)', 'batch_bray_r2': 'batch R2 (Bray-Curtis)',
           'biovar_aitch_r2': 'condition R2 (Aitchison)', 'biovar_bray_r2': 'condition R2 (Bray-Curtis)',
           'FDR': 'FDR', 'power': 'sensitivity', 'rf_auc': 'RF AUC', 'rf_weighted_f1': 'RF weighted F1'}
ORDER = ['nobc', 'combat_seq', 'limma', 'MMUPHin', 'ConQuR', 'ConQuR_libsize', 'ruv3', 'bmc', 'bat', 'rbe', 'svd',
         'harmony', 'percentile_norm', 'debiasm']

rows = []
for name, ev, tag in SETUPS:
    for it in range(1, N + 1):
        p = f'{ev}/out_{tag}_iter_{it}/global_benchmarking_stats_nocovar_batch_0.csv'
        if not os.path.exists(p):
            continue
        g = pd.read_csv(p, index_col=0).apply(pd.to_numeric, errors='coerce')
        for m in g.columns:
            rows.append(dict(setup=name, iter=it, method=m, **{k: g.at[k, m] for k in METRICS if k in g.index}))
    print(name, flush=True)

d = pd.DataFrame(rows)
out = os.path.join(B, 'ablation_summary')
os.makedirs(out, exist_ok=True)
setup_order = [s[0] for s in SETUPS]
mean = d.groupby(['setup', 'method'])[list(METRICS)].mean()
n_it = d.groupby('setup')['iter'].nunique().reindex(setup_order)
mean.reset_index().to_csv(f'{out}/ablation_summary_long.csv', index=False)
tables = {}
for k, label in METRICS.items():
    t = mean[k].unstack('method').reindex(index=setup_order, columns=[m for m in ORDER if m in mean.index.get_level_values(1)])
    t.insert(0, 'n_iterations', n_it)
    t.round(4).to_csv(f'{out}/ablation_summary_{k}.csv')
    tables[label] = t.round(4)
try:
    with pd.ExcelWriter(f'{out}/ablation_summary.xlsx') as xw:
        for label, t in tables.items():
            t.to_excel(xw, sheet_name=label[:31])
except ImportError:
    print('openpyxl not available; wrote CSVs only')
print(out)
