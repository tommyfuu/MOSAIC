"""Aggregate per-iteration integration runtimes on the main simulation grid (mean + 95% CI, median, quartiles per
method and setting) for supp_runtime_plot.py. Runtimes come from the 'runtime' row of each iteration's
global_benchmarking_stats_nocovar_batch_0.csv (seconds; nobc has none). One call per (datatype, relation, OR).
Usage: python3 runtime_aggregate.py <count|relab> <no|yes> <1|1.25|1.5> <out_dir> [n_iter]
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # repository root
from paths import SIM_OUTPUT_ROOT

dt, rel, OR, out = sys.argv[1:5]
N = int(sys.argv[5]) if len(sys.argv) > 5 else 1000
GRID = [(0.25, 0), (0.5, 0), (0.75, 0), (1, 0), (0.25, 0.25), (0.5, 0.25), (0.75, 0.25), (0.25, 0.5), (0.5, 0.5), (0.25, 0.75)]
os.makedirs(out, exist_ok=True)


def fmt(v):
    return str(int(v)) if float(v).is_integer() else str(v)


base = f'{SIM_OUTPUT_ROOT}/082924check/simulation_data_eval_new_{dt}_{rel}relation_102023'
rows = []
for c, b in GRID:
    setting = f'{fmt(c)}_{fmt(b)}'
    for it in range(1, N + 1):
        p = f'{base}/out_{OR}_{fmt(c)}_{fmt(b)}_iter_{it}/global_benchmarking_stats_nocovar_batch_0.csv'
        if not os.path.exists(p):
            continue
        g = pd.read_csv(p, index_col=0)
        if 'runtime' not in g.index:
            continue
        for m, v in pd.to_numeric(g.loc['runtime'], errors='coerce').items():
            if np.isfinite(v):
                rows.append(dict(setting=setting, it=it, method=m, runtime=v))
d = pd.DataFrame(rows)
summ = []
for (s, m), g in d.groupby(['setting', 'method'], sort=False):
    x = g.runtime.values
    ci = stats.t.ppf(0.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.nan
    summ.append(dict(setting=s, method=m, mean=x.mean(), ci=ci, median=np.median(x), q25=np.percentile(x, 25),
                     q75=np.percentile(x, 75), n=len(x)))
pd.DataFrame(summ).to_csv(f'{out}/runtime_{dt}_{rel}_OR{OR}.csv', index=False)
print('done', dt, rel, OR, d.groupby('setting').it.nunique().to_dict())
