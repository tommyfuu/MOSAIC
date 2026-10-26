"""Aggregate per-iteration simulation stats (nocovar file) into mean + 95% CI per method and setting, for
fig2_style_plots.py. One call per (datatype, relation, OR); condsweep only exists at OR 1.25.
Usage: python3 fig2_style_aggregate.py <count|relab> <no|yes> <1|1.25|1.5> <out_dir> [n_iter]
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import SIM_OUTPUT_ROOT

dt, rel, OR, out = sys.argv[1:5]
N = int(sys.argv[5]) if len(sys.argv) > 5 else 1000
KEYS = ['batch_aitch_r2', 'batch_bray_r2', 'biovar_aitch_r2', 'biovar_bray_r2', 'FDR', 'power', 'rf_auc', 'rf_weighted_f1']
GRID = [(0.25, 0), (0.5, 0), (0.75, 0), (1, 0), (0.25, 0.25), (0.5, 0.25), (0.75, 0.25), (0.25, 0.5), (0.5, 0.5), (0.25, 0.75)]
SWEEP = [0.005, 0.01, 0.02, 0.05, 0.1]


def fmt(v):
    return str(int(v)) if float(v).is_integer() else str(v)


def collect(base, tags):
    rows = []
    for setting, tag in tags:
        for it in range(1, N + 1):
            p = f'{base}/{tag}_iter_{it}/global_benchmarking_stats_nocovar_batch_0.csv'
            if not os.path.exists(p):
                continue
            g = pd.read_csv(p, index_col=0).apply(pd.to_numeric, errors='coerce')
            for m in g.columns:
                r = dict(setting=setting, it=it, method=m, **{k: g.at[k, m] for k in KEYS if k in g.index})
                r['shan_batch'] = float(g.at['batch_shannon_pval', m] < 0.05)
                r['shan_cond'] = float(g.at['biovar_shannon_pval', m] < 0.05)
                rows.append(r)
    return pd.DataFrame(rows)


def summarize(d):
    out = []
    for (s, m), g in d.groupby(['setting', 'method'], sort=False):
        for k in KEYS + ['shan_batch', 'shan_cond']:
            if k not in g or g[k].isna().all():
                continue
            x = g[k].dropna().values
            if k.startswith('shan'):
                p = x.mean(); ci = 1.96 * np.sqrt(p * (1 - p) / len(x))
            else:
                ci = stats.t.ppf(0.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.nan
            out.append(dict(setting=s, method=m, metric=k, mean=x.mean(), ci=ci, n=len(x)))
    return pd.DataFrame(out)


grid_base = f'{SIM_OUTPUT_ROOT}/082924check/simulation_data_eval_new_{dt}_{rel}relation_102023'
g = collect(grid_base, [(f'{fmt(c)}_{fmt(b)}', f'out_{OR}_{fmt(c)}_{fmt(b)}') for c, b in GRID])
summarize(g).to_csv(f'{out}/grid_{dt}_{rel}_OR{OR}.csv', index=False)
sw_base = f'{SIM_OUTPUT_ROOT}/scenarios_simulation_data_052626/eval_{dt}_condsweep_{rel}_relation'
if os.path.exists(f'{sw_base}/out_{OR}_0.005_0.5_iter_1'):
    sw_base = f'{SIM_OUTPUT_ROOT}/scenarios_simulation_data_052626/eval_{dt}_condsweep_{rel}_relation'
    s = collect(sw_base, [(fmt(c), f'out_{OR}_{fmt(c)}_0.5') for c in SWEEP])
    summarize(s).to_csv(f'{out}/condsweep_{dt}_{rel}_OR{OR}.csv', index=False)
print('done', dt, rel, OR, g.groupby('setting').it.nunique().to_dict())
