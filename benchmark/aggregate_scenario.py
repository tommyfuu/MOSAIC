"""Summarize the evaluation of a simulation scenario over replicates (last step of Snakefile_simulation).

Reads <eval-dir>/out_<or>_<cond>_<batch>_iter_<i>/global_benchmarking_stats_nocovar_batch_0.csv for every setting
and replicate (falls back to the regular file when the leakage-free one is missing; the two differ only in the
random forest metrics) and writes
  --values   one row per setting x replicate x method, all metrics; Shannon p-values are also turned into rejection
             indicators (shannon_batch_rej, shannon_cond_rej: p < 0.05)
  --summary  one row per setting x method x metric: mean, 95% CI half-width (t-based; binomial for rejection
             proportions), SD, median, 2.5% and 97.5% percentiles, number of replicates

Usage: python3 aggregate_scenario.py --eval-dir <dir> --or-l 1 --cond-l 0.5 --batch-l 0.5 --n-iter 1000 \
           --values values.csv --summary summary.csv
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats


def vals(x):
    return [int(float(v)) if float(v).is_integer() else float(v) for v in str(x).split(',')]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--eval-dir', required=True)
    ap.add_argument('--or-l', default='1')
    ap.add_argument('--cond-l', default='0.5')
    ap.add_argument('--batch-l', default='0.5')
    ap.add_argument('--n-iter', type=int, required=True)
    ap.add_argument('--values', required=True)
    ap.add_argument('--summary', required=True)
    a = ap.parse_args()

    rows, missing = [], 0
    for o in vals(a.or_l):
        for c in vals(a.cond_l):
            for b in vals(a.batch_l):
                if c + b > 1:
                    continue
                for it in range(1, a.n_iter + 1):
                    d = f'{a.eval_dir}/out_{o}_{c}_{b}_iter_{it}'
                    f = next((p for p in (f'{d}/global_benchmarking_stats_nocovar_batch_0.csv',
                                          f'{d}/global_benchmarking_stats_batch_0.csv') if os.path.exists(p)), None)
                    if f is None:
                        missing += 1
                        continue
                    g = pd.read_csv(f, index_col=0).apply(pd.to_numeric, errors='coerce')
                    for m in g.columns:
                        r = dict(odds_ratio=o, cond_effect=c, batch_effect=b, iteration=it, method=m, **g[m].to_dict())
                        for key, new in (('batch_shannon_pval', 'shannon_batch_rej'), ('biovar_shannon_pval', 'shannon_cond_rej')):
                            if key in g.index and pd.notna(g.at[key, m]):
                                r[new] = float(g.at[key, m] < 0.05)
                        rows.append(r)
    if not rows:
        raise SystemExit(f'no evaluation files found under {a.eval_dir}')
    v = pd.DataFrame(rows)
    v.to_csv(a.values, index=False)

    keys = ['odds_ratio', 'cond_effect', 'batch_effect', 'method']
    metrics = [c for c in v.columns if c not in keys + ['iteration']]
    out = []
    for k, grp in v.groupby(keys, sort=False):
        for met in metrics:
            x = grp[met].dropna().to_numpy(dtype=float)
            if x.size == 0:
                continue
            if met.endswith('_rej'):
                p = x.mean()
                ci = 1.96 * np.sqrt(p * (1 - p) / x.size)
            else:
                ci = stats.t.ppf(0.975, x.size - 1) * x.std(ddof=1) / np.sqrt(x.size) if x.size > 1 else np.nan
            out.append(dict(zip(keys, k), metric=met, mean=x.mean(), ci95=ci, sd=x.std(ddof=1) if x.size > 1 else np.nan,
                            median=np.median(x), p2_5=np.percentile(x, 2.5), p97_5=np.percentile(x, 97.5), n=x.size))
    pd.DataFrame(out).to_csv(a.summary, index=False)
    print(f'{len(v)} rows from {v.iteration.nunique()} replicates; {missing} setting-replicates missing -> {a.summary}')


if __name__ == '__main__':
    main()
