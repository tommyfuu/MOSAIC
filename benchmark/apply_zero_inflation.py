"""Raise the zero fraction of simulated count tables in place (zero-inflation scenarios of Snakefile_simulation).

For every count table of one replicate in --data-dir, sets the smallest non-zero counts to zero (ties broken at
random, seeded by replicate and target level) until the overall zero fraction reaches --target; a sample's largest
count is never zeroed, so no sample becomes empty. The matching relative abundance table is rewritten from the
new counts. Same procedure and seeding as generate_zeroinf_ablation.py. Running it again is a no-op.

Usage: python3 apply_zero_inflation.py --data-dir <dir> --iter <i> --target 0.9
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd


def add_zeros(X, target, rng):
    X = X.copy()
    n_extra = int(round(target * X.size)) - int((X == 0).sum())
    if n_extra <= 0:
        return X
    protected = np.zeros_like(X, dtype=bool)
    protected[np.arange(X.shape[0]), X.argmax(axis=1)] = True
    rows, cols = np.where((X > 0) & ~protected)
    order = np.lexsort((rng.random(rows.size), X[rows, cols]))  # smallest counts first, random among ties
    pick = order[:n_extra]
    X[rows[pick], cols[pick]] = 0
    return X


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--data-dir', required=True)
    ap.add_argument('--iter', type=int, required=True)
    ap.add_argument('--target', type=float, required=True, help='target overall zero fraction, e.g. 0.9')
    a = ap.parse_args()
    level = int(round(a.target * 100))
    files = sorted(glob.glob(os.path.join(a.data_dir, f'ibd_150_count_*_iter_{a.iter}.csv')))
    if not files:
        raise SystemExit(f'no count tables for replicate {a.iter} in {a.data_dir}')
    for f in files:
        counts = pd.read_csv(f)
        rng = np.random.default_rng([a.iter, level])
        Z = add_zeros(counts.values, a.target, rng)
        zc = pd.DataFrame(Z, columns=counts.columns)
        zc.to_csv(f, index=False)
        zc.div(zc.sum(axis=1), axis=0).to_csv(f.replace('ibd_150_count_', 'ibd_150_relab_'), index=False)
        print(f'{os.path.basename(f)}: zero fraction {(Z == 0).mean():.4f}', flush=True)


if __name__ == '__main__':
    main()
