"""Zero-inflation ablation: raise the zero fraction of existing simulated count tables.

Starts from the main-grid data (count, no relation, OR=1, cond=0.5, batch=0.5; zero fraction ~0.70) and sets the
smallest non-zero counts to zero (ties broken at random, seeded by iteration and level) until the overall zero
fraction reaches each target level. Each sample's largest count is never zeroed, so no sample becomes empty.
Writes the same file set as the other ablation scenarios to
<SIM_OUTPUT_ROOT>/scenarios_simulation_data_052626/count_zeroinf<level>/.

Usage: python3 generate_zeroinf_ablation.py -i <iter>
"""
import argparse
import os
import shutil
import sys

import numpy as np
import pandas as pd

sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import SIM_OUTPUT_ROOT

LEVELS = [75, 80, 85, 90]
TAG = '1_0.5_0.5'
SRC = os.path.join(SIM_OUTPUT_ROOT, 'simulation_data_MIDAS_1000_norelation_102023')
DST = os.path.join(SIM_OUTPUT_ROOT, 'scenarios_simulation_data_052626')


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


ap = argparse.ArgumentParser()
ap.add_argument('-i', '--iter', type=int, required=True)
it = ap.parse_args().iter
base = f'ibd_150_{{}}_{TAG}_iter_{it}'
counts = pd.read_csv(os.path.join(SRC, base.format('count') + '.csv'))  # no sample-ID column; rows follow the metadata
for level in LEVELS:
    out = os.path.join(DST, f'count_zeroinf{level}')
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng([it, level])
    Z = add_zeros(counts.values, level / 100, rng)
    zc = pd.DataFrame(Z, columns=counts.columns)
    zc.to_csv(os.path.join(out, base.format('count') + '.csv'), index=False)
    zc.div(zc.sum(axis=1), axis=0).to_csv(os.path.join(out, base.format('relab') + '.csv'), index=False)
    for kind, ext in [('meta', 'csv'), ('id_batch', 'txt'), ('id_cond', 'txt')]:
        shutil.copy2(os.path.join(SRC, base.format(kind) + '.' + ext), os.path.join(out, base.format(kind) + '.' + ext))
    print(f'iter {it} level {level}: zero fraction {(Z == 0).mean():.4f}', flush=True)
