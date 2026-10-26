#!/usr/bin/env python3
"""
Patch ablation FDR/power: re-runs only the OLS diff_abund_test with the
correct binarizing_agent='cond_1' (was defaulting to 'H', giving all p=1).
Updates _summary.csv, _summary_nocovar.csv, and the FDR/power rows of
global_benchmarking_stats_batch_0.csv / global_benchmarking_stats_nocovar_batch_0.csv.
All other metrics (RF, PERMANOVA, alpha diversity) are left untouched.

Usage:
    python3 patch_ablation_fdr.py -i <iter_num>
"""

import argparse, os
import pandas as pd
import numpy as np
from scipy import stats as scipy_stats
from statsmodels.stats.multitest import multipletests

BASE      = '/athena/linglab/scratch/chf4012/simulation_outputs/scenarios_simulation_data_052626'
SCENARIOS = ['n150_no_relation', 'halftaxa', 'batchfrac40', 'condfrac40']
AGENT     = 'cond_1'
OR, COND, BATCH = 1, 0.5, 0.5

# methods whose corrected data is not raw counts (skip rel-abund normalization)
NO_RELAB = {'limma', 'harmony', 'debiasm'}


def pvals_vectorized(df, bio_var_arr):
    """
    Compute OLS p-values for y ~ const + x for every column of df at once.
    Uses Pearson r -> t -> p (mathematically identical to statsmodels OLS t-test
    for a single predictor with intercept).
    df shape: (n_samples, n_taxa)
    """
    X = df.values.astype(float)
    y = bio_var_arr.astype(float)
    n = len(y)

    Xc = X - X.mean(axis=0)
    yc = y - y.mean()

    num   = (Xc * yc[:, None]).sum(axis=0)
    denom = np.sqrt((Xc ** 2).sum(axis=0)) * np.sqrt((yc ** 2).sum())

    with np.errstate(divide='ignore', invalid='ignore'):
        r = np.where(denom == 0, 0.0, num / denom)

    with np.errstate(divide='ignore', invalid='ignore'):
        t = r * np.sqrt(n - 2) / np.sqrt(np.clip(1 - r ** 2, 1e-15, None))

    p = 2 * scipy_stats.t.sf(np.abs(t), df=n - 2)
    return np.where(np.isnan(p), 1.0, p)


def rerun_diff_abund(corrected_csv, meta_csv, gt_path, method):
    """
    Re-run diff_abund_test with correct binarizing agent.
    Returns (fdr_val, power_val, out_df) or None if data is unusable.
    """
    if not os.path.exists(corrected_csv) or not os.path.exists(meta_csv):
        return None

    X = pd.read_csv(corrected_csv, index_col=0)
    meta = pd.read_csv(meta_csv, index_col=0)
    meta.index = X.index

    df = X.copy()

    # relative-abundance normalization for count methods that output non-negative data
    if method not in NO_RELAB and df.values.min() >= 0:
        rs = df.sum(axis=1)
        if (rs == 0).any():
            return None
        df = df.div(rs, axis=0)

    # drop NaN/Inf — OLS handles scale, so no upper-bound check
    if df.isnull().values.any() or np.isinf(df.values).any():
        return None

    bio_var_arr = np.array([1 if c == AGENT else 0 for c in meta['cond'].values])
    taxa_names  = df.columns

    pvals     = pvals_vectorized(df, bio_var_arr)
    fdr_pvals = multipletests(pvals, method='fdr_bh')[1]

    # build output in taxa × samples layout (same as original evaluate code)
    out = X[taxa_names].T.copy()
    out['p_value']     = pvals
    out['FDR_p_value'] = fdr_pvals

    taxa_gt = None
    if os.path.exists(gt_path):
        with open(gt_path) as f:
            line = f.readline().strip()
            if line:
                taxa_gt = set(int(t) - 1 for t in line.split())

    if taxa_gt is not None:
        out['ground_truth'] = [1 if i in taxa_gt else 0
                               for i in range(len(taxa_names))]
        TP = int(((out['ground_truth'] == 1) & (out['FDR_p_value'] < 0.05)).sum())
        FP = int(((out['ground_truth'] == 0) & (out['FDR_p_value'] < 0.05)).sum())
        FN = int(((out['ground_truth'] == 1) & (out['FDR_p_value'] >= 0.05)).sum())
        fdr_val   = FP / (TP + FP) if (TP + FP) > 0 else 0.0
        power_val = TP / (TP + FN) if (TP + FN) > 0 else 0.0
    else:
        fdr_val, power_val = 0.0, 0.0

    return fdr_val, power_val, out


def patch_global_stats(eval_dir, fdr_power_dict, suffix=''):
    """Update only FDR and power rows in the existing global stats CSV."""
    gs_path = os.path.join(eval_dir, f'global_benchmarking_stats{suffix}_batch_0.csv')
    if not os.path.exists(gs_path):
        return
    gs = pd.read_csv(gs_path, index_col=0)
    for method, (fdr_val, power_val) in fdr_power_dict.items():
        if method in gs.columns:
            gs.loc['FDR',   method] = fdr_val
            gs.loc['power', method] = power_val
    gs.to_csv(gs_path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--iter', type=int, required=True)
    args = parser.parse_args()
    iter_num = args.iter

    for scenario in SCENARIOS:
        data_dir   = f'{BASE}/count_{scenario}'
        output_dir = f'{BASE}/output_count_{scenario}'
        eval_dir   = f'{BASE}/eval_count_{scenario}/out_{OR}_{COND}_{BATCH}_iter_{iter_num}'

        if not os.path.isdir(eval_dir):
            print(f"[skip] {scenario} iter {iter_num}: eval dir missing", flush=True)
            continue

        gt_path  = f'{data_dir}/ibd_150_id_cond_{OR}_{COND}_{BATCH}_iter_{iter_num}.txt'
        meta_csv = f'{data_dir}/ibd_150_meta_{OR}_{COND}_{BATCH}_iter_{iter_num}.csv'
        out_stem = (f'{output_dir}/out_{OR}_{COND}_{BATCH}_iter_{iter_num}'
                    f'/ibd_{OR}_{COND}_{BATCH}_iter_{iter_num}')

        corrected_map = {
            'nobc':            f'{data_dir}/ibd_150_count_{OR}_{COND}_{BATCH}_iter_{iter_num}.csv',
            'harmony':         out_stem + '_harmony_adjusted_count.csv',
            'percentile_norm': out_stem + '_percentile_norm.csv',
            'debiasm':         out_stem + '_debiasm_adjusted_count.csv',
        }

        methods_in_dir = [d for d in os.listdir(eval_dir)
                          if os.path.isdir(os.path.join(eval_dir, d))]

        fdr_power = {}

        for method in methods_in_dir:
            method_dir      = os.path.join(eval_dir, method)
            summary_path    = os.path.join(method_dir, f'{method}_summary.csv')
            summary_nc_path = os.path.join(method_dir, f'{method}_summary_nocovar.csv')
            diff_path       = os.path.join(method_dir, f'{method}_diff_abund_test.csv')

            if not os.path.exists(summary_path):
                continue

            corrected_csv = corrected_map.get(method, out_stem + f'_{method}.csv')

            result = rerun_diff_abund(corrected_csv, meta_csv, gt_path, method)
            if result is None:
                print(f"  [{scenario}] {method}: skipped (missing/bad data)", flush=True)
                continue

            fdr_val, power_val, out_df = result
            fdr_power[method] = (fdr_val, power_val)

            out_df.to_csv(diff_path)

            for spath in [summary_path, summary_nc_path]:
                if os.path.exists(spath):
                    s = pd.read_csv(spath, index_col=0)
                    s.loc['FDR',   '0'] = fdr_val
                    s.loc['power', '0'] = power_val
                    s.to_csv(spath)

            print(f"  [{scenario}] {method}: FDR={fdr_val:.3f}, power={power_val:.3f}", flush=True)

        patch_global_stats(eval_dir, fdr_power, suffix='')
        patch_global_stats(eval_dir, fdr_power, suffix='_nocovar')
        print(f"  [{scenario}] global stats patched", flush=True)

    print(f"DONE iter {iter_num}", flush=True)


if __name__ == '__main__':
    main()
