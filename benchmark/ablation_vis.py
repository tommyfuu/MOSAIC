"""Alternative visualizations of the ablation scenarios vs the original setting (count or relab, no relation, OR=1,
cond 0.5, batch 0.5), using global_benchmarking_stats_nocovar_batch_0.csv over iterations 1..N.

Plots per metric include:
  original_style_<metric>.pdf   points per scenario, as in evaluate.py -o 8
  reference_ci_<metric>.pdf     one panel per method; original mean as a dashed line with its 95% CI band,
                                each ablation as a point with its 95% CI
  difference_<metric>.pdf       ablation minus original, with 95% CI (Welch), per method; line at 0
  absolute_reference_spread_<metric>.pdf   as absolute_reference, but bars = 2.5-97.5 percentile across iterations
  highlight_{ci,spread}_<metric>.pdf       as absolute_reference with plain points; only the best method per setup
                                           (BEST_RULE) gets a bar (95% CI of the mean, or the 2.5-97.5 percentile spread)

Usage: python3 ablation_vis.py -i 50 [-d relab] [--setups n150,zeroinf90] [--out DIR] [--reuse]
       python3 ablation_vis.py -i 1000 --root <Snakefile_simulation output_root> --setups n150,zeroinf90 --out DIR
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42  # embed TrueType so text stays editable in Illustrator
matplotlib.rcParams['ps.fonttype'] = 42
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import SIM_OUTPUT_ROOT

ap = argparse.ArgumentParser()
ap.add_argument('-i', '--iterations', type=int, default=50)
ap.add_argument('--out', default=None)
ap.add_argument('-d', '--datatype', default='count', choices=['count', 'relab'])
ap.add_argument('--setups', default=None, help='comma-separated ablation scenarios to compare with the original (default: n150,halftaxa,batchfrac40,condfrac40,zeroinf90)')
ap.add_argument('--reuse', action='store_true', help='read OUT/values_per_iteration.csv instead of the per-iteration stats files')
ap.add_argument('--root', default=None, help='scenario folder of Snakefile_simulation: compare <root>/<setup>/eval_<datatype> with <root>/baseline/eval_<datatype>')
args = ap.parse_args()
N = args.iterations
B = os.path.join(SIM_OUTPUT_ROOT, 'scenarios_simulation_data_052626')
DT = args.datatype
OUT = args.out or os.path.join(B, f'ablation_vis_test_{N}iters' + ('' if DT == 'count' else f'_{DT}'))
os.makedirs(OUT, exist_ok=True)

if args.root:  # Snakefile_simulation layout
    SETUPS = {'original': f'{args.root}/baseline/eval_{DT}',
              **{s: f'{args.root}/{s}/eval_{DT}' for s in (args.setups or 'n150,halftaxa,zeroinf90,batchfrac10,condfrac10').split(',')}}
else:
    SETUPS = {'original': os.path.join(SIM_OUTPUT_ROOT, '082924check', f'simulation_data_eval_new_{DT}_norelation_102023'),
              **{s: (f'{B}/eval_{DT}_n150_no_relation' if s == 'n150' else f'{B}/eval_{DT}_{s}')
                 for s in (args.setups or 'n150,halftaxa,batchfrac40,condfrac40,zeroinf90').split(',')}}
ABLATIONS = [s for s in SETUPS if s != 'original']
# succinct method set (same as evaluate.py --succinct: no bat/ruv3/svd/rbe/bmc)
METHODS = (['nobc', 'combat_seq', 'limma', 'MMUPHin', 'ConQuR', 'ConQuR_libsize', 'harmony', 'percentile_norm', 'debiasm']
           if DT == 'count' else
           ['nobc', 'combat', 'limma', 'MMUPHin', 'ConQuR_rel', 'harmony', 'percentile_norm', 'debiasm'])
METRICS = {'power': 'Sensitivity', 'FDR': 'FDR', 'batch_aitch_r2': 'PERMANOVA batch R2 (Aitchison)',
           'biovar_aitch_r2': 'PERMANOVA condition R2 (Aitchison)', 'batch_bray_r2': 'PERMANOVA batch R2 (Bray-Curtis)',
           'biovar_bray_r2': 'PERMANOVA condition R2 (Bray-Curtis)', 'rf_auc': 'RF AUC',
           'rf_weighted_f1': 'RF weighted F1',
           'batch_shannon_rej': 'Shannon diversity vs batch: proportion p < 0.05',
           'biovar_shannon_rej': 'Shannon diversity vs condition: proportion p < 0.05'}
SHANNON_P = {'batch_shannon_rej': 'batch_shannon_pval', 'biovar_shannon_rej': 'biovar_shannon_pval'}
METHOD_COLORS = {'nobc': 'black', 'combat_seq': 'gold', 'combat': 'gold', 'limma': 'blue', 'MMUPHin': 'green', 'ConQuR': 'orange',
                 'ConQuR_rel': 'orange', 'ConQuR_libsize': 'purple', 'bmc': 'coral', 'harmony': 'turquoise', 'percentile_norm': 'gray', 'debiasm': 'magenta'}
METHOD_MARKERS = {'nobc': 'o', 'combat_seq': 'D', 'combat': 'D', 'limma': '^', 'MMUPHin': '<', 'ConQuR': '>', 'ConQuR_rel': '>',
                  'ConQuR_libsize': 's',
                  'bmc': 'H', 'harmony': 'P', 'percentile_norm': '*', 'debiasm': 'd'}
ABL_COLORS = dict(zip(ABLATIONS, ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#666666']))
ABL_MARKERS = dict(zip(ABLATIONS, ['o', 's', '^', 'v', 'D', 'P', 'X', 'h']))


def setup_legend(ax, extra=()):
    from matplotlib.lines import Line2D
    h = list(extra) + [Line2D([], [], color='black', marker=ABL_MARKERS[s], ls='none', ms=6, label=s) for s in ABLATIONS]
    ax.legend(handles=h, bbox_to_anchor=(1.02, 0.5), loc='center left', frameon=False, title='setup')


def color_ticks(ax):
    for t in ax.get_xticklabels():
        t.set_color(METHOD_COLORS.get(t.get_text(), 'black') if t.get_text() not in ('combat_seq', 'combat') else 'darkgoldenrod')

rows = []
for s, ev in ({} if args.reuse else SETUPS).items():
    for it in range(1, N + 1):
        p = f'{ev}/out_1_0.5_0.5_iter_{it}/global_benchmarking_stats_nocovar_batch_0.csv'
        if os.path.exists(p):
            g = pd.read_csv(p, index_col=0).apply(pd.to_numeric, errors='coerce')
            for m in METHODS:
                if m in g.columns:
                    rows.append(dict(setup=s, iter=it, method=m, **{k: g.at[k, m] for k in list(METRICS) + list(SHANNON_P.values()) if k in g.index}))
d = pd.read_csv(f'{OUT}/values_per_iteration.csv') if args.reuse else pd.DataFrame(rows)
for k, pcol in SHANNON_P.items():  # alpha diversity: per-iteration rejection at p < 0.05, averaged = rejection proportion
    if pcol in d.columns:
        d[k] = (d[pcol] < 0.05).astype(float).where(d[pcol].notna())
METRICS = {k: v for k, v in METRICS.items() if k in d.columns}  # relab evals have no Aitchison rows
d.to_csv(f'{OUT}/values_per_iteration.csv', index=False)


def mean_ci(x):
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    if len(x) < 2:  # a single replicate: mean only, no CI
        return (x.mean() if len(x) else np.nan), np.nan, len(x)
    return x.mean(), stats.t.ppf(0.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x)), len(x)


summ = []
for (s, m), g in d.groupby(['setup', 'method']):
    for k in METRICS:
        mu, hw, n = mean_ci(g[k])
        lo, hi = np.nanpercentile(g[k], [2.5, 97.5]) if g[k].notna().any() else (np.nan, np.nan)
        summ.append(dict(setup=s, method=m, metric=k, mean=mu, ci=hw, n=n, sd=np.nanstd(g[k], ddof=1), p2_5=lo, p97_5=hi))
summ = pd.DataFrame(summ).set_index(['setup', 'method', 'metric'])
summ.reset_index().to_csv(f'{OUT}/summary_mean_ci.csv', index=False)


def get(s, m, k):
    try:
        r = summ.loc[(s, m, k)]
        return r['mean'], r['ci'], r['sd'], r['n']
    except KeyError:
        return np.nan, np.nan, np.nan, 0


def get_spread(s, m, k):
    try:
        r = summ.loc[(s, m, k)]
        return r['p2_5'], r['p97_5']
    except KeyError:
        return np.nan, np.nan


# best method per setup, used for the highlight plots; nobc is never chosen
BEST_RULE = {'FDR': 'closest to 0.05', 'power': 'highest', 'rf_auc': 'highest', 'rf_weighted_f1': 'highest',
             'batch_shannon_rej': 'lowest', 'biovar_shannon_rej': 'highest', 'batch_aitch_r2': 'lowest',
             'batch_bray_r2': 'lowest', 'biovar_aitch_r2': 'closest to nobc', 'biovar_bray_r2': 'closest to nobc'}


def best_method(s, k):
    mus = {m: get(s, m, k)[0] for m in METHODS if m != 'nobc' and not np.isnan(get(s, m, k)[0])}
    rule = BEST_RULE[k]
    if rule == 'highest':
        return max(mus, key=mus.get)
    if rule == 'lowest':
        return min(mus, key=mus.get)
    target = 0.05 if rule == 'closest to 0.05' else get(s, 'nobc', k)[0]
    return min(mus, key=lambda m: abs(mus[m] - target))


def absolute_reference(k, label, bars, highlight, fname):
    # x = setup; each method's original mean as a thin dashed line in its color; at each ablation a thin dotted vertical
    # line from that reference to the ablation mean (= the change). bars: 'ci' (95% CI of the mean) or 'spread'
    # (2.5-97.5 percentile across iterations); highlight: bars only for the best method per setup
    fig, ax = plt.subplots(figsize=(9, 5))
    setups = ['original'] + ABLATIONS
    off = np.linspace(-0.3, 0.3, len(METHODS))
    best = {s: best_method(s, k) for s in setups} if highlight else {}
    for i, m in enumerate(METHODS):
        c = METHOD_COLORS[m]
        mu0 = get('original', m, k)[0]
        ax.axhline(mu0, color=c, ls=(0, (4, 3)), lw=0.6, alpha=0.55, zorder=1)
        for j, s in enumerate(setups):
            x = j + off[i]
            mu, ci, _, _ = get(s, m, k)
            if s != 'original' and not np.isnan(mu):
                ax.vlines(x, mu0, mu, color=c, ls=':', lw=0.9, alpha=0.8, zorder=2)
            if bars == 'spread':
                lo, hi = get_spread(s, m, k)
                yerr = [[max(mu - lo, 0)], [max(hi - mu, 0)]]
            else:
                yerr = ci
            mec = 'black' if m in ('combat_seq', 'combat') else c
            if highlight and best[s] != m:
                ax.plot(x, mu, METHOD_MARKERS[m], color=c, mec=mec, mew=0.4, ms=3.5, zorder=3, label=m if j == 0 else None)
            elif highlight:
                ax.errorbar(x, mu, yerr=yerr, fmt=METHOD_MARKERS[m], color=c, mec='black', mew=0.8, ms=6, capsize=3,
                            elinewidth=1.6, capthick=1.2, zorder=4, label=m if j == 0 else None)
                ax.annotate(m, (x, mu + (yerr[1][0] if bars == 'spread' else yerr)), textcoords='offset points', xytext=(0, 4),
                            ha='center', fontsize=6.5, color=c if m not in ('combat_seq', 'combat') else 'darkgoldenrod')
            else:
                ax.errorbar(x, mu, yerr=yerr, fmt=METHOD_MARKERS[m], color=c, mec=mec, mew=0.4, ms=3.5, capsize=1,
                            elinewidth=0.7, capthick=0.7, zorder=3, label=m if j == 0 else None)
    for j in range(len(setups) - 1):
        ax.axvline(j + 0.5, color='0.85', lw=0.4, zorder=0)
    ax.axvspan(-0.5, 0.5, color='0.94', lw=0, zorder=0)
    if k in ('FDR', 'batch_shannon_rej'):
        ax.axhline(0.05, color='red', ls=':', lw=0.8)
    ax.set_xticks(range(len(setups)), setups, rotation=30, ha='right', fontsize=9)
    ax.tick_params(axis='y', labelsize=9)
    ax.set_xlim(-0.5, len(setups) - 0.5)
    ax.set_ylabel(label, fontsize=10)
    what = 'mean ± 95% CI' if bars == 'ci' else 'mean, bars = 2.5-97.5 percentile across iterations'
    if highlight:
        what = f'bars only for the best method per setup ({BEST_RULE[k]}): ' + ('95% CI of the mean' if bars == 'ci' else '2.5-97.5 percentile across iterations')
        ax.set_title(f'{label} vs original setting (dashed: original mean; dotted: change)\n{what}; {N} iterations', fontsize=9)
    else:
        ax.set_title(f'{label} vs original setting (dashed: original mean; dotted: change; {what}, {N} iterations)', fontsize=9)
    ax.spines[['top', 'right']].set_visible(False)
    for sp in ['left', 'bottom']:
        ax.spines[sp].set_linewidth(0.6)
    h, l = ax.get_legend_handles_labels()
    order = sorted(range(len(l)), key=lambda i: METHODS.index(l[i]))  # keep METHODS order when highlighted entries are errorbars
    ax.legend([h[i] for i in order], [l[i] for i in order], bbox_to_anchor=(1.01, 0.5), loc='center left', frameon=False,
              fontsize=8, markerscale=1.2)
    fig.tight_layout(); fig.savefig(f'{OUT}/{fname}_{k}.pdf', bbox_inches='tight'); plt.close(fig)


for k, label in METRICS.items():
    if all(np.isnan(get(s, m, k)[0]) for s in SETUPS for m in METHODS):  # e.g. Aitchison metrics on relative abundance data
        print(f'skipping {k}: no values')
        continue
    # 1. original style: x = setup, one marker per method (no dodge), means only
    fig, ax = plt.subplots(figsize=(7, 5))
    for m in METHODS:
        ax.scatter(list(SETUPS), [get(s, m, k)[0] for s in SETUPS], color=METHOD_COLORS[m], marker=METHOD_MARKERS[m], label=m, s=40)
    ax.set_title(label); ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(axis='x', rotation=90)
    if k in ('FDR', 'batch_shannon_rej'):
        ax.axhline(0.05, color='red', ls='--', lw=1)
    ax.legend(bbox_to_anchor=(1.02, 0.5), loc='center left', fontsize=8, frameon=False)
    fig.savefig(f'{OUT}/original_style_{k}.pdf', bbox_inches='tight'); plt.close(fig)

    # 2. reference line + CI: one panel per method; shared y axis
    fig, axes = plt.subplots(2, 5, figsize=(16, 6.5), sharey=True)
    for ax, m in zip(axes.flat, METHODS):
        mu0, ci0, _, _ = get('original', m, k)
        ax.axhspan(mu0 - ci0, mu0 + ci0, color='lightgray', alpha=0.6, lw=0)
        ax.axhline(mu0, color='black', ls='--', lw=1)
        for j, s in enumerate(ABLATIONS):
            mu, ci, _, _ = get(s, m, k)
            ax.errorbar(j, mu, yerr=ci, fmt='o', color=ABL_COLORS[s], capsize=3, ms=5)
        ax.set_xticks(range(len(ABLATIONS)), ABLATIONS, rotation=90, fontsize=8)
        ax.set_title(m, fontsize=10, color=METHOD_COLORS[m] if m not in ('combat_seq', 'combat') else 'darkgoldenrod')
        ax.spines[['top', 'right']].set_visible(False)
        if k in ('FDR', 'batch_shannon_rej'):
            ax.axhline(0.05, color='red', ls=':', lw=1)
    for ax in axes.flat[len(METHODS):]:
        ax.set_visible(False)
    fig.suptitle(f'{label}: ablations (points, 95% CI) vs original (dashed line, gray = 95% CI); {N} iterations', fontsize=11)
    fig.tight_layout(); fig.savefig(f'{OUT}/reference_ci_{k}.pdf', bbox_inches='tight'); plt.close(fig)

    # 2a. absolute by setup: x = setup (original first), one line per method (main-figure colors/markers), 95% CI
    fig, ax = plt.subplots(figsize=(9, 5))
    setups = ['original'] + ABLATIONS
    off = np.linspace(-0.2, 0.2, len(METHODS))
    for i, m in enumerate(METHODS):
        xs = np.arange(len(setups)) + off[i]
        mus = [get(s, m, k)[0] for s in setups]
        cis = [get(s, m, k)[1] for s in setups]
        ax.plot(xs, mus, color=METHOD_COLORS[m], lw=1, alpha=0.6)
        ax.errorbar(xs, mus, yerr=cis, fmt=METHOD_MARKERS[m], color=METHOD_COLORS[m],
                    mec='black' if m in ('combat_seq', 'combat') else METHOD_COLORS[m], mew=0.5, capsize=2, ms=5, label=m)
    ax.axvspan(-0.5, 0.5, color='lightgray', alpha=0.4, lw=0)
    if k in ('FDR', 'batch_shannon_rej'):
        ax.axhline(0.05, color='red', ls=':', lw=1)
    ax.set_xticks(range(len(setups)), setups, rotation=45, ha='right')
    ax.set_xlim(-0.5, len(setups) - 0.5)
    ax.set_ylabel(label)
    ax.set_title(f'{label} by setup (mean, 95% CI; original shaded); {N} iterations', fontsize=11)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(bbox_to_anchor=(1.02, 0.5), loc='center left', frameon=False, fontsize=9)
    fig.tight_layout(); fig.savefig(f'{OUT}/absolute_by_setup_{k}.pdf', bbox_inches='tight'); plt.close(fig)

    # 2c. absolute with references (main style), plus spread and highlight variants
    absolute_reference(k, label, 'ci', False, 'absolute_reference')
    absolute_reference(k, label, 'spread', False, 'absolute_reference_spread')
    absolute_reference(k, label, 'ci', True, 'highlight_ci')
    absolute_reference(k, label, 'spread', True, 'highlight_spread')

    # 2b. single panel: per method, original as a short flat line (+ 95% CI band), ablations dodged around it
    fig, ax = plt.subplots(figsize=(12, 4.8))
    w = 0.14
    half = (len(ABLATIONS) - 1) / 2 * w + 0.1
    for i, m in enumerate(METHODS):
        mu0, ci0, _, _ = get('original', m, k)
        ax.fill_between([i - half, i + half], mu0 - ci0, mu0 + ci0, color='lightgray', alpha=0.7, lw=0)
        ax.plot([i - half, i + half], [mu0, mu0], color='black', ls='--', lw=1.2, label='original (95% CI)' if i == 0 else None)
    for j, s in enumerate(ABLATIONS):
        for i, m in enumerate(METHODS):
            mu, ci, _, _ = get(s, m, k)
            ax.errorbar(i + (j - (len(ABLATIONS) - 1) / 2) * w, mu, yerr=ci, fmt=ABL_MARKERS[s], color=METHOD_COLORS[m],
                        mec='black' if m in ('combat_seq', 'combat') else METHOD_COLORS[m], mew=0.5, capsize=2, ms=5)
    for i in range(len(METHODS) - 1):
        ax.axvline(i + 0.5, color='lightgray', lw=0.5)
    if k in ('FDR', 'batch_shannon_rej'):
        ax.axhline(0.05, color='red', ls=':', lw=1)
    ax.set_xticks(range(len(METHODS)), METHODS, rotation=45, ha='right')
    ax.set_ylabel(label)
    ax.set_title(f'{label}: original (dashed line, gray = 95% CI) vs ablations (mean, 95% CI); {N} iterations', fontsize=11)
    ax.spines[['top', 'right']].set_visible(False)
    from matplotlib.lines import Line2D
    setup_legend(ax, [Line2D([], [], color='black', ls='--', label='original (95% CI)')])
    color_ticks(ax)
    fig.tight_layout(); fig.savefig(f'{OUT}/reference_ci_single_{k}.pdf', bbox_inches='tight'); plt.close(fig)

    # 3. difference from original: x = method, dodged by ablation
    fig, ax = plt.subplots(figsize=(11, 4.5))
    w = 0.14
    for j, s in enumerate(ABLATIONS):
        xs, ys, es = [], [], []
        for i, m in enumerate(METHODS):
            mu0, _, sd0, n0 = get('original', m, k)
            mu, _, sd, n = get(s, m, k)
            if n0 > 1 and n > 1:
                se = np.sqrt(sd0 ** 2 / n0 + sd ** 2 / n)
                df = se ** 4 / ((sd0 ** 2 / n0) ** 2 / (n0 - 1) + (sd ** 2 / n) ** 2 / (n - 1)) if se > 0 else n0 + n - 2
                ax.errorbar(i + (j - (len(ABLATIONS) - 1) / 2) * w, mu - mu0, yerr=stats.t.ppf(0.975, df) * se, fmt=ABL_MARKERS[s],
                            color=METHOD_COLORS[m], mec='black' if m in ('combat_seq', 'combat') else METHOD_COLORS[m], mew=0.5, capsize=2, ms=5)
    ax.axhline(0, color='black', lw=1)
    for i in range(len(METHODS) - 1):
        ax.axvline(i + 0.5, color='lightgray', lw=0.5)
    ax.set_xticks(range(len(METHODS)), METHODS, rotation=45, ha='right')
    ax.set_ylabel(f'{label}\n(ablation - original)')
    ax.set_title(f'Change from original setting, mean difference with 95% CI; {N} iterations', fontsize=11)
    ax.spines[['top', 'right']].set_visible(False)
    setup_legend(ax)
    color_ticks(ax)
    fig.tight_layout(); fig.savefig(f'{OUT}/difference_{k}.pdf', bbox_inches='tight'); plt.close(fig)

print(OUT)
print(d.groupby('setup')['iter'].nunique().to_string())
