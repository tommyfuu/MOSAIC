"""Supplemental Figure 5: runtime of the integration strategies (Illustrator-editable, Fig 2 colors and fonts).
  a-c. simulated data at OR 1, 1.25, 1.5: count / relative abundance x batch-library size confounding,
       median runtime over 1000 iterations (bars: interquartile range) per setting, from runtime_aggregate.py
       (median rather than mean: a few iterations ran on loaded nodes and inflate the mean)
  d.   real-world datasets (single run each), from the 'runtime' row of global_benchmarking_stats_<dataset>.csv
Runtimes are in seconds on a log scale; the uncorrected data (nobc) has no runtime. Text is embedded Arial (TrueType).
Usage: python3 supp_runtime_plot.py <runtime_agg_dir> <out_dir>
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, NullLocator

sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import REPO_ROOT

AGG, OUT = sys.argv[1:3]
os.makedirs(OUT, exist_ok=True)
FONT_DIR = '/athena/linglab/scratch/chf4012/tools/fonts_arial_dl/x'  # Arial from Microsoft core fonts (arial32.exe)
for f in ('Arial.TTF', 'Arialbd.TTF'):
    fm.fontManager.addfont(os.path.join(FONT_DIR, f))
matplotlib.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42, 'font.family': 'Arial', 'axes.unicode_minus': False,
                            'font.size': 10, 'axes.linewidth': 0.8, 'xtick.major.width': 0.8, 'ytick.major.width': 0.8})

COLORS = {'combat_seq': '#faec28', 'combat': '#faec28', 'limma': '#343c97', 'MMUPHin': '#00a24e',
          'ConQuR': '#f7af26', 'ConQuR_rel': '#f7af26', 'ConQuR_libsize': '#902e8f', 'harmony': '#00b9d9',
          'percentile_norm': '#949698', 'debiasm': '#ec2390'}
MARKERS = {'combat': 'v', 'combat_seq': 'D', 'limma': '^', 'MMUPHin': '<', 'ConQuR': '>', 'ConQuR_rel': 'p',
           'ConQuR_libsize': 's', 'harmony': '+', 'percentile_norm': '*', 'debiasm': 'd'}
LABELS = {'combat_seq': 'ComBat_seq', 'combat': 'ComBat', 'limma': 'limma', 'MMUPHin': 'MMUPHin', 'ConQuR': 'ConQuR',
          'ConQuR_rel': 'ConQuR_rel', 'ConQuR_libsize': 'ConQuR_libsize', 'harmony': 'Harmony',
          'percentile_norm': 'Percentile Normalization', 'debiasm': 'DEBIAS-M'}
METHODS = {'count': ['combat_seq', 'limma', 'MMUPHin', 'ConQuR', 'ConQuR_libsize', 'harmony', 'percentile_norm', 'debiasm'],
           'relab': ['combat', 'limma', 'MMUPHin', 'ConQuR_rel', 'harmony', 'percentile_norm', 'debiasm']}
GRID = ['0.25_0', '0.5_0', '0.75_0', '1_0', '0.25_0.25', '0.5_0.25', '0.75_0.25', '0.25_0.5', '0.5_0.5', '0.25_0.75']
RW = {'count': [('autism_2_microbiomeHD', 'autism (2 studies)', 'o', '#d38fc6'), ('cdi_3_microbiomeHD_clean', 'CDI (3 studies)', 's', 'black')],
      'relab': [('ibd_3_CMD', 'IBD (3 studies)', 'D', '#f7941d'), ('crc_11_CMD', 'CRC (11 studies)', 'o', '#3a53a4')]}
YTICKS = [0.01, 0.1, 1, 10, 100, 1000]


def style_y(ax, lo, hi, show):
    ax.set_yscale('log', base=2)
    ax.set_ylim(lo, hi)
    ax.yaxis.set_major_locator(FixedLocator([t for t in YTICKS if lo <= t <= hi]))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_yticklabels([f'{t:g}' for t in YTICKS if lo <= t <= hi] if show else [])


sim = {(dt, rel, OR): pd.read_csv(f'{AGG}/runtime_{dt}_{rel}_OR{OR}.csv', dtype={'setting': str})
       for dt in ('count', 'relab') for rel in ('yes', 'no') for OR in ('1', '1.25', '1.5')}
allv = pd.concat(sim.values())
lo, hi = 2 ** np.floor(np.log2(allv['q25'].min() * 0.8)), 2 ** np.ceil(np.log2(allv['q75'].max() * 1.25))

fig = plt.figure(figsize=(8.5, 16))
L, R, CW, GAPX = 0.14, 0.75, 0.29, 0.03  # left edge, legend start, column width, column gap
top, BH, RH, GAPR, GAPB = 0.975, 0.20, 0.078, 0.022, 0.03  # block height, row height, gaps
for bi, OR in enumerate(('1', '1.25', '1.5')):
    y_block = top - bi * (BH + GAPB)
    fig.text(0.02, y_block, f'{"abc"[bi]}. runtime on simulated data, odds ratio = {OR}', fontsize=13, va='top')
    for ri, rel in enumerate(('yes', 'no')):
        for ci, dt in enumerate(('count', 'relab')):
            y0 = y_block - 0.022 - (ri + 1) * RH - ri * GAPR
            ax = fig.add_axes([L + ci * (CW + GAPX), y0, CW, RH])
            s = sim[(dt, rel, OR)]
            for m in METHODS[dt]:
                r = s[s.method == m].set_index('setting').reindex(GRID)
                if r['median'].isna().all():
                    continue
                ax.errorbar(range(len(GRID)), r['median'], yerr=[r['median'] - r['q25'], r['q75'] - r['median']], color=COLORS[m], marker=MARKERS[m], ms=2.5,
                            mew=0.6, lw=1, capsize=1.2, elinewidth=0.7,
                            mec='black' if m in ('combat_seq', 'combat') else COLORS[m])
            style_y(ax, lo, hi, ci == 0)
            ax.set_xlim(-0.5, len(GRID) - 0.5)
            ax.spines[['top', 'right']].set_visible(False)
            for j in range(len(GRID)):
                ax.axvline(j, color='0.9', lw=0.5, zorder=0)
            ax.set_title(f'{"count" if dt == "count" else "relative abundance"}, '
                         f'{"library size confounded" if rel == "yes" else "no library size confounding"}',
                         fontsize=9.5, pad=3)
            last = bi == 2 and ri == 1
            ax.set_xticks(range(len(GRID)), [''] * len(GRID))
            if last:
                for j, sset in enumerate(GRID):
                    c, b = sset.split('_')
                    ax.text(j, -0.06, f'{float(c):.2f}', transform=ax.get_xaxis_transform(), rotation=90, ha='center', va='top', fontsize=9)
                    ax.text(j, -0.42, f'{float(b):.2f}', transform=ax.get_xaxis_transform(), rotation=90, ha='center', va='top', fontsize=9)
                if ci == 0:
                    ax.text(-0.03, -0.24, 'condition effect', transform=ax.transAxes, ha='right', va='center', fontsize=9.5)
                    ax.text(-0.03, -0.60, 'batch effect', transform=ax.transAxes, ha='right', va='center', fontsize=9.5)
            if ci == 0 and ri == 0:
                ax.set_ylabel('runtime (s)', fontsize=9.5)
                ax.yaxis.set_label_coords(-0.16, -0.14)

# d. real-world datasets
y_d = top - 3 * (BH + GAPB) - 0.035
fig.text(0.02, y_d, 'd. runtime on real-world datasets', fontsize=13, va='top')
for ci, dt in enumerate(('count', 'relab')):
    ax = fig.add_axes([L + ci * (CW + GAPX), y_d - 0.03 - 0.115, CW, 0.115])
    vals = []
    for ds, lab, mk, col in RW[dt]:
        g = pd.read_csv(f'{REPO_ROOT}/trial_rw/outputs/{ds}/global_benchmarking_stats_{ds}.csv', index_col=0)
        rt = pd.to_numeric(g.loc['runtime'], errors='coerce').reindex(METHODS[dt])
        vals += list(rt.dropna())
        ax.plot(range(len(METHODS[dt])), rt.values, marker=mk, color=col, lw=1.2, ms=4, label=lab)
    style_y(ax, 2 ** np.floor(np.log2(min(vals) * 0.8)), 2 ** np.ceil(np.log2(max(vals) * 1.25)), True)
    ax.set_xticks(range(len(METHODS[dt])), [LABELS[m] for m in METHODS[dt]], rotation=90, fontsize=9)
    ax.set_xlim(-0.5, len(METHODS[dt]) - 0.5)
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_title(f'{"count" if dt == "count" else "relative abundance"} datasets', fontsize=9.5, pad=3)
    if ci == 0:
        ax.set_ylabel('runtime (s)', fontsize=9.5)
    if ci == 1:
        rw_handles = [Line2D([], [], color=col, marker=mk, lw=1.2, ms=4, label=lab) for d_ in ('count', 'relab') for _, lab, mk, col in RW[d_]]
        fig.legend(handles=rw_handles, loc='upper left', bbox_to_anchor=(L + 2 * CW + GAPX + 0.015, y_d - 0.03), fontsize=9,
                   frameon=True, edgecolor='0.8').get_frame().set_linewidth(0.8)

# shared legend for a-c
ms = ['combat_seq', 'combat', 'limma', 'MMUPHin', 'ConQuR', 'ConQuR_rel', 'ConQuR_libsize', 'harmony', 'percentile_norm', 'debiasm']
handles = [Line2D([], [], color=COLORS[m], marker=MARKERS[m], ms=5, lw=1.2,
                  mec='black' if m in ('combat_seq', 'combat') else COLORS[m], label=LABELS[m]) for m in ms]
leg = fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(L + 2 * CW + GAPX + 0.015, top - 0.03), fontsize=9,
                 frameon=True, edgecolor='0.8', handlelength=2.0, borderpad=0.6, labelspacing=0.5)
leg.get_frame().set_linewidth(0.8)
fig.savefig(f'{OUT}/SuppFig5_runtime.pdf')
print(f'{OUT}/SuppFig5_runtime.pdf')
