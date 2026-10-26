"""Figure 3: ablations vs the baseline simulation, in the Figure 2 layout (fig2_style_plots.py).

Reads the summary written by ablation_vis.py (summary_mean_ci.csv: setup, method, metric, mean, ci) and plots the
baseline ("original": count, OR 1, condition 0.5, batch 0.5, no library-size confounding) next to each ablation.

Usage: python3 fig3_ablation_plots.py <summary_mean_ci.csv> <out_dir> [--setups n150,halftaxa,zeroinf90,batchfrac10,condfrac10]
       [--datatype count|relab]
"""
import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

import fig2_style_plots as F

XLABELS = {'original': 'Original', 'n150': 'n = 150', 'halftaxa': 'Half taxa', 'zeroinf90': '90% zeros',
           'batchfrac10': '10% batch', 'condfrac10': '10% condition', 'batchfrac20': '20% batch',
           'condfrac20': '20% condition', 'batchfrac40': '40% batch', 'condfrac40': '40% condition',
           'zeroinf75': '75% zeros', 'zeroinf80': '80% zeros', 'zeroinf85': '85% zeros'}
LABELLED = {'biovar_aitch_r2', 'biovar_bray_r2', 'shan_batch', 'shan_cond', 'rf_auc', 'rf_weighted_f1'}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('summary')
    ap.add_argument('out_dir')
    ap.add_argument('--setups', default='n150,halftaxa,zeroinf90,batchfrac10,condfrac10')
    ap.add_argument('--name', default='Figure3_ablations.pdf')
    ap.add_argument('--datatype', default='count', choices=['count', 'relab'])
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    setups = ['original'] + a.setups.split(',')
    summ = pd.read_csv(a.summary).rename(columns={'setup': 'setting'})
    summ['metric'] = summ['metric'].replace({'batch_shannon_rej': 'shan_batch', 'biovar_shannon_rej': 'shan_cond'})
    summ = summ[summ.setting.isin(setups)]
    methods = [m for m in F.METHODS[a.datatype] if m in set(summ.method)]

    fig = plt.figure(figsize=(F.W, F.H))
    axes = {}
    for key in F.AX:
        if summ[summ.metric == key].empty:
            continue
        ax = F.draw(fig, key, summ, methods, setups, False)
        ax.axvspan(-0.5, 0.5, color='0.92', zorder=-1, lw=0)  # baseline column
        if key in LABELLED:
            ax.set_xticks(range(len(setups)), [XLABELS.get(s, s) for s in setups], rotation=35, ha='right',
                          rotation_mode='anchor', fontsize=9.5)
        else:
            ax.set_xticks(range(len(setups)), [''] * len(setups))
        axes[key] = ax
    # each panel keeps its own y axis, since the ablations shift some metrics only slightly; the right panel of
    # each pair moves right to make room for its tick labels
    for _, q in F.PAIRS:
        if q in axes:
            b = axes[q].get_position()
            axes[q].set_position([b.x0 + 0.022, b.y0, b.width - 0.022, b.height])
    for k in ('FDR', 'shan_batch', 'shan_cond'):
        if k in axes:
            axes[k].set_ylim(0, 1.0)

    fig.text(73 / 1173.28, 1 - 10 / 714.88, 'a. Beta diversity', fontsize=14, va='top')
    fig.text(617 / 1173.28, 1 - 10 / 714.88, 'b. Differential abundance analysis', fontsize=14, va='top')
    fig.text(617 / 1173.28, 1 - 210 / 714.88, 'c. Alpha diversity', fontsize=14, va='top')
    fig.text(617 / 1173.28, 1 - 459 / 714.88, 'd. Machine learning', fontsize=14, va='top')
    handles = [Line2D([], [], color=F.COLORS[m], marker=F.MARKERS[m], ms=5, lw=1.2,
                      mec='black' if m in ('combat_seq', 'combat') else F.COLORS[m], label=F.LABELS[m]) for m in methods]
    leg = fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(598 / 1800, 1 - 745 / 1098), fontsize=11,
                     frameon=True, edgecolor='0.8', handlelength=2.2, borderpad=0.6, labelspacing=0.55)
    leg.get_frame().set_linewidth(0.8)
    items = [(None, 'Data type: ', 'count data' if a.datatype == 'count' else 'relative abundance data'),
             ('Odds ratio between batch and', 'conditions: ', '1'),
             ('Confounding between batch', 'effect and library sizes: ', 'No')]
    F.boxed_text(fig, F.frac(173, 520, 772, 980), 'Simulation Dataset Parameters', items)
    fig.savefig(os.path.join(a.out_dir, a.name))
    plt.close(fig)
    print(os.path.join(a.out_dir, a.name))


if __name__ == '__main__':
    main()
