"""Fig 2-style composite simulation figures (Illustrator-editable), one per data type x odds ratio x
batch-library size confounding, from the aggregates written by fig2_style_aggregate.py.

Layout, colors, font sizes and page size follow Fig2_v10 (Illustrator, 1173 x 715 pt, Arial):
  a. Beta diversity      batch / condition PERMANOVA R2 (Aitchison, Bray-Curtis)
  b. Differential abundance analysis   FDR, sensitivity
  c. Alpha diversity     Shannon ~ batch / condition rejection proportion (p < 0.05)
  d. Machine learning    random forest AUROC / weighted F1 (condition-effect sweep at OR 1.25; main grid otherwise)
Methods are the succinct set (no bat/ruv3/svd/rbe/bmc); points are means over 1000 iterations with 95% CIs,
spread slightly within each setting; RF metrics are leakage-free (nocovar). Text is embedded as TrueType
(Arial) so it stays editable in Illustrator.

Usage: python3 fig2_style_plots.py <agg_dir> <out_dir> [--main]
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
from matplotlib.patches import FancyBboxPatch


FONT_DIR = '/athena/linglab/scratch/chf4012/tools/fonts_arial_dl/x'  # Arial from Microsoft core fonts (arial32.exe)
for f in ('Arial.TTF', 'Arialbd.TTF'):
    fm.fontManager.addfont(os.path.join(FONT_DIR, f))
matplotlib.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42, 'font.family': 'Arial', 'axes.unicode_minus': False,
                            'font.size': 10, 'axes.linewidth': 0.8, 'xtick.major.width': 0.8, 'ytick.major.width': 0.8})

W, H = 1173.28 / 72, 714.88 / 72  # Fig2_v10 page size in inches
COLORS = {'nobc': '#231f20', 'combat_seq': '#faec28', 'combat': '#faec28', 'limma': '#343c97', 'MMUPHin': '#00a24e',
          'ConQuR': '#f7af26', 'ConQuR_rel': '#f7af26', 'ConQuR_libsize': '#902e8f', 'harmony': '#00b9d9',
          'percentile_norm': '#949698', 'debiasm': '#ec2390'}
MARKERS = {'nobc': 'o', 'combat': 'v', 'combat_seq': 'D', 'limma': '^', 'MMUPHin': '<', 'ConQuR': '>', 'ConQuR_rel': 'p',
           'ConQuR_libsize': 's', 'harmony': '+', 'percentile_norm': '*', 'debiasm': 'd'}
LABELS = {'nobc': 'original', 'combat_seq': 'ComBat_seq', 'combat': 'ComBat', 'limma': 'limma', 'MMUPHin': 'MMUPHin',
          'ConQuR': 'ConQuR', 'ConQuR_rel': 'ConQuR_rel', 'ConQuR_libsize': 'ConQuR_libsize', 'harmony': 'Harmony',
          'percentile_norm': 'Percentile Normalization', 'debiasm': 'DEBIAS-M'}
METHODS = {'count': ['nobc', 'combat_seq', 'limma', 'MMUPHin', 'ConQuR', 'ConQuR_libsize', 'harmony', 'percentile_norm', 'debiasm'],
           'relab': ['nobc', 'combat', 'limma', 'MMUPHin', 'ConQuR_rel', 'harmony', 'percentile_norm', 'debiasm']}
GRID = ['0.25_0', '0.5_0', '0.75_0', '1_0', '0.25_0.25', '0.5_0.25', '0.75_0.25', '0.25_0.5', '0.5_0.5', '0.25_0.75']
SWEEP = ['0.005', '0.01', '0.02', '0.05', '0.1']


def frac(x0, x1, y0, y1):
    """Axes rectangle from Fig2_v10 pixel coordinates (1800 x 1098 render)."""
    return [x0 / 1800, 1 - y1 / 1098, (x1 - x0) / 1800, (y1 - y0) / 1098]


AX = {'batch_aitch_r2': frac(127, 515, 70, 312), 'batch_bray_r2': frac(537, 910, 70, 312),
      'biovar_aitch_r2': frac(127, 515, 375, 618), 'biovar_bray_r2': frac(537, 910, 375, 618),
      'FDR': frac(970, 1362, 68, 312), 'power': frac(1368, 1757, 68, 312),
      'shan_batch': frac(970, 1362, 372, 620), 'shan_cond': frac(1368, 1757, 372, 620),
      'rf_auc': frac(970, 1362, 750, 997), 'rf_weighted_f1': frac(1368, 1757, 750, 997)}
TITLES = {'batch_aitch_r2': 'batch PERMANOVA R² (Aitchison)', 'batch_bray_r2': 'batch PERMANOVA R² (Bray-Curtis)',
          'biovar_aitch_r2': 'condition PERMANOVA R² (Aitchison)', 'biovar_bray_r2': 'condition PERMANOVA R² (Bray-Curtis)',
          'FDR': 'False discovery rate', 'power': 'Sensitivity',
          'shan_batch': 'rejection proportion (Shannon ~ batches, p-value<0.05)',
          'shan_cond': 'rejection proportion (Shannon ~ conditions, p-value<0.05)',
          'rf_auc': 'random forest AUROC scores', 'rf_weighted_f1': 'random forest weighted F1 scores'}
PAIRS = [('batch_aitch_r2', 'batch_bray_r2'), ('biovar_aitch_r2', 'biovar_bray_r2'), ('FDR', 'power'),
         ('shan_batch', 'shan_cond'), ('rf_auc', 'rf_weighted_f1')]
LEFT = {a for a, _ in PAIRS}


def setting_labels(fig, ax, settings, sweep, row_names):
    """Two rows of rotated tick values below the axis: condition effect (upper) and batch effect (lower).
    row_names: None, or 'figleft' / 'axleft' for where the row names go."""
    ax.set_xticks(range(len(settings)), [''] * len(settings))
    y_c, y_b = (-0.05, -0.25) if sweep else (-0.05, -0.20)
    for i, s in enumerate(settings):
        c, b = (s, '0.5') if sweep else s.split('_')
        cl = c if sweep else f'{float(c):.2f}'
        bl = b if sweep else f'{float(b):.2f}'
        ax.text(i, y_c, cl, transform=ax.get_xaxis_transform(), rotation=90, ha='center', va='top', fontsize=10)
        ax.text(i, y_b, bl, transform=ax.get_xaxis_transform(), rotation=90, ha='center', va='top', fontsize=10)
    if row_names:
        yc, yb = (y_c - 0.075, y_b - 0.05) if sweep else (y_c - 0.06, y_b - 0.06)
        if row_names == 'figleft':
            for lab, yy in (('condition effect', yc), ('batch effect', yb)):
                ya = ax.transAxes.transform((0, yy))[1]
                fig.text(10 / 1800 + (0.02 if lab == 'batch effect' else 0), fig.transFigure.inverted().transform((0, ya))[1],
                         lab, ha='left', va='center', fontsize=10.5)
        else:
            ax.text(-0.02, yc, 'condition effect', transform=ax.transAxes, ha='right', va='center', fontsize=10.5)
            ax.text(-0.02, yb, 'batch effect', transform=ax.transAxes, ha='right', va='center', fontsize=10.5)


def draw(fig, key, summ, methods, settings, sweep):
    ax = fig.add_axes(AX[key])
    n = len(methods)
    off = np.linspace(-0.22, 0.22, n)
    sub = summ[summ.metric == key]
    for i, m in enumerate(methods):
        r = sub[sub.method == m].set_index('setting').reindex(settings)
        if r['mean'].isna().all():
            continue
        x = np.arange(len(settings)) + off[i]
        ax.errorbar(x, r['mean'], yerr=r['ci'], color=COLORS[m], marker=MARKERS[m], ms=2.5, mew=0.6, lw=1,
                    capsize=1.5, elinewidth=0.8, label=LABELS[m],
                    mec='black' if m in ('combat_seq', 'combat') else COLORS[m])
    for j in range(len(settings)):
        ax.axvline(j, color='0.88', lw=0.6, zorder=0)
    ax.set_xlim(-0.6, len(settings) - 0.4)
    if key in ('FDR', 'shan_batch'):
        ax.axhline(0.05, color='#ed1c24', ls=(0, (4, 2)), lw=1.4)
    ax.set_title(TITLES[key], fontsize=10 if key.startswith('shan') else 10.5, pad=4)
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelsize=10)
    return ax


def boxed_text(fig, rect, title, items):
    """items: list of (first_line or None, second_line_label, bold_value); laid out like Fig2_v10."""
    x0, y0, w, h = rect
    fig.patches.append(FancyBboxPatch((x0, y0), w, h, boxstyle='round,pad=0,rounding_size=0.004', transform=fig.transFigure,
                                      fc='white', ec='0.8', lw=0.8))
    fig.text(x0 + 0.012, y0 + h - 0.022, title, fontsize=12, va='top')
    y = y0 + h - 0.062
    lh = 0.024
    for first, lab, val in items:
        if first:
            fig.text(x0 + 0.018, y, first, fontsize=11, va='baseline')
            y -= lh
        ind = 0.018 + (0.016 if first else 0)
        t = fig.text(x0 + ind, y, lab, fontsize=11, va='baseline')
        fig.canvas.draw()
        bb = t.get_window_extent().transformed(fig.transFigure.inverted())
        fig.text(bb.x1 + 0.002, y, val, fontsize=11, fontweight='bold', va='baseline')
        y -= lh + 0.006


def make(dt, rel, OR, ml_sweep=False, name=None):
    """ml_sweep=False: panel d = RF on the same settings as a-c (supplementary figures);
    ml_sweep=True: panel d = condition-effect sweep at this OR/relation (main Figure 2)."""
    grid = pd.read_csv(f'{AGG}/grid_{dt}_{rel}_OR{OR}.csv', dtype={'setting': str})
    cs_file = f'{AGG}/condsweep_{dt}_{rel}_OR{OR}.csv'
    sweep = ml_sweep and os.path.exists(cs_file)
    if ml_sweep and not sweep:
        print(f'missing {cs_file}; skipping'); return
    cs = pd.read_csv(cs_file, dtype={'setting': str}) if sweep else None
    methods = [m for m in METHODS[dt] if m in set(grid.method)]
    fig = plt.figure(figsize=(W, H))
    axes = {}
    for key in AX:
        ml = key in ('rf_auc', 'rf_weighted_f1')
        src, settings = (cs, SWEEP) if (ml and sweep) else (grid, GRID)
        if src[src.metric == key].empty:  # Aitchison is not computed for relative abundance data
            ax = fig.add_axes(AX[key]); ax.axis('off')
            ax.text(0.5, 0.5, f'{TITLES[key]}\nnot computed for relative abundance data', ha='center', va='center',
                    fontsize=10, color='0.4', transform=ax.transAxes)
            continue
        axes[key] = draw(fig, key, src, methods, settings, ml and sweep)
        row_names = {'biovar_aitch_r2': 'figleft', 'rf_auc': 'axleft'}.get(key)
        if key == 'biovar_bray_r2' and 'biovar_aitch_r2' not in axes:  # relab: no Aitchison panels
            row_names = 'axleft'
        if key in ('biovar_aitch_r2', 'biovar_bray_r2', 'shan_batch', 'shan_cond', 'rf_auc', 'rf_weighted_f1'):
            setting_labels(fig, axes[key], settings, ml and sweep, row_names)
        else:
            axes[key].set_xticks(range(len(settings)), [''] * len(settings))
    for a, b in PAIRS:  # shared y within each pair, right panel without tick labels
        if a in axes and b in axes:
            lo = min(axes[a].get_ylim()[0], axes[b].get_ylim()[0]); hi = max(axes[a].get_ylim()[1], axes[b].get_ylim()[1])
            axes[a].set_ylim(lo, hi); axes[b].set_ylim(lo, hi); axes[b].tick_params(labelleft=False)
    for k in ('FDR', 'power', 'shan_batch', 'shan_cond'):
        if k in axes:
            axes[k].set_ylim(0, 1.0 if k != 'power' else max(0.65, axes[k].get_ylim()[1]))
    fig.text(73 / 1173.28, 1 - 10 / 714.88, 'a. Beta diversity', fontsize=14, va='top')
    fig.text(617 / 1173.28, 1 - 10 / 714.88, 'b. Differential abundance analysis', fontsize=14, va='top')
    fig.text(617 / 1173.28, 1 - 210 / 714.88, 'c. Alpha diversity', fontsize=14, va='top')
    fig.text(617 / 1173.28, 1 - 459 / 714.88, 'd. Machine learning', fontsize=14, va='top')
    handles = [Line2D([], [], color=COLORS[m], marker=MARKERS[m], ms=5, lw=1.2,
                      mec='black' if m in ('combat_seq', 'combat') else COLORS[m], label=LABELS[m]) for m in methods]
    leg = fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(598 / 1800, 1 - 745 / 1098), fontsize=11,
                     frameon=True, edgecolor='0.8', handlelength=2.2, borderpad=0.6, labelspacing=0.55)
    leg.get_frame().set_linewidth(0.8)
    items = [(None, 'Data type: ', 'count data' if dt == 'count' else 'relative abundance data'),
             ('Odds ratio between batch and', 'conditions: ', OR),
             ('Confounding between batch', 'effect and library sizes: ', 'Yes' if rel == 'yes' else 'No')]
    boxed_text(fig, frac(173, 520, 772, 980), 'Simulation Dataset Parameters', items)
    name = name or f'Fig2style_{dt}_OR{OR}_{"libsize_confounded" if rel == "yes" else "no_libsize_confounding"}.pdf'
    fig.savefig(f'{OUT}/{name}')
    plt.close(fig)
    print(name)


if __name__ == '__main__':
    AGG, OUT = sys.argv[1:3]
    os.makedirs(OUT, exist_ok=True)
    MAIN = '--main' in sys.argv  # main Figure 2 only: count, OR 1, no confounding, panel d from the OR 1 condition sweep
    if MAIN:
        make('count', 'no', '1', ml_sweep=True, name='Figure2_count_OR1_no_libsize_confounding_condsweepML.pdf')
    for dt in ([] if MAIN else ('count', 'relab')):
        for rel in ('no', 'yes'):
            for OR in ('1', '1.25', '1.5'):
                if os.path.exists(f'{AGG}/grid_{dt}_{rel}_OR{OR}.csv'):
                    make(dt, rel, OR)
