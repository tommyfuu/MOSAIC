#!/usr/bin/env python3
"""Summarize the condsweep ablation (OR=1.25, cond_effect in {0.005,0.01,0.02,0.05,0.1},
batch_effect=0.5) across metrics and methods, and compare the spread across methods to
the original {0,0.25,0.5,0.75,1} grid.

Outputs: condsweep_ablation_report.txt, condsweep_ablation_results.csv,
condsweep_ablation_spread_<metric>.png
"""
import os
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import REPO_ROOT, SIM_OUTPUT_ROOT, BENCH_DIR, DATA_DIR

SIM_BASE = SIM_OUTPUT_ROOT
SCEN_BASE = f"{SIM_BASE}/scenarios_simulation_data_052626"

COND_SWEEP = [0.005, 0.01, 0.02, 0.05, 0.1]
OR = 1.25
BATCH = 0.5
N_ITERS = 1000

ALL_COMBOS = [
    ("count", "no_relation"),
    ("count", "yes_relation"),
    ("relab", "no_relation"),
    ("relab", "yes_relation"),
]
COMBO_NAMES = {
    "count_no": ("count", "no_relation"),
    "count_yes": ("count", "yes_relation"),
    "relab_no": ("relab", "no_relation"),
    "relab_yes": ("relab", "yes_relation"),
}

# Full metric set: RF/ML, diff-abundance, PERMANOVA R^2, alpha-diversity.
METRICS = [
    "rf_auc", "rf_weighted_f1",
    "FDR", "power",
    "batch_aitch_r2", "biovar_aitch_r2",
    "batch_bray_r2", "biovar_bray_r2",
    "batch_shannon_pval", "biovar_shannon_pval",
]
PRIMARY_METRIC = "rf_auc"  # used for the original-grid comparison/verdict

# original main-grid combo (cond=0.5, batch=0.5, OR=1) for the count/no-relation
# comparison baseline.
ORIG_EVAL_DIR = f"{SIM_BASE}/082924check/simulation_data_eval_new_count_norelation_102023"
ORIG_COMBO = "out_1_0.5_0.5"


def load_combo_stats(eval_dir_path, datatype, relation):
    """Read global_benchmarking_stats_batch_0.csv across all cond_effect
    values and iterations for one (datatype, relation) combo. Returns one
    row per (cond_effect, method) with a mean column per metric."""
    rows = []
    for cond in COND_SWEEP:
        combo_name = f"out_{OR}_{cond}_{BATCH}"
        vals_by_method_metric = {}  # method -> metric -> [values]
        n_found = 0
        for it in range(1, N_ITERS + 1):
            p = os.path.join(eval_dir_path, f"{combo_name}_iter_{it}",
                              "global_benchmarking_stats_batch_0.csv")
            if not os.path.exists(p):
                continue
            n_found += 1
            df = pd.read_csv(p, index_col=0)
            df.columns = [c.strip() for c in df.columns]
            for method in df.columns:
                for metric in METRICS:
                    if metric not in df.index:
                        continue
                    v = df.loc[metric, method]
                    try:
                        v = float(v)
                    except (ValueError, TypeError):
                        continue
                    vals_by_method_metric.setdefault(method, {}).setdefault(metric, []).append(v)

        for method, metric_vals in vals_by_method_metric.items():
            row = {
                "datatype": datatype, "relation": relation,
                "cond_effect": cond, "method": method,
                "n_iters_found": n_found,
            }
            has_any = False
            for metric in METRICS:
                vals = metric_vals.get(metric, [])
                if vals:
                    row[f"{metric}_mean"] = np.mean(vals)
                    row[f"{metric}_n"] = len(vals)
                    has_any = True
                else:
                    row[f"{metric}_mean"] = np.nan
                    row[f"{metric}_n"] = 0
            if has_any:
                rows.append(row)
    return pd.DataFrame(rows)


def load_original_grid_spread(metric):
    """Spread across methods for the original main grid's cond=0.5,batch=0.5,
    OR=1 combo (count/no-relation), all 1000 iterations, for one metric."""
    vals_by_method = {}
    n_found = 0
    for it in range(1, 1001):
        p = os.path.join(ORIG_EVAL_DIR, f"{ORIG_COMBO}_iter_{it}",
                          "global_benchmarking_stats_batch_0.csv")
        if not os.path.exists(p):
            continue
        n_found += 1
        df = pd.read_csv(p, index_col=0)
        df.columns = [c.strip() for c in df.columns]
        if metric not in df.index:
            continue
        for method in df.columns:
            v = df.loc[metric, method]
            try:
                v = float(v)
            except (ValueError, TypeError):
                continue
            vals_by_method.setdefault(method, []).append(v)
    means = {m: np.mean(v) for m, v in vals_by_method.items() if v}
    return means, n_found


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--combo", choices=list(COMBO_NAMES.keys()), default=None,
                         help="run a single combo only (e.g. count_no) instead of all 4, so "
                              "each variation's report/plots can be made as soon as its own "
                              "evaluate data is ready")
    args = parser.parse_args()

    if args.combo:
        COMBOS = [COMBO_NAMES[args.combo]]
        suffix = f"_{args.combo}"
    else:
        COMBOS = ALL_COMBOS
        suffix = ""

    all_rows = []
    report_lines = []
    report_lines.append("CONDSWEEP ABLATION -- FULL METRIC ANALYSIS")
    report_lines.append("=" * 60)
    report_lines.append(f"Grid: OR={OR}, cond_effect in {COND_SWEEP}, batch={BATCH}")
    report_lines.append(f"Iterations: 1-{N_ITERS}")
    report_lines.append(f"Metrics: {', '.join(METRICS)}")
    report_lines.append("")

    for datatype, relation in COMBOS:
        eval_dir = f"{SCEN_BASE}/eval_{datatype}_condsweep_{relation}"
        df = load_combo_stats(eval_dir, datatype, relation)
        if df.empty:
            report_lines.append(f"[{datatype}/{relation}] NO DATA FOUND at {eval_dir}")
            continue
        all_rows.append(df)

        report_lines.append(f"--- {datatype} / {relation} ---")
        for cond in COND_SWEEP:
            sub = df[df["cond_effect"] == cond]
            if sub.empty:
                report_lines.append(f"  cond_effect={cond}: no data")
                continue
            n = sub[f"{PRIMARY_METRIC}_n"].max()
            report_lines.append(f"  cond_effect={cond:<5} n_iters~{n:<3}")
            for metric in METRICS:
                col = f"{metric}_mean"
                vals = sub[col].dropna()
                if vals.empty:
                    continue
                spread = vals.max() - vals.min()
                report_lines.append(
                    f"      {metric:<20} range=[{vals.min():.4f}, {vals.max():.4f}] "
                    f"spread={spread:.4f} std_across_methods={vals.std():.4f}"
                )
        report_lines.append("")

    if all_rows:
        full_df = pd.concat(all_rows, ignore_index=True)
        full_df.to_csv(os.path.join(BENCH_DIR, f"condsweep_ablation_results{suffix}.csv"), index=False)
    else:
        full_df = pd.DataFrame()
        report_lines.append("WARNING: no condsweep data found at all -- pipeline may not have finished yet.")

    # compare RF AUC spread to the original {0,0.25,0.5,0.75,1} grid
    orig_means, orig_n = load_original_grid_spread(PRIMARY_METRIC)
    if orig_means:
        orig_spread = max(orig_means.values()) - min(orig_means.values())
        orig_std = np.std(list(orig_means.values()))
        report_lines.append(f"--- ORIGINAL GRID (count/no-relation, OR=1, cond=0.5, batch=0.5) [{PRIMARY_METRIC}] ---")
        report_lines.append(f"  n_iters={orig_n} range=[{min(orig_means.values()):.4f}, {max(orig_means.values()):.4f}] "
                             f"spread={orig_spread:.4f} std_across_methods={orig_std:.4f}")
        report_lines.append("")

        if not full_df.empty:
            report_lines.append(f"--- VERDICT ({PRIMARY_METRIC} spread vs original grid) ---")
            for datatype, relation in COMBOS:
                sub = full_df[(full_df["datatype"] == datatype) & (full_df["relation"] == relation)]
                if sub.empty:
                    continue
                best_cond = None
                best_spread = -1
                col = f"{PRIMARY_METRIC}_mean"
                for cond in COND_SWEEP:
                    s = sub[sub["cond_effect"] == cond][col].dropna()
                    if s.empty:
                        continue
                    sp = s.max() - s.min()
                    if sp > best_spread:
                        best_spread = sp
                        best_cond = cond
                if best_cond is not None:
                    verdict = "BETTER differentiation" if best_spread > orig_spread else "NOT better"
                    report_lines.append(
                        f"  {datatype}/{relation}: best spread at cond_effect={best_cond} "
                        f"(spread={best_spread:.4f} vs original {orig_spread:.4f}) -> {verdict}"
                    )
            report_lines.append("")

    report_text = "\n".join(report_lines)
    print(report_text)
    with open(os.path.join(BENCH_DIR, f"condsweep_ablation_report{suffix}.txt"), "w") as f:
        f.write(report_text + "\n")

    # one PNG per metric: 2x2 panels for all combos, 1 panel with --combo
    if not full_df.empty:
        nrows, ncols = (1, 1) if args.combo else (2, 2)
        for metric in METRICS:
            col = f"{metric}_mean"
            if col not in full_df.columns or full_df[col].dropna().empty:
                continue
            fig, axes = plt.subplots(nrows, ncols, figsize=(7, 5) if args.combo else (14, 10), sharex=True)
            axes_flat = [axes] if args.combo else axes.flat
            for ax, (datatype, relation) in zip(axes_flat, COMBOS):
                sub = full_df[(full_df["datatype"] == datatype) & (full_df["relation"] == relation)]
                if sub.empty:
                    ax.set_title(f"{datatype}/{relation} (no data)")
                    continue
                for method in sorted(sub["method"].unique()):
                    m = sub[sub["method"] == method].sort_values("cond_effect")
                    ax.plot(m["cond_effect"], m[col], marker="o", label=method, linewidth=1.5)
                ax.set_title(f"{datatype} / {relation}")
                ax.set_xlabel("cond_effect")
                ax.set_ylabel(metric)
                ax.set_xscale("log")
            handles, labels = axes_flat[0].get_legend_handles_labels()
            fig.legend(handles, labels, loc="center right", bbox_to_anchor=(1.15, 0.5))
            fig.suptitle(f"Condsweep ablation: {metric} vs cond_effect (OR={OR}, batch={BATCH}, n_iters={N_ITERS})")
            fig.tight_layout()
            out_path = os.path.join(BENCH_DIR, f"condsweep_ablation_spread_{metric}{suffix}.png")
            fig.savefig(out_path, dpi=150, bbox_inches="tight")
            plt.close(fig)
            print(f"Plot saved to {out_path}")


if __name__ == "__main__":
    main()
