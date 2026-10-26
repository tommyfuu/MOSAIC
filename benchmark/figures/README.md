# Figure and table scripts

These scripts turn the per-replicate evaluation results into the figures and tables of the manuscript. Run them
from the repository root inside the `mic_data_integration` environment. The aggregation scripts read the
evaluation folders under `sim_output_root` (`paths.yml`) used for the manuscript runs; `ablation_vis.py --root`
reads the folders written by `Snakefile_simulation`.

| Output | Script(s) |
|---|---|
| Figure 2 | `fig2_style_aggregate.py count no 1 <agg>`, then `fig2_style_plots.py <agg> <out> --main` |
| Supplemental Figures 1-4 | `fig2_style_aggregate.py <count\|relab> <no\|yes> <1\|1.25\|1.5> <agg>` for each setting, `fig2_style_plots.py <agg> <pages>`, then `supp_figs_1to4.py <pages> <out>` (needs PyMuPDF) |
| Figure 3 | `ablation_vis.py -i 1000 --setups n150,halftaxa,zeroinf90,batchfrac10,condfrac10 --out <dir>` (writes `summary_mean_ci.csv`), then `fig3_ablation_plots.py <dir>/summary_mean_ci.csv <out>` |
| Supplemental Figure 5 | `runtime_aggregate.py <count\|relab> <no\|yes> <1\|1.25\|1.5> <dir>` for each setting, then `supp_runtime_plot.py <dir> <out>` |
| Ablation summary table | `ablation_summary_table.py` |

Figures 4 and 5 (real-world data) are made by `benchmark/evaluate_alt_rw_plots.py --succinct` and the PCoA plots
written by `benchmark/evaluate.py -o 4`. A worked example from simulated data to each figure is in
`docs/vignette.html`.
