# MOSAIC <img src="https://github.com/tommyfuu/MOSAIC/blob/main/logo/MOSAIC_logo.jpeg" width="25" height="25"> : a pipeline for MicrobiOme Studies Analytical Integration and Correction
Repository for the manuscript entitled __MOSAIC: MicrobiOme Studies Analytical Integration and Correction pipeline__ by _Fu et al._. 

<p align="center">
  <img src="https://github.com/tommyfuu/MOSAIC/blob/main/logo/MOSAIC_logo.jpeg" width="500" height="500">
</p>



(Logo generated with [Google Gemini](https://gemini.google.com/app))



_Developed and maintained by_:  [Chenlian (Tom) Fu](fcl200089@outlook.com)\
_Advised by_: Dr. Wodan Ling\
_Affiliations_: Weill Cornell Medicine, Memorial Sloan Kettering Cancer Center


This repository holds the code for the __MOSAIC__ data integration pipeline described in the manuscript. __MOSAIC__ is a three-part [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline. It preprocesses microbiome data from multiple sources (or simulates data), applies the data integration methods chosen by the user, and reports four types of evaluation metrics. The repository also contains code to reproduce the manuscript figures.

## 1. Environment set-up

The pipeline uses Python 3.8(.17) and R 4.3 in a conda environment. First install the conda packages from the environment file:

```
conda env create -f env.yml
conda activate mic_data_integration
```

Then open `R` and install MIDASim, FDboost and mixOmics individually, since they are not available through conda. Installation instructions:

- [MIDASim](https://github.com/mengyu-he/MIDASim) 
- [FDboost](https://github.com/boost-R/FDboost)
- [mixOmics](https://bioconductor.org/packages/release/bioc/html/mixOmics.html)

**Optional: DEBIAS-M** (`debiasm`) is supported as a Python batch correction method. It is not installed by default because it requires PyTorch (~2.4 GB of dependencies). To use it, install it after activating the environment:

```bash
pip install DEBIAS-M
```

Then add `"debiasm"` to `used_Python_methods` in your `config.yml`. If it is not listed there, it is never imported.


**Paths.** The only machine-specific setting is `sim_output_root` in `paths.yml` (where
simulation data and results are written); the repository location is detected
automatically. Python scripts read it through `paths.py` and R scripts through `paths.R`;
setting the environment variable `SIM_OUTPUT_ROOT` overrides the file. SLURM scripts write
their logs to `logs/slurm/` and assume they are submitted from `benchmark/` (or `pubtator/`
for the PubTator job). Change `--mail-user` in the `#SBATCH` headers to your own address.
The simulation config `config_simulate_052626.yml` also has a `sim_output_root` entry (used
by `evaluate.py -o 1/5` and the Snakemake interface); keep it the same as `paths.yml`.

## 2. Workflow used for the manuscript (SLURM scripts in `benchmark/`)

The simulation study, the ablations and the real-world analysis were run with the SLURM
scripts in `benchmark/` (run them from inside `benchmark/`). The Snakemake interface in
section 2a runs a single real-world dataset end to end, and `Snakefile_simulation` (section 3)
runs any simulation design, including all the ablations. `benchmark/evaluate.py`
is the Python entry point; `-o` selects the step:

| `-o` | Step |
|---|---|
| 1 / 2 | main simulation grid: Python methods / evaluation |
| 3 / 4 | real-world: Python methods / evaluation (`--config <yml>`) |
| 5 | main simulation grid: line plots |
| 6 / 7 | ablation scenarios incl. condsweep: Python methods / evaluation |
| 8 / 9 | ablation plots / condsweep plots |
| 13 / 14 / 15 | zero-inflation and named ablation scenarios: Python methods / evaluation / plots |
| 16 / 17 | any simulation scenario, given its folders (`--data-dir`, `--output-dir`, `--eval-dir`) and grid (`--or-l`, `--cond-l`, `--batch-l`): Python methods / evaluation; used by `Snakefile_simulation` (section 3) |

`--rerun-rf` evaluates RF on the batch-corrected data obtained without covariates
(leakage-free RF, saved as `*_nocovar`), and `--succinct` drops bat/ruv3/svd/rbe/bmc from plots (rbe is limma's removeBatchEffect with the standard design).
Check the `#SBATCH --array` line of each script before submitting; some still list the
iterations of their last rerun. The simulation R scripts (`methods_benchmarking_sim_052626*.R`) run whichever
scenario blocks are uncommented at the bottom of the file.

**Main simulation grid** (OR x cond_effect x batch_effect, count and relab, with and
without library-size confounding, 1000 iterations; config `config_simulate_052626.yml`):

1. data: `slurm_bash_scripts/step0_run_simulate.sh` (`generate_data_MIDASim.R`)
2. R methods: `step1_methods_run_sim_batch_og.sh` (`methods_benchmarking_sim_052626_og.R`)
3. Python methods: `step1_python_methods_run_sim.sh` (`evaluate.py -o 1`)
4. evaluation: `step2_evaluate_run_sim_no.sh`, `step2_evaluate_run_sim_yes.sh` (`-o 2`, add `--rerun-rf` for nocovar RF)
5. plots: `step3_plot_*.sh` (`-o 5`)

**Ablations** (n=150, half taxa, 40% batch, 40% condition; count):

1. data: `step0_run_simulate.sh` (`generate_data_MIDASim_052626.R`)
2. R methods: `step1_methods_run_sim_batch_ablations.sh` and `step1_methods_run_sim_batch_n150.sh`
3. Python methods: `slurm_bash_scripts/step1_python_methods_run_sim_ablations.sh` (`-o 6`)
4. evaluation: `slurm_bash_scripts/step2_evaluate_run_sim_ablation_nocovar_rf.sh` (`-o 7`)
5. plots: `evaluate.py -o 8 -i 1000`

**Condsweep ablation** (OR=1.25, cond_effect in {0.005, 0.01, 0.02, 0.05, 0.1},
batch_effect=0.5; count/relab x no/yes confounding). Each step has a script for
iterations 1-50 and a `_51-1000` script:

1. data: `step0_run_simulate_condsweep*.sh`
2. R methods: `step1_methods_run_sim_condsweep_{countno,countyes,relabno,relabyes}*.sh`
3. Python methods: `step1b_python_methods_run_condsweep*.sh` (`-o 6`, can run alongside step 2)
4. evaluation: `step2_evaluate_condsweep_{countno,countyes,relabno,relabyes}*.sh` (`-o 7`)
5. plots: `step9_visualize_condsweep_*.sh` (`-o 9`), summary across metrics: `python3 analyze_condsweep_results.py`

**Real-world datasets** (`autism_2_microbiomeHD`, `cdi_3_microbiomeHD_clean`, `ibd_3_CMD`,
`crc_11_CMD`; see `configs_README.md` for the config of each):

1. preprocessing: scripts in `data/` (`step0_microbiomeHD_precleaning.py`, `step1_clean_and_prune.R`,
   `step2_preprocessing_summarystats.py`; `step2b_filter_cdi_clean.py` for the CDI-vs-healthy subset,
   `build_crc_11_CMD.R` for crc_11_CMD), output in `data/cleaned_data/`
2. R methods: `methods_benchmarking_fixed_052626.R`, which reads the config from the first line of
   the `Snakefile` (see `slurm_bash_scripts/step1_methods_run_cdi_clean.sh`)
3. Python methods, run from the repo root (config paths are relative to it):
   `python3 benchmark/evaluate.py -o 3 --config <config.yml>`
4. evaluation: `python3 benchmark/evaluate.py -o 4 --config <config.yml>` (add `--rerun-rf` for nocovar RF)
5. figures: `python3 benchmark/evaluate_alt_rw_plots.py [--succinct]` (comparison plots across methods for
   Fig 3/4 b-d and real-world runtime, written to `trial_rw/rw_data_alt_plots[_succinct]/`); PCoA plots
   are written by step 4 to `trial_rw/outputs/<dataset>/`
6. literature check: `pubtator/` (section 4)

`benchmark/one_off_scripts/` and `benchmark/archive/` contain one-off fix/monitoring
scripts and superseded versions; they are not needed to rerun the pipeline.

## 2a. Running the pipeline on real-world datasets (legacy Snakemake interface)

__MOSAIC__ supports datasets from two widely used databases, [CuratedMetagenomicData](https://waldronlab.io/curatedMetagenomicData/) (CMD) and [MicrobiomeHD](https://zenodo.org/records/1146764). See `config.yml` for a MicrobiomeHD example and `./additional_ymls/config_IBD.yml` for a CMD example. For other data, make a folder with one subfolder per dataset to be integrated. Each subfolder needs a file ending in `count_data.csv` (OTU table: taxa as rows/index, samples as columns) and a file ending in `meta_data.csv` (one row per sample). Then edit `config.yml`; the `dataset_name` must not contain `'microbiomeHD'` or `'CMD'`.

Entries in the yml file:
- `rules_to_run`: 'all_rw', 'preprocess', 'integrate_rw' or 'evaluation_rw' for real-world datasets. 'all_rw' runs preprocessing, data integration and evaluation in order.
- `data_source`: 'microbiomeHD', 'CMD', or 'custom'.
- `src`: where the source data is stored (MicrobiomeHD and custom data). Preprocessed data is written to the `cleaned_data` subfolder of `src`.
- `dataset_name`: name of the dataset, ending with its source: 'microbiomeHD', 'CMD', or 'custom'.
- `datatype`: 'count' or 'relab'.
- `RW_LIBSIZE_THRESHOLD`: in preprocessing, samples with library size below this quantile (e.g. 0.05 = 5th percentile) are removed.
- `RW_RELAB_THRESHOLD`: in preprocessing, taxa with relative abundance below this quantile are removed.
- `COVAR_L`: if not empty, samples with NA in any of these covariates are removed.
- `post_integration_outputs`: directory for the batch-corrected data.
- `used_R_methods`: R integration methods to run. Options for count data: `["combat_seq", "limma", "MMUPHin", 'ConQuR', 'ConQuR_libsize']`; for relab data: `["combat", "limma", "MMUPHin", 'ConQuR_rel']`.
- `used_Python_methods`: Python integration methods to run, for count or relab data: `['harmony', "percentile_norm"]`. `"debiasm"` is also supported but needs `pip install DEBIAS-M` first (see Environment set-up above).
- `batch_ref`: reference batch, needed by some methods.
- `covar`: covariates whose biological signal integration should preserve. The first element should be the condition variable of interest, such as `'disease'`.
- `condition_value`: positive value for the condition of interest.
- `evaluation_outputs`: directory for evaluation outputs and summary statistics.
- `binarizing_agent`: the value of the condition variable (e.g. `'disease'`) used to binarize samples for all evaluations: `disease==binarizing_agent` vs all other samples.

Then set the first line of the `Snakefile` to `configfile: "./config.yml"` (or the path of your own config). It must stay on the first line. Run the pipeline from the command line:

```
# run from the directory containing the Snakefile
conda activate mic_data_integration
snakemake -c1 --latency-wait 5
```

This runs the steps set in `rules_to_run`, provided their inputs exist.

## 3. Simulation workflow (`Snakefile_simulation`)

All simulation designs used in the manuscript are built into one Snakemake workflow, configured by
`config_simulation_scenarios.yml`. Each **scenario** is a baseline design with some parameters changed:

| Scenario | Design (changes from the baseline) | Manuscript |
|---|---|---|
| `main_grid_no_libsize_confounding` | OR 1/1.25/1.5 x condition effect 0.25-1 x batch effect 0-0.75 | Fig 2a-c, Supp Figs 2 and 4 |
| `main_grid_libsize_confounded` | same grid, batch membership depends on library size | Supp Figs 1 and 3 |
| `condsweep` | condition effect 0.005-0.1, batch effect 0.5 | Fig 2d |
| `baseline` | 450 samples, all 301 taxa batch-affected, ~70% zeros, balanced batches and conditions, OR 1, condition and batch effect 0.5 | Fig 3 reference |
| `n150` | 150 samples (template library sizes, no resampling) | Fig 3 |
| `halftaxa` | batch effect on 151 of 301 taxa | Fig 3 |
| `zeroinf75` ... `zeroinf90` | smallest non-zero counts set to zero until 75-90% of entries are zero | Fig 3 (90%) |
| `batchfrac40/20/10` | 40/20/10% of samples in batch 1 | Fig 3 (10%) |
| `condfrac40/20/10` | 40/20/10% of samples in condition 1 | Fig 3 (10%) |

Parameters available to any scenario: `n_samples`, `resample_libsize`, `p_cond`, `p_batch`, `or_l`,
`cond_effect_val_l`, `batch_effect_val_l`, `batch_libsize_related`, `batch_taxa_fraction`, `zero_fraction`
(see the comments in the config). To add a design, add an entry under `scenarios:` and list it in
`run_scenarios`.

```bash
conda activate mic_data_integration
# from the repository root
snakemake -s Snakefile_simulation -n                                   # dry run: what would be run
snakemake -s Snakefile_simulation --cores 4 --config run_scenarios=[zeroinf90] iterations=5
bash run_simulation_workflow.sh                                        # full study on SLURM (one job per step)
bash run_simulation_workflow.sh --config run_scenarios=[n150,halftaxa] datatypes=[count,relab]
```

For every replicate the workflow simulates the data (`benchmark/simulate_scenario.R`, seed = replicate
number; `benchmark/apply_zero_inflation.py` for zero-inflation scenarios), applies the R strategies with and
without the condition covariate (`benchmark/run_sim_methods.R`) and the Python strategies (`evaluate.py -o 16`),
evaluates them (`evaluate.py -o 17`, then `--rerun-rf` for the leakage-free random forest), and finally
summarizes each scenario over replicates (`benchmark/aggregate_scenario.py`). Outputs, under `output_root`
(default `<sim_output_root>/mosaic_scenarios`):

```
<scenario>/data/                     simulated count and relative abundance tables, metadata, true condition/batch taxa
<scenario>/output_<datatype>[_nocovar]/   integrated data
<scenario>/eval_<datatype>/          per-replicate evaluation (global_benchmarking_stats[_nocovar]_batch_0.csv, ...)
<scenario>/summary_<datatype>.csv    mean, 95% CI, SD, median, 2.5/97.5 percentiles per setting x method x metric
<scenario>/values_<datatype>.csv     all metrics per replicate
ablation_figures_<datatype>/         when `baseline` and ablations are run: each ablation vs the baseline per
                                     metric and method, and Figure 3 (benchmark/figures/)
```

Steps leave flag files in `<scenario>/.flags/`, so an interrupted run resumes where it stopped. The
simulation (`benchmark/sim_lib/midas_sim_lib.R`) and integration (`benchmark/sim_lib/sim_methods_lib.R`) code is
shared with the scripts used for the manuscript, and reproduces their output exactly for the same seed and
software versions (see `benchmark/tests/test_simulation_workflow.sh`).

The older single-iteration Snakemake rules for simulations in `Snakefile` (`rules_to_run: all_simulated`) are
kept for reference only; use `Snakefile_simulation` instead.

## 3a. Manuscript figures (`benchmark/figures/`)

The scripts that build Figures 2-3 and Supplemental Figures 1-5 from the evaluation results are in
`benchmark/figures/`; `benchmark/figures/README.md` lists which script makes which figure. The website in `docs/`
(GitHub Pages) has a how-to-run guide and a vignette that follows one simulated dataset to each figure.

## 4. PubTator literature check (`pubtator/`)

As a check on the real-world differential abundance results, each significant taxon is searched in
PubTator3 together with five clinical key terms per disease. A taxon is "strongly relevant" if it has
at least 1 publication for all five terms, or at least 10 publications for three of them; the number of
strongly relevant taxa per method gives Supplemental Table 3, and the per-taxon counts Supplemental
Table 4. See `pubtator/README.md` for the steps.

## License

The MOSAIC code is licensed under the GNU General Public License v3.0 (`LICENSE`). The bundled ConQuR code
(`ConQuR/`) is also GPL-3.0. The real-world data in `data/` keep the licenses of their sources; the
MicrobiomeHD-derived files are CC BY-NC 4.0 (non-commercial). See `data/README.md`.

## 5. Acknowledgements
This project started as a class project in Dr. Wesley Tansey's class Foundation of Data Science at Memorial Sloan Kettering and Weill Cornell, and would not have been possible without the generous support of Dr. Wodan Ling, Jiuyao Lu, Dr. Quaid Morris, and Dr. Wesley Tansey. The funding support for my PhD comes from the Tri-institutional Program of Computational Biology and Medicine.

If you have questions regarding this benchmarking project or other inquiries, please reach out to me at ___chf4012@med.cornell.edu___. I hope you have a great day!

Cheers,\
Tom Fu

