# microbiome batch correction benchmarking project
Code and partial data repository for the manuscript entitled __A comprehensive benchmark of batch correction methods in microbiome data__ by _Fu et al._. 

_Developed and maintained by_:  [Chenlian (Tom) Fu](fcl200089@outlook.com)\
_Advised by_: Dr. Wodan Ling, Dr. Quaid Morris\
_Affiliations_: Weill Cornell Medicine, Memorial Sloan Kettering Cancer Center


This github repository stores the code for benchmarking microbiome batch correction methods, enabling readers to reproduce the analyses done in the paper as well as leveraging their own data along with the provided code to do similar analyses. Our batch correction method evaluation consists of three steps:
- Data generation/collection
- Running batch correction methods on datasets to generate batch corrected data
- Evaluation of the methods


## Table of content
- 0. [Environment set up](https://github.com/tommyfuu/mic_bc_benchmark/blob/main/README.md#0-environment-set-up)
- 1. [(Simulated) Data generation/ (real-world data) collection](https://github.com/tommyfuu/mic_bc_benchmark/blob/main/README.md#1-data-generationcollection)
    - 1.1 [Simulation data generation](https://github.com/tommyfuu/mic_bc_benchmark/blob/main/README.md#111-simulation-data-generation)
    - 1.2 [Real world microbiome data collection and cleaning](https://github.com/tommyfuu/mic_bc_benchmark/blob/main/README.md#12-real-world-microbiome-data-collection-and-cleaning)
- 2. [Running/benchmarking methods](https://github.com/tommyfuu/mic_bc_benchmark/blob/main/README.md#2-benchmarking-methods)
    - 2.1 [For methods implemented in R](https://github.com/tommyfuu/mic_bc_benchmark/blob/main/README.md#21-for-methods-implemented-in-r)
        - 2.1.1 [Running these methods for simulation](https://github.com/tommyfuu/mic_bc_benchmark/blob/main/README.md#211-running-these-methods-for-simulation)
        - 2.1.2 [Running these methods for real-world datasets](https://github.com/tommyfuu/mic_bc_benchmark/blob/main/README.md#212-running-these-methods-for-real-world-datasets)
    - 2.2 [For methods implemented in Python](https://github.com/tommyfuu/mic_bc_benchmark/blob/main/README.md#22-for-methods-implemented-in-python)
        - 2.2.1 [Running these methods for simulation](https://github.com/tommyfuu/mic_bc_benchmark/blob/main/README.md#221-simulation)
        - 2.2.2 [Running these methods for real-world datasets](https://github.com/tommyfuu/mic_bc_benchmark/blob/main/README.md#222-real-world-dataset)
- 3. [Evaluation](https://github.com/tommyfuu/mic_bc_benchmark/blob/main/README.md#3-evaluation)
    - 3.1 [Generating summary statistics for each simulated dataset or for each real-world dataset](https://github.com/tommyfuu/mic_bc_benchmark/blob/main/README.md#31-generating-summary-statistics-for-each-simulated-dataset-or-for-each-real-world-dataset)
        - 3.1.1 [Summary stats for simulated datasets](https://github.com/tommyfuu/mic_bc_benchmark/blob/main/README.md#311-simulated-datasets)
        - 3.1.2 [Summary stats for real-world datasets](https://github.com/tommyfuu/mic_bc_benchmark/blob/main/README.md#312-real-world-dataset)
    - 3.2 [Generating intuitive visualizations based on the summary statistics](https://github.com/tommyfuu/mic_bc_benchmark/blob/main/README.md#32-generating-intuitive-visualizations-based-on-the-summary-statistics)

## 0. Environment set up

You can potentially set up the environment using conda by executing the following command:

```
conda env create -f environment.yml
```

The above way to set up the environment might lead to deprecated versions of files. To create a clean environment with all the necessary libraries, you can also manually install the following packages in a virtual environment with Python 3.8(.17) and R 4.3:

- Python packages: pandas, numpy, matplotlib, scikit-learn, seaborn, scipy, skbio, rpy2, statsmodels
- R packages: phyloseq, bindata, MIDAS, tibble, xtable, sva, limma, vegan, MMUPHin, FDboost, doParallel, dplyr, readr, mixOmics, parallel, ade4, compositions, cqrReg, fastDummies, coda.base


## 1. Data generation/collection

### 1.1.1 Simulation data generation

A comprehensive evaluation requires stringent simulation and corresponding analyses. Here, we employ [MIDAS](https://pubmed.ncbi.nlm.nih.gov/36993431/), an intuitive microbiome data simulator that recapitulates microbiome data structure without distributional assumptions while allowing users to manipulate variables of interest, in our case batch effect and conditional effect. We further edited codes to enable the simulation to include varying odds ratio (between biological signals and batch in each taxa), existing and non-existing relationship between batch effect and library size, as well as the generation of both count and relative abundance data.

To reproduce one iteration of our data using MIDAS, one can do the following in command line to generate minimally viable simulated dataset

```
cd benchmark # enter the benchmark directory
R # to initiate the R environment
source('./generate_data_MIDAS.R') # load generation file
# the odds ratio/conditional (biological) effect/batch effect you wish to generate the data for
or_l = c(1.25)
cond_effect_val_l = c(0.5)
batch_effect_val_l = c(0.5)
output_root = './trial/simulate' # please make sure the output_root exists
# run functions
scaled_slurm_midas_data_generation(output_root, otu_original, n, or_l, cond_effect_val_l, batch_effect_val_l, iter=GLOBAL_ITER, batch_libsize_related = FALSE, libsize_l=sampled_libsize_l)
```

To explain the arguments above:
- `otu_original` is real-world microbiome data we generate the simulation upon (default is a cleaned 150-sample 301-taxa ibd cohort from the Human Microbiome Project, loaded by `ibd_150.Rdata` in the same directory).
- `n` is the number of samples in the simulated dataset (default=450).
- `or_l, cond_effect_val_l, batch_effect_val_l`: three parameter lists to generate data for.
- `iter` is simulation iteration for data generation (default=1 without user selection).
- `batch_libsize_related`, a boolean decides whether the simulated dataset has an established relationship between batch effect and library, if `TRUE` then we turn the odds ratio into binary correlation between the occurence likehood of batch effect on a taxa and library size, otherwise the relationship is random.
- `libsize_l`, a list of library_sizes for samples defined, only relevant if `batch_libsize_related == TRUE`.

This code chunk will generate a count and its corresponding relative abundance datasets for each paratemer combination in the three lists defined above, along with the metadata for samples (which ones in which batch) as well as metadata for taxa (which ones are ground truth perturbed biomarkers) for later experiments.

Note that you may always replace the input parameter list with a complete list so that you generate data for all parameter combinations, which might take a bit longer:

```
or_l = c(1, 1.25, 1.5)
cond_effect_val_l = c(0, 0.25, 0.5, 0.75, 1)
batch_effect_val_l = c(0, 0.25, 0.5, 0.75, 1)
```

To run the simulation script in scale, in the folder `benchmark/slurm_bash_scripts` there is a bash script called `step0_run_simulate_sim.sh`, which one can revise to generate their own slurm bash scripts based on and move back into the `benchmark` folder for running the `generate_data_MIDAS.R` script in scale with `slurm`. Note that for running MIDAS data generation in slurm, you might need to uncomment the commented lines at the bottom of the `generate_data_MIDAS.R` script.

### 1.2 Real world microbiome data collection and cleaning

A comprehensive evaluation also requires well-collected and well-preprocessed real-world microbiome data to ensure that the theoretically excelling batch correction methods from simulation experiments work in practice to enabled biomarker and trend discovery in batch corrected datasets while removing batch effect.

To this end, we collect and clean both count and relative abundance real world datasets from two databases, [MicrobiomeHD](https://zenodo.org/records/569601) and [CuratedMetagenomicsData](https://waldronlab.io/curatedMetagenomicData/) (CMD). The real world data preparation can be done by running the following code in the command line in the data directory:

- pre-clean the two microbiomeHD datasets, moving them from different batches to the same files
```
ipython # enter interactive python mode
> from step0_microbiomeHD_precleaning import * 
> overall_path = '/athena/linglab/scratch/chf4012'
> # autism 2 microbiomeHD
> output_dir_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/intermediate_autism_2_microbiomeHD/autism_2_microbiomeHD'
> address_directory = overall_path+'/mic_bc_benchmark/data/autism_2_microbiomeHD'
> data_mat, meta_data = load_data_microbiomeHD(address_directory, output_dir_path)
> # cdi 3 microbiomeHD
> output_dir_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/intermediate_cdi_3_microbiomeHD/cdi_3_microbiomeHD'
> address_directory = overall_path+'/mic_bc_benchmark/data/cdi_3_microbiomeHD'
> data_mat, meta_data = load_data_microbiomeHD(address_directory, output_dir_path
```

- loads the four datasets into phyloseq objects, clean and prune samples and taxa below 0.05% abundance, and save to a folder
```
R # enter intractive R mode
> source('step1_clean_and_prune.R')
> overall_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/'
> autism_phyloseq_obj <- load_phyloseq_from_merged_microbiomeHD(paste0(overall_path, 'intermediate_autism_2_microbiomeHD/autism_2_microbiomeHD_count_data.csv'), paste0(overall_path, 'intermediate_autism_2_microbiomeHD/autism_2_microbiomeHD_meta_data.csv'))
> clean_prune_save_phyloseq(autism_phyloseq_obj, 'autism_2_microbiomeHD', 0.05, 0.05, save = TRUE, save_to = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/pruned_autism_2_microbiomeHD')

> cdi_phyloseq_obj <- load_phyloseq_from_merged_microbiomeHD(paste0(overall_path, 'intermediate_cdi_3_microbiomeHD/cdi_3_microbiomeHD_count_data.csv'), paste0(overall_path, 'intermediate_cdi_3_microbiomeHD/cdi_3_microbiomeHD_meta_data.csv'))
> clean_prune_save_phyloseq(cdi_phyloseq_obj, 'cdi_2_microbiomeHD', 0.05, 0.05, save = TRUE, save_to = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/pruned_cdi_3_microbiomeHD')

> ibd_phyloseq_obj <- load_phyloseq_from_merged_CMD(c("HMP_2019_ibdmdb", "LiJ_2014", "NielsenHB_2014"), c("IBD", 'healthy'))
> clean_prune_save_phyloseq(ibd_phyloseq_obj[[1]], ibd_phyloseq_obj[[2]], 0.05, 0.05, save = TRUE, save_to = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/pruned_ibd_3_CMD')

> crc_phyloseq_obj <- load_phyloseq_from_merged_CMD(c("FengQ_2015", "HanniganGD_2017", "ThomasAM_2018a", "YachidaS_2019", "ZellerG_2014"), c("adenoma", "healthy"))
> clean_prune_save_phyloseq(crc_phyloseq_obj[[1]], crc_phyloseq_obj[[2]], 0.05, 0.05, save = TRUE, save_to = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/pruned_crc_8_CMD')
```

- loads the phyloseq formatted data from step1 to save to a directory in a format easily digestible by later steps of the pipeline, along with generating a complete confounding check table
```
ipython
> from step2_preprocessing_summarystats import *
> overall_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data'
> # autism_2_microbiomeHD
> data_mat, meta_data = preprocess_data_phyloseq(f'{overall_path}/pruned_autism_2_microbiomeHD', f'{overall_path}/cleaned_data/autism_2_microbiomeHD/autism_2_microbiomeHD', id = 'Sam_id', covar_l = [], relab = False)
> data_mat, meta_data = load_results_from_benchmarked_methods(f'{overall_path}/cleaned_data/autism_2_microbiomeHD/autism_2_microbiomeHD_count_data.csv', f'{overall_path}/cleaned_data/autism_2_microbiomeHD/autism_2_microbiomeHD_meta_data.csv')
> check_complete_confounding(meta_data, 'Dataset', 'DiseaseState', f'{overall_path}/cleaned_data/autism_2_microbiomeHD/autism_2_microbiomeHD')

> # cdi_3_microbiomeHD
> data_mat, meta_data = preprocess_data_phyloseq(f'{overall_path}/pruned_cdi_3_microbiomeHD', f'{overall_path}/cleaned_data/cdi_3_microbiomeHD/cdi_3_microbiomeHD', id = 'Sam_id', covar_l = [], relab = False)
> data_mat, meta_data = load_results_from_benchmarked_methods(f'{overall_path}/cleaned_data/cdi_3_microbiomeHD/cdi_3_microbiomeHD_count_data.csv', f'{overall_path}/cleaned_data/cdi_3_microbiomeHD/cdi_3_microbiomeHD_meta_data.csv')
> check_complete_confounding(meta_data, 'Dataset', 'DiseaseState', f'{overall_path}/cleaned_data/cdi_3_microbiomeHD/cdi_3_microbiomeHD')

> # ibd_3_CMD
> data_mat, meta_data = preprocess_data_phyloseq(f'{overall_path}/pruned_ibd_3_CMD', f'{overall_path}/cleaned_data/ibd_3_CMD/ibd_3_CMD', id = 'Sam_id', covar_l = ['disease', 'gender', 'age_category'])
> data_mat, meta_data = load_results_from_benchmarked_methods(f'{overall_path}/cleaned_data/ibd_3_CMD/ibd_3_CMD_count_data.csv', f'{overall_path}/cleaned_data/ibd_3_CMD/ibd_3_CMD_meta_data.csv')
> check_complete_confounding(meta_data, "study_name", "disease", f'{overall_path}/cleaned_data/ibd_3_CMD/ibd_3_CMD')

> # crc_8_CMD
> data_mat, meta_data = preprocess_data_phyloseq(f'{overall_path}/pruned_crc_8_CMD', f'{overall_path}/cleaned_data/crc_8_CMD/crc_8_CMD', id = 'Sam_id', covar_l = [])
> data_mat, meta_data = load_results_from_benchmarked_methods(f'{overall_path}/cleaned_data/crc_8_CMD/crc_8_CMD_count_data.csv', f'{overall_path}/cleaned_data/crc_8_CMD/crc_8_CMD_meta_data.csv')
> check_complete_confounding(meta_data, "study_name", "disease", f'{overall_path}/cleaned_data/crc_8_CMD/crc_8_CMD')
```

These three code chunks result in the four subfolders in `data/cleaned_data`, each containing three files: the count (relative abundance) OTU data, the corresponding metadata, and the complete confounding checking table (number of samples in each biological condition in each cleaned dataset).

Note that while the crc dataset started with 8 batches, after pre-cleaning, there are actually 5 datasets (batches) that participate in batch correction.

## 2. Benchmarking methods

### 2.1 For methods implemented in R

For methods implemented in R including `combat (combat/combat-seq), limma, MMUPHin, ConQuR (ConQuR/ConQuR_libsize/ConQuR_rel)`, they can be run on a dataset along the preprocessing steps a dataset needs prior to running each of these methods using the scripts `/benchmark/methods_benchmarking.R` for the real-world datasets and `/benchmark/methods_benchmarking_sim.R` for simulation datasets. Note that before running, please edit line 27 of `methods_benchmarking_sim.R` and line 29 of `methods_benchmarking.R` to ensure that you link to the correct folder path for the folder `ConQuR` under the main source folder.

#### 2.1.1 Running these methods for simulation
To run the methods in scale in slurm, relevant parallelization and iterative running has been set up in the script `methods_benchmarking_sim.R`. For example, one can run the following to run all data on iteration `1` of the simulated dataset for running methods and generating batch corrected datasets for the minimally viable parameter set:
```
source('./methods_benchmarking_sim.R')
overall_path = './trial/simulate'
output_dir = './trial/simulation_data_output_count_norelation' # pls make sure this dir exists
or_l = c(1.25)
cond_effect_val_l = c(0.5)
batch_effect_val_l = c(0.5)
method_l = c("combat_seq", "limma", "MMUPHin", 'ConQuR', 'ConQuR_libsize')
scaled_slurm_methods_bencharking(output_dir, overall_path, method_l, or_l, cond_effect_val_l, batch_effect_val_l, GLOBAL_ITER, count = TRUE)
```

Note that the script has been written in a way that makes it easy to call and run on different iterations of input, as long as you uncomment the last commented code chunk in the `methods_benchmarking_sim.R`. For example, running `Rscript methods_benchmarking_sim.R 10` will run iteration 10. To scale run simulation in slurm, one can reference file `/benchmark/slurm_bash_scripts/step1_methods_run_sim_batch.sh` and move the revised file back to `/benchmark` before running on slurm in scale.

#### 2.1.2 Running these methods for real-world datasets
Similarly to running these methdos for simulation, you can run the methods for the real-world dataset by calling the script like below after running `R` in the `benchmark` folder:

```
source('./methods_benchmarking.R')
# autism 2 microbiomeHD
current_root = '../data/cleaned_data/autism_2_microbiomeHD/autism_2_microbiomeHD'
output_root = './trial/rw_methods_out/autism_2_microbiomeHD/autism_2_microbiomeHD'
run_methods(paste0(current_root, "_count_data.csv"),
    paste0(current_root, "_meta_data.csv"),
    output_root,
    dataset = "Dataset",
    batch_ref = 'asd_son',
    covar = c("DiseaseState"),
    count = TRUE,
    used_methods =  c("combat_seq", "limma", "MMUPHin", 'ConQuR', 'ConQuR_libsize')
)
```

To simplify the running of methods on each of the four real-world dataset, one can simply run the following in command line.
```
Rscript methods_benchmarking.R 1 # 1 for autism_2_microbiomeHD, 2 for cdi_3_microbiomeHD, 3 for ibd_3_CMD, 4 for crc_8_CMD
```

And of course, you can run this in scale by referencing `/benchmark/slurm_bash_scripts/step1_methods_run_rw.sh` and move the revised file back to `/benchmark` for slurm runs.

### 2.2 For methods implemented in Python

#### 2.2.1 Simulation

For methods implemented in Python including `harmony, percentile_normalization`, they can be run on a dataset along the preprocessing steps a dataset needs prior to running each of these methods using the script `/benchmark/evaluate.py -o 1`. Similarly to before, one can run this option in both an iterative python kernel (iPython) or through command line.

In iPython, once you enter the `benchmark` directory, you may do the following for simulation as a minimum viable example:
```
from evaluate import *
overall_path = './trial' # pls make sure this dir exists
or_l = [1.25]
cond_effect_val_l = [0.5]
batch_effect_val_l = [0.5]
GLOBAL_DATATYPE = 'count'
related = 'no'
binarizing_agent_biovar = 'cond_1'
iterative_methods_running_evaluate(run_or_evaluate = 'run', datatype = GLOBAL_DATATYPE, iter = 1, or_l = or_l, cond_effect_val_l = cond_effect_val_l, batch_effect_val_l = batch_effect_val_l, 
                    address_XY_dir_path = './trial/simulate', 
                    output_dir_path = f'./trial/simulation_data_output_{GLOBAL_DATATYPE}_{related}relation', 
                    eval_dir_path = overall_path+f"/simulation_data_eval_{GLOBAL_DATATYPE}_{related}relation",
                    methods_list = methods_list_dict[GLOBAL_DATATYPE], 
                    binarizing_agent_biovar = binarizing_agent_biovar)

```

To explain the parameters, we first note that `iterative_methods_running_evaluate` is a generic wrapper function, where arguments
- `run_or_evaluate`, a `str` type argument, should be either `run` or `evaluate`. In this case, you should use `run` because we are running the two batch correction methods implemented in python; for the later evaluation section, you use `evaluate`
- `datatype`, a `str` type argument, should be either `count` or `relab`. 
- `iter` is simulation iteration for data generation (default=1 without user selection).
- `or_l, cond_effect_val_l, batch_effect_val_l`: three parameter lists to generate data for.
- `address_XY_dir_path`, a `str` type argument denoting the path of the directory where the original (uncorrected) count/relab dataset and the corresponding metadata (Y) should be.
- `output_dir_path`, a `str` type argument denoting the path of the directory where the batch corrected datasets after running the methods and the related statistics should be.
- `eval_dir_path`, a `str` type argument denoting the path of the directory where the evaluation results after analyzing the characteristics of each corrected dataset should be. (Need to be defined but not super relevant for this module)
- `methods_list`, a `list` of `str` denoting the methods whose batch corrected datasets are then evaluated. For example, it can be `["nobc", "harmony", "combat_seq", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]` for count datasets, and `["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR_rel", "percentile_norm"]` for relab sets. (Need to be defined but not super relevant for this module)
- `binarizing_agent_biovar`, a `str` type argument denoting the name of the column in the input metadata for the column of the binary biological condition variable. For the simulated dataset in our setup, it should be `cond_1`.

To make it easier for streamlining, we also defined the following:
- `related`, a `str` type argument, should be either `yes` or `no` denoting whether confounding exists or not between library size and the simulated batch effect.

Alternatively, we provide a command line option, which might require you to go in and change certain codes in `evaluate.py`, but allows you a more streamlined running experience:

1. Replace the parameter combination set up (line 1119-1122).
2. Ensure that the methods run and the paths where the input files are taken and output files generated are correct (lines 1132-1136)
3. Run `python3 /benchmark/evaluate.py -o 1 -d relab -r yes -i 1 -a cond_1 -p /athena/linglab/scratch/chf4012`. You can easily replace `relab` with `count` (the default) or `yes` with `no` (default), 1 with other iteration numbers, and the path with the overall path where all your files are generated. For details, please check lines 38-43 of the file `evaluate.py`.

For this alternative option, you can scale run with slurm by referencing the file `/benchmark/slurm_bash_scripts/step1_python_methods_run_sim.sh`, revising it, and moving it to `benchmark` folder to slurm run.

#### 2.2.2 Real-world dataset

For the iPython option,
```
vars_use = ["Dataset"]
IDCol = 'Sam_id'
overall_path = './trial' # pls make sure this dir exists
address_X = '../data/cleaned_data/autism_2_microbiomeHD/autism_2_microbiomeHD_count_data.csv'
address_Y = '../data/cleaned_data/autism_2_microbiomeHD/autism_2_microbiomeHD_meta_data.csv'
data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
res_h, meta_data_h = generate_harmony_results(data_mat, meta_data, IDCol, vars_use, overall_path+"/rw_methods_out/autism_2_microbiomeHD/autism_2_microbiomeHD_harmony")
percentile_norm(address_X, address_Y, "DiseaseState", "ASD", "comma", overall_path+"/rw_methods_out/autism_2_microbiomeHD/autism_2_microbiomeHD")
```

One can scale run the methods on simulated dataset by referencing the lines containing `-o 1` in the slurm file `/benchmark/slurm_bash_scripts/evaluate_run_sim.sh` and scale run in slurm.


For the ease of running, one can use option `4` in the `evaluate.py` script once they have checked to ensure all the variables and file paths are correct in lines 1196-1236.


## 3. evaluation

Comprehensive evaluation is implemented in the script `/benchmark/evaluate.py`, with partial visualizations for the real-world datasets implemented in `/benchmark/evaluate_alt_rw_plots.py`. 

Similarly to the methods benchmarking, for simulation, it can be wrapped and run in command line using options `2, 3` in `/benchmark/evaluate.py` with the `iterative_methods_running_evaluate` function with `run_or_evaluate = 'evaluate'`. On the other hand, for real-world dataset, it can be wrapped and run in command line using 
options `4, 5` in `/benchmark/evaluate.py` with the `iterative_methods_running_evaluate` function with `run_or_evaluate = 'evaluate'`. For reproducibility, the line plots for real-world evaluations are implemented in `/benchmark/evaluate_alt_rw_plots.py`.

The evaluation section is divided into two sections: __3.1__ generating summary statistics for each simulated dataset or for each real-world dataset; and __3.2__ generating intuitive visualizations based on the summary statistics from __3.1__.

### 3.1 Generating summary statistics for each simulated dataset or for each real-world dataset

#### 3.1.1 Simulated datasets
To aid your intuitive understanding, you can do the following for a minimally viable run in iPython. This will generate a folder of files similar to `benchmark/trial/simulation_data_eval_count_norelation/out_1.25_0.5_0.5_iter_1`, which has output subfolders consisting of evaluation results from each method, a csv file with global summary stats, and a pdf figure for visualizing batch effect removal across all methods.

```
from evaluate import *
overall_path = './trial' # pls make sure this dir exists
or_l = [1.25]
cond_effect_val_l = [0.5]
batch_effect_val_l = [0.5]
GLOBAL_DATATYPE = 'count'
related = 'no'
binarizing_agent_biovar = 'cond_1'
iterative_methods_running_evaluate(run_or_evaluate = 'evaluate', datatype = GLOBAL_DATATYPE, iter = 1, or_l = or_l, cond_effect_val_l = cond_effect_val_l, batch_effect_val_l = batch_effect_val_l, 
                    address_XY_dir_path = './trial/simulate', 
                    output_dir_path = f'./trial/simulation_data_output_{GLOBAL_DATATYPE}_{related}relation', 
                    eval_dir_path = overall_path+f"/simulation_data_eval_{GLOBAL_DATATYPE}_{related}relation",
                    methods_list = methods_list_dict[GLOBAL_DATATYPE], 
                    binarizing_agent_biovar = binarizing_agent_biovar)

```

Alternatively, we provide a command line option, which might require you to go in and change certain codes in `evaluate.py`, but allows you a more streamlined running experience:

1. Replace the parameter combination set up (line 1119-1122).
2. Ensure that the methods run and the paths where the input files are taken and output files generated are correct (lines 1140-1144)
3. Run `python3 /benchmark/evaluate.py -o 2 -d relab -r yes -i 1 -a cond_1 -p /athena/linglab/scratch/chf4012`. You can easily replace `relab` with `count` (the default) or `yes` with `no` (default), 1 with other iteration numbers, and the path with the overall path where all your files are generated. For details, please check lines 38-43 of the file `evaluate.py`.

For this alternative option, you can scale run with slurm by referencing the file `/benchmark/slurm_bash_scripts/step2_evaluate_run_sim.sh`, revising it, and moving it to `benchmark` folder to slurm run.

#### 3.1.2 Real-world dataset

For the iPython option, you can do the following for a minimally viable run in iPython. This will generate a folder of files similar to `benchmark/trial/autism_2_microbiomeHD_eval`, which has output subfolders consisting of evaluation results from each method, a csv file with global summary stats, and a pdf figure for visualizing batch effect removal across all methods.
```
# autism 2 microbiomeHD
output_dir_path = './trial/autism_2_microbiomeHD_eval'
address_directory = '../data/cleaned_data/autism_2_microbiomeHD'
vars_use = ["Dataset"]
IDCol = 'Sam_id'

address_Y = "../data/cleaned_data/autism_2_microbiomeHD/autism_2_microbiomeHD_meta_data.csv"

# nobc
address_X = "../data/cleaned_data/autism_2_microbiomeHD/autism_2_microbiomeHD_count_data.csv"
data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
### nobc
Evaluate(data_mat, meta_data, 'Dataset', output_dir_path + '/output_autism_2_microbiomeHD_nobc/autism_2_microbiomeHD_nobc', "DiseaseState", 30, [], 'Sam_id', method = 'nobc')

### harmony
address_X = "./trial/rw_methods_out/autism_2_microbiomeHD/autism_2_microbiomeHD_harmony_adjusted_count.csv"
res_h, meta_data_h = load_results_from_benchmarked_methods(address_X, address_Y)
Evaluate(res_h, meta_data_h, 'Dataset', output_dir_path + '/output_autism_2_microbiomeHD_harmony/autism_2_microbiomeHD_harmony', "DiseaseState", 30, [], 'Sam_id', method = 'harmony')

# benchmarking other methods: 
### combat (combat_seq)
address_X = "./trial/rw_methods_out/autism_2_microbiomeHD/autism_2_microbiomeHD_combat_seq.csv"
data_mat_combat, meta_data_combat = load_results_from_benchmarked_methods(address_X, address_Y)
Evaluate(data_mat_combat, meta_data_combat, 'Dataset', output_dir_path + '/output_autism_2_microbiomeHD_combat_seq/autism_2_microbiomeHD_combat_seq', "DiseaseState", 30, [], 'Sam_id', method = 'combat_seq')

### limma
address_X = "./trial/rw_methods_out/autism_2_microbiomeHD/autism_2_microbiomeHD_limma.csv"
data_mat_limma, meta_data_limma = load_results_from_benchmarked_methods(address_X, address_Y)
Evaluate(data_mat_limma, meta_data_limma, 'Dataset', output_dir_path + '/output_autism_2_microbiomeHD_limma/autism_2_microbiomeHD_limma', "DiseaseState", 30, [], 'Sam_id', method = 'limma')

### MMUPHin
address_X = "./trial/rw_methods_out/autism_2_microbiomeHD/autism_2_microbiomeHD_MMUPHin.csv"
data_mat_mmuphin, meta_data_mmuphin = load_results_from_benchmarked_methods(address_X, address_Y)
Evaluate(data_mat_mmuphin, meta_data_mmuphin, 'Dataset', output_dir_path + '/output_autism_2_microbiomeHD_MMUPHin/autism_2_microbiomeHD_MMUPHin', "DiseaseState", 30, [], 'Sam_id', method = 'MMUPHin')

### ConQuR
address_X = "./trial/rw_methods_out/autism_2_microbiomeHD/autism_2_microbiomeHD_ConQuR.csv"
data_mat_conqur, meta_data_conqur = load_results_from_benchmarked_methods(address_X, address_Y)
Evaluate(data_mat_conqur, meta_data_conqur, 'Dataset', output_dir_path + '/output_autism_2_microbiomeHD_ConQuR/autism_2_microbiomeHD_ConQuR', "DiseaseState", 30, [], 'Sam_id', method = 'ConQuR')

### ConQuR_libsize
address_X = "./trial/rw_methods_out/autism_2_microbiomeHD/autism_2_microbiomeHD_ConQuR_libsize.csv"
data_mat_conqur_libsize, meta_data_conqur_libsize = load_results_from_benchmarked_methods(address_X, address_Y)
Evaluate(data_mat_conqur_libsize, meta_data_conqur_libsize, 'Dataset', output_dir_path + '/output_autism_2_microbiomeHD_ConQuR_libsize/autism_2_microbiomeHD_ConQuR_libsize', "DiseaseState", 30, [], 'Sam_id', method = 'ConQuR_libsize')

### percentile_norm
address_X = "./trial/rw_methods_out/autism_2_microbiomeHD/autism_2_microbiomeHD_percentile_norm.csv"
data_mat_percentile_norm, meta_data_percentile_norm = load_results_from_benchmarked_methods(address_X, address_Y)
Evaluate(data_mat_percentile_norm, meta_data_percentile_norm, 'Dataset', output_dir_path + '/output_autism_2_microbiomeHD_percentile_norm/autism_2_microbiomeHD_percentile_norm', "DiseaseState", 30, [], 'Sam_id', method = 'percentile_norm')

### global evaluation
input_frame_path = "../data/cleaned_data/autism_2_microbiomeHD/autism_2_microbiomeHD_count_data.csv"
bio_var = "DiseaseState"
dataset_name = "autism_2_microbiomeHD"
methods_list = ["nobc", "harmony", "combat_seq", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
global_eval_dataframe(input_frame_path, bio_var, dataset_name, methods_list, "./trial/rw_methods_out", './trial/autism_2_microbiomeHD_eval', simulate = False)

### multi-method plot
df_l = [data_mat, res_h, data_mat_combat, data_mat_limma, data_mat_mmuphin, data_mat_conqur, data_mat_conqur_libsize, data_mat_percentile_norm]
methods = ["nobc", "harmony", "combat_seq", "limma", "MMUPhin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
meta_data_l = [meta_data, meta_data_h, meta_data_combat, meta_data_limma, meta_data_mmuphin, meta_data_conqur, meta_data_conqur_libsize, meta_data_percentile_norm]
plot_PCOA_multiple('autism_2_microbiomeHD', df_l, methods, meta_data_l, used_var="Dataset", output_root= output_dir_path + '/')

```

### 3.2 Generating intuitive visualizations based on the summary statistics

After generating summary statistics over 1000 iterations in simulations, or across all real-world datasets, it is important to visualize in an intuitive manner to allow the users to understand what methods are better in what scenarios. Different from the other sections, this section requires users to either have run all iterations of simulations, or to have run all the real-world datasets of interest. This is because these summary plots would not make sense without all the iterations/diversity of real-world datasets.

This section is made easy with the options `3` of `evaluate.py` for simulation dataset and directly running `python3 evaluate_alt_rw_plots` and for real-world datasets respectively. These generate line plots to visualize the statistics of different batch corrected methods in various simulation settings and in various real-world datasets. If you seek to run these codes, please double check the lines 1146-1193 of `evaluate.py` and lines 299-311 of `evaluate_alt_rw_plots` carefully to ensure the path of the files are correct. Note that there is a certain level of customness to the codes to ensure the generated figures are good-looking.

For example, one can run the following code for simulation:
```
python3 evaluate.py -o 3 -d relab -r yes
```

And the following code for the real-world data:
```
python3 evaluate_alt_rw_plots.py
```


## 4. Additional notes

This project started as a class project in Dr. Wesley Tansey's class Foundation of Data Science at Memorial Sloan Kettering and Weill Cornell, and would not have been possible without the generous support of Dr. Wodan Ling, Jiuyao Lu, Dr. Quaid Morris, and Dr. Wesley Tansey. The funding support for my PhD comes from the Tri-institutional Program of Computational Biology and Medicine.

If you have questions regarding this benchmarking project or other inquries, please reach out to me at ___chf4012@med.cornell.edu___. I hope you have a great day!

Cheers,\
Tom Fu

