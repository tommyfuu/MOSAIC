# microbiome batch correction benchmarking project
Code and partial data repository for the manuscript entitled __A comprehensive benchmark of batch correction methods in microbiome data__ by _Fu et al.. 

_Developed and maintained by_:  [Chenlian (Tom) Fu](fcl200089@outlook.com)\
_Advised by_: Dr. Wodan Ling, Dr. Quaid Morris\
_Affiliations_: Weill Cornell Medicine, Memorial Sloan Kettering Cancer Center

## 0. Environment set up

You can potentially set up the environment using conda by executing the following command:

```
conda env create -f /benchmark/final_env.yml
```

Or, you can manually install the following packages:

- Python packages: pandas, numpy, matplotlib, scikit-learn, seaborn, scipy, skbio, rpy2, statsmodels
- R packages: phyloseq, bindata, MIDAS, tibble, xtable, sva, limma, vegan, MMUPHin, FDboost, doParallel, dylyr, readr, mixOmics, parallel, ade4, compositions

## 1. Functionalities

This github repository stores the code for benchmarking microbiome batch correction methods, enabling readers to reproduce the analyses done in the paper as well as leveraging their own data along with the provided code to do similar analyses. Our batch correction method evaluation consists of three steps:

### 1.1 Data generation/collection

#### 1.1.1 Simulation data generation

A comprehensive evaluation requires stringent simulation and corresponding analyses. Here, we employs [MIDAS](https://pubmed.ncbi.nlm.nih.gov/36993431/), an intuitive microbiome data simulator that recapitulates microbiome data structure without distributional assumptions while allowing users to manipulate variables of interest, in our case batch effect and conditional effect. We further edited codes to enable the simulation to include varying odds ratio (between biological signals and batch in each taxa), existing and non-existing relationship between batch effect and library size, as well as the generation of both count and relative abundance data.

To reproduce one iteration of our data using MIDAS, one can do the following in command line

```
cd benchmark # enter the benchmark directory
R # to initiate the R environment
source('./generate_data_MIDAS.R') # load generation file
# the odds ratio/conditional (biological) effect/batch effect you wish to generate the data for
or_l = c(1, 1.25, 1.5)
cond_effect_val_l = c(0, 0.25, 0.5, 0.75, 1)
batch_effect_val_l = c(0, 0.25, 0.5, 0.75, 1)
output_root = './trial' # please make sure the output_root exists
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

#### 1.2 Real world microbiome data collection and cleaning

A comprehensive evaluation also requires well-collected and well-preprocessed real-world microbiome data to ensure that the theoretically excelling batch correction methods from simulation experiments work in practice to enabled biomarker and trend discovery in batch corrected datasets while removing batch effect.

To this end, we collect and clean both count and relative abundance real world datasets from two databases, [MicrobiomeHD](https://zenodo.org/records/569601) and [CuratedMetagenomicsData](https://waldronlab.io/curatedMetagenomicData/) (CMD). The real world data preparation can be done by running the following code in the command line:

```
cd data # enter the data directory
# pre-clean the two microbiomeHD datasets, moving them from different batches to the same files
python3 step0_microbiomeHD_precleaning.py 
# (1) loads the four datasets into phyloseq objects (2) clean and prune samples and taxa below 0.05% abundance (3) save to a folder
Rscript step1_clean_and_prune.R 
# loads the phyloseq formatted data from step1 to save to a directory in a format easily digestible by later steps of the pipeline, along with generating a complete confounding check table
python3 step2_preprocessing_summarystats.py
```

These three commands result in the four subfolders in `data/cleaned_data`, each containing three files: the count (relative abundance) OTU data, the corresponding metadata, and the complete confounding checking table (number of samples in each biological condition in each cleaned dataset).

Note that while the crc dataset started with 8 batches, after pre-cleaning, there are actually 5 datasets (batches) that participate in batch correction.

(1) Code harmonicMic

 - add beta diversity to objective function # DONE -> to be improved, add local
 - add alpha diversity to objective function  # DONE -> to be improved, add local
 - account for the fact that the counts/data might be unevenly distributed between different batches and covariate # DONE
 - added an additional weight to the Y-update to make dominant features be affected more # DONE
     In HarmonicMic, Yk is not only influenced by the updated cluster membership Rki, but also by the ratio of a respective feature’s count sum to a reference value
 - revert back to count state # DONE

(2) evaluate pipeline

SINGLE METHOD EVALUATION:

 - Visualizing (PCA) batch-corrected dataset in terms of batches # DONE
 - Visualizing (PCA) batch-corrected dataset in terms of biological variables of interest # DONE
 - Evaluating whether the biological significance is retained in the top 2 principal components # DONE
 - Visualizing alpha diversity (Shannon) in terms of batches # DONE
 - Visualizing alpha diversity (Shannon) in terms of biological variables of interest # DONE
 - Visualizing beta diversity (bray-curtis) in terms of batches # DONE
 - Visualizing beta diversity (bray-curtis) in terms of biological variables of interest # DONE

MULTI-METHOD BENCHMARKING:

 - Visualizing distance between batches (bray-curtis and Aitchson) before bc and after using different methods # TODO
    - reference: https://www.nature.com/articles/s41467-022-33071-9#:~:text=Numerically%2C%20although%20ConQuR%20did%20not,from%205.66%25%20to%200.10%25.
- Visualizing/tablizing p-values among batches before bc and after using different methods # TODO
- Visualizing/tablizing p-values for alpha diversity among batches before bc and after using different methods # TODO
- Visualizing/tablizing p-values for beta diversity among batches before bc and after using different methods # TODO
- OMNIBUS testing for each meta analysis real datasets
    - Figure 3 of https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02753-4
- Identify biomarkers/consensus biomarkers from each meta analysis real datasets to demonstrate the power of harmonicMic or one of the methods
- Benchmark runtime 
    - kind of done for the benchmarked methods r script, to be added for Percentile-Norm and harmony-related methods

(3) Benchmarked methods for running
- R script for # DONE-ish
    - MMUPHin
    - combat
    - limma
    - ConQuR
- Python script for Percentile-Normalization # DONE

(3) Simulation
- Generate simulation experiments with 2, 3, 4, 5, 8, 10 batches using SparseDOSSA (ongoing)
- Run batch correction methods (1) and (3) on them and benchmark with (2)

(4) Real datasets
- code to preprocess microbiomeHD data format # DONE
- code to preprocess curatedMetagenomics data format # DONE
- fetch 2/3/4/5/6-sized metadatasets # ongoing
- run on all methods  # ongoing
    - individual method benchmark
    - cross-method benchmark

(5) biological insights/discussion
 - biological variable info retained in real world data? # todo
 - biomarker identification # todo


## 1. Usage

### 1a. HarmonicMic usage

### 1b. Data benchmarking

#### 1b.i Simulation data

#### 1b.ii Real data

#### 1b.iii Methods benchmarked

#### 1b.iv PCA-based statistical testing and visualization for batch effect correction

#### 1b.v Alpha-/beta- diversity based statisical testing and visualization for batch effect correction

#### 1b.vi Evaluating whether biological significance is retained

## 2. Replicating the results

Codes for replicating the results can be found in the results directory.


## 3. Additional notes

The name of this method, harmonicMic, is a combination of [harmony](https://www.nature.com/articles/s41592-019-0619-0), the tried-and-true single-cell batch correction method I revise from, and [microbiome data](https://www.niehs.nih.gov/health/topics/science/microbiome/index.cfm); as well as inspired by the time I spent at the [One Night Stanza a cappella group](https://www.instagram.com/stanza.gram/?hl=en).

Cheers,\
Tom Fu


To save the environment to an yml file, we can do
```
conda env export -p /home/fuc/anaconda3/envs/harmonicMic_env > environment.yml
```
