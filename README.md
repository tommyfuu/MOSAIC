# MOSAIC <img src="https://github.com/tommyfuu/MOSAIC/blob/main/logo/MOSAIC_logo.jpeg" width="25" height="25"> : a pipeline for MicrobiOme Studies Analytical Integration and Correction
Repository for the manuscript entitled __MOSAIC: MicrobiOme Studies Analytical Integration and Correction pipeline__ by _Fu et al._. 

<p align="center">
  <img src="https://github.com/tommyfuu/MOSAIC/blob/main/logo/MOSAIC_logo.jpeg" width="500" height="500">
</p>



(Logo generated with [Google Gemini](https://gemini.google.com/app))



_Developed and maintained by_:  [Chenlian (Tom) Fu](fcl200089@outlook.com)\
_Advised by_: Dr. Wodan Ling\
_Affiliations_: Weill Cornell Medicine, Memorial Sloan Kettering Cancer Center


This github repository stores the code for the __MOSAIC__ data integration pipeline as described in the manuscript. __MOSAIC__ is a three-part [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline that takes in and preprocesses input microbiome data from multiple sources (or generating the data in the case of simulation), apply multiple microbiome data integration methods of users' choosing, and generates evaluative reports based on four types of metrics important for microbiome analyses. This repository also stores some additional codes for reproducing the figures as presented in the manuscript.

The usage of the pipeline is quite simple and user-friendly:

## 1. Environment set-up

To create a clean environment with all the necessary libraries, you can manually install the following packages in a virtual environment (conda recommended) with Python 3.8(.17) and R 4.3:

Firstly, please install all conda-installable packages via our environment yaml file in the command line with:

```
conda env create -f env.yml
conda activate mic_data_integration
```

After this, open `R` in command line, and then indivudally install the three packages MIDASim, FDboost, and mixOmics from the respective sites linked above. This is necessary as these packages are not yet in the conda ecosystem. Instructions for installing these packages are as follows:

- [MIDASim](https://github.com/mengyu-he/MIDASim) 
- [FDboost](https://github.com/boost-R/FDboost)
- [mixOmics](https://bioconductor.org/packages/release/bioc/html/mixOmics.html)


## 2. Running the pipeline on real-world datasets

__MOSAIC__ has built-in support for microbiome datasets coming from two highly used databases, namely [CuratedMetagenomicData](https://waldronlab.io/curatedMetagenomicData/) (referred to as CMD from now on) and [MicrobiomeHD](https://zenodo.org/records/1146764). To use any of these two databases, please refer to `config.yml` for a MicrobiomeHD example or `./additional_ymls/config_IBD.yml` for a CMD example. To use data from any other databases or custom-made datasets, please have a folder in which there are subfolders of datasets to be integrated; in each subfolder, there should be a (taxon rows * sample columns with taxons as index and sample names as column names) file that ends with `count_data.csv` and another ending with `meta_data.csv` (where each row is a sample with their meta data), serving as OTU_table and meta data respectively. The user can then revise the `config.yml` accordingly, making sure that the `dataset_name` item does not have `'microbiomeHD'` nor `'CMD'` in it.

Brief explanations on the entries in the yml file:
- `rules_to_run`: currently supports 'all_rw', 'preprocess', 'integrate', 'evaluation' for real-world datasets. While the rest are pretty self-apparent in meanings, 'all_rw' simply means to run the whole pipeline in the order of preprocessing, data integration, and evaluation.
- `data_source`: please specify if you are using 'microbiomeHD', 'CMD', or 'custom' for real-world datasets.
- `src`: path where the source real-world data is stored for MicrobiomeHD and custom data, as well as where all the preprocessed data will be stored after that step under the `cleaned_data` subfolder of `src`.
- `dataset_name`: the name of your dataset. Please end your dataset name with the database source your data is from: 'microbiomeHD', 'CMD', or 'custom'.
- `datatype`: currently supports microbiome data in 'count' or 'relab' forms.
- `RW_LIBSIZE_THRESHOLD`: the threshold for filtering out samples with a less than this number * 100 th quantile of library size in preprocessing.
- `RW_RELAB_THRESHOLD`: the threshold for filtering out taxa with a less than this number * 100 th quantile of relative abundance in preprocessing.
- `COVAR_L`: if not empty, this list of covariates will be checked against the metadata to remove samples that have any NAs for each of the included covariates.
- `post_integration_outputs`: the path of directory where the data after data integration/batch correction method is applied is saved.
- `used_R_methods`: the list of implemented R methods for data integration; all options for count data include `["combat_seq", "limma", "MMUPHin", 'ConQuR', 'ConQuR_libsize']`, all options for relab data include `["combat", "limma", "MMUPHin", 'ConQuR_rel']`.
- `used_Python_methods`: : the list of implemented Python methods for data integration; all optionsfor both count and relab data include `['harmony', "percentile_norm"]`.
- `batch_ref`: reference batch name in data integration necessary for some methods.
- `covar`: a list of covariates to be considered for biological signal retainment in data integration. The first element of the list should be the condition variable of interest, such as `'disease'`.
- `condition_value`: positive value for the condition of interest.
- `evaluation_outputs`: the path of directory where the evaluation outputs and summary statistics will be saved.
- `binarizing_agent`: for the given condition variable, such as `'disease'`, please specify the condition value based on which all the samples will be binarized into patients with this `disease==binarizing_agent` vs samples that do not have this agent condition for all evaluations.

After this yml file is configured, please make sure the line `configfile: "./config.yml"` in the `Snakefile` is changed into the path of your custom configuration yml file. Then you can run the pipeline with Snakemake in command line:

```
# make sure you are at the path where the Snakefile is located
conda activate mic_data_integration
snakemake -c1 --latency-wait 5
```

If you wish to run each of the three steps separately, you may run:

```
snakemake -R preprocess -c1 --latency-wait 5 # preprocess can be replaced with integrate or evaluation
```


#### Acknowledgements
This project started as a class project in Dr. Wesley Tansey's class Foundation of Data Science at Memorial Sloan Kettering and Weill Cornell, and would not have been possible without the generous support of Dr. Wodan Ling, Jiuyao Lu, Dr. Quaid Morris, and Dr. Wesley Tansey. The funding support for my PhD comes from the Tri-institutional Program of Computational Biology and Medicine.

If you have questions regarding this benchmarking project or other inquries, please reach out to me at ___chf4012@med.cornell.edu___. I hope you have a great day!

Cheers,\
Tom Fu

