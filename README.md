# harmonicMic🎤
Harmony batch correction and meta-analysis method revised to be used for microbiome data.

_Developed and maintained by_:  [Chenlian (Tom) Fu](tfu@g.hmc.edu)\
_Advised by_: Dr. Wodan Ling, Dr. Wesley Tansey, Dr. Quaid Morris\
_Affiliations_: Weill Cornell Medicine, Memorial Sloan Kettering Cancer Center

## 0. TODO

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
