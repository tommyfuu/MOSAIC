# harmonicMic🎤
Harmony batch correction and meta-analysis method revised to be used for microbiome data.

_Developed and maintained by_:  [Chenlian (Tom) Fu](tfu@g.hmc.edu)\
_Advised by_: Dr. Wesley Tansey, Dr. Quaid Morris\
_Collaborators (TB data source)_: Dr. Michael Glickman\
_Affiliations_: Memorial Sloan Kettering Cancer Center, Weill Cornell Medicine

## 0. TODO

(1) Write pipeline to combine data in microbiomeHD format to the curatedMetagenomicData format
 - find overlapping species in the RDP/datasetID.otu_table.100.denovo.rdp_assigned files across all data, concatenate
 - combine metadata and label useful information (Dataset, Sex, + biological information labelled)
 - preprocessing (standard scaler normalization)
 - pre-batch-correction evaluation
    - PCA visualization, with the most important biological information colored, and batches styled # DONE
    - statistically evaluate batch effect with the kw # DONE
    - Evaluate sex bias and bias induced by other covariates # DONE
    - Evaluate alpha and beta diversity

(2) Code harmonicMicPy


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
