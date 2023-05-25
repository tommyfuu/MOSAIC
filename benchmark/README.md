harmonypy
=========

[![Latest PyPI Version][pb]][pypi] [![PyPI Downloads][db]][pypi] [![DOI](https://zenodo.org/badge/229105533.svg)](https://zenodo.org/badge/latestdoi/229105533)

[pb]: https://img.shields.io/pypi/v/harmonypy.svg
[pypi]: https://pypi.org/project/harmonypy/

[db]: https://img.shields.io/pypi/dm/harmonypy?label=pypi%20downloads

Harmony is an algorithm for integrating multiple high-dimensional datasets.

harmonypy is a port of the [harmony] R package by [Ilya Korsunsky].

Example
-------

<p align="center">
  <img src="https://i.imgur.com/lqReopf.gif">
</p>

This animation shows the Harmony alignment of three single-cell RNA-seq datasets from different donors.

[→ How to make this animation.](https://slowkow.com/notes/harmony-animation/)

Installation
------------

This package has been tested with Python 3.7.

Use [pip] to install:

```bash
pip install harmonypy
```

Usage
-----

Here is a brief example using the data that comes with the R package:

```python
# Load data
import pandas as pd

meta_data = pd.read_csv("data/meta.tsv.gz", sep = "\t")
vars_use = ['dataset']

# meta_data
#
#                  cell_id dataset  nGene  percent_mito cell_type
# 0    half_TGAAATTGGTCTAG    half   3664      0.017722    jurkat
# 1    half_GCGATATGCTGATG    half   3858      0.029228      t293
# 2    half_ATTTCTCTCACTAG    half   4049      0.015966    jurkat
# 3    half_CGTAACGACGAGAG    half   3443      0.020379    jurkat
# 4    half_ACGCCTTGTTTACC    half   2813      0.024774      t293
# ..                   ...     ...    ...           ...       ...
# 295  t293_TTACGTACGACACT    t293   4152      0.033997      t293
# 296  t293_TAGAATTGTTGGTG    t293   3097      0.021769      t293
# 297  t293_CGGATAACACCACA    t293   3157      0.020411      t293
# 298  t293_GGTACTGAGTCGAT    t293   2685      0.027846      t293
# 299  t293_ACGCTGCTTCTTAC    t293   3513      0.021240      t293

data_mat = pd.read_csv("data/pcs.tsv.gz", sep = "\t")
data_mat = np.array(data_mat)

# data_mat[:5,:5]
#
# array([[ 0.0071695 , -0.00552724, -0.0036281 , -0.00798025,  0.00028931],
#        [-0.011333  ,  0.00022233, -0.00073589, -0.00192452,  0.0032624 ],
#        [ 0.0091214 , -0.00940727, -0.00106816, -0.0042749 , -0.00029096],
#        [ 0.00866286, -0.00514987, -0.0008989 , -0.00821785, -0.00126997],
#        [-0.00953977,  0.00222714, -0.00374373, -0.00028554,  0.00063737]])

# meta_data.shape # 300 cells, 5 variables
# (300, 5)
#
# data_mat.shape  # 300 cells, 20 PCs
# (300, 20)

# Run Harmony
import harmonypy as hm
ho = hm.run_harmony(data_mat, meta_data, vars_use)

# Write the adjusted PCs to a new file.
res = pd.DataFrame(ho.Z_corr)
res.columns = ['X{}'.format(i + 1) for i in range(res.shape[1])]
res.to_csv("data/adj.tsv.gz", sep = "\t", index = False)
```

[harmony]: https://github.com/immunogenomics/harmony
[Ilya Korsunsky]: https://github.com/ilyakorsunsky
[pip]: https://pip.readthedocs.io/

import numpy as np
import pandas as pd
data_mat = pd.read_csv("/home/fuc/HRZE_TB/tom_organized_codes/batch_correction_PCA/1021_microbiome_batchcorrection/microbiome_merged_intersect_1023.csv", index_col="Unnamed: 0")
data_mat = np.array(data_mat)
meta_data = pd.read_csv("/home/fuc/HRZE_TB/tom_organized_codes/batch_correction_PCA/1021_microbiome_batchcorrection/intersect_metadata_1023.csv")
vars_use = ["Dataset", "Sex"]
from harmony import run_harmony
ho = run_harmony(data_mat, meta_data, vars_use)


# the gist of new ideas

1. in preprocessing, get rid of features and samples that does not satisfy at least one of the following conditions
  - features with all zeroes
  - features where the sum of counts are below 0.01% compared to the total sum of all counts (include those features with all zeroes)
  - samples with all zeroes

2. instead of using PCA, we maintain the structure of the data and use the original representation of the data instead and do harmony on that

In the objective function, multiply the maximum diversity clustering term with a coefficient mu which represents the ratio of that sample/feature's count value to the whole. 
  - This might come out to be an extremely small value or even zero. So we standardize the ratio by the maximum count for any feature in that sample
     which is to say, the more abundant the feature, the more we want to maximize the diversity of that feature in the particular sample

3. In the objective function, multiply everything by an indicator function which gives zero for a feature/cell when the corresponding feature is one of
  - features where the sum of counts are below 0.01% compared to the total sum of all counts in at least one of all batches (include all zeroes)
  - features where the sum of counts are below 0.01% compared to the total sum of all counts in at least one of all covariates (include all zeroes)
