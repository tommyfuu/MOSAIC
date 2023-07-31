#!/usr/bin/env python
"""
This script converts an input OTU table with cases and controls into
percentiles of control samples.
"""
__author__ = "Sean Gibbons and Claire Duvallet; edited by Tom Fu to fit our datatypes"
__copyright__ = "Copyright 2017"
__credits__ = ["Sean Gibbons; Claire Duvallet; Eric Alm"]
__reference__ = "PLoS Computational Biology DOI: https://doi.org/10.1371/journal.pcbi.1006102"
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Sean Gibbons"
__email__ = "sgibbons@isbscience.org"

import numpy as np
import scipy.stats as sp
import pandas as pd
import argparse
import time

# Passing through \n doesn't work...
seps = {'tab': '\t', 'newline': '\n', 'comma': ','}

## Read data
print('Loading data...')

def percentile_norm(input_file, meta_file, diseaseName, diseaseCase, delimiter = 'comma', output_root = '', simulate = False):
    if not simulate:
        df = pd.read_csv(input_file, sep=seps[delimiter], header=0, index_col=0)
    else:
        df = pd.read_csv(input_file, sep=seps[delimiter])

    # documenting time elapsed
    import time
    start = time.time()


    #replace zeros with random draw from uniform(0, 10**-9)
    df = df.replace(0.0,np.nan)
    df_rand = pd.DataFrame(np.random.uniform(0.0,10**-9,size=(df.shape[0],df.shape[1])),index=df.index,columns=df.columns)
    df[pd.isnull(df)] = df_rand[pd.isnull(df)]

    # Get numpy array
    x = df.values

    # get case info
    metadata = pd.read_csv(meta_file)
    sample_ids = list(df.index)
    disease_column = list(metadata[diseaseName])
    case_list = [sample for idx, sample in enumerate(sample_ids) if disease_column[idx] == diseaseCase]
    control_list = [sample for idx, sample in enumerate(sample_ids) if disease_column[idx] != diseaseCase]
    # # Read case and control samples as lists
    # with open(diseaseCase, 'r') as f:
    #     case_list = f.read().rstrip().split(seps[args.sample_d])
    # with open(args.control, 'r') as f:
    #     control_list = f.read().rstrip().split(seps[args.sample_d])
    # print('Loading data complete.')



    # Get control and case indices
    control_indices = [df.index.get_loc(i) for i in control_list]
    case_indices = [df.index.get_loc(i) for i in case_list]

    print("control indices")
    print(control_indices)
    print("case indices")
    print(case_indices)
    all_samples = control_list + case_list
    all_indices = control_indices + case_indices

    ## Normalize control and case samples to percentiles of control distribution
    print('Running percentile-normalization...')
    import time

    start = time.time()

    norm_x = np.array(
        [
            [sp.percentileofscore(x[control_indices, i], x[j, i], kind='mean')
                for j in all_indices]
        for i in range(x.shape[1])
        ]).T
    print('Running percentile-normalization complete.')

    print("elapsed time:")
    end = time.time()
    print(end - start)

    # write the time elapsed to a separate file
    with open(output_root + "_percentile_norm_elapsed_time.txt", "w") as f:
        f.write(str(end - start))

    ## Put back into dataframe and write to file
    norm_df = pd.DataFrame(data=norm_x, columns=df.columns, index=all_samples)
    norm_df.to_csv(output_root+ "_percentile_norm.csv")

    print('Percentile-normalized data written to {}'.format(output_root))
    return

# overall_path = "/athena/linglab/scratch/chf4012/mic_bc_benchmark"
# percentile_norm(overall_path+"/benchmark/benchmarked_data/autism_2_microbiomeHD_count_data.csv", overall_path+"/benchmark/benchmarked_data/autism_2_microbiomeHD_meta_data.csv", "DiseaseState", "ASD", "comma", overall_path+"/benchmark/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD")
