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

## Input arguments
parser = argparse.ArgumentParser(description='Script to convert case control '
    + 'OTU tables into percentiles of control samples.')
parser.add_argument('-i', help='input OTU table text file (rows = samples, '
    + ' columns = OTUs; default format = tab-delimited)', required=True)
# parser.add_argument('-case', help='input case sample list', required=True)
# parser.add_argument('-control', help='input control sample list', required=True)
parser.add_argument('-meta', help='input path to metadata file, containing a case-vs-control column', required=True)
parser.add_argument('-diseaseName', help='name of the column representing the batch in the metadata', required=True)
parser.add_argument('-case', help='name of the case (disease) scenario ', required=True)
parser.add_argument('-otu-d', help='OTU table field delimiter [default: '
    + '%(default)s]', default='comma', choices=['tab', 'newline', 'comma'])
parser.add_argument('-sample-d', help='sample list delimiters [default: '
    + '%(default)s]', default='tab', choices=['tab', 'newline', 'comma'])
parser.add_argument('-o', help='output file name [default: %(default)s]',
    default='out_percentile_norm.txt')
args = parser.parse_args()

# Passing through \n doesn't work...
seps = {'tab': '\t', 'newline': '\n', 'comma': ','}

## Read data
print('Loading data...')
df = pd.read_csv(args.i, sep=seps[args.otu_d], header=0, index_col=0)

#replace zeros with random draw from uniform(0, 10**-9)
df = df.replace(0.0,np.nan)
df_rand = pd.DataFrame(np.random.uniform(0.0,10**-9,size=(df.shape[0],df.shape[1])),index=df.index,columns=df.columns)
df[pd.isnull(df)] = df_rand[pd.isnull(df)]

# Get numpy array
x = df.values

# get case info
metadata = pd.read_csv(args.meta)
sample_ids = list(df.index)
disease_column = list(metadata[args.diseaseName])
case_list = [sample for idx, sample in enumerate(sample_ids) if disease_column[idx] == args.case]
control_list = [sample for idx, sample in enumerate(sample_ids) if disease_column[idx] != args.case]
# # Read case and control samples as lists
# with open(args.case, 'r') as f:
#     case_list = f.read().rstrip().split(seps[args.sample_d])
# with open(args.control, 'r') as f:
#     control_list = f.read().rstrip().split(seps[args.sample_d])
# print('Loading data complete.')



# Get control and case indices
control_indices = [df.index.get_loc(i) for i in control_list]
case_indices = [df.index.get_loc(i) for i in case_list]

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

print("elapsed time:", args.o)
end = time.time()
print(end - start)

## Put back into dataframe and write to file
norm_df = pd.DataFrame(data=norm_x, columns=df.columns, index=all_samples)
norm_df.to_csv(args.o)

print('Percentile-normalized data written to {}'.format(args.o))

# python percentile_norm.py -i /home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/Glickman_count_data.csv -meta /home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/Glickman_meta_data.csv -diseaseName Visit -case Day_0 -o /home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/Glickman/Glickman_percentile_norm.csv
# python percentile_norm.py -i /home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/autism_2_microbiomeHD_count_data.csv -meta /home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/autism_2_microbiomeHD_meta_data.csv -diseaseName DiseaseState -case ASD -o /home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD_percentile_norm.csv
# python percentile_norm.py -i /home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/cdi_3_microbiomeHD_count_data.csv -meta /home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/cdi_3_microbiomeHD_meta_data.csv -diseaseName DiseaseState -case CDI -o /home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD_percentile_norm.csv
# python percentile_norm.py -i /home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/adenoma_5_CMD_count_data.csv -meta /home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/adenoma_5_CMD_meta_data.csv -diseaseName disease -case adenoma -o /home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/adenoma_5_CMD/adenoma_5_CMD_percentile_norm.csv
# python percentile_norm.py -i /home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/ibd_3_CMD_count_data.csv -meta /home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/ibd_3_CMD_meta_data.csv -diseaseName disease -case IBD -o /home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/ibd_3_CMD/ibd_3_CMD_percentile_norm.csv
# python percentile_norm.py -i /home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/CRC_8_CMD_count_data.csv -meta /home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/CRC_8_CMD_meta_data.csv -diseaseName disease -case CRC -o /home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/CRC_8_CMD/CRC_8_CMD_percentile_norm.csv
# python percentile_norm.py -i /home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/T2D_10_CMD_count_data.csv -meta /home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/T2D_10_CMD_meta_data.csv -diseaseName disease -case T2D -o /home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/T2D_10_CMD/T2D_10_CMD_percentile_norm.csv

# python percentile_norm.py -i /Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/ibdmdb_interval_0.0_count_data.csv -meta /Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/ibdmdb_interval_0.0_meta_data.csv -diseaseName disease -case IBD -o /Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_results/ibdmdb_interval_0.0/ibdmdb_interval_0.0_percentile_norm.csv
# python percentile_norm.py -i /Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/hanninganGD_noBoston_count_data.csv -meta /Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/hanninganGD_noBoston_meta_data.csv -diseaseName disease -case adenoma -o /Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_results/hanninganGD_noBoston/hanninganGD_noBoston_percentile_norm.csv

# bin_corr_val_l = [0.3]
# cond_effect_val_l = [0, 0.099, 0.899]
# batch_effect_val_l = [0, 0.099, 0.899]
# num_iters = 1
# IDCol = 'subjectid_text'
# methods_list = ["nobc", "combat_seq", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "harmony", "Percentile_norm"]
# for bin_corr_val in bin_corr_val_l:
#     for cond_effect_val in cond_effect_val_l:
#         for batch_effect_val in batch_effect_val_l:
#             for iter in list(range(1, num_iters+1)):
# python percentile_norm.py -i /Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/hanninganGD_noBoston_count_data.csv -meta /Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/hanninganGD_noBoston_meta_data.csv -diseaseName disease -case adenoma -o /Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_results/hanninganGD_noBoston/hanninganGD_noBoston_percentile_norm.csv
# python percentile_norm.py -i /Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/hanninganGD_noBoston_count_data.csv -meta /Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/hanninganGD_noBoston_meta_data.csv -diseaseName disease -case adenoma -o /Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_results/hanninganGD_noBoston/hanninganGD_noBoston_percentile_norm.csv

