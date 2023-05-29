import pandas as pd
import numpy as np
import os

def stratify_data(otu_path, sample_table_path, var_to_stratify, interval, out_dir, out_postfix):
    """
    Stratify data based on a variable in the sample table. The variable must be
    numeric. The data is stratified into intervals of size interval. The
    stratified data is saved in a new file.

    Parameters
    ----------
    otu_path : str
        Path to the OTU table.
    sample_table_path : str
        Path to the sample table.
    var_to_stratify : str
        Name of the variable to stratify the data by.
    interval : float
        Size of the interval to stratify the data by.

    Returns
    -------
    None
    """
    # Load the data
    otu = pd.read_csv(otu_path, index_col=0)
    sample_table = pd.read_csv(sample_table_path, index_col=0)
    print(sample_table)
    # Check that the variable to stratify by is numeric
    if not np.issubdtype(sample_table[var_to_stratify].dtype, np.number):
        raise ValueError('The variable to stratify by must be numeric.')

    # Create a new column in the sample table that contains the interval that
    # each sample belongs to
    sample_table['interval'] = np.floor(sample_table[var_to_stratify] / interval)

    # Group the samples by interval
    grouped = sample_table.groupby('interval')

    # Create a new OTU table for each interval
    for interval, group in grouped:
        # Get the samples in the interval
        samples = group.index

        # Get the OTUs in the interval
        otus = otu[samples]

        if not os.path.exists(out_dir+out_postfix+'_interval_'+ str(interval)):
            os.makedirs(out_dir+out_postfix+'_interval_'+ str(interval))

        # Save the OTU table
        otus.to_csv(out_dir+out_postfix+'_interval_'+ str(interval)+'/otu_table+'+out_postfix+'_interval_' + str(interval) + '.tsv', sep='\t')

        # Save the sample table
        group.to_csv(out_dir+out_postfix+'_interval_'+ str(interval)+'/sample_table_'+out_postfix+'_interval_' + str(interval) + '.tsv', sep='\t')

    return


otu_path = '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/data/ibdmdb/otu_table_HMP_2019_ibdmdb.csv'
sample_table_path = '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/data/ibdmdb/sample_table_HMP_2019_ibdmdb.csv'
var_to_stratify = 'days_from_first_collection'
out_dir = '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/data/'
out_postfix = 'ibdmdb'
stratify_data(otu_path, sample_table_path, var_to_stratify, 30, out_dir, out_postfix)