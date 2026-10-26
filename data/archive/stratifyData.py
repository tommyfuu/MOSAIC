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
        otus.to_csv(out_dir+out_postfix+'_interval_'+ str(interval)+'/otu_table_'+out_postfix+'_interval_' + str(interval) + '.csv')

        # Save the sample table
        group.to_csv(out_dir+out_postfix+'_interval_'+ str(interval)+'/sample_table_'+out_postfix+'_interval_' + str(interval) + '.csv')

    return


# otu_path = '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/data/ibdmdb/otu_table_HMP_2019_ibdmdb.csv'
# sample_table_path = '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/data/ibdmdb/sample_table_HMP_2019_ibdmdb.csv'
# var_to_stratify = 'days_from_first_collection'
# out_dir = '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/data/'
# out_postfix = 'ibdmdb'
# stratify_data(otu_path, sample_table_path, var_to_stratify, 30, out_dir, out_postfix)

def remove_certain_var_values(otu_path, sample_table_path, var_to_remove, var_values_to_remove, out_dir, out_postfix):
    """
    Remove certain values of a variable from the sample table. The data is saved
    in a new file.

    Parameters
    ----------
    otu_path : str
        Path to the OTU table.
    sample_table_path : str
        Path to the sample table.
    var_to_remove : str
        Name of the variable to remove values from.
    var_values_to_remove : list
        List of values to remove from the variable.
    out_dir : str
        Path to the directory to save the new data in.
    out_postfix : str
        Postfix to add to the file names of the new data.

    Returns
    -------
    None
    """
    # Load the data
    otu = pd.read_csv(otu_path, index_col=0)
    sample_table = pd.read_csv(sample_table_path, index_col=0)

    # Remove the samples with the values to remove
    sample_table = sample_table[~sample_table[var_to_remove].isin(var_values_to_remove)]

    # Get the samples that remain
    samples = sample_table.index

    # Get the OTUs that remain
    otus = otu[samples]

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Save the OTU table
    otus.to_csv(out_dir+'otu_table_'+out_postfix+'.csv')

    # Save the sample table
    sample_table.to_csv(out_dir+'sample_table_'+out_postfix+'.csv')

    return

otu_path = '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/data/hanninganGD/otu_table_HanniganGD_2017.csv'
sample_table_path = '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/data/hanninganGD/sample_table_HanniganGD_2017.csv'
var_to_stratify = 'days_from_first_collection'
out_dir = '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/data/hanninganGD_noBoston/'
out_postfix = 'HanniganGD_2017_filtered'
remove_certain_var_values(otu_path, sample_table_path, 'location', ['Boston'], out_dir, out_postfix)