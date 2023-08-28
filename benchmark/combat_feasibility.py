import pandas as pd
import numpy as np


def combat_feasibility(counts_path, meta_path, batch_var='study_name', sam_id = 'Sam_id'):
    if sam_id == 'Sam_id':
        counts_df = pd.read_csv(counts_path, index_col=0)
    else:
        counts_df = pd.read_csv(counts_path)
    print(counts_df.shape)
    print("____________")
    # get unique batches
    meta = pd.read_csv(meta_path)
    batches = meta[batch_var].unique()

    # for each unique batch, find columns (taxa) with 0 variance
    for batch in batches:
        print(batch)
        # get counts for batch
        batch_meta = meta[meta[batch_var] == batch]
        if sam_id == 'subjectid_text':
            batch_indices = [int(x.split('_')[-1])-1 for x in batch_meta[sam_id]]
            batch_counts = counts_df.iloc[batch_indices, :]
        else:
            batch_counts = counts_df.loc[list(batch_meta[sam_id]), :]
        # get taxa with 0 variance
        zero_var_taxa = batch_counts.columns[batch_counts.var() == 0]
        print(zero_var_taxa.shape)

        # find samples (rows) with 0 variance
        zero_var_samples = batch_counts.index[batch_counts.var(axis=1) == 0]
        print(zero_var_samples.shape)


if __name__ == '__main__':
    print("ibd, filled with NAs in combat output")
    meta_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark/benchmarked_data/ibd_3_CMD_meta_data.csv'
    counts_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark/benchmarked_data/ibd_3_CMD_count_data.csv'
    combat_feasibility(counts_path, meta_path)
    print("____________")
    print("CRC, filled with NAs in combat output")
    meta_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark/benchmarked_data/CRC_8_CMD_meta_data.csv'
    counts_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark/benchmarked_data/CRC_8_CMD_count_data.csv'
    combat_feasibility(counts_path, meta_path)
    print("____________")
    print("yesrelation, simulation_ibd_150_1_0_0.899_iter_2, filled with NAs in combat output")
    meta_path = '/athena/linglab/scratch/chf4012/simulation_data_MIDAS_small_yesrelation_080723/ibd_150_meta_1_0_0.899_iter_2.csv'
    counts_path = '/athena/linglab/scratch/chf4012/simulation_data_MIDAS_small_yesrelation_080723/ibd_150_relab_1_0_0.899_iter_2.csv'
    combat_feasibility(counts_path, meta_path, batch_var='batchid', sam_id='subjectid_text')
    print("____________")
    print("norelation, simulation_ibd_150_1_0_0.499_iter_4, filled with NAs in combat output")
    meta_path = '/athena/linglab/scratch/chf4012/simulation_data_MIDAS_small_norelation_080723/ibd_150_meta_1_0_0.499_iter_4.csv'
    counts_path = '/athena/linglab/scratch/chf4012/simulation_data_MIDAS_small_norelation_080723/ibd_150_relab_1_0_0.499_iter_4.csv'
    combat_feasibility(counts_path, meta_path, batch_var='batchid', sam_id='subjectid_text')
    print("____________")
    print("yesrelation, simulation_ibd_150_1_0_0.299_iter_2, no NAs in combat output")
    meta_path = '/athena/linglab/scratch/chf4012/simulation_data_MIDAS_small_yesrelation_080723/ibd_150_meta_1_0_0.299_iter_2.csv'
    counts_path = '/athena/linglab/scratch/chf4012/simulation_data_MIDAS_small_yesrelation_080723/ibd_150_relab_1_0_0.299_iter_2.csv'
    combat_feasibility(counts_path, meta_path, batch_var='batchid', sam_id='subjectid_text')
    print("____________")
    print("norelation, simulation_ibd_150_1_0_0.099_iter_5,no NAs in combat output")
    meta_path = '/athena/linglab/scratch/chf4012/simulation_data_MIDAS_small_norelation_080723/ibd_150_meta_1_0_0.099_iter_5.csv'
    counts_path = '/athena/linglab/scratch/chf4012/simulation_data_MIDAS_small_norelation_080723/ibd_150_relab_1_0_0.099_iter_5.csv'
    combat_feasibility(counts_path, meta_path, batch_var='batchid', sam_id='subjectid_text')