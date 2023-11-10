import os
import pandas as pd
from step0_microbiomeHD_precleaning import check_complete_confounding, preprocess

def preprocess_data_phyloseq(address_directory, output_root = False, id = 'Sam_id', covar_l = []):
    ### CuratedMetagenomicsDataset provides way more metadata in a congestible manner
    cur_dir_names = os.listdir(address_directory)
    address_X = address_directory + '/'+ [result for result in cur_dir_names if "otu_table_" in result][0]
    address_Y = address_directory + '/'+ [result for result in cur_dir_names if "sample_table_" in result][0]
    data_mat = pd.read_csv(address_X, index_col=0)
    meta_data = pd.read_csv(address_Y, index_col=0)
    meta_data[id] = list(meta_data.index)
    data_mat, meta_data = preprocess(data_mat.T, meta_data, id, covar_l)
    
    # TODO: ensure that the sample ids are correctly aligned in metadata and count_table
    data_mat_ids = list(data_mat.index)
    meta_data_ids = list(meta_data.index)
    intersection_ids = list(set(meta_data_ids).intersection(data_mat_ids))

    # drop rows where indexes are not overlapping
    data_mat_non_intersecting = [id for id in data_mat_ids if id not in intersection_ids]
    data_mat = data_mat.drop(data_mat_non_intersecting)
    meta_data_non_intersecting = [id for id in meta_data_ids if id not in intersection_ids]
    meta_data = meta_data.drop(meta_data_non_intersecting)
    data_mat = data_mat.reindex(intersection_ids)
    meta_data = meta_data.reindex(intersection_ids)

    # drop taxa where all values are 0
    data_mat = data_mat.loc[:, (data_mat != 0).any(axis=0)]
    # data_mat = data_mat.loc[(data_mat != 0).any(axis=1), :]

    # convert data_mat to relative abundance for each sample
    data_mat = data_mat.div(data_mat.sum(axis=1), axis=0) # divide by row sum

    # save stuff if needed
    if output_root != False:
        data_mat.to_csv(output_root+"_count_data.csv") # this is actually not count, but relative abundance
        meta_data.to_csv(output_root+"_meta_data.csv", index=False)
    return data_mat, meta_data

def load_results_from_benchmarked_methods(address_X, address_Y):
    data_mat = pd.read_csv(address_X, index_col=0)
    meta_data = pd.read_csv(address_Y)

    # check for missing values and fill out with mean
    if data_mat.isnull().values.any():
        print("fill out nans")
        data_mat = data_mat.fillna(data_mat.mean())
        # remove columns (taxa) filled with NAs
        data_mat = data_mat.dropna(axis=1, how='all')

    return data_mat, meta_data

overall_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data'
# autism_2_microbiomeHD
data_mat, meta_data = preprocess_data_phyloseq(f'{overall_path}/pruned_autism_2_microbiomeHD', f'{overall_path}/cleaned_data/autism_2_microbiomeHD/autism_2_microbiomeHD', id = 'Sam_id', covar_l = [])
data_mat, meta_data = load_results_from_benchmarked_methods(f'{overall_path}/cleaned_data/autism_2_microbiomeHD/autism_2_microbiomeHD_count_data.csv', f'{overall_path}/cleaned_data/autism_2_microbiomeHD/autism_2_microbiomeHD_meta_data.csv')
check_complete_confounding(meta_data, 'Dataset', 'DiseaseState', f'{overall_path}/cleaned_data/autism_2_microbiomeHD/autism_2_microbiomeHD')

# cdi_3_microbiomeHD
data_mat, meta_data = preprocess_data_phyloseq(f'{overall_path}/pruned_cdi_3_microbiomeHD', f'{overall_path}/cleaned_data/cdi_3_microbiomeHD/cdi_3_microbiomeHD', id = 'Sam_id', covar_l = [])
data_mat, meta_data = load_results_from_benchmarked_methods(f'{overall_path}/cleaned_data/cdi_3_microbiomeHD/cdi_3_microbiomeHD_count_data.csv', f'{overall_path}/cleaned_data/cdi_3_microbiomeHD/cdi_3_microbiomeHD_meta_data.csv')
check_complete_confounding(meta_data, 'Dataset', 'DiseaseState', f'{overall_path}/cleaned_data/cdi_3_microbiomeHD/cdi_3_microbiomeHD')

# ibd_3_CMD
data_mat, meta_data = preprocess_data_phyloseq(f'{overall_path}/pruned_ibd_3_CMD', f'{overall_path}/cleaned_data/ibd_3_CMD/ibd_3_CMD', id = 'Sam_id', covar_l = [])
data_mat, meta_data = load_results_from_benchmarked_methods(f'{overall_path}/cleaned_data/ibd_3_CMD/ibd_3_CMD_count_data.csv', f'{overall_path}/cleaned_data/ibd_3_CMD/ibd_3_CMD_meta_data.csv')
check_complete_confounding(meta_data, "study_name", "disease", f'{overall_path}/cleaned_data/ibd_3_CMD/ibd_3_CMD')

# crc_8_CMD
data_mat, meta_data = preprocess_data_phyloseq(f'{overall_path}/pruned_crc_8_CMD', f'{overall_path}/cleaned_data/crc_8_CMD/crc_8_CMD', id = 'Sam_id', covar_l = [])
data_mat, meta_data = load_results_from_benchmarked_methods(f'{overall_path}/cleaned_data/crc_8_CMD/crc_8_CMD_count_data.csv', f'{overall_path}/cleaned_data/crc_8_CMD/crc_8_CMD_meta_data.csv')
check_complete_confounding(meta_data, "study_name", "disease", f'{overall_path}/cleaned_data/crc_8_CMD/crc_8_CMD')