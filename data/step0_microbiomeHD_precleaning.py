import pandas as pd
import numpy as np
import os
import sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import REPO_ROOT, SIM_OUTPUT_ROOT, BENCH_DIR, DATA_DIR

def preprocess(data_mat, meta_data, IDCol, covar_l = []):
    # remove samples with all zeros
    data_mat = data_mat.loc[~(data_mat==0).all(axis=1)]
    kept_samples = data_mat.index
    meta_data = meta_data[meta_data[IDCol].isin(kept_samples)]
    # remove features with all zeros
    col_names = list(data_mat.columns)
    col_sums = data_mat.sum(axis = 1)
    removable_feature_names = [col_names[index] for index, col_sum in enumerate(col_sums) if col_sum==0]
    data_mat.drop(removable_feature_names, axis=1, inplace=True)
    # for each covar used, remove samples with missing values
    for covar in covar_l:
        meta_data = meta_data[meta_data[covar].notna()]
    data_mat = data_mat.loc[meta_data[IDCol]]
    # after cleaning samples, remove features with all zeros again
    col_names = list(data_mat.columns)
    col_sums = data_mat.sum(axis = 1)
    removable_feature_names = [col_names[index] for index, col_sum in enumerate(col_sums) if col_sum==0]
    data_mat.drop(removable_feature_names, axis=1, inplace=True)
    return data_mat, meta_data

def load_data_microbiomeHD(address_directory, output_root = False, id = 'Sam_id', covar_l = []):
    # MicrobiomeHD metadata differ across studies, so by default (covar_l=[]) only the batch
    # (Dataset) and DiseaseState are loaded.
    subdir_names = os.listdir(address_directory)
    subdir_names = [result for result in subdir_names if "results" in result]
    count_data_l = []
    intersect_taxa = []
    metadata_l = []
    for subdir in subdir_names:
        # 1. OTU table
        subdir_path = address_directory + '/' + subdir
        current_RDP_names = os.listdir(subdir_path + '/RDP')
        current_dbotu_count_data_path = [result for result in current_RDP_names if "100.denovo.rdp_assigned" in result][0]
        current_dbotu_count_data = pd.read_csv(subdir_path+'/RDP/'+current_dbotu_count_data_path, delimiter='\t', index_col=0)
        current_taxa = list(current_dbotu_count_data.index)
        # truncate taxonomy strings to genus level
        current_taxa = [";".join(taxa.split(';')[:-2]) for taxa in current_taxa]
        print(current_dbotu_count_data.shape)
        current_dbotu_count_data.index = current_taxa
        print(current_dbotu_count_data.shape)
        # sum OTUs assigned to the same genus
        current_dbotu_count_data = current_dbotu_count_data.groupby(level=0).sum()
        
        count_data_l.append(current_dbotu_count_data)
        
        if intersect_taxa == []:
            intersect_taxa = current_taxa
        else:
            intersect_taxa = list(set(intersect_taxa).intersection(current_taxa))

        # 2. metadata
        current_files_names = os.listdir(subdir_path)
        current_metadata_path = subdir_path + '/' + [result for result in current_files_names if "metadata" in result][0]
        print("metadata path", current_metadata_path)
        current_metadata = pd.read_csv(current_metadata_path, delimiter='\t', index_col=0, encoding='ISO-8859-1')['DiseaseState'].to_frame()
        current_metadata['Dataset'] = ["_".join(subdir.split("_")[:-1])]*current_metadata.shape[0]
        print(current_metadata.shape)

        if covar_l != []:
            current_covars = pd.read_csv(current_metadata_path, delimiter='\t', index_col=0, encoding='ISO-8859-1')[covar_l]
            current_metadata = pd.concat([current_metadata, current_covars], axis=1)
        metadata_l.append(current_metadata)

    # keep taxa present in every study
    intersect_count_data_l = [count_data[count_data.index.isin(intersect_taxa)] for count_data in count_data_l]

    combined_countdf = pd.concat(intersect_count_data_l, axis=1)
    combined_countdf = combined_countdf.dropna().T
    combined_metadata = pd.concat(metadata_l)
    combined_metadata[id] = list(combined_metadata.index) # the default IDCol for microbiomeHD will be Sam_id

    data_mat, meta_data = preprocess(combined_countdf, combined_metadata, id, covar_l)

    # align sample IDs between metadata and count table
    data_mat_ids = list(data_mat.index)
    meta_data_ids = list(meta_data.index)
    intersection_ids = list(set(meta_data_ids).intersection(data_mat_ids))

    data_mat_non_intersecting = [id for id in data_mat_ids if id not in intersection_ids]
    data_mat = data_mat.drop(data_mat_non_intersecting)
    meta_data_non_intersecting = [id for id in meta_data_ids if id not in intersection_ids]
    meta_data = meta_data.drop(meta_data_non_intersecting)
    data_mat = data_mat.reindex(intersection_ids)
    meta_data = meta_data.reindex(intersection_ids)

    if output_root != False:
        data_mat.to_csv(output_root+"_count_data.csv")
        meta_data.to_csv(output_root+"_meta_data.csv", index=False)
    return data_mat, meta_data


def load_data_custom(address_directory, output_root = False, id = 'Sam_id', covar_l = []):
    # MicrobiomeHD metadata differ across studies, so by default (covar_l=[]) only the batch
    # (Dataset) and DiseaseState are loaded.
    subdir_names = os.listdir(address_directory)
    count_data_l = []
    intersect_taxa = []
    metadata_l = []
    for subdir in subdir_names:
        # 1. OTU table
        subdir_path = address_directory + '/' + subdir
        current_files_names = os.listdir(subdir_path)
        current_dbotu_count_data_path = [result for result in current_files_names if "count_data.csv" in result][0]
        current_dbotu_count_data = pd.read_csv(subdir_path+'/'+current_dbotu_count_data_path, index_col=0)
        print(current_dbotu_count_data.shape)
                
        count_data_l.append(current_dbotu_count_data)
        
        if intersect_taxa == []:
            intersect_taxa = current_taxa
        else:
            intersect_taxa = list(set(intersect_taxa).intersection(current_taxa))

        # 2. metadata
        current_metadata_path = [result for result in current_files_names if "meta_data.csv" in result][0]
        current_metadata = pd.read_csv(subdir_path+'/'+current_metadata_path, index_col=0)

        if covar_l != []:
            current_covars = pd.read_csv(current_metadata_path, delimiter='\t', index_col=0, encoding='ISO-8859-1')[covar_l]
            current_metadata = pd.concat([current_metadata, current_covars], axis=1)
        metadata_l.append(current_metadata)

    # keep taxa present in every study
    intersect_count_data_l = [count_data[count_data.index.isin(intersect_taxa)] for count_data in count_data_l]

    combined_countdf = pd.concat(intersect_count_data_l, axis=1)
    combined_countdf = combined_countdf.dropna().T
    combined_metadata = pd.concat(metadata_l)
    combined_metadata[id] = list(combined_metadata.index) # the default IDCol for microbiomeHD will be Sam_id

    data_mat, meta_data = preprocess(combined_countdf, combined_metadata, id, covar_l)

    # align sample IDs between metadata and count table
    data_mat_ids = list(data_mat.index)
    meta_data_ids = list(meta_data.index)
    intersection_ids = list(set(meta_data_ids).intersection(data_mat_ids))

    data_mat_non_intersecting = [id for id in data_mat_ids if id not in intersection_ids]
    data_mat = data_mat.drop(data_mat_non_intersecting)
    meta_data_non_intersecting = [id for id in meta_data_ids if id not in intersection_ids]
    meta_data = meta_data.drop(meta_data_non_intersecting)
    data_mat = data_mat.reindex(intersection_ids)
    meta_data = meta_data.reindex(intersection_ids)

    if output_root != False:
        data_mat.to_csv(output_root+"_count_data.csv")
        meta_data.to_csv(output_root+"_meta_data.csv", index=False)
    return data_mat, meta_data



# overall_path = os.path.dirname(REPO_ROOT)
# # autism 2 microbiomeHD
# output_dir_path = os.path.join(REPO_ROOT, 'data/intermediate_autism_2_microbiomeHD/autism_2_microbiomeHD')
# address_directory = overall_path+'/mic_bc_benchmark/data/autism_2_microbiomeHD'
# data_mat, meta_data = load_data_microbiomeHD(address_directory, output_dir_path)

# # cdi 3 microbiomeHD
# output_dir_path = os.path.join(REPO_ROOT, 'data/intermediate_cdi_3_microbiomeHD/cdi_3_microbiomeHD')
# address_directory = overall_path+'/mic_bc_benchmark/data/cdi_3_microbiomeHD'
# data_mat, meta_data = load_data_microbiomeHD(address_directory, output_dir_path)

import argparse
parser = argparse.ArgumentParser(description='Preprocess microbiomeHD datasets')
parser.add_argument('-d', '--dataset_name', type=str, default='autism_2_microbiomeHD', help='Name of the dataset')
parser.add_argument('-s', '--source_dir', type=str, default=os.path.join(REPO_ROOT, 'data/'), help='Directory where the dataset is stored')

args = parser.parse_args()
dataset_name = args.dataset_name
source_dir = args.source_dir

source_dir_path = f'{source_dir}/{dataset_name}'
output_dir_path = f'{source_dir}/intermediate_{dataset_name}/{dataset_name}'

if not os.path.exists(output_dir_path):
    os.makedirs(output_dir_path)

if 'microbiomeHD' in dataset_name:
    data_mat, meta_data = load_data_microbiomeHD(source_dir_path, output_dir_path)
else:
    data_mat, meta_data = load_data_custom(source_dir_path, output_dir_path)
# data_mat, meta_data = load_data_microbiomeHD(source_dir_path, output_dir_path)
