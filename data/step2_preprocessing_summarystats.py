import os, sys
import pandas as pd
import numpy as np

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

def preprocess_data_phyloseq(address_directory, output_root = False, id = 'Sam_id', covar_l = [], relab = True):
    ### CuratedMetagenomicsDataset provides way more metadata in a congestible manner
    cur_dir_names = os.listdir(address_directory)
    address_X = address_directory + '/'+ [result for result in cur_dir_names if "otu_table_" in result][0]
    address_Y = address_directory + '/'+ [result for result in cur_dir_names if "sample_table_" in result][0]
    data_mat = pd.read_csv(address_X, index_col=0)
    meta_data = pd.read_csv(address_Y, index_col=0)
    meta_data[id] = list(meta_data.index)
    data_mat, meta_data = preprocess(data_mat.T, meta_data, id, covar_l)
    
    if not relab:
        # get taxa names
        taxa_df = pd.read_csv(address_directory + '/'+ [result for result in cur_dir_names if "tax_table_" in result][0], index_col=0)
        if 'ta6' in taxa_df.columns:
            taxa_names = list(taxa_df['ta6'])
        elif 'genus' in taxa_df.columns:
            taxa_names = list(taxa_df['genus'])
        else:
            taxa_names = list(taxa_df['Genus'])
        data_mat.columns = taxa_names
        
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
    if relab:
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

def check_complete_confounding(meta_data, batch_var, bio_var, output_root = ''):
    # make a pandas dataframe where rows are batches whereas columns are the bio_var options
    # each entry is the number of samples in that batch with that bio_var option
    # if there is a batch with only one bio_var option, then it is a complete confounder
    print(meta_data)
    # get the list of batches
    batch_l = list(meta_data[batch_var])
    # batch_l = [x for x in batch_l if str(x) != 'nan']
    batch_l = list(np.unique(batch_l))

    # get the list of bio_var options
    bio_var_l = list(meta_data[bio_var])
    # bio_var_l = [x for x in bio_var_l if str(x) != 'nan']
    bio_var_l = list(np.unique(bio_var_l))

    # generate a dataframe
    df = pd.DataFrame(columns=bio_var_l, index=batch_l)
    for batch in batch_l:
        for bio_var_val in bio_var_l:
            df.loc[batch, bio_var_val] = len(meta_data.loc[(meta_data[batch_var]==batch) & (meta_data[bio_var]==bio_var_val)])
    
    if output_root != '':
        df.to_csv(output_root+"_complete_confounding.csv")
    print(df)
    # check if there is a batch with only one bio_var option
    for batch in batch_l:
        if len(df.loc[batch].unique())==1:
            print("batch", batch, "is a complete confounder")
    return



# 
# overall_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data'
# # autism_2_microbiomeHD
# data_mat, meta_data = preprocess_data_phyloseq(f'{overall_path}/pruned_autism_2_microbiomeHD', f'{overall_path}/cleaned_data/autism_2_microbiomeHD/autism_2_microbiomeHD', id = 'Sam_id', covar_l = [], relab = False)
# data_mat, meta_data = load_results_from_benchmarked_methods(f'{overall_path}/cleaned_data/autism_2_microbiomeHD/autism_2_microbiomeHD_count_data.csv', f'{overall_path}/cleaned_data/autism_2_microbiomeHD/autism_2_microbiomeHD_meta_data.csv')
# check_complete_confounding(meta_data, 'Dataset', 'DiseaseState', f'{overall_path}/cleaned_data/autism_2_microbiomeHD/autism_2_microbiomeHD')

# # cdi_3_microbiomeHD
# data_mat, meta_data = preprocess_data_phyloseq(f'{overall_path}/pruned_cdi_3_microbiomeHD', f'{overall_path}/cleaned_data/cdi_3_microbiomeHD/cdi_3_microbiomeHD', id = 'Sam_id', covar_l = [], relab = False)
# data_mat, meta_data = load_results_from_benchmarked_methods(f'{overall_path}/cleaned_data/cdi_3_microbiomeHD/cdi_3_microbiomeHD_count_data.csv', f'{overall_path}/cleaned_data/cdi_3_microbiomeHD/cdi_3_microbiomeHD_meta_data.csv')
# check_complete_confounding(meta_data, 'Dataset', 'DiseaseState', f'{overall_path}/cleaned_data/cdi_3_microbiomeHD/cdi_3_microbiomeHD')

# # ibd_3_CMD
# data_mat, meta_data = preprocess_data_phyloseq(f'{overall_path}/pruned_ibd_3_CMD', f'{overall_path}/cleaned_data/ibd_3_CMD/ibd_3_CMD', id = 'Sam_id', covar_l = ['disease', 'gender', 'age_category'])
# data_mat, meta_data = load_results_from_benchmarked_methods(f'{overall_path}/cleaned_data/ibd_3_CMD/ibd_3_CMD_count_data.csv', f'{overall_path}/cleaned_data/ibd_3_CMD/ibd_3_CMD_meta_data.csv')
# check_complete_confounding(meta_data, "study_name", "disease", f'{overall_path}/cleaned_data/ibd_3_CMD/ibd_3_CMD')

# # crc_8_CMD
# data_mat, meta_data = preprocess_data_phyloseq(f'{overall_path}/pruned_crc_8_CMD', f'{overall_path}/cleaned_data/crc_8_CMD/crc_8_CMD', id = 'Sam_id', covar_l = [])
# data_mat, meta_data = load_results_from_benchmarked_methods(f'{overall_path}/cleaned_data/crc_8_CMD/crc_8_CMD_count_data.csv', f'{overall_path}/cleaned_data/crc_8_CMD/crc_8_CMD_meta_data.csv')
# check_complete_confounding(meta_data, "study_name", "disease", f'{overall_path}/cleaned_data/crc_8_CMD/crc_8_CMD')

import argparse
parser = argparse.ArgumentParser(description='Calculating summary stats for real-world microbiome datasets')
parser.add_argument('-d', '--dataset_name', type=str, default='autism_2_microbiomeHD', help='Name of the dataset')
parser.add_argument('-s', '--source_dir', type=str, default='/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/', help='Directory where the dataset is stored')
parser.add_argument('-b', '--batch_var', type=str, default='Dataset', help='Variable to calculate summary stats for')
parser.add_argument('-c', '--condition_var', type=str, default='DiseaseState', help='Variable to calculate summary stats for')

args = parser.parse_args()
dataset_name = args.dataset_name
source_dir = args.source_dir
batch_var = args.batch_var
condition_var = args.condition_var

os.chdir(sys.path[0])

# load covar_l from yaml file
import yaml
with open(f'../config.yml') as file:
    covar_l = yaml.load(file, Loader=yaml.FullLoader)['COVAR_L']

# mkdir cleaned_data if not exists
if not os.path.exists(f'{source_dir}/cleaned_data/{dataset_name}'):
    os.makedirs(f'{source_dir}/cleaned_data/{dataset_name}')

data_mat, meta_data = preprocess_data_phyloseq(f'{source_dir}/pruned_{dataset_name}', f'{source_dir}/cleaned_data/{dataset_name}/{dataset_name}', id = 'Sam_id', covar_l = covar_l, relab = False)
data_mat, meta_data = load_results_from_benchmarked_methods(f'{source_dir}/cleaned_data/{dataset_name}/{dataset_name}_count_data.csv', f'{source_dir}/cleaned_data/{dataset_name}/{dataset_name}_meta_data.csv')
check_complete_confounding(meta_data, batch_var, condition_var, f'{source_dir}/cleaned_data/{dataset_name}/{dataset_name}')
