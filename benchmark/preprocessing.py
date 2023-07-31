# harmonicMic - A data alignment algorithm dedicated to microbiome data.
# Copyright (C) 2022  Chenlian (Tom) Fu <chf4012@med.cornell.edu; tfu@g.hmc.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import pandas as pd
import numpy as np
import os

def load_data_microbiomeHD(address_directory, output_root = False, id = 'Sam_id', covar_l = []):
    ### note that due to the complexity of metadata, the current microbiomeHD loading does 
    ### not take into account the covariates other than batches and diseaseStates
    ### so default vars_use will just be Dataset
    subdir_names = os.listdir(address_directory)
    subdir_names = [result for result in subdir_names if "results" in result]
    count_data_l = []
    intersect_taxa = []
    metadata_l = []
    for subdir in subdir_names:
        ## 1. get the otu table part
        # get the otu file
        subdir_path = address_directory + '/' + subdir
        current_RDP_names = os.listdir(subdir_path + '/RDP')
        current_dbotu_count_data_path = [result for result in current_RDP_names if "100.denovo.rdp_assigned" in result][0]
        current_dbotu_count_data = pd.read_csv(subdir_path+'/RDP/'+current_dbotu_count_data_path, delimiter='\t', index_col=0)
        current_taxa = list(current_dbotu_count_data.index)
        # set index with genus level
        current_taxa = [";".join(taxa.split(';')[:-2]) for taxa in current_taxa]
        print(current_dbotu_count_data.shape)
        current_dbotu_count_data.index = current_taxa
        print(current_dbotu_count_data.shape)
        # remove duplicate indexes by summing up rows with the same index
        current_dbotu_count_data = current_dbotu_count_data.groupby(level=0).sum()
        
        # save dataframe and feature list
        count_data_l.append(current_dbotu_count_data)
        
        if intersect_taxa == []:
            intersect_taxa = current_taxa
        else:
            intersect_taxa = list(set(intersect_taxa).intersection(current_taxa))

        ## 2. get the metadata
        # get the metadata file
        current_files_names = os.listdir(subdir_path)
        current_metadata_path = subdir_path + '/' + [result for result in current_files_names if "metadata" in result][0]
        print("metadata path", current_metadata_path)
        current_metadata = pd.read_csv(current_metadata_path, delimiter='\t', index_col=0, encoding='ISO-8859-1')['DiseaseState'].to_frame()
        current_metadata['Dataset'] = ["_".join(subdir.split("_")[:-1])]*current_metadata.shape[0]
        print(current_metadata.shape)
        metadata_l.append(current_metadata)

    # intersect count data list
    intersect_count_data_l = [count_data[count_data.index.isin(intersect_taxa)] for count_data in count_data_l]

    # generate results
    combined_countdf = pd.concat(intersect_count_data_l, axis=1)
    combined_countdf = combined_countdf.dropna().T
    combined_metadata = pd.concat(metadata_l)
    combined_metadata[id] = list(combined_metadata.index) # the default IDCol for microbiomeHD will be Sam_id

    data_mat, meta_data = preprocess(combined_countdf, combined_metadata, id, covar_l)

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

    # save stuff if needed
    if output_root != False:
        data_mat.to_csv(output_root+"_count_data.csv")
        meta_data.to_csv(output_root+"_meta_data.csv", index=False)
    return data_mat, meta_data

def load_data_CMD(address_directory, output_root = False, id = 'Sam_id', covar_l = []):
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

def load_data(address_X, address_Y, IDCol, index_col = False, output_root = False, covar_l = []):
    if index_col != False:
        data_mat = pd.read_csv(address_X, index_col=index_col)
    else:
        data_mat = pd.read_csv(address_X)
    meta_data = pd.read_csv(address_Y)
    data_mat.index = list(meta_data[IDCol])
    data_mat, meta_data = preprocess(data_mat, meta_data, IDCol, covar_l)

    # save stuff if needed
    if output_root != False:
        data_mat.to_csv(output_root+"_count_data.csv")
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

def load_data_simulation(address, output_root = False):

    return

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

# output_root = '/home/fuc/harmonicMic/harmonypy/harmonypy/percentile_norm_data/ibd_3_CMD'
def save_data_percentile_norm(data_mat, meta_data, output_root, bio_var, bio_var_disease):
    data_mat.to_csv(output_root+"_percentile_norm.txt", sep='\t')
    meta_data_disease = meta_data[meta_data[bio_var] == bio_var_disease]
    meta_data_control = meta_data[meta_data[bio_var] != bio_var_disease]
    samples_disease = list(meta_data_disease.index)
    samples_control = list(meta_data_control.index)
    with open(output_root+'_case.txt', 'w') as f:
        for line in samples_disease:
            f.write(f"{line}\t")
    with open(output_root+'_control.txt', 'w') as f:
        for line in samples_control:
            f.write(f"{line}\t")
    return 