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

def load_data_microbiomeHD(address_directory):
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
        current_metadata = pd.read_csv(current_metadata_path, delimiter='\t', index_col=0)['DiseaseState'].to_frame()
        current_metadata['Dataset'] = ["_".join(subdir.split("_")[:-1])]*current_metadata.shape[0]
        metadata_l.append(current_metadata)

    # intersect count data list
    intersect_count_data_l = [count_data[count_data.index.isin(intersect_taxa)] for count_data in count_data_l]

    # generate results
    combined_countdf = pd.concat(intersect_count_data_l, axis=1)
    combined_metadata = pd.concat(metadata_l)
    combined_metadata['Sam_id'] = list(combined_metadata.index) # the default IDCol for microbiomeHD will be Sam_id
    return combined_countdf.dropna().T, combined_metadata

def load_data(address_X, address_Y, IDCol, index_col = False):
    if index_col != False:
        data_mat = pd.read_csv(address_X, index_col=index_col)
    else:
        data_mat = pd.read_csv(address_X)
    meta_data = pd.read_csv(address_Y)
    data_mat, meta_data = preprocess(data_mat, meta_data, IDCol)

    return data_mat, meta_data

def preprocess(data_mat, meta_data, IDCol):
    # remove samples with all zeros
    data_mat = data_mat.loc[~(data_mat==0).all(axis=1)]
    kept_samples = data_mat.index
    meta_data = meta_data[meta_data[IDCol].isin(kept_samples)]
    print(data_mat.shape, meta_data.shape)
    # remove features with all zeros
    col_names = list(data_mat.columns)
    col_sums = data_mat.sum(axis = 1)
    print(len(col_sums))
    removable_feature_names = [col_names[index] for index, col_sum in enumerate(col_sums) if col_sum==0]
    data_mat.drop(removable_feature_names, axis=1, inplace=True)
    print(data_mat.shape, meta_data.shape)

    # print(data_mat.shape)
    # # remove features where the sum of counts are below 0.01% compared to the total sum of all counts (include those features with all zeroes)
    # # Arumugam, Manimozhiyan, Jeroen Raes, Eric Pelletier, Denis Le Paslier, Takuji Yamada, Daniel R Mende, Gabriel R Fernandes, et al. 2011. “Enterotypes of the Human Gut Microbiome.” Nature 473 (7346). Nature Publishing Group: 174.
    # col_names = list(data_mat.columns)
    # col_sums = data_mat.sum(axis = 0)
    # total = data_mat.to_numpy().sum()
    # col_sums_ratios = [col_sum/total for col_sum in col_sums]
    # removable_feature_names = [col_names[index] for index, ratio in enumerate(col_sums_ratios) if ratio<0.0001]
    # data_mat.drop(removable_feature_names, axis=1, inplace=True)
    return data_mat, meta_data


