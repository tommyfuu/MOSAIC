import pandas as pd
import numpy as np
import os

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

        # get covariates if exists
        if covar_l != []:
            current_covars = pd.read_csv(current_metadata_path, delimiter='\t', index_col=0, encoding='ISO-8859-1')[covar_l]
            current_metadata = pd.concat([current_metadata, current_covars], axis=1)
        metadata_l.append(current_metadata)

    # intersect count data list
    intersect_count_data_l = [count_data[count_data.index.isin(intersect_taxa)] for count_data in count_data_l]

    # generate results
    combined_countdf = pd.concat(intersect_count_data_l, axis=1)
    combined_countdf = combined_countdf.dropna().T
    combined_metadata = pd.concat(metadata_l)
    combined_metadata[id] = list(combined_metadata.index) # the default IDCol for microbiomeHD will be Sam_id

    data_mat, meta_data = preprocess(combined_countdf, combined_metadata, id, covar_l)

    # ensure that the sample ids are correctly aligned in metadata and count_table
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


# overall_path = '/athena/linglab/scratch/chf4012'
# # autism 2 microbiomeHD
# output_dir_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/intermediate_autism_2_microbiomeHD/autism_2_microbiomeHD'
# address_directory = overall_path+'/mic_bc_benchmark/data/autism_2_microbiomeHD'
# data_mat, meta_data = load_data_microbiomeHD(address_directory, output_dir_path)

# # cdi 3 microbiomeHD
# output_dir_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/intermediate_cdi_3_microbiomeHD/cdi_3_microbiomeHD'
# address_directory = overall_path+'/mic_bc_benchmark/data/cdi_3_microbiomeHD'
# data_mat, meta_data = load_data_microbiomeHD(address_directory, output_dir_path)

