import os
import pandas as pd
import numpy as np

# Build a CDI-vs-Healthy-only variant of cdi_3_microbiomeHD (drops nonCDI and postFMT_CDI
# samples, which are complete confounders of study/batch -- see cdi_3_microbiomeHD_complete_confounding.csv:
# nonCDI is 100% cdi_schubert and postFMT_CDI is 100% cdi_youngster).
#
# NOTE: does not import step2_preprocessing_summarystats.py -- that module has unguarded
# top-level side effects (it regenerates crc_11_CMD and a default dataset on import).


def preprocess(data_mat, meta_data, IDCol, covar_l=[]):
    data_mat = data_mat.loc[~(data_mat == 0).all(axis=1)]
    kept_samples = data_mat.index
    meta_data = meta_data[meta_data[IDCol].isin(kept_samples)]
    col_names = list(data_mat.columns)
    col_sums = data_mat.sum(axis=1)
    removable_feature_names = [col_names[index] for index, col_sum in enumerate(col_sums) if col_sum == 0]
    data_mat.drop(removable_feature_names, axis=1, inplace=True)
    for covar in covar_l:
        meta_data = meta_data[meta_data[covar].notna()]
    data_mat = data_mat.loc[meta_data[IDCol]]
    col_names = list(data_mat.columns)
    col_sums = data_mat.sum(axis=1)
    removable_feature_names = [col_names[index] for index, col_sum in enumerate(col_sums) if col_sum == 0]
    data_mat.drop(removable_feature_names, axis=1, inplace=True)
    return data_mat, meta_data


def check_complete_confounding(meta_data, batch_var, bio_var, output_root=''):
    batch_l = list(np.unique(list(meta_data[batch_var])))
    bio_var_l = list(np.unique(list(meta_data[bio_var])))
    df = pd.DataFrame(columns=bio_var_l, index=batch_l)
    for batch in batch_l:
        for bio_var_val in bio_var_l:
            df.loc[batch, bio_var_val] = len(meta_data.loc[(meta_data[batch_var] == batch) & (meta_data[bio_var] == bio_var_val)])
    if output_root != '':
        df.to_csv(output_root + "_complete_confounding.csv")
    print(df)
    for batch in batch_l:
        if len(df.loc[batch].unique()) == 1:
            print("batch", batch, "is a complete confounder")
    return


overall_path = os.path.dirname(os.path.abspath(__file__)) + "/.."
src_root = f"{overall_path}/data/cleaned_data/cdi_3_microbiomeHD"
dst_root = f"{overall_path}/data/cleaned_data/cdi_3_microbiomeHD_clean"
os.makedirs(dst_root, exist_ok=True)

count_data = pd.read_csv(f"{src_root}/cdi_3_microbiomeHD_count_data.csv", index_col=0)
meta_data = pd.read_csv(f"{src_root}/cdi_3_microbiomeHD_meta_data.csv")

print("Before filtering:")
print(meta_data['DiseaseState'].value_counts())

meta_data = meta_data[meta_data['DiseaseState'].isin(['CDI', 'H'])]
count_data = count_data.loc[meta_data['Sam_id']]

count_data, meta_data = preprocess(count_data, meta_data, 'Sam_id', covar_l=[])

print("\nAfter filtering to CDI/H and dropping newly-all-zero taxa:")
print(meta_data['DiseaseState'].value_counts())
print(f"{count_data.shape[1]} taxa retained")

count_data.to_csv(f"{dst_root}/cdi_3_microbiomeHD_clean_count_data.csv")
meta_data.to_csv(f"{dst_root}/cdi_3_microbiomeHD_clean_meta_data.csv", index=False)
check_complete_confounding(meta_data, 'Dataset', 'DiseaseState', f"{dst_root}/cdi_3_microbiomeHD_clean")
