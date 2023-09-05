import pandas as pd
import numpy as np
import os


def consensus_biomarker(methods_l, paths, dataset):
    print("______________________________________________________")
    print('Dataset: ', dataset)
    method_biomarker_dict = {}
    all_unique_taxa = []
    for idx, path in enumerate(paths):
        files_in_path = os.listdir(path)
        current_file_path = [file for file in files_in_path if 'diff_abund_test.csv' in file][0]
        df = pd.read_csv(path+'/'+current_file_path)
        sig_taxa = list(df[df['FDR_p_value']<0.05]['Unnamed: 0'])
        print(methods_l[idx], len(sig_taxa))
        method_biomarker_dict[methods_l[idx]] = sig_taxa
        all_unique_taxa += sig_taxa
    all_unique_taxa = list(set(all_unique_taxa))

    # turn into a presence/absence matrix
    method_biomarker_dict_presence = {}
    for method in method_biomarker_dict.keys():
        for taxon in all_unique_taxa:
            if method not in method_biomarker_dict_presence.keys():
                method_biomarker_dict_presence[method] = []
            if taxon not in method_biomarker_dict[method]:
                method_biomarker_dict_presence[method].append(0)
            else:
                method_biomarker_dict_presence[method].append(1)

    # turn dict into dataframe which is a presence/absence matrix where rows are taxa and columns are methods
    df = pd.DataFrame.from_dict(method_biomarker_dict_presence, orient='index')
    df = df.transpose()
    df.index = all_unique_taxa
    print(df)
    # sort the rows by sum, descending
    df['sum'] = df.sum(axis=1)
    df = df.sort_values(by=['sum'], ascending=False)
    df.to_csv('/athena/linglab/scratch/chf4012/mic_bc_benchmark/outputs/consensus_taxa/consensus_diff_abn_'+dataset+'.csv')
    return df

methods_l = ['combat_seq', 'ConQuR', 'ConQuR_libsize', 'harmony', 'limma', 'MMUPHin', 'nobc', 'percentile_norm']
general_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/outputs/'
dataset = 'autism_2_microbiomeHD'
paths = [general_path+dataset+'/output_'+dataset+'_'+method+'/' for method in methods_l]
df = consensus_biomarker(methods_l, paths, dataset)

dataset = 'cdi_3_microbiomeHD'
paths = [general_path+dataset+'/output_'+dataset+'_'+method+'/' for method in methods_l]
df = consensus_biomarker(methods_l, paths, dataset)

methods_l = ['combat', 'ConQuR_rel', 'harmony', 'limma', 'MMUPHin', 'nobc', 'percentile_norm']

dataset = 'ibd_3_CMD'
paths = [general_path+dataset+'/output_'+dataset+'_'+method+'/' for method in methods_l]
df = consensus_biomarker(methods_l, paths, dataset)

dataset = 'CRC_8_CMD'
paths = [general_path+dataset+'/output_'+dataset+'_'+method+'/' for method in methods_l]
df = consensus_biomarker(methods_l, paths, dataset)
        