import pandas as pd

dataset_names = ['autism_2_microbiomeHD', 'cdi_3_microbiomeHD', 'ibd_3_CMD', 'crc_8_CMD']
# dataset_names = ['autism_2_microbiomeHD', 'cdi_3_microbiomeHD']
dir_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/outputs'
dataset_methods_dict = {'autism_2_microbiomeHD': ["nobc", "harmony", "combat_seq", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"],
                     'cdi_3_microbiomeHD': ["nobc", "harmony", "combat_seq", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"], 
                     'ibd_3_CMD': ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR_rel", "percentile_norm"], 
                     'crc_8_CMD': ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR_rel", "percentile_norm"]}
# get differential abundance results
for dataset_name in dataset_names:
    current_dataset_method_to_taxa_dict = {}
    for method in dataset_methods_dict[dataset_name]:
        print(dataset_name, method)
        # read the runtime file
        runtime_df = pd.read_csv(f'{dir_path}/{dataset_name}/output_{dataset_name}_{method}/{dataset_name}_{method}_diff_abund_test.csv', index_col=0)
        # add to dict with taxa (index) with FDR_p_value < 0.05
        runtime_df_sig = runtime_df[runtime_df['FDR_p_value'] < 0.05]
        print(runtime_df_sig.shape)
        taxa_l = list(runtime_df_sig.index)
        current_dataset_method_to_taxa_dict[method] = [taxa_l]
        # clean up taxa names, for the sake of consistency, get genus level data
        if dataset_name in ['autism_2_microbiomeHD', 'cdi_3_microbiomeHD']:
            delimiter = ';'
        else:
            delimiter = '|'
        genus_idx = -1
        # taxa_l = [taxa.split(delimiter)[genus_idx].split('g__')[-1] for taxa in taxa_l]
        # taxa_l = [taxa.replace('_', ' ') for taxa in taxa_l]
        if dataset_name in ['autism_2_microbiomeHD', 'cdi_3_microbiomeHD']:
            taxa_l = [taxa.split(delimiter)[genus_idx].split('g__')[-1] for taxa in taxa_l]
        else:
            taxa_l = [taxa.split(delimiter)[genus_idx].split('g__')[-1].split('__')[1] for taxa in taxa_l]
        current_dataset_method_to_taxa_dict[method].append(taxa_l)
        current_dataset_method_to_taxa_dict[method] = dict(zip(current_dataset_method_to_taxa_dict[method][0], current_dataset_method_to_taxa_dict[method][1]))
    
    # get the universal dict of original taxa name to genus name
    original_taxa_to_genus_dict = {}
    for method in current_dataset_method_to_taxa_dict:
        original_taxa_to_genus_dict.update(current_dataset_method_to_taxa_dict[method])

    # get the union of taxa, form pd dataframe, with columns showing which methods have the taxa
    taxa_union = set()
    for method in current_dataset_method_to_taxa_dict:
        taxa_union.update(set(current_dataset_method_to_taxa_dict[method].keys()))
    taxa_union = list(taxa_union)
    taxa_union.sort()
    taxa_union_df = pd.DataFrame(index=taxa_union)
    # add a column for cleaned genus names
    taxa_union_df['cleaned_taxa'] = [original_taxa_to_genus_dict[taxa] for taxa in taxa_union_df.index]
    # add columns for each method
    for method in current_dataset_method_to_taxa_dict:
        taxa_union_df[method] = taxa_union_df.index.isin(current_dataset_method_to_taxa_dict[method])
    # write the taxa union to a file, this will serve as a backbone for pubtator search
    taxa_union_df.to_csv(f'{dir_path}/{dataset_name}_taxa_union_wilcoxon.csv')

