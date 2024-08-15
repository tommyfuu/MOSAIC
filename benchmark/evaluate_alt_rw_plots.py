import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib as mpl
import os
import yaml

def alt_visualize_simulation_stats(output_root, output_dir_l, datasets, methods, highlighted_method, ylim_dict, simulate = False, sim_num_iters = 1000, dimensions = (7, 5), taxa_gt = None, count_l = [True, True, False, False], 
    marker_dict = {'autism_2_microbiomeHD': 'o', 'cdi_3_microbiomeHD': 's', 'ibd_3_CMD': 'd', 'crc_8_CMD': 'H'},
    ds_colors_dict = {'autism_2_microbiomeHD': 'violet', 'cdi_3_microbiomeHD': 'black', 'ibd_3_CMD': 'orange', 'crc_8_CMD': 'blue'},
    postfix = '.png', demonstrate = False):
    '''visualize the PERMANOVA batch R2 (Bray/Aitch), PERMANOVA condition R2 (Bray/Aitch), ROC-AUC and FDR/sensitivity'''
    # global set up
    print("initializing")
    global_methods_batch_bray_r2_l_dict = {}
    global_methods_batch_aitch_r2_l_dict = {}
    global_methods_biovar_bray_r2_l_dict = {}
    global_methods_biovar_aitch_r2_l_dict = {}
    global_methods_batch_shannon_pval_l_dict = {}
    global_methods_biovar_shannon_pval_l_dict = {}
    global_methods_rf_auc_l_dict = {}
    global_methods_rf_f1_l_dict = {}
    global_methods_rf_precision_l_dict = {}
    global_methods_rf_recall_l_dict = {}
    global_methods_runtime_l_dict = {}

    if taxa_gt is not None:
        global_methods_FDR_r2_l_dict = {}
        global_methods_sensitivity_r2_l_dict = {}

    # loop through datasets to get the stats
    def get_stats(output_dir, taxa_gt, idx):
        methods_batch_aitch_r2_dict = {}
        methods_biovar_aitch_r2_dict = {}
        # get the global stats path
        files = os.listdir(output_dir)
        global_stats_path = [file for file in files if "global_benchmarking_stats" in file][0]

        # read the global stats df
        global_stats_df = pd.read_csv(output_dir+"/"+global_stats_path, index_col=0)
        # get the stats of interest and visualize
        if count_l[idx]:
            batch_aitch_r2 = global_stats_df.loc["batch_aitch_r2"]
            methods_batch_aitch_r2_dict = batch_aitch_r2.to_dict()
            biovar_aitch_r2 = global_stats_df.loc["biovar_aitch_r2"]
            methods_biovar_aitch_r2_dict = biovar_aitch_r2.to_dict()
        
        batch_bray_r2 = global_stats_df.loc["batch_bray_r2"]  
        methods_batch_bray_r2_dict = batch_bray_r2.to_dict()
        biovar_bray_r2 = global_stats_df.loc["biovar_bray_r2"]
        methods_biovar_bray_r2_dict = biovar_bray_r2.to_dict()
        batch_shannon_pval = global_stats_df.loc["batch_shannon_pval"]
        methods_batch_shannon_pval_dict = batch_shannon_pval.to_dict()
        biovar_shannon_pval = global_stats_df.loc["biovar_shannon_pval"]
        rf_auc = global_stats_df.loc["rf_auc"]
        methods_rf_auc_dict = rf_auc.to_dict()
        rf_f1 = global_stats_df.loc["rf_weighted_f1"]
        methods_rf_f1_dict = rf_f1.to_dict()
        rf_precision = global_stats_df.loc["rf_weighted_precision"]
        methods_rf_precision_dict = rf_precision.to_dict()
        rf_recall = global_stats_df.loc["rf_weighted_recall"]
        methods_rf_recall_dict = rf_recall.to_dict()
        methods_biovar_shannon_pval_dict = biovar_shannon_pval.to_dict()
        runtime = global_stats_df.loc["runtime"]
        methods_runtime = runtime.to_dict()

        if taxa_gt is not None:
            FDR_r2 = global_stats_df.loc["FDR"]
            methods_FDR_r2_dict = FDR_r2.to_dict()
            sensitivity_r2 = global_stats_df.loc["power"]
            methods_sensitivity_r2_dict = sensitivity_r2.to_dict()

        if taxa_gt is not None:
            return methods_batch_aitch_r2_dict, methods_batch_bray_r2_dict, methods_biovar_aitch_r2_dict, methods_biovar_bray_r2_dict, methods_FDR_r2_dict, methods_sensitivity_r2_dict, methods_batch_shannon_pval_dict, methods_biovar_shannon_pval_dict, methods_rf_auc_dict, methods_rf_f1_dict, methods_rf_precision_dict, methods_rf_recall_dict, methods_runtime
        else:
            return methods_batch_aitch_r2_dict, methods_batch_bray_r2_dict, methods_biovar_aitch_r2_dict, methods_biovar_bray_r2_dict, methods_batch_shannon_pval_dict, methods_biovar_shannon_pval_dict, methods_rf_auc_dict, methods_rf_f1_dict, methods_rf_precision_dict, methods_rf_recall_dict, methods_runtime
    print("extracting data")
    for idx, output_dir in enumerate(output_dir_l):
        print(idx)
        if not simulate:   
            if taxa_gt is not None:
                methods_batch_aitch_r2_dict, methods_batch_bray_r2_dict, methods_biovar_aitch_r2_dict, methods_biovar_bray_r2_dict, methods_FDR_r2_dict, methods_sensitivity_r2_dict, methods_batch_shannon_pval_dict, methods_biovar_shannon_pval_dict, methods_rf_auc_dict, methods_rf_f1_dict, methods_rf_precision_dict, methods_rf_recall_dict, methods_runtime = get_stats(output_dir, taxa_gt, idx)    
            else:
                methods_batch_aitch_r2_dict, methods_batch_bray_r2_dict, methods_biovar_aitch_r2_dict, methods_biovar_bray_r2_dict, methods_batch_shannon_pval_dict, methods_biovar_shannon_pval_dict, methods_rf_auc_dict, methods_rf_f1_dict, methods_rf_precision_dict, methods_rf_recall_dict, methods_runtime = get_stats(output_dir, taxa_gt, idx)    
            ## append to global dict
            for method in methods:
                if method not in global_methods_batch_bray_r2_l_dict.keys():
                    if count_l[idx]:
                        global_methods_batch_aitch_r2_l_dict[method] = []
                        global_methods_biovar_aitch_r2_l_dict[method] = []
                    if taxa_gt is not None:
                        global_methods_FDR_r2_l_dict[method] = []
                        global_methods_sensitivity_r2_l_dict[method] = []
                    global_methods_batch_bray_r2_l_dict[method] = []
                    global_methods_biovar_bray_r2_l_dict[method] = []
                    global_methods_batch_shannon_pval_l_dict[method] = []
                    global_methods_biovar_shannon_pval_l_dict[method] = []
                    global_methods_rf_auc_l_dict[method] = []
                    global_methods_rf_f1_l_dict[method] = []
                    global_methods_rf_precision_l_dict[method] = []
                    global_methods_rf_recall_l_dict[method] = []
                    global_methods_runtime_l_dict[method] = []
                                
                if count_l[idx]:
                    global_methods_batch_aitch_r2_l_dict[method].append(methods_batch_aitch_r2_dict[method])
                    global_methods_biovar_aitch_r2_l_dict[method].append(methods_biovar_aitch_r2_dict[method])
                if taxa_gt is not None:
                    global_methods_FDR_r2_l_dict[method].append(methods_FDR_r2_dict[method])
                    global_methods_sensitivity_r2_l_dict[method].append(methods_sensitivity_r2_dict[method])
                global_methods_batch_bray_r2_l_dict[method].append(methods_batch_bray_r2_dict[method])
                global_methods_biovar_bray_r2_l_dict[method].append(methods_biovar_bray_r2_dict[method])
                global_methods_batch_shannon_pval_l_dict[method].append(methods_batch_shannon_pval_dict[method])
                global_methods_biovar_shannon_pval_l_dict[method].append(methods_biovar_shannon_pval_dict[method])
                global_methods_rf_auc_l_dict[method].append(methods_rf_auc_dict[method])
                global_methods_rf_f1_l_dict[method].append(methods_rf_f1_dict[method])
                global_methods_rf_precision_l_dict[method].append(methods_rf_precision_dict[method])
                global_methods_rf_recall_l_dict[method].append(methods_rf_recall_dict[method])
                global_methods_runtime_l_dict[method].append(methods_runtime[method])

        else:
            # in this case, for each dataset, we have a preset number of iterations
            cross_iter_batch_aitch_r2_dict = {}
            cross_iter_batch_bray_r2_dict = {}
            cross_iter_biovar_aitch_r2_dict = {}
            cross_iter_biovar_bray_r2_dict = {}
            cross_iter_FDR_r2_dict = {}
            cross_iter_sensitivity_r2_dict = {}
            cross_iter_batch_shannon_pval_dict = {}
            cross_iter_biovar_shannon_pval_dict = {}
            cross_iter_rf_auc_dict = {}
            cross_iter_rf_f1_dict = {}
            cross_iter_rf_precision_dict = {}
            cross_iter_rf_recall_dict = {}
            cross_iter_runtime_dict = {}
            for iter in range(sim_num_iters):
                if taxa_gt is not None:
                    methods_batch_aitch_r2_dict, methods_batch_bray_r2_dict, methods_biovar_aitch_r2_dict, methods_biovar_bray_r2_dict, methods_FDR_r2_dict, methods_sensitivity_r2_dict, methods_batch_shannon_pval_dict, methods_biovar_shannon_pval_dict, methods_rf_auc_dict, methods_rf_f1_dict, methods_rf_precision_dict, methods_rf_recall_dict, methods_runtime = get_stats(output_dir+f"_iter_{iter+1}", taxa_gt, idx)    
                else:
                    methods_batch_aitch_r2_dict, methods_batch_bray_r2_dict, methods_biovar_aitch_r2_dict, methods_biovar_bray_r2_dict, methods_rf_auc_dict, methods_rf_f1_dict, methods_rf_precision_dict, methods_rf_recall_dict = get_stats(output_dir+f"_iter_{iter+1}", taxa_gt, idx)    
                cross_iter_batch_aitch_r2_dict[iter] = methods_batch_aitch_r2_dict
                cross_iter_batch_bray_r2_dict[iter] = methods_batch_bray_r2_dict
                cross_iter_biovar_aitch_r2_dict[iter] = methods_biovar_aitch_r2_dict
                cross_iter_biovar_bray_r2_dict[iter] = methods_biovar_bray_r2_dict
                cross_iter_batch_shannon_pval_dict[iter] = methods_batch_shannon_pval_dict
                cross_iter_biovar_shannon_pval_dict[iter] = methods_biovar_shannon_pval_dict
                cross_iter_rf_auc_dict[iter] = methods_rf_auc_dict
                cross_iter_rf_f1_dict[iter] = methods_rf_f1_dict
                cross_iter_rf_precision_dict[iter] = methods_rf_precision_dict
                cross_iter_rf_recall_dict[iter] = methods_rf_recall_dict
                cross_iter_runtime_dict[iter] = methods_runtime
                if taxa_gt is not None:
                    cross_iter_FDR_r2_dict[iter] = methods_FDR_r2_dict
                    cross_iter_sensitivity_r2_dict[iter] = methods_sensitivity_r2_dict

            # calculate mean across iterations and append to global dict
            for method in methods:
                if method not in global_methods_batch_bray_r2_l_dict.keys():
                    if count_l[idx]:
                        global_methods_batch_aitch_r2_l_dict[method] = []
                        global_methods_biovar_aitch_r2_l_dict[method] = []
                    global_methods_batch_bray_r2_l_dict[method] = []
                    global_methods_biovar_bray_r2_l_dict[method] = []
                    global_methods_batch_shannon_pval_l_dict[method] = []
                    global_methods_biovar_shannon_pval_l_dict[method] = []
                    global_methods_rf_auc_l_dict[method] = []
                    global_methods_rf_f1_l_dict[method] = []
                    global_methods_rf_precision_l_dict[method] = []
                    global_methods_rf_recall_l_dict[method] = []
                    global_methods_runtime_l_dict[method] = []
                    if taxa_gt is not None:
                        global_methods_FDR_r2_l_dict[method] = []
                        global_methods_sensitivity_r2_l_dict[method] = []
                if count_l[idx]:
                    global_methods_batch_aitch_r2_l_dict[method].append(np.mean([cross_iter_batch_aitch_r2_dict[iter][method] for iter in range(sim_num_iters)]))
                    global_methods_biovar_aitch_r2_l_dict[method].append(np.mean([cross_iter_biovar_aitch_r2_dict[iter][method] for iter in range(sim_num_iters)]))
                global_methods_batch_bray_r2_l_dict[method].append(np.mean([cross_iter_batch_bray_r2_dict[iter][method] for iter in range(sim_num_iters)]))
                global_methods_biovar_bray_r2_l_dict[method].append(np.mean([cross_iter_biovar_bray_r2_dict[iter][method] for iter in range(sim_num_iters)]))
                global_methods_batch_shannon_pval_l_dict[method].append(np.mean([cross_iter_batch_shannon_pval_dict[iter][method] for iter in range(sim_num_iters)]))
                global_methods_biovar_shannon_pval_l_dict[method].append(np.mean([cross_iter_biovar_shannon_pval_dict[iter][method] for iter in range(sim_num_iters)]))
                global_methods_rf_auc_l_dict[method].append(np.mean([cross_iter_rf_auc_dict[iter][method] for iter in range(sim_num_iters)]))
                global_methods_rf_f1_l_dict[method].append(np.mean([cross_iter_rf_f1_dict[iter][method] for iter in range(sim_num_iters)]))
                global_methods_rf_precision_l_dict[method].append(np.mean([cross_iter_rf_precision_dict[iter][method] for iter in range(sim_num_iters)]))
                global_methods_rf_recall_l_dict[method].append(np.mean([cross_iter_rf_recall_dict[iter][method] for iter in range(sim_num_iters)]))
                global_methods_runtime_l_dict[method].append(np.mean([cross_iter_runtime_dict[iter][method] for iter in range(sim_num_iters)]))
                if taxa_gt is not None:
                    global_methods_FDR_r2_l_dict[method].append(np.mean([cross_iter_FDR_r2_dict[iter][method] for iter in range(sim_num_iters)]))
                    global_methods_sensitivity_r2_l_dict[method].append(np.mean([cross_iter_sensitivity_r2_dict[iter][method] for iter in range(sim_num_iters)]))

    print("plotting")
    def plot_stats(stats_summary_name, stats_name_l, stats_dict_1, stats_dict_2 = {}, ylim_l = [], postfix = '.png', pvalline = False):
        mpl.rcParams['pdf.fonttype'] = 42 # ensure exported pdf has edited text
        linestyle = '-'
        ## plot the dictionaries in two matplotlib subplots as line plots
        plt.clf()
        if stats_dict_2 != {}:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(dimensions[0]*2, dimensions[1]))
        else:
            print("only one plot")
            fig, ax1 = plt.subplots(1, 1, figsize=(dimensions[0], dimensions[1]))

        # log2 case
        if 'FDR' in stats_name_l or 'runtime' in stats_name_l:
            # deep copy dictionaries
            import copy
            stats_dict_1_clone = copy.deepcopy(stats_dict_1)
            stats_dict_2_clone = copy.deepcopy(stats_dict_2)
            for method in stats_dict_1.keys():
                # very small added to avoid log(0)
                stats_dict_1[method] = [np.log2(x+1e-4) for x in stats_dict_1[method]]
                if stats_dict_2 != {}:
                    stats_dict_2[method] = [np.log2(x+1e-4) for x in stats_dict_2[method]]

        # now use methods as x axis and stats as y axis, with two lines corresponding to the two datasets
        ax1.tick_params(axis='both', which='major', labelsize=14)
        if stats_dict_2 != {}:
            ax2.tick_params(axis='both', which='major', labelsize=14)
        for dataset in datasets:
            color = ds_colors_dict[dataset]
            marker = marker_dict[dataset]
            current_ds_stats_dict_1 = {method: stats_dict_1[method][datasets.index(dataset)] for method in methods}
            if stats_dict_2 != {}:
                current_ds_stats_dict_2 = {method: stats_dict_2[method][datasets.index(dataset)] for method in methods}

            # plot the line for each dataset
            ax1.plot(current_ds_stats_dict_1.keys(), current_ds_stats_dict_1.values(),label=dataset, color=color, marker=marker, linestyle=linestyle)
            if stats_dict_2 != {}:
                ax2.plot(current_ds_stats_dict_2.keys(), current_ds_stats_dict_2.values(),label=dataset, color=color, marker=marker, linestyle=linestyle)
            
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)
        ax1.set_title(stats_name_l[0])
        ax1.set_xticks(methods)
        # set xticks to be verticle
        for tick in ax1.get_xticklabels():
            tick.set_rotation(90)
        
        if 'runtime' in stats_name_l:
            ax1.set_yticks([np.log2(1e-4), np.log2(10), np.log2(25), np.log2(50), np.log2(100), np.log2(200), np.log2(350), np.log2(800), np.log2(2000)], ["~0", "10", "25", "50", "100", "200", "350", "800", "2000"])

        if stats_dict_2 != {}:
            ax2.spines["top"].set_visible(False)
            ax2.spines["right"].set_visible(False)
            ax2.set_title(stats_name_l[1])
            ax2.set_xticks(methods)
            # set xticks to be verticle
            for tick in ax2.get_xticklabels():
                tick.set_rotation(90)

        plt.subplots_adjust(right=0.8)
        if stats_dict_2 != {}:
            ax2.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        else:
            ax1.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)

        # optionally set ylim
        if ylim_l != [] and 'rw' in output_root:
            ax1.set_ylim(ylim_l[0])
            if stats_dict_2 != {}:
                ax2.set_ylim(ylim_l[1])
        
        # optionally add 0.05 significance line
        if pvalline:
            ax1.axhline(y=0.05, color='r', linestyle='--')
            if stats_dict_2 != {}:
                ax2.axhline(y=0.05, color='r', linestyle='--')

        # find parent dir of output_dir_l[0]
        plt.savefig(output_root+"_"+stats_summary_name+postfix, bbox_inches="tight")
        plt.clf()
        plt.close()

    # plot
    if taxa_gt is not None:
        plot_stats('FDR_sensitivity', ["FDR", "Sensitivity"], global_methods_FDR_r2_l_dict, global_methods_sensitivity_r2_l_dict, ylim_l = ylim_dict['FDR_sensitivity'], postfix=postfix)
    
    plot_stats('runtime', ["runtime"], global_methods_runtime_l_dict, ylim_l = ylim_dict['runtime'], postfix=postfix)
    plot_stats('auc and weighted f1', ["auc", "weighted f1"], global_methods_rf_auc_l_dict, global_methods_rf_f1_l_dict, ylim_l = ylim_dict['auc_f1'], postfix=postfix)
    plot_stats('weighted precision and weighted recall', ["weighted precision", "weighted recall"], global_methods_rf_precision_l_dict, global_methods_rf_recall_l_dict, ylim_l = ylim_dict['precision_recall'], postfix=postfix)
    plot_stats('shannon_pval', ["PERMANOVA batch Shannon pval", "PERMANOVA biovar Shannon pval"], global_methods_batch_shannon_pval_l_dict, global_methods_biovar_shannon_pval_l_dict, ylim_l = ylim_dict['shannon_pval'], postfix=postfix, pvalline=True)

    if not demonstrate:
        if count_l[0]:
            plot_stats('PERMANOVA_batch_R2', ["PERMANOVA batch R2 (Aitchinson)", "PERMANOVA batch R2 (Bray-Curtis)"], global_methods_batch_aitch_r2_l_dict, global_methods_batch_bray_r2_l_dict, ylim_l = ylim_dict['batch_r2'], postfix=postfix)
            plot_stats('PERMANOVA_biovar_R2', ["PERMANOVA biovar R2 (Aitchinson)", "PERMANOVA biovar R2 (Bray-Curtis)"], global_methods_biovar_aitch_r2_l_dict, global_methods_biovar_bray_r2_l_dict, ylim_l = ylim_dict['biovar_r2'], postfix=postfix)
        else:
            plot_stats('PERMANOVA_batch_R2', ["PERMANOVA batch R2 (Bray-Curtis)"], global_methods_batch_bray_r2_l_dict, ylim_l = ylim_dict['batch_r2'], postfix=postfix)
            plot_stats('PERMANOVA_biovar_R2', ["PERMANOVA biovar R2 (Bray-Curtis)"], global_methods_biovar_bray_r2_l_dict, ylim_l = ylim_dict['biovar_r2'], postfix=postfix)
    else:
        if count_l[0]:
            plot_stats('PERMANOVA_batch_R2', ["PERMANOVA batch R2 (Aitchinson)", "PERMANOVA batch R2 (Bray-Curtis)"], global_methods_batch_aitch_r2_l_dict, global_methods_batch_bray_r2_l_dict, postfix=postfix)
            plot_stats('PERMANOVA_biovar_R2', ["PERMANOVA biovar R2 (Aitchinson)", "PERMANOVA biovar R2 (Bray-Curtis)"], global_methods_biovar_aitch_r2_l_dict, global_methods_biovar_bray_r2_l_dict, postfix=postfix)
        else:
            plot_stats('PERMANOVA_batch_R2', ["PERMANOVA batch R2 (Bray-Curtis)"], global_methods_batch_bray_r2_l_dict, postfix=postfix)
            plot_stats('PERMANOVA_biovar_R2', ["PERMANOVA biovar R2 (Bray-Curtis)"], global_methods_biovar_bray_r2_l_dict, postfix=postfix)        
    return

# ## VISUALIZE LINE PLOTS FOR 2 COUNT-TYPE RW DATASETS and 2 RELAB-TYPE RW DATASETS
# #############################################################################
# output_dir_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/outputs'
# methods = ["nobc", "harmony", "combat_seq", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
# datasets = ["autism_2_microbiomeHD", "cdi_3_microbiomeHD"]
# output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
# ylim_dict = {'batch_r2': [[0, 0.35], [0, 0.35]], 'biovar_r2': [[0, 0.2], [0, 0.2]], 'shannon_pval': [[0, 1], [0, 1]], 'auc_f1': [[0.2, 0.9], [0.2, 0.9]], 'runtime': [[np.log2(5e-5), np.log2(500)]], 'precision_recall': [[0, 1], [0, 1]]}
# alt_visualize_simulation_stats('/athena/linglab/scratch/chf4012/mic_bc_benchmark/outputs/rw_data_alt_plots/count_rw', output_dir_l, datasets, methods, highlighted_method = "ConQuR",  simulate = False, count_l = [True, True], postfix = '.pdf', ylim_dict = ylim_dict)

# output_dir_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/outputs'
# methods = ["nobc", "combat", "harmony", "limma", "MMUPHin", "ConQuR_rel", "percentile_norm"]
# datasets = ["ibd_3_CMD", "crc_8_CMD"]
# output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
# ylim_dict = {'batch_r2': [[0, 0.14], [0, 0.14]], 'biovar_r2': [[0, 0.030], [0, 0.030]], 'shannon_pval': [[0, 0.1], [0, 1]], 'auc_f1': [[0.4, 1], [0.4, 1]], 'runtime': [[np.log2(5e-5), np.log2(2500)]], 'precision_recall': [[0, 1], [0, 1]]}
# alt_visualize_simulation_stats('/athena/linglab/scratch/chf4012/mic_bc_benchmark/outputs/rw_data_alt_plots/relab_rw', output_dir_l, datasets, methods, highlighted_method = "ConQuR_rel", count_l = [False, False], simulate = False, postfix = '.pdf', ylim_dict = ylim_dict)

current_path = os.path.dirname(os.path.abspath(__file__))
print("current_path", current_path)
overall_path = current_path + "/../.."
with open(f'{current_path}/../config.yml') as file:
    config_data = yaml.load(file, Loader=yaml.FullLoader)
# ## VISUALIZE LINE PLOTS FOR 2 COUNT-TYPE RW DATASETS and 2 RELAB-TYPE RW DATASETS
##############################################################################
# create dir if not exists
if not os.path.exists(f'{config_data["evaluation_outputs"]}/rw_data_plots'):
    os.makedirs(f'{config_data["evaluation_outputs"]}/rw_data_plots')
if not os.path.exists(f'{config_data["evaluation_outputs"]}/rw_data_plots/count_rw'):
    os.makedirs(f'{config_data["evaluation_outputs"]}/rw_data_plots/count_rw')
if not os.path.exists(f'{config_data["evaluation_outputs"]}/rw_data_plots/relab_rw'):
    os.makedirs(f'{config_data["evaluation_outputs"]}/rw_data_plots/relab_rw')
methods = ["nobc", "harmony", "combat_seq", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
datasets = ["autism_2_microbiomeHD", "cdi_3_microbiomeHD"]
output_dir_l = [config_data['evaluation_outputs']+'/'+dataset for dataset in datasets]
ylim_dict = {'batch_r2': [[0, 0.35], [0, 0.35]], 'biovar_r2': [[0, 0.2], [0, 0.2]], 'shannon_pval': [[0, 1], [0, 1]], 'auc_f1': [[0.2, 0.9], [0.2, 0.9]], 'runtime': [[np.log2(5e-5), np.log2(500)]], 'precision_recall': [[0, 1], [0, 1]]}
alt_visualize_simulation_stats(f'{config_data["evaluation_outputs"]}/rw_data_plots/count_rw', output_dir_l, datasets, methods, highlighted_method = "ConQuR",  simulate = False, count_l = [True, True], postfix = '.pdf', ylim_dict = ylim_dict)

methods = ["nobc", "combat", "harmony", "limma", "MMUPHin", "ConQuR_rel", "percentile_norm"]
datasets = ["ibd_3_CMD", "crc_8_CMD"]
output_dir_l = [config_data['evaluation_outputs']+'/'+dataset for dataset in datasets]
alt_visualize_simulation_stats(f'{config_data["evaluation_outputs"]}/rw_data_plots/relab_rw', output_dir_l, datasets, methods, highlighted_method = "ConQuR_rel", count_l = [False, False], simulate = False, postfix = '.pdf', ylim_dict = ylim_dict)
