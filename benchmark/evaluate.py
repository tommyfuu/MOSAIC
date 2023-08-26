from preprocessing import *
from percentile_norm_func import *
from harmony import run_harmony
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
from sklearn.decomposition import PCA
import seaborn as sns
from sklearn.preprocessing import StandardScaler
import itertools
from scipy import stats
from skbio.diversity import alpha_diversity
import os
import time
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, KFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, f1_score, precision_score, recall_score, RocCurveDisplay, roc_curve
from sklearn.preprocessing import OneHotEncoder
import time
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.packages import importr
from statsmodels.stats.multitest import multipletests


# data_mat, meta_data = load_data(address_X, address_Y, IDCol, index_col, PCA_first = False)
def generate_harmony_results(data_mat, meta_data, IDCol, vars_use, output_root):
    start_time = time.time()
    ho = run_harmony(data_mat, meta_data, vars_use)
    res = pd.DataFrame(ho.Z_corr)
    elapsed = time.time() - start_time
    print("time")
    print(elapsed)
    # give the index back
    res.index = data_mat.columns
    res.columns = list(meta_data[IDCol])
    
    with open(output_root+"_elapsed_time.txt", "w") as text_file:
        print(str(elapsed), "seconds", file=text_file)
        print("\n", file=text_file)

    res.T.to_csv(output_root+"_adjusted_count.csv")
    return res.T, meta_data


def run_eval(batch_corrected_df, meta_data, batch_var, output_root, bio_var = False, n_pc=30, covar = False, IDCol = None):
    a = Evaluate(batch_corrected_df, meta_data, batch_var, output_root, bio_var, n_pc, covar, IDCol)
    return

def plot_PCOA_multiple(dataset_name, batch_corrected_df_l, methods, meta_data_l, used_var, output_root, datatype = 'count'):
    # prepare metadata
    r = robjects.r
    numpy2ri.activate()
    pandas2ri.activate()
    grdevices = importr('grDevices')
    r.source('./PERMANOVA_supporting.R')
    # r_used_var = meta_data[used_var]
    # r_bio_var = meta_data[bio_var]

    # initial graphics
    
    r.pdf(output_root+dataset_name+"_multi_PCOA_both_batch.pdf", len(batch_corrected_df_l)*6, height=12)
    if datatype == 'count':
        r.par(mfrow = robjects.IntVector([2, len(batch_corrected_df_l)]))
    else:
        r.par(mfrow = robjects.IntVector([1, len(batch_corrected_df_l)]))

    # plot the subplots in order
    if datatype == 'count':
        for idx, batch_corrected_df in enumerate(batch_corrected_df_l):
            data = np.array(batch_corrected_df)
            data = np.where(data<np.percentile(data.flatten(), 0.01), 0, data)
            data = data+np.abs(np.min(data))
            r_used_var = meta_data_l[idx][used_var]
            r.Plot_single_PCoA(data, r_used_var, dissimilarity="Aitch", bc_method = methods[idx])
    
    for idx, batch_corrected_df in enumerate(batch_corrected_df_l):
        data = np.array(batch_corrected_df)
        data = np.where(data<np.percentile(data.flatten(), 0.01), 0, data)
        data = data+np.abs(np.min(data))
        r_used_var = meta_data_l[idx][used_var]
        r.Plot_single_PCoA(data, r_used_var, dissimilarity="Bray", bc_method = methods[idx])
    grdevices.dev_off()
    return

class Evaluate(object):
    def __init__(
            self, batch_corrected_df, meta_data, batch_var, output_root, bio_var = False, n_pc=30, covar_l = [], IDCol = None, test_percent = 0.2, poslabel = '', pipeline = 'default', taxa_gt = None, datatype = "count"
    ):
        self.batch_corrected_df = batch_corrected_df
        self.meta_data = meta_data
        self.batch_var = batch_var
        self.output_root = output_root
        self.bio_var = bio_var # var retaining the biological info
        self.n_pc = n_pc
        self.covar_l = covar_l
        self.rng = np.random.RandomState(100)
        self.IDCol = IDCol
        self.test_percent = test_percent
        self.poslabel = poslabel
        self.datatype = datatype
        self.taxa_gt = taxa_gt # a list of taxonomy indices for ground truth positive taxa for differential abundance
        self.pipeline = pipeline # "scaled" or "default"
        ## note: for pipelines
        # 1. scaled pipeline generates:
        # 1.1 _summary.csv: contains summary stats for biovar/batch aitch/bray r2 and pval + biovar/batch shannon pval
        # 1.2  _diff_abund_test.csv: contains diff abund test results for biovar/batch concatenated to taxa level data
        # 2. default pipeline additionally generates,
        # 2.1 4* PCoA_Aitch/Bray.pdf: PCoA plots for Aitch/Bray dissimilarity
        # 2.2 5* _rf_fold_k_roc.pdf: ROC plots for concatenated 5F CV prediction of biovar states
        # 2.3 2* _shannon_df_batch/biovar.pdf: Shannon diversity indices with groupings for batch/biovar
        self.short_summary_dict = {}
        # make directory if not exists
        directory_path = "/".join(output_root.split("/")[:-1])
        if not os.path.exists(directory_path):
            os.makedirs(directory_path)
        # functions executed
        self.diff_abund_test()
        # self.alpha_diversity_and_tests(self.batch_var) # has to happen before standard scaler
        self.alpha_diversity_and_tests_for_batches_and_biovars()
        self.predict_difference_RF()
        self.calculate_R_sq()
        ## TODO: add FP stuff for simulation

        # convert summary dict to df
        self.summary_df = pd.DataFrame.from_dict(self.short_summary_dict, orient='index')
        print(self.summary_df)
        self.summary_df.to_csv(self.output_root+"_summary.csv")

    def diff_abund_test(self):
        ### if two categories then wilcox test, three or more then kruskal-wallis test
        df = self.batch_corrected_df
        meta_data = self.meta_data
        ### if we know the ground truth, so we can calculate the stats such as FDR, TPR, FPR, sensitivity, etc
        if self.taxa_gt is not None:
            gt_dict = {}
        # use the metadata to get indices of samples in each bio_var category
        bio_var_l = list(meta_data[self.bio_var].values)
        # get the indices for each category
        bio_var_indices = {}
        print("unique biovar variables", np.unique(bio_var_l))
        for var in np.unique(bio_var_l):
            bio_var_indices[var] = np.where(np.array(bio_var_l) == var)[0].tolist()
        p_values_by_taxa = {}
        for idx, taxa in enumerate(df.columns):
            current_taxa_values_by_bio_var = {}
            current_taxa_values = df[taxa].values.tolist()
            for var in np.unique(bio_var_l):
                current_taxa_values_by_bio_var[var] = [current_taxa_values[i] for i in bio_var_indices[var]]
            
            # conduct statistical test
            if len(np.unique(bio_var_l)) == 2:
                # wilcox test (unpaired)
                test_type = 'wilcox'
                p = stats.ranksums(current_taxa_values_by_bio_var[np.unique(bio_var_l)[0]], current_taxa_values_by_bio_var[np.unique(bio_var_l)[1]])[1]
            else:
                # kruskal-wallis test
                test_type = 'kruskal'
                try:
                    p = stats.kruskal(*current_taxa_values_by_bio_var.values())[1]
                except ValueError:
                    p = 1
            p_values_by_taxa[taxa] = p

            ### if we know the ground truth, so we can calculate the stats such as FDR, TPR, FPR, sensitivity, etc
            if self.taxa_gt is not None:
                gt_dict[taxa] = 1 if idx in self.taxa_gt else 0
    

        # transpose the batch corrected dataframe and save the p-values alongside
        df = df.T
        df[test_type+"_p_value"] = [p_values_by_taxa[taxa] for taxa in df.index]
        # add FDR corrected p-values
        df["FDR_p_value"] = multipletests(df[test_type+"_p_value"].values, method='fdr_bh')[1]
        # add ground truth
        if self.taxa_gt is not None:
            df["ground_truth"] = [gt_dict[taxa] if taxa in gt_dict.keys() else 0 for taxa in df.index]
            # calculate stats
            TP = len(df[(df["ground_truth"] == 1) & (df["FDR_p_value"] < 0.05)])
            print(df["ground_truth"])
            FP = len(df[(df["ground_truth"] == 0) & (df["FDR_p_value"] < 0.05)])
            TN = len(df[(df["ground_truth"] == 0) & (df["FDR_p_value"] >= 0.05)])
            FN = len(df[(df["ground_truth"] == 1) & (df["FDR_p_value"] >= 0.05)])
            print("TP, FP, TN, FN", TP, FP, TN, FN)

            # calculate stats
            if TP+FP == 0:
                self.short_summary_dict["FDR"] = 0
            else:
                self.short_summary_dict["FDR"] = FP/(FP+TP)
            if TP+FN == 0:
                self.short_summary_dict["power"] = 0
            else:
                self.short_summary_dict["power"] = TP/(TP+FN)


        # save the dataframe
        df.to_csv(self.output_root+"_diff_abund_test.csv")
        return

    def calculate_R_sq(self):
        data = np.array(self.batch_corrected_df)
        data = np.where(data<np.percentile(data.flatten(), 0.01), 0, data)
        data = data+np.abs(np.min(data))

        # attempting rpy2
        #Must be activated
        pandas2ri.activate()
        numpy2ri.activate()

        # import r packages/functions
        r = robjects.r
        r.source('./PERMANOVA_supporting.R')
        r_batch = self.meta_data[self.batch_var]

        if self.covar_l != False:
            r_covariate = self.meta_data[self.covar_l]
            PERMANOVA_R2_results = r.PERMANOVA_R2(data, r_batch, r_covariate, 'batch', self.covar_l)
            idx_l = ['batch'] + self.covar_l
            aitchinson_permanova_df = pd.DataFrame(PERMANOVA_R2_results.rx2("tab_rel"), columns=["standard", "sqrt.dist=T", "add=T", 'Pr(>F)'], index=idx_l)
            bray_curtis_permanova_df = pd.DataFrame(PERMANOVA_R2_results.rx2("tab_count"), columns=["standard", "sqrt.dist=T", "add=T", 'Pr(>F)'], index=idx_l)
        else:
            PERMANOVA_R2_results = r.PERMANOVA_R2(data, r_batch, robjects.NULL, 'batch')
            aitchinson_permanova_df = pd.DataFrame(PERMANOVA_R2_results.rx2("tab_rel"), columns=["standard", "sqrt.dist=T", "add=T", 'Pr(>F)'], index=['batch'])
            bray_curtis_permanova_df = pd.DataFrame(PERMANOVA_R2_results.rx2("tab_count"), columns=["standard", "sqrt.dist=T", "add=T", 'Pr(>F)'], index=['batch'])

        print("______________batch_PERMANOVA_R2_results______________")
        print(aitchinson_permanova_df)
        print(bray_curtis_permanova_df)
        print("______________batch_PERMANOVA_R2_results______________")
        self.short_summary_dict["batch_aitch_r2"] = aitchinson_permanova_df.loc["batch", "standard"]
        self.short_summary_dict["batch_bray_r2"] = bray_curtis_permanova_df.loc["batch", "standard"]
        self.short_summary_dict["batch_aitch_pval"] = aitchinson_permanova_df.loc["batch", 'Pr(>F)']
        self.short_summary_dict["batch_bray_pval"] = bray_curtis_permanova_df.loc["batch", 'Pr(>F)']
        # calculate the equivalent stats but for biological variable
        r_bio_var = self.meta_data[self.bio_var]
        if self.covar_l != []:
            r_covariate = self.meta_data[self.covar_l]
            PERMANOVA_R2_results = r.PERMANOVA_R2(data, r_bio_var, r_covariate, 'biovar', self.covar_l)
            idx_l = ['batch'] + self.covar_l
            aitchinson_permanova_df = pd.DataFrame(PERMANOVA_R2_results.rx2("tab_rel"), columns=["standard", "sqrt.dist=T", "add=T", 'Pr(>F)'], index=idx_l)
            bray_curtis_permanova_df = pd.DataFrame(PERMANOVA_R2_results.rx2("tab_count"), columns=["standard", "sqrt.dist=T", "add=T", 'Pr(>F)'], index=idx_l)
        else:
            PERMANOVA_R2_results = r.PERMANOVA_R2(data, r_bio_var, robjects.NULL, 'biovar')
            aitchinson_permanova_df = pd.DataFrame(PERMANOVA_R2_results.rx2("tab_rel"), columns=["standard", "sqrt.dist=T", "add=T", 'Pr(>F)'], index=['batch'])
            bray_curtis_permanova_df = pd.DataFrame(PERMANOVA_R2_results.rx2("tab_count"), columns=["standard", "sqrt.dist=T", "add=T", 'Pr(>F)'], index=['batch'])
        print("______________biovar_PERMANOVA_R2_results______________")
        print(bray_curtis_permanova_df)
        print(aitchinson_permanova_df)
        # bray_anova = PERMANOVA_R2_results[0]
        print("______________biovar_PERMANOVA_R2_results______________")
        self.short_summary_dict["biovar_aitch_r2"] = aitchinson_permanova_df.loc["batch", "standard"]
        self.short_summary_dict["biovar_bray_r2"] = bray_curtis_permanova_df.loc["batch", "standard"]
        self.short_summary_dict["biovar_aitch_pval"] = aitchinson_permanova_df.loc["batch", 'Pr(>F)']
        self.short_summary_dict["biovar_bray_pval"] = bray_curtis_permanova_df.loc["batch", 'Pr(>F)']

        # try plotting stuff
        if self.pipeline == "default":
            if self.datatype == "count":
                r.Plot_PCoA(self.output_root+'_batch', data, r_batch, dissimilarity="Aitch", main="Aitchinson")
            r.Plot_PCoA(self.output_root+'_batch', data, r_batch, dissimilarity="Bray", main="Bray-Curtis")
            r.Plot_PCoA(self.output_root+'_biovar', data, r_bio_var, dissimilarity="Aitch", main="Aitchinson")
            r.Plot_PCoA(self.output_root+'_biovar', data, r_bio_var, dissimilarity="Bray", main="Bray-Curtis")
        return

    
    def alpha_diversity_and_tests(self, test_var):
        print("______________alpha_diversity_and_tests______________")
        
        data = np.array(self.batch_corrected_df)
        data = np.where(data<np.percentile(data.flatten(), 0.01), 0, data)
        data = data+np.abs(np.min(data))
        ids = list(self.meta_data[self.IDCol])
        shannon_div = alpha_diversity('shannon', data, ids)
        

        # visualize and save
        shannon_df = shannon_div.to_frame()
        shannon_df.columns = ["shannon"]
        shannon_df[test_var] = list(self.meta_data[test_var])
        if self.pipeline == "default":
            shannon_fig = shannon_df.boxplot(column='shannon', by=test_var)
            shannon_fig.figure.savefig(self.output_root+"_"+test_var+"_alpha_shannon_boxplot.png")

        # get list of unique test_var options
        test_var_col_l = list(self.meta_data[test_var])
        test_var_l = list(np.unique(self.meta_data[test_var].dropna()))

        # conduct statistical tests
        ## 1. alpha diversity
        alpha_kruskal_data = []
        for var in test_var_l:
            removed_indices = [i for i, e in enumerate(test_var_col_l) if e != var]
            removed_indices = [var for i, var in enumerate(list(self.meta_data[self.IDCol])) if i in removed_indices]
            current_data_condition = list(shannon_div.drop(index = removed_indices))
            alpha_kruskal_data.append(current_data_condition)

        # calculate global kw p_val for alpha diversity
        if len(alpha_kruskal_data) == 2:
            # wilcox test (unpaired)
            shannon_global_pval = stats.ranksums(alpha_kruskal_data[0], alpha_kruskal_data[1])[1]
        else:
            shannon_global_pval = stats.kruskal(*alpha_kruskal_data)[1]

        return shannon_global_pval, shannon_df
    
    def alpha_diversity_and_tests_for_batches_and_biovars(self):
        print("______________alpha_diversity_and_tests_for_batches_and_biovars______________")
        shannon_global_pval_batch, shannon_df_batch = self.alpha_diversity_and_tests(self.batch_var)
        shannon_global_pval_biovar, shannon_df_biovar = self.alpha_diversity_and_tests(self.bio_var)

        # save to csv
        if self.pipeline == 'default':
            shannon_df_batch.to_csv(self.output_root+"_shannon_df_batch.csv")
            shannon_df_biovar.to_csv(self.output_root+"_shannon_df_biovar.csv")

        # save to summary dict
        self.short_summary_dict["batches_shannon_pval"] = shannon_global_pval_batch
        self.short_summary_dict["biovar_shannon_pval"] = shannon_global_pval_biovar
        return
    
    def predict_difference_RF(self):
        # get the data
        used_x = self.batch_corrected_df.copy() # note that standard scaler is already conducted
        used_y = list(self.meta_data[self.bio_var])
        sampld_ids = list(self.meta_data[self.IDCol])
        # Creating the Training and Test set from data
        # X_train, X_test, y_train, y_test = train_test_split(used_x, used_y, test_size = self.test_percent, random_state = self.rng)
        kf = KFold(n_splits=5, random_state=self.rng, shuffle = True)# Define the split - into n_splits folds 
        kf.get_n_splits(used_x, used_y) # returns the number of splitting iterations in the cross-validator

        # train random forest classifier and evaluate
        used_x.reset_index(drop=True, inplace=True)
                
        if self.pipeline == 'default':
            tprs = []
            base_fpr = np.linspace(0, 1, 101)
            plt.figure(figsize=(5, 5))
            plt.axes().set_aspect('equal', 'datalim')
        
        y_test_concat = []
        y_pred_concat = []
        y_pred_prob_concat = []
        IDs_concat = []
        for i, (train_index, test_index) in enumerate(kf.split(used_x, used_y)):
            # get training and testing data
            X_train = used_x.loc[train_index]
            y_train = [used_y[idx] for idx in train_index]
            X_test = used_x.loc[test_index]
            y_test = [used_y[idx] for idx in test_index]
            
            # fit rf model and evaluate
            clf = RandomForestClassifier(n_estimators=100, max_depth=2, random_state=0)
            clf.fit(X_train, y_train)
            y_pred = clf.predict(X_test)
            y_pred_prob = clf.predict_proba(X_test)

            # save to concatenated list
            y_test_concat.extend(y_test)
            y_pred_concat.extend(y_pred)
            y_pred_prob_concat.extend(y_pred_prob)
            IDs_concat.extend([sampld_ids[idx] for idx in test_index])
        average_acc = sum([y_pred_concat[i] == y_test_concat[i] for i in range(len(y_pred_concat))])/len(y_test_concat)

        macro_precision = precision_score(y_test_concat, y_pred_concat, average = 'macro')
        weighted_precision = precision_score(y_test_concat, y_pred_concat, average = 'weighted')
        macro_recall = recall_score(y_test_concat, y_pred_concat, average = 'macro')
        weighted_recall = recall_score(y_test_concat, y_pred_concat, average = 'weighted')
        macro_f1 = f1_score(y_test_concat, y_pred_concat, average = 'macro')
        weighted_f1 = f1_score(y_test_concat, y_pred_concat, average = 'weighted')
        
        enc = OneHotEncoder(handle_unknown='ignore')
        enc.fit(np.array(y_test_concat).reshape(-1, 1))
        y_test_oh = enc.transform(np.array(y_test_concat).reshape(-1, 1)).toarray()
        y_pred_oh = enc.transform(np.array(y_pred_concat).reshape(-1, 1)).toarray()
        auc = roc_auc_score(y_test_oh,  y_pred_oh)
        # find the most common element in y_test_concat
        most_common_element = max(set(y_test_concat), key = y_test_concat.count)
        baseline_likelihood = y_test_concat.count(most_common_element)/len(y_test_concat)

        self.short_summary_dict["rf_baseline_likelihood"] = baseline_likelihood
        self.short_summary_dict["rf_average_acc"] = average_acc
        self.short_summary_dict["rf_macro_precision"] = macro_precision
        self.short_summary_dict["weighted_precision"] = weighted_precision
        self.short_summary_dict["macro_recall"] = macro_recall
        self.short_summary_dict["weighted_recall"] = weighted_recall
        self.short_summary_dict["macro_f1"] = macro_f1
        self.short_summary_dict["weighted_f1"] = weighted_f1
        self.short_summary_dict["auc"] = auc

        # for each fold, also plot the roc plot
        # one hot encode y_test and y_pred
        y_test_zeroone = np.where(np.array(y_test) == most_common_element, 1, 0)
        y_pred_zeroone = np.where(np.array(y_pred) == most_common_element, 1, 0)
            
        if self.pipeline == 'default':
            RocCurveDisplay.from_predictions(
                y_test_zeroone, y_pred_zeroone,
                name=f"{self.bio_var} classification",
                color="darkorange",
            )
            plt.plot([0, 1], [0, 1], "k--", label="chance level (AUC = 0.5)")
            plt.axis("square")
            plt.xlabel("False Positive Rate")
            plt.ylabel("True Positive Rate")
            plt.title(f"{self.bio_var} classification concatenated over 5 folds")
            plt.legend()
            plt.savefig(self.output_root+"_"+self.bio_var+"_rf_5fold_roc.png")

        eval_df = pd.DataFrame({"Sam_ID": IDs_concat, "y_true": y_test_concat, "y_pred": y_pred_concat})
        eval_df.to_csv(self.output_root+"_"+self.bio_var+"_rf_evaluate.csv")
        return
        

def global_eval_dataframe(input_frame_path, bio_var, dataset_name, methods_list, methods_output_dir, output_dir_path = ".", taxa_gt = None, simulate = False, datatype = "count"):
    # fetch dimension (number of samples and number of features)
    data_mat = pd.read_csv(input_frame_path, index_col=0)

    # rest of the stats
    method_dict = {}
    for method in methods_list:
        print(method)
        current_method_dict = {}

        # methods_output_dir = overall_path+'/simulation_data_output_small_072623/'
        # output_dir_path = overall_path+"/simulation_data_eval_small_072623/out_"+str(odds_ratio)+"_"+ str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter))


        if not simulate:
            current_dir_path = output_dir_path + "/output_" + dataset_name+"_"+method 
        else:
            current_dir_path = output_dir_path + '/' + method
        print(current_dir_path)
        current_dir_file_names = os.listdir(current_dir_path)

        # fetch stats available in the summary file: batch/biovar PERMANOVA R2 in both Aitchinson and Bray-curtis + Shannon pval
        aitch_method_perm_path = [result for result in current_dir_file_names if "_summary" in result][0]
        summary_df = pd.read_csv(current_dir_path+'/'+aitch_method_perm_path)
        summary_dict = dict(zip(summary_df["Unnamed: 0"], summary_df["0"]))
        aitch_r2_batch = summary_dict["batch_aitch_r2"]
        aitch_r2_biovar = summary_dict["biovar_aitch_r2"]
        bray_r2_batch = summary_dict["batch_bray_r2"]
        bray_r2_biovar = summary_dict["biovar_bray_r2"]
        shannon_pval_batch = summary_dict["batches_shannon_pval"]
        shannon_pval_biovar = summary_dict["biovar_shannon_pval"]
        aitch_pval_batch = summary_dict["batch_aitch_pval"]
        aitch_pval_biovar = summary_dict["biovar_aitch_pval"]
        bray_pval_batch = summary_dict["batch_bray_pval"]
        bray_pval_biovar = summary_dict["biovar_bray_pval"]
        # rf stuff
        baseline_likelihood = summary_dict["rf_baseline_likelihood"]
        average_acc = summary_dict["rf_average_acc"]
        macro_precision = summary_dict["rf_macro_precision"]
        weighted_precision = summary_dict["weighted_precision"]
        macro_recall = summary_dict["macro_recall"]
        weighted_recall = summary_dict["weighted_recall"]
        macro_f1 = summary_dict["macro_f1"]
        weighted_f1 = summary_dict["weighted_f1"]
        auc = summary_dict["auc"]

        if taxa_gt is not None:
            FDR = summary_dict["FDR"]
            power = summary_dict["power"]

        if datatype == 'count':
            current_method_dict.update({"batch_aitch_r2": aitch_r2_batch})
            current_method_dict.update({"biovar_aitch_r2": aitch_r2_biovar})
            current_method_dict.update({"aitch_pval_batch": aitch_pval_batch})
            current_method_dict.update({"aitch_pval_biovar": aitch_pval_biovar})
        current_method_dict.update({"batch_bray_r2": bray_r2_batch})
        current_method_dict.update({"biovar_bray_r2": bray_r2_biovar})
        current_method_dict.update({"shannon_pval": shannon_pval_batch})
        current_method_dict.update({"shannon_pval": shannon_pval_biovar})
        current_method_dict.update({"bray_pval_batch": bray_pval_batch})
        current_method_dict.update({"bray_pval_biovar": bray_pval_biovar})
        current_method_dict.update({"rf_baseline_likelihood": baseline_likelihood})
        current_method_dict.update({"rf_average_acc": average_acc})
        current_method_dict.update({"rf_macro_precision": macro_precision})
        current_method_dict.update({"rf_weighted_precision": weighted_precision})
        current_method_dict.update({"rf_macro_recall": macro_recall})
        current_method_dict.update({"rf_weighted_recall": weighted_recall})
        current_method_dict.update({"rf_macro_f1": macro_f1})
        current_method_dict.update({"rf_weighted_f1": weighted_f1})
        current_method_dict.update({"rf_auc": auc})
        
        if taxa_gt is not None:
            current_method_dict.update({"FDR": FDR})
            current_method_dict.update({"power": power})

        method_dict[method] = current_method_dict

    # fetch time spent running
    if not simulate:
        benchmarked_results_dir = methods_output_dir + "/" + dataset_name
    else:
        benchmarked_results_dir = methods_output_dir
    benchmarked_results_dir_files = os.listdir(benchmarked_results_dir)
    print(benchmarked_results_dir_files)

    current_runtime_kw_path = [result for result in benchmarked_results_dir_files if "runtime.txt" in result][0]
    with open(benchmarked_results_dir+'/'+current_runtime_kw_path) as f:
        lines = f.readlines()

    # print(method_dict)
    for line in lines:
        if "combat" in line and "combat" in method_dict.keys():
            method_dict["combat"]["runtime"] = float(line.split(" ")[-2])
        # if "combat_seq" in line:
        #     method_dict["combat_seq"]["runtime"] = float(line.split(" ")[-2])
        if "limma" in line and "limma" in method_dict.keys():
            method_dict["limma"]["runtime"] = float(line.split(" ")[-2])
        if "MMUPHin" in line and "MMUPHin" in method_dict.keys():
            method_dict["MMUPHin"]["runtime"] = float(line.split(" ")[-2])
        if "ConQuR" in line and "Tune" not in line and "ConQuR" in method_dict.keys():
            method_dict["ConQuR"]["runtime"] = float(line.split(" ")[-2])
        if "ConQuR_libsize" in line and "ConQuR_libsize" in method_dict.keys():
            method_dict["ConQuR_libsize"]["runtime"] = float(line.split(" ")[-2])
        if "Tune_ConQuR runtime" in line and "Tune_ConQuR" in method_dict.keys():
            method_dict["Tune_ConQuR"]["runtime"] = float(line.split(" ")[-2])
        if "Tune_ConQuR_libsize" in line and "Tune_ConQuR_libsize" in method_dict.keys():
            method_dict["Tune_ConQuR_libsize"]["runtime"] = float(line.split(" ")[-2])
        if "ConQuR_rel runtime" in line and "Tune" not in line and "ConQuR_rel" in method_dict.keys():
            method_dict["ConQuR_rel"]["runtime"] = float(line.split(" ")[-2])
        if "Tune_ConQuR_rel" in line and "Tune_ConQuR_rel" in method_dict.keys():
            method_dict["Tune_ConQuR_rel"]["runtime"] = float(line.split(" ")[-2])

    if not simulate:
        benchmarked_data_harmony_dir = benchmarked_results_dir
    else:
        benchmarked_data_harmony_dir = output_dir_path + "/harmony"
    benchmarked_data_harmony_dir_files = os.listdir(benchmarked_data_harmony_dir)
    current_runtime_kw_paths = [result for result in benchmarked_data_harmony_dir_files if "elapsed_time.txt" in result]
    print("runtime files")
    print(benchmarked_data_harmony_dir)
    print(benchmarked_data_harmony_dir_files)
    print(current_runtime_kw_paths)
    for time_file in current_runtime_kw_paths:
        time_file = benchmarked_data_harmony_dir + "/"+time_file
        if "harmony_elapsed_time" in time_file:
            with open(time_file) as f:
                lines = f.readlines()
            method_dict["harmony"]["runtime"] = float(lines[0].split(" ")[-2])

    if not simulate:
        benchmarked_data_harmony_dir = benchmarked_results_dir
    else:
        benchmarked_data_harmony_dir = output_dir_path + "/percentile_norm"
    for time_file in current_runtime_kw_paths:
        time_file = benchmarked_data_harmony_dir + "/"+time_file
        if "percentile_norm_elapsed_time" in time_file:
            with open(time_file) as f:
                lines = f.readlines()
            method_dict["percentile_norm"]["runtime"] = float(lines[0].split(" ")[0])
    if 'nobc' in method_dict:
        method_dict['nobc']['runtime'] = 'NA'

    ## TODO: add FP stuff for simulation after those calculations are done
        
    # print into a pandas df where rows are datasets and columns and different stats
    results_df = pd.DataFrame.from_dict(method_dict, orient ='index') 
    results_df.T.to_csv(output_dir_path+"/global_benchmarking_stats_"+dataset_name+".csv")
    return pd.DataFrame.from_dict(method_dict, orient ='index') 


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
    print(df)
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
    
def visualize_simulation_stats(output_root, output_dir_l, datasets, methods, highlighted_method, simulate = False, sim_num_iters = 5, dimensions = (5, 5), taxa_gt = None, line = True, count_l = [True, True, False, False], marker_dict = {'nobc': 'o', "combat": 'v', "limma": '^', "MMUPHin": '<', "ConQuR": '>', "ConQuR_libsize": 's', "ConQuR_rel": 'p', "harmony":"+", "percentile_norm": "*"}):
    '''visualize the PERMANOVA batch R2 (Bray/Aitch), PERMANOVA condition R2 (Bray/Aitch), ROC-AUC and FDR/sensitivity'''
    # global set up
    global_methods_batch_bray_r2_l_dict = {}
    global_methods_batch_aitch_r2_l_dict = {}
    global_methods_biovar_bray_r2_l_dict = {}
    global_methods_biovar_aitch_r2_l_dict = {}
    if taxa_gt is not None:
        global_methods_FDR_r2_l_dict = {}
        global_methods_sensitivity_r2_l_dict = {}
    if line:
        linestyle = '-'
    else:
        linestyle = 'None'

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
        
        if taxa_gt is not None:
            FDR_r2 = global_stats_df.loc["FDR"]
            methods_FDR_r2_dict = FDR_r2.to_dict()
            sensitivity_r2 = global_stats_df.loc["power"]
            methods_sensitivity_r2_dict = sensitivity_r2.to_dict()
        if taxa_gt is not None:
            return methods_batch_aitch_r2_dict, methods_batch_bray_r2_dict, methods_biovar_aitch_r2_dict, methods_biovar_bray_r2_dict, methods_FDR_r2_dict, methods_sensitivity_r2_dict
        else:
            return methods_batch_aitch_r2_dict, methods_batch_bray_r2_dict, methods_biovar_aitch_r2_dict, methods_biovar_bray_r2_dict
    for idx, output_dir in enumerate(output_dir_l):
        if not simulate:   
            # # get the global stats path
            # files = os.listdir(output_dir)
            # global_stats_path = [file for file in files if "global_benchmarking_stats" in file][0]

            # # read the global stats df
            # global_stats_df = pd.read_csv(output_dir+"/"+global_stats_path, index_col=0)

            # # get the stats of interest and visualize
            # if count_l[idx]:
            #     batch_aitch_r2 = global_stats_df.loc["batch_aitch_r2"]
            #     methods_batch_aitch_r2_dict = batch_aitch_r2.to_dict()
            #     biovar_aitch_r2 = global_stats_df.loc["biovar_aitch_r2"]
            #     methods_biovar_aitch_r2_dict = biovar_aitch_r2.to_dict()
            #     if taxa_gt is not None:
            #         FDR_r2 = global_stats_df.loc["FDR"]
            #         methods_FDR_r2_dict = FDR_r2.to_dict()
            #         sensitivity_r2 = global_stats_df.loc["power"]
            #         methods_sensitivity_r2_dict = sensitivity_r2.to_dict()

            #     ## append to global dict
            #     for method in methods:
            #         if method not in global_methods_batch_aitch_r2_l_dict.keys():
            #             global_methods_batch_aitch_r2_l_dict[method] = []
            #         global_methods_batch_aitch_r2_l_dict[method].append(methods_batch_aitch_r2_dict[method])
            #         if method not in global_methods_biovar_aitch_r2_l_dict.keys():
            #             global_methods_biovar_aitch_r2_l_dict[method] = []
            #         global_methods_biovar_aitch_r2_l_dict[method].append(methods_biovar_aitch_r2_dict[method])
            #         if taxa_gt is not None:
            #             if method not in global_methods_FDR_r2_l_dict.keys():
            #                 global_methods_FDR_r2_l_dict[method] = []
            #             global_methods_FDR_r2_l_dict[method].append(methods_FDR_r2_dict[method])
            #             if method not in global_methods_sensitivity_r2_l_dict.keys():
            #                 global_methods_sensitivity_r2_l_dict[method] = []
            #             global_methods_sensitivity_r2_l_dict[method].append(methods_sensitivity_r2_dict[method])
            
            # batch_bray_r2 = global_stats_df.loc["batch_bray_r2"]  
            # methods_batch_bray_r2_dict = batch_bray_r2.to_dict()
            # biovar_bray_r2 = global_stats_df.loc["biovar_bray_r2"]
            # methods_biovar_bray_r2_dict = biovar_bray_r2.to_dict()
            # if taxa_gt is not None:
            #     FDR_r2 = global_stats_df.loc["FDR"]
            #     methods_FDR_r2_dict = FDR_r2.to_dict()
            #     sensitivity_r2 = global_stats_df.loc["power"]
            #     methods_sensitivity_r2_dict = sensitivity_r2.to_dict()

            if taxa_gt is not None:
                methods_batch_aitch_r2_dict, methods_batch_bray_r2_dict, methods_biovar_aitch_r2_dict, methods_biovar_bray_r2_dict, methods_FDR_r2_dict, methods_sensitivity_r2_dict = get_stats(output_dir, taxa_gt, idx)    
            else:
                methods_batch_aitch_r2_dict, methods_batch_bray_r2_dict, methods_biovar_aitch_r2_dict, methods_biovar_bray_r2_dict = get_stats(output_dir, taxa_gt, idx)    
            ## append to global dict
            for method in methods:
                if count_l[idx]:
                    if method not in global_methods_batch_aitch_r2_l_dict.keys():
                        global_methods_batch_aitch_r2_l_dict[method] = []
                    global_methods_batch_aitch_r2_l_dict[method].append(methods_batch_aitch_r2_dict[method])
                    if method not in global_methods_biovar_aitch_r2_l_dict.keys():
                        global_methods_biovar_aitch_r2_l_dict[method] = []
                    global_methods_biovar_aitch_r2_l_dict[method].append(methods_biovar_aitch_r2_dict[method])
                if method not in global_methods_batch_bray_r2_l_dict.keys():
                    global_methods_batch_bray_r2_l_dict[method] = []
                global_methods_batch_bray_r2_l_dict[method].append(methods_batch_bray_r2_dict[method])
                if method not in global_methods_biovar_bray_r2_l_dict.keys():
                    global_methods_biovar_bray_r2_l_dict[method] = []
                global_methods_biovar_bray_r2_l_dict[method].append(methods_biovar_bray_r2_dict[method])
                if taxa_gt is not None:
                    if method not in global_methods_FDR_r2_l_dict.keys():
                        global_methods_FDR_r2_l_dict[method] = []
                    global_methods_FDR_r2_l_dict[method].append(methods_FDR_r2_dict[method])
                    if method not in global_methods_sensitivity_r2_l_dict.keys():
                        global_methods_sensitivity_r2_l_dict[method] = []
                    global_methods_sensitivity_r2_l_dict[method].append(methods_sensitivity_r2_dict[method])

        else:
            # in this case, for each dataset, we have a preset number of iterations
            cross_iter_batch_aitch_r2_dict = {}
            cross_iter_batch_bray_r2_dict = {}
            cross_iter_biovar_aitch_r2_dict = {}
            cross_iter_biovar_bray_r2_dict = {}
            cross_iter_FDR_r2_dict = {}
            cross_iter_sensitivity_r2_dict = {}
            for iter in range(sim_num_iters):
                if taxa_gt is not None:
                    methods_batch_aitch_r2_dict, methods_batch_bray_r2_dict, methods_biovar_aitch_r2_dict, methods_biovar_bray_r2_dict, methods_FDR_r2_dict, methods_sensitivity_r2_dict = get_stats(output_dir+f"_iter_{iter+1}", taxa_gt, idx)    
                else:
                    methods_batch_aitch_r2_dict, methods_batch_bray_r2_dict, methods_biovar_aitch_r2_dict, methods_biovar_bray_r2_dict = get_stats(output_dir+f"_iter_{iter+1}", taxa_gt, idx)    
                cross_iter_batch_aitch_r2_dict[iter] = methods_batch_aitch_r2_dict
                cross_iter_batch_bray_r2_dict[iter] = methods_batch_bray_r2_dict
                cross_iter_biovar_aitch_r2_dict[iter] = methods_biovar_aitch_r2_dict
                cross_iter_biovar_bray_r2_dict[iter] = methods_biovar_bray_r2_dict
                if taxa_gt is not None:
                    cross_iter_FDR_r2_dict[iter] = methods_FDR_r2_dict
                    cross_iter_sensitivity_r2_dict[iter] = methods_sensitivity_r2_dict

            # calculate mean across iterations and append to global dict
            for method in methods:
                if count_l[idx]:
                    if method not in global_methods_batch_aitch_r2_l_dict.keys():
                        global_methods_batch_aitch_r2_l_dict[method] = []
                    global_methods_batch_aitch_r2_l_dict[method].append(np.mean([cross_iter_batch_aitch_r2_dict[iter][method] for iter in range(sim_num_iters)]))
                    if method not in global_methods_biovar_aitch_r2_l_dict.keys():
                        global_methods_biovar_aitch_r2_l_dict[method] = []
                    global_methods_biovar_aitch_r2_l_dict[method].append(np.mean([cross_iter_biovar_aitch_r2_dict[iter][method] for iter in range(sim_num_iters)]))
                if method not in global_methods_batch_bray_r2_l_dict.keys():
                    global_methods_batch_bray_r2_l_dict[method] = []
                global_methods_batch_bray_r2_l_dict[method].append(np.mean([cross_iter_batch_bray_r2_dict[iter][method] for iter in range(sim_num_iters)]))
                if method not in global_methods_biovar_bray_r2_l_dict.keys():
                    global_methods_biovar_bray_r2_l_dict[method] = []
                global_methods_biovar_bray_r2_l_dict[method].append(np.mean([cross_iter_biovar_bray_r2_dict[iter][method] for iter in range(sim_num_iters)]))
                if taxa_gt is not None:
                    if method not in global_methods_FDR_r2_l_dict.keys():
                        global_methods_FDR_r2_l_dict[method] = []
                    global_methods_FDR_r2_l_dict[method].append(np.mean([cross_iter_FDR_r2_dict[iter][method] for iter in range(sim_num_iters)]))
                    if method not in global_methods_sensitivity_r2_l_dict.keys():
                        global_methods_sensitivity_r2_l_dict[method] = []
                    global_methods_sensitivity_r2_l_dict[method].append(np.mean([cross_iter_sensitivity_r2_dict[iter][method] for iter in range(sim_num_iters)]))
                    
    print("global_methods_batch_aitch_r2_l_dict")
    print(global_methods_batch_aitch_r2_l_dict)
    print("global_methods_batch_bray_r2_l_dict")
    print(global_methods_batch_bray_r2_l_dict)

    def plot_stats(stats_summary_name, stats_name_l, stats_dict_1, stats_dict_2):
        ## plot the dictionaries in two matplotlib subplots as line plots
        plt.clf()
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(dimensions[0]*2, dimensions[1]))
        for method in methods:
            if method != highlighted_method:
                color = 'grey'
                alpha = 0.6
            else:
                color = 'red'
                alpha = 1
            if count_l[0]:
                ax1.plot(datasets, stats_dict_1[method], label=method, color=color, alpha=alpha, marker=marker_dict[method], linestyle=linestyle)
                ax1.set_title(stats_name_l[0])
                ax1.set_xticks(datasets)
                # set xticks to be verticle
                for tick in ax1.get_xticklabels():
                    tick.set_rotation(60)

            ax2.plot(datasets, stats_dict_2[method], label=method, color=color, alpha=alpha, marker=marker_dict[method], linestyle=linestyle)
            ax2.set_title(stats_name_l[1])
            ax2.set_xticks(datasets)
            # set xticks to be verticle
            for tick in ax2.get_xticklabels():
                tick.set_rotation(60)


        plt.subplots_adjust(right=0.8)
        ax2.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        # find parent dir of output_dir_l[0]
        parent_dir = os.path.dirname(output_dir_l[0])
        plt.savefig(parent_dir+"/"+output_root+"_"+stats_summary_name+".png", bbox_inches="tight")
        plt.clf()
        plt.close()
    plot_stats('PERMANOVA_batch_R2', ["PERMANOVA batch R2 (Aitchinson)", "PERMANOVA batch R2 (Bray-Curtis)"], global_methods_batch_aitch_r2_l_dict, global_methods_batch_bray_r2_l_dict)
    plot_stats('PERMANOVA_biovar_R2', ["PERMANOVA biovar R2 (Aitchinson)", "PERMANOVA biovar R2 (Bray-Curtis)"], global_methods_biovar_aitch_r2_l_dict, global_methods_biovar_bray_r2_l_dict)
    if taxa_gt is not None:
        plot_stats('FDR_sensitivity', ["FDR", "Sensitivity"], global_methods_FDR_r2_l_dict, global_methods_sensitivity_r2_l_dict)


    return

overall_path = '/athena/linglab/scratch/chf4012'

## GENERATE RW DATA FOR RUNNING R METHODS + RUN ON HARMONY/PERCENTILE_NORM
# # autism 2 microbiomeHD
# ################################################################################
# ### STEP 1. GENERATE DATA FROM DATABASE
# address_directory = overall_path+'/mic_bc_benchmark/data/autism_2_microbiomeHD'
# output_root = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/autism_2_microbiomeHD"
# data_mat, meta_data = load_data_microbiomeHD(address_directory)
# vars_use = ["Dataset"]
# IDCol = 'Sam_id'
# check_complete_confounding(meta_data, "Dataset", "DiseaseState", output_root)
# ### STEP 2. RUNNING METHODS (FOR THE TWO RUNNING IN PYTHON)
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/autism_2_microbiomeHD_count_data.csv"
# address_Y = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/autism_2_microbiomeHD_meta_data.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# res_h, meta_data_h = generate_harmony_results(data_mat, meta_data, IDCol, vars_use, overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/autism_2_microbiomeHD/"+"autism_2_microbiomeHD_harmony")
# percentile_norm(overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/autism_2_microbiomeHD_count_data.csv", overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/autism_2_microbiomeHD_meta_data.csv", "DiseaseState", "ASD", "comma", overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD")


# # cdi 3 microbiomeHD
# ################################################################################
# ### STEP 1. GENERATE DATA FROM DATABASE
# address_directory = overall_path+'/mic_bc_benchmark/data/cdi_3_microbiomeHD'
# output_root = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/cdi_3_microbiomeHD"
# data_mat, meta_data = load_data_microbiomeHD(address_directory, output_root)
# vars_use = ["Dataset"]
# IDCol = 'Sam_id'
# check_complete_confounding(meta_data, 'Dataset', "DiseaseState", output_root)
# ### STEP 2. RUNNING METHODS (FOR THE TWO RUNNING IN PYTHON)
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/cdi_3_microbiomeHD_count_data.csv"
# address_Y = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/cdi_3_microbiomeHD_meta_data.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# res_h, meta_data_h = generate_harmony_results(data_mat, meta_data, IDCol, vars_use, overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/cdi_3_microbiomeHD/"+"cdi_3_microbiomeHD_harmony")
# percentile_norm(overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/cdi_3_microbiomeHD_count_data.csv", overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/cdi_3_microbiomeHD_meta_data.csv", "DiseaseState", "CDI", "comma", overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD")


# # ibd_3_CMD
# ################################################################################
# ### STEP 1. GENERATE DATA FROM DATABASE
# address_directory = overall_path+'/mic_bc_benchmark/data/ibd_3_CMD'
# output_root = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/ibd_3_CMD"
# data_mat, meta_data = load_data_CMD(address_directory, output_root, covar_l=['age', 'gender'] )
# print("ibd_3_CMD loaded")
# vars_use = ["study_name"]
# IDCol = 'Sam_id'
# check_complete_confounding(meta_data, "study_name", "disease", output_root)
# ### STEP 2. RUNNING METHODS (FOR THE TWO RUNNING IN PYTHON)
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/ibd_3_CMD_count_data.csv"
# address_Y = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/ibd_3_CMD_meta_data.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# res_h, meta_data_h = generate_harmony_results(data_mat, meta_data, IDCol, vars_use, overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/ibd_3_CMD/"+"ibd_3_CMD_harmony")
# percentile_norm(overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/ibd_3_CMD_count_data.csv", overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/ibd_3_CMD_meta_data.csv", "disease", "IBD", "comma", overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/ibd_3_CMD/ibd_3_CMD")

# # CRC_8_CMD
# ################################################################################
# ### STEP 1. GENERATE DATA FROM DATABASE
# address_directory = overall_path+'/mic_bc_benchmark/data/CRC_8_CMD'
# output_root = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/CRC_8_CMD"
# data_mat, meta_data = load_data_CMD(address_directory, output_root, covar_l=['age', 'gender'] )
# print("CRC_8_CMD loaded")
# vars_use = ["study_name"]
# IDCol = 'Sam_id'
# check_complete_confounding(meta_data, "study_name", "disease", output_root)
# ### STEP 2. RUNNING METHODS (FOR THE TWO RUNNING IN PYTHON)
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/CRC_8_CMD_count_data.csv"
# address_Y = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/CRC_8_CMD_meta_data.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# res_h, meta_data_h = generate_harmony_results(data_mat, meta_data, IDCol, vars_use, overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/CRC_8_CMD/"+"CRC_8_CMD_harmony")
# percentile_norm(overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/CRC_8_CMD_count_data.csv", overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/CRC_8_CMD_meta_data.csv", "disease", "CRC", "comma", overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/CRC_8_CMD/CRC_8_CMD")





## simulation evaluation - MIDAS
## null data
# ## STEP 1. GENERATE DATA FROM DATABASE
# or_l = [1, 1.25, 1.5]
# cond_effect_val_l = [0, 0.099, 0.299, 0.499, 0.699, 0.899]
# batch_effect_val_l = [0, 0.099, 0.299, 0.499, 0.699, 0.899]
# num_iters = 5
# IDCol = 'subjectid_text'
# methods_list = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
# for odds_ratio in or_l:
#     for cond_effect_val in cond_effect_val_l:
#         for batch_effect_val in batch_effect_val_l:
#             if cond_effect_val + batch_effect_val <= 1:
#                 for iter in list(range(1, num_iters+1)):
#                     print("odds_ratio", odds_ratio)
#                     print("cond_effect_val", cond_effect_val)
#                     print("batch_effect_val", batch_effect_val)
#                     print("iter", iter)
                    # ## STEP 1. GENERATE DATA FROM DATABASE
                    # ### already done
                    # ### STEP 2. RUNNING METHODS (FOR THE TWO RUNNING IN PYTHON)
                    # address_X = overall_path+'/simulation_data_MIDAS_small_norelation_080723/ibd_150_count_' + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val) + '_iter_' + str(iter) + '.csv'
                    # address_Y = overall_path+'/simulation_data_MIDAS_small_norelation_080723/ibd_150_meta_' + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val) + '_iter_' + str(iter) + '.csv'
                    # data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
                    # output_root = overall_path+"/simulation_data_output_small_norelation_080723/out_"+str(odds_ratio)+"_"+ str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter) 
                    # if not os.path.exists(output_root):
                    #     os.makedirs(output_root)
                    # res_h, meta_data_h = generate_harmony_results(data_mat, meta_data, IDCol, ["batchid"],  output_root+"/ibd_"+str(odds_ratio)+ "_" +str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter) +"_harmony")
                    # percentile_norm(address_X, address_Y, "cond", "cond_1", "comma", output_root+"/ibd_"+str(odds_ratio)+ "_"+str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter), simulate=True)
                    # # THE ONE BELOW IS CORRECT BUT UNUSED
                    # # res_h, meta_data_h = generate_harmony_results(data_mat, meta_data, IDCol, ["batchid"],  output_root+"/ibd_"+str(odds_ratio)+str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter) +"_harmony")
                    # # percentile_norm(address_X, address_Y, "cond", "cond_1", "comma", output_root+"/ibd_"+str(odds_ratio)+str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter), simulate=True)

                    # address_X = overall_path+'/simulation_data_MIDAS_small_yesrelation_080723/ibd_150_count_' + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val) + '_iter_' + str(iter) + '.csv'
                    # address_Y = overall_path+'/simulation_data_MIDAS_small_yesrelation_080723/ibd_150_meta_' + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val) + '_iter_' + str(iter) + '.csv'
                    # data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
                    # output_root = overall_path+"/simulation_data_output_small_yesrelation_080723/out_"+str(odds_ratio)+"_"+ str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter) 
                    # if not os.path.exists(output_root):
                    #     os.makedirs(output_root)
                    # res_h, meta_data_h = generate_harmony_results(data_mat, meta_data, IDCol, ["batchid"],  output_root+"/ibd_"+str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter) +"_harmony")
                    # percentile_norm(address_X, address_Y, "cond", "cond_1", "comma", output_root+"/ibd_"+str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter), simulate=True)
                    # # THE ONE BELOW IS CORRECT BUT UNUSED
                    # # res_h, meta_data_h = generate_harmony_results(data_mat, meta_data, IDCol, ["batchid"],  output_root+"/ibd_"+str(odds_ratio)+str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter) +"_harmony")
                    # # percentile_norm(address_X, address_Y, "cond", "cond_1", "comma", output_root+"/ibd_"+str(odds_ratio)+str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter), simulate=True)

                    
                    # ### STEP 3. EVALUATE ALL THE METHODS
                    # if iter == 1:
                    #     pipeline = "default"
                    # else:
                    #     pipeline = "scaled"

                    # df_l = []
                    # metadata_l = []

                    # address_Y = overall_path+'/simulation_data_MIDAS_small_yesrelation_080723/ibd_150_meta_' + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val) + '_iter_' + str(iter) + '.csv'
                    # output_root = overall_path+'/simulation_data_eval_small_yesrelation_080723'
                    # # check if already exists
                    # output_dir = output_root+'/out_' + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val)+ '_iter_' + str(iter)
                    # if not os.path.exists(output_dir):
                    #     os.makedirs(output_dir)
                    # for method in methods_list:   
                    #     print("___________________________________________________________")
                        # print("odds_ratio", odds_ratio)
                        # print("cond_effect_val", cond_effect_val)
                        # print("batch_effect_val", batch_effect_val)
                        # print("iter", iter)
                    #     print("method", method)                     
                    #     # check if already exists
                    #     if not os.path.exists(output_dir+'/'+method):
                    #         os.makedirs(output_dir+'/'+method)
                        
                    #     if method == 'harmony':
                    #         # address_X = overall_path+"/simulation_data_output_small_yesrelation_080723/out_" + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val)+ '_iter_' + str(iter)+"/ibd_"+ str(odds_ratio)+ "_"+ str(cond_effect_val)+ "_"+ str(batch_effect_val)+ '_iter_'+ str(iter)+ "_harmony_adjusted_count.csv"
                    #         address_X = overall_path+"/simulation_data_output_small_yesrelation_080723/out_" + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val)+ '_iter_' + str(iter)+"/ibd_"+ str(cond_effect_val)+ "_"+ str(batch_effect_val)+ '_iter_'+ str(iter)+ "_harmony_adjusted_count.csv"
                    #     elif method == 'percentile_norm':
                    #         address_X = overall_path+"/simulation_data_output_small_yesrelation_080723/out_" + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val)+ '_iter_' + str(iter)+"/ibd_"+ str(cond_effect_val)+ "_"+ str(batch_effect_val)+ '_iter_'+ str(iter)+ "_percentile_norm.csv"
                    #     elif method == 'nobc':
                    #         address_X = overall_path+"/simulation_data_MIDAS_small_yesrelation_080723/ibd_150_count_"+ str(odds_ratio)+ "_"+ str(cond_effect_val)+ "_"+ str(batch_effect_val)+ '_iter_'+ str(iter)+ ".csv"
                    #     else:
                    #         address_X = overall_path+"/simulation_data_output_small_yesrelation_080723/out_" + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val)+ '_iter_' + str(iter)+"/ibd_"+ str(odds_ratio)+ "_"+ str(cond_effect_val)+ "_"+ str(batch_effect_val)+ '_iter_'+ str(iter)+ "_"+method+".csv"
                    #     data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
                    #     df_l.append(data_mat)
                    #     metadata_l.append(meta_data)

                    #     if not os.path.exists(output_dir+'/'+method+'/' + method +'_summary.csv'):
                    #         Evaluate(data_mat, meta_data, 'batchid', output_dir+'/'+method+'/'+method,
                    #                 "cond", 30, IDCol=IDCol, pipeline = pipeline, taxa_gt = list(range(150)))
                    #     else:
                    #         print("already exists", output_dir+'/'+method+'/out_' + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val) + '_iter_' + str(iter) + '_summary.csv')
                    # # TO FIX TOMORROW:
                    # input_frame_path = overall_path+"/simulation_data_MIDAS_small_yesrelation_080723/ibd_150_count_"+ str(odds_ratio)+ "_"+ str(cond_effect_val)+ "_"+ str(batch_effect_val)+ '_iter_'+ str(iter)+ ".csv"
                    # bio_var = "cond"
                    # dataset_name = "batch_0"
                    # global_eval_dataframe(input_frame_path, bio_var, dataset_name, methods_list, overall_path+"/simulation_data_output_small_yesrelation_080723/out_" + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val)+ '_iter_' + str(iter)+"/",
                    #  output_dir_path = overall_path+"/simulation_data_eval_small_yesrelation_080723/out_"+str(odds_ratio)+"_"+ str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter), simulate = True, taxa_gt = True)

                    # if not os.path.exists(overall_path+"/simulation_data_eval_small_yesrelation_080723/out_"+str(odds_ratio)+"_"+ str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter)+"/ibd_150_"+str(odds_ratio) + "_" + str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter)+'_multi_PCOA_both_batch.pdf'):
                    #     plot_PCOA_multiple("ibd_150_"+str(odds_ratio) + "_" + str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter), df_l, methods_list, metadata_l, used_var="batchid", output_root= overall_path+"/simulation_data_eval_small_yesrelation_080723/out_"+str(odds_ratio)+"_"+ str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter)+"/")

                    # continue

                    # ### STEP 3. EVALUATE ALL THE METHODS
                    # if iter == 1:
                    #     pipeline = "default"
                    # else:
                    #     pipeline = "scaled"

                    # df_l = []
                    # metadata_l = []

                    # address_Y = overall_path+'/simulation_data_MIDAS_small_norelation_080723/ibd_150_meta_' + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val) + '_iter_' + str(iter) + '.csv'
                    # output_root = overall_path+'/simulation_data_eval_small_norelation_080723'
                    # # check if already exists
                    # output_dir = output_root+'/out_' + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val)+ '_iter_' + str(iter)
                    # if not os.path.exists(output_dir):
                    #     os.makedirs(output_dir)
                    # for method in methods_list:   
                    #     print("___________________________________________________________")
                    #     print("odds_ratio", odds_ratio)
                    #     print("cond_effect_val", cond_effect_val)
                    #     print("batch_effect_val", batch_effect_val)
                    #     print("iter", iter)
                    #     print("method", method)                     
                    #     # check if already exists
                    #     if not os.path.exists(output_dir+'/'+method):
                    #         os.makedirs(output_dir+'/'+method)
                        
                    #     if method == 'harmony':
                    #         address_X = overall_path+"/simulation_data_output_small_norelation_080723/out_" + str(odds_ratio) + '_'+  str(cond_effect_val) + '_' + str(batch_effect_val)+ '_iter_' + str(iter)+"/ibd_"+ str(odds_ratio)+ str(cond_effect_val)+ "_"+ str(batch_effect_val)+ '_iter_'+ str(iter)+ "_harmony_adjusted_count.csv"
                    #         # address_X = overall_path+"/simulation_data_output_small_norelation_080723/out_" + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val)+ '_iter_' + str(iter)+"/ibd_"+ str(cond_effect_val)+ "_"+ str(batch_effect_val)+ '_iter_'+ str(iter)+ "_harmony_adjusted_count.csv"
                    #     elif method == 'percentile_norm':
                    #         address_X = overall_path+"/simulation_data_output_small_norelation_080723/out_" + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val)+ '_iter_' + str(iter)+"/ibd_"+ str(odds_ratio)+ str(cond_effect_val)+"_"+ str(batch_effect_val)+ '_iter_'+ str(iter)+ "_percentile_norm.csv"
                    #         # address_X = overall_path+"/simulation_data_output_small_norelation_080723/out_" + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val)+ '_iter_' + str(iter)+"/ibd_"+ str(odds_ratio)+ str(cond_effect_val)+ "_"+ str(batch_effect_val)+ '_iter_'+ str(iter)+ "_percentile_norm.csv"
                    #     elif method == 'nobc':
                    #         address_X = overall_path+"/simulation_data_MIDAS_small_norelation_080723/ibd_150_count_"+ str(odds_ratio)+ "_"+ str(cond_effect_val)+ "_"+ str(batch_effect_val)+ '_iter_'+ str(iter)+ ".csv"
                    #     else:
                    #         address_X = overall_path+"/simulation_data_output_small_norelation_080723/out_" + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val)+ '_iter_' + str(iter)+"/ibd_"+ str(odds_ratio)+ "_"+ str(cond_effect_val)+ "_"+ str(batch_effect_val)+ '_iter_'+ str(iter)+ "_"+method+".csv"
                    #     data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
                    #     df_l.append(data_mat)
                    #     metadata_l.append(meta_data)

                    #     if not os.path.exists(output_dir+'/'+method+'/' + method +'_summary.csv'):
                    #         Evaluate(data_mat, meta_data, 'batchid', output_dir+'/'+method+'/'+method,
                    #                 "cond", 30, IDCol=IDCol, pipeline = pipeline, taxa_gt = list(range(150)))
                    #     else:
                    #         print("already exists", output_dir+'/'+method+'/out_' + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val) + '_iter_' + str(iter) + '_summary.csv')
                    # # TO FIX TOMORROW:
                    # input_frame_path = overall_path+"/simulation_data_MIDAS_small_norelation_080723/ibd_150_count_"+ str(odds_ratio)+ "_"+ str(cond_effect_val)+ "_"+ str(batch_effect_val)+ '_iter_'+ str(iter)+ ".csv"
                    # bio_var = "cond"
                    # dataset_name = "batch_0"
                    # global_eval_dataframe(input_frame_path, bio_var, dataset_name, methods_list, overall_path+"/simulation_data_output_small_norelation_080723/out_" + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val)+ '_iter_' + str(iter)+"/",
                    #  output_dir_path = overall_path+"/simulation_data_eval_small_norelation_080723/out_"+str(odds_ratio)+"_"+ str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter), simulate = True, taxa_gt = True)

                    # if not os.path.exists(overall_path+"/simulation_data_eval_small_norelation_080723/out_"+str(odds_ratio)+"_"+ str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter)+"/ibd_150_"+str(odds_ratio) + "_" + str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter)+'_multi_PCOA_both_batch.pdf'):
                    #     plot_PCOA_multiple("ibd_150_"+str(odds_ratio) + "_" + str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter), df_l, methods_list, metadata_l, used_var="batchid", output_root= overall_path+"/simulation_data_eval_small_norelation_080723/out_"+str(odds_ratio)+"_"+ str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter)+"/")

                    # ### STEP 4. visualize simulation stats
                    # visualize_simulation_stats()
                    # continue


################################################################################
# # autism 2 microbiomeHD
# output_dir_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/outputs/autism_2_microbiomeHD'
# address_directory = overall_path+'/mic_bc_benchmark/data/autism_2_microbiomeHD'
# output_root = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/autism_2_microbiomeHD"
# vars_use = ["Dataset"]
# IDCol = 'Sam_id'

# address_Y = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/autism_2_microbiomeHD_meta_data.csv"

# # data_mat, meta_data = load_data_microbiomeHD(address_directory)
# # nobc
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/autism_2_microbiomeHD_count_data.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# ### nobc
# # Evaluate(data_mat, meta_data, 'Dataset', output_dir_path + '/output_autism_2_microbiomeHD_nobc/autism_2_microbiomeHD_nobc', "DiseaseState", 30, [], 'Sam_id')

# ### harmony
# # res_h, meta_data_h = generate_harmony_results(data_mat, meta_data, IDCol, vars_use, overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/autism_2_microbiomeHD/"+"harmony")
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD_harmony_adjusted_count.csv"
# res_h, meta_data_h = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(res_h, meta_data_h, 'Dataset', output_dir_path + '/output_autism_2_microbiomeHD_harmony/autism_2_microbiomeHD_nob', "DiseaseState", 30, [], 'Sam_id')

# # benchmarking other methods: 
# ### combat (combat_seq)
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD_combat.csv"
# data_mat_combat, meta_data_combat = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(data_mat_combat, meta_data_combat, 'Dataset', output_dir_path + '/output_autism_2_microbiomeHD_combat/autism_2_microbiomeHD_combat', "DiseaseState", 30, [], 'Sam_id')

# ### limma
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD_limma.csv"
# data_mat_limma, meta_data_limma = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(data_mat_limma, meta_data_limma, 'Dataset', output_dir_path + '/output_autism_2_microbiomeHD_limma/autism_2_microbiomeHD_limma', "DiseaseState", 30, [], 'Sam_id')

# ### MMUPHin
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD_MMUPHin.csv"
# data_mat_mmuphin, meta_data_mmuphin = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(data_mat_mmuphin, meta_data_mmuphin, 'Dataset', output_dir_path + '/output_autism_2_microbiomeHD_MMUPHin/autism_2_microbiomeHD_MMUPHin', "DiseaseState", 30, [], 'Sam_id')

# ### ConQuR
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD_ConQuR.csv"
# data_mat_conqur, meta_data_conqur = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(data_mat_conqur, meta_data_conqur, 'Dataset', output_dir_path + '/output_autism_2_microbiomeHD_ConQuR/autism_2_microbiomeHD_ConQuR', "DiseaseState", 30, [], 'Sam_id')

# ### ConQuR_libsize
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD_ConQuR_libsize.csv"
# data_mat_conqur_libsize, meta_data_conqur_libsize = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(data_mat_conqur_libsize, meta_data_conqur_libsize, 'Dataset', output_dir_path + '/output_autism_2_microbiomeHD_ConQuR_libsize/autism_2_microbiomeHD_ConQuR_libsize', "DiseaseState", 30, [], 'Sam_id')

# ### percentile_norm
# # percentile_norm(overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/autism_2_microbiomeHD_count_data.csv", overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/autism_2_microbiomeHD_meta_data.csv", "DiseaseState", "ASD", "comma", overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD")
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD_percentile_norm.csv"
# data_mat_percentile_norm, meta_data_percentile_norm = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(data_mat_percentile_norm, meta_data_percentile_norm, 'Dataset', output_dir_path + '/output_autism_2_microbiomeHD_percentile_norm/autism_2_microbiomeHD_percentile_norm', "DiseaseState", 30, [], 'Sam_id')

# ### global evaluation
# input_frame_path = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/autism_2_microbiomeHD_count_data.csv"
# bio_var = "DiseaseState"
# dataset_name = "autism_2_microbiomeHD"
# # methods_list = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "Tune_ConQuR", "Tune_ConQuR_libsize", "percentile_norm"]
# methods_list = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
# global_eval_dataframe(input_frame_path, bio_var, dataset_name, methods_list, overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results", output_dir_path, simulate = False)

# ### multi-method plot
# # df_l = [data_mat, res_h, data_mat_combat, data_mat_limma, data_mat_mmuphin, data_mat_conqur, data_mat_conqur_libsize, data_mat_conqur_tune, data_mat_conqur_tune_libsize, data_mat_percentile_norm]
# # methods = ["nobc", "harmony", "combat", "limma", "MMUPhin", "ConQuR", "ConQuR_libsize", "Tune_ConQuR", "Tune_ConQuR_libsize", "percentile_norm"]
# # meta_data_l = [meta_data, meta_data_h, meta_data_combat, meta_data_limma, meta_data_mmuphin, meta_data_conqur, meta_data_conqur_libsize, meta_data_conqur_tune, meta_data_conqur_tune_libsize, meta_data_percentile_norm]
# df_l = [data_mat, res_h, data_mat_combat, data_mat_limma, data_mat_mmuphin, data_mat_conqur, data_mat_conqur_libsize, data_mat_percentile_norm]
# methods = ["nobc", "harmony", "combat", "limma", "MMUPhin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
# meta_data_l = [meta_data, meta_data_h, meta_data_combat, meta_data_limma, meta_data_mmuphin, meta_data_conqur, meta_data_conqur_libsize, meta_data_percentile_norm]
# plot_PCOA_multiple('autism_2_microbiomeHD', df_l, methods, meta_data_l, used_var="Dataset", output_root= output_dir_path + '/')

################################################################################
# # cdi 3 microbiomeHD
# output_dir_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/outputs/cdi_3_microbiomeHD'
# address_directory = overall_path+'/mic_bc_benchmark/data/cdi_3_microbiomeHD'
# output_root = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/cdi_3_microbiomeHD"
# vars_use = ["Dataset"]
# IDCol = 'Sam_id'

# address_Y = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/cdi_3_microbiomeHD_meta_data.csv"

# ### nobc
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/cdi_3_microbiomeHD_count_data.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(data_mat, meta_data, 'Dataset', output_dir_path + '/output_cdi_3_microbiomeHD_nobc/cdi_3_microbiomeHD_nobc', "DiseaseState", 30, [], 'Sam_id')

# ### harmony
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD_harmony_adjusted_count.csv"
# res_h, meta_data_h = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(res_h, meta_data_h, 'Dataset', output_dir_path + '/output_cdi_3_microbiomeHD_harmony/cdi_3_microbiomeHD_harmony', "DiseaseState", 30, [], 'Sam_id')

# ### combat
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD_combat.csv"
# data_mat_combat, meta_data_combat = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(data_mat_combat, meta_data_combat, 'Dataset', output_dir_path + '/output_cdi_3_microbiomeHD_combat/cdi_3_microbiomeHD_combat', "DiseaseState", 30, [], 'Sam_id')

# ### limma
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD_limma.csv"
# data_mat_limma, meta_data_limma = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(data_mat_limma, meta_data_limma, 'Dataset', output_dir_path + '/output_cdi_3_microbiomeHD_limma/cdi_3_microbiomeHD_limma', "DiseaseState", 30, [], 'Sam_id')

# ### MMUPHin
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD_MMUPHin.csv"
# data_mat_mmuphin, meta_data_mmuphin = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(data_mat_mmuphin, meta_data_mmuphin, 'Dataset', output_dir_path + '/output_cdi_3_microbiomeHD_MMUPHin/cdi_3_microbiomeHD_MMUPHin', "DiseaseState", 30, [], 'Sam_id')

# ### ConQuR
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD_ConQuR.csv"
# data_mat_conqur, meta_data_conqur = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(data_mat_conqur, meta_data_conqur, 'Dataset', output_dir_path + '/output_cdi_3_microbiomeHD_ConQuR/cdi_3_microbiomeHD_ConQuR', "DiseaseState", 30, [], 'Sam_id')

# ### ConQuR_libsize
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD_ConQuR_libsize.csv"
# data_mat_conqur_libsize, meta_data_conqur_libsize = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(data_mat_conqur_libsize, meta_data_conqur_libsize, 'Dataset', output_dir_path + '/output_cdi_3_microbiomeHD_ConQuR_libsize/cdi_3_microbiomeHD_ConQuR_libsize', "DiseaseState", 30, [], 'Sam_id')

# ### percentile_norm
# # percentile_norm(overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/cdi_3_microbiomeHD_count_data.csv", overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/cdi_3_microbiomeHD_meta_data.csv", "DiseaseState", "CDI", "comma", overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD")
# # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD_percentile_norm.csv"
# # data_mat_percentile_norm, meta_data_percentile_norm = load_results_from_benchmarked_methods(address_X, address_Y)
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD_percentile_norm.csv"
# data_mat_percentile_norm, meta_data_percentile_norm = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(data_mat_percentile_norm, meta_data_percentile_norm, 'Dataset', output_dir_path + '/output_cdi_3_microbiomeHD_percentile_norm/cdi_3_microbiomeHD_percentile_norm', "DiseaseState", 30, [], 'Sam_id')

# #### global evaluation
# input_frame_path = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/cdi_3_microbiomeHD_count_data.csv"
# bio_var = "DiseaseState"
# dataset_name = "cdi_3_microbiomeHD"
# # methods_list = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "Tune_ConQuR", "Tune_ConQuR_libsize", "percentile_norm"]
# methods_list = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
# global_eval_dataframe(input_frame_path, bio_var, dataset_name, methods_list, overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results", output_dir_path, simulate = False)

# #### multi-method plot
# # df_l = [data_mat, res_h, data_mat_combat, data_mat_limma, data_mat_mmuphin, data_mat_conqur, data_mat_conqur_libsize, data_mat_conqur_tune, data_mat_conqur_tune_libsize, data_mat_percentile_norm]
# # methods = ["nobc", "harmony", "combat", "limma", "MMUPhin", "ConQuR", "ConQuR_libsize", "Tune_ConQuR", "Tune_ConQuR_libsize", "percentile_norm"]
# # meta_data_l = [meta_data, meta_data_h, meta_data_combat, meta_data_limma, meta_data_mmuphin, meta_data_conqur, meta_data_conqur_libsize, meta_data_conqur_tune, meta_data_conqur_tune_libsize, meta_data_percentile_norm]
# df_l = [data_mat, res_h, data_mat_combat, data_mat_limma, data_mat_mmuphin, data_mat_conqur, data_mat_conqur_libsize, data_mat_percentile_norm]
# methods = ["nobc", "harmony", "combat", "limma", "MMUPhin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
# meta_data_l = [meta_data, meta_data_h, meta_data_combat, meta_data_limma, meta_data_mmuphin, meta_data_conqur, meta_data_conqur_libsize, meta_data_percentile_norm]
# plot_PCOA_multiple('cdi_3_microbiomeHD', df_l, methods, meta_data_l, used_var="Dataset", output_root= output_dir_path + '/')


##############################################################################
# # ibd_3_CMD
# output_dir_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/outputs/ibd_3_CMD'
# address_directory = overall_path+'/mic_bc_benchmark/data/ibd_3_CMD'
# output_root = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/ibd_3_CMD"
# vars_use = ["study_name"]
# IDCol = 'Sam_id'

# # ### nobc
# address_Y = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/ibd_3_CMD_meta_data.csv"

# # ### nobc
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/ibd_3_CMD_count_data.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(data_mat, meta_data, "study_name", output_dir_path + '/output_ibd_3_CMD_nobc/ibd_3_CMD_nobc', "disease", 30, ['gender', 'age'], 'Sam_id', datatype = 'relab')

# # ### harmony
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/ibd_3_CMD/ibd_3_CMD_harmony_adjusted_count.csv"
# res_h, meta_data_h = load_results_from_benchmarked_methods(address_X, address_Y)
# # # res_h, meta_data_h = generate_harmony_results(data_mat, meta_data, IDCol, vars_use, overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/ibd_3_CMD/"+"harmony")
# # Evaluate(res_h, meta_data_h, "study_name", output_dir_path + '/output_ibd_3_CMD_harmony/ibd_3_CMD_harmony', "disease", 30, ['gender', 'age'], 'Sam_id', datatype = 'relab')

# # # benchmarking other methods: 
# address_Y = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/ibd_3_CMD_meta_data.csv"

# # ### limma
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/ibd_3_CMD/ibd_3_CMD_limma.csv"
# data_mat_limma, meta_data_limma = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(data_mat_limma, meta_data_limma,  "study_name", output_dir_path + '/output_ibd_3_CMD_limma/ibd_3_CMD_limma', "disease", 30, ['gender', 'age'], 'Sam_id', datatype = 'relab')

# # ### MMUPHin
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/ibd_3_CMD/ibd_3_CMD_MMUPHin.csv"
# data_mat_mmuphin, meta_data_mmuphin = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(data_mat_mmuphin, meta_data_mmuphin,  "study_name", output_dir_path + '/output_ibd_3_CMD_MMUPHin/ibd_3_CMD_MMUPHin', "disease", 30, ['gender', 'age'], 'Sam_id', datatype = 'relab')

# # ### ConQuR_rel
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/ibd_3_CMD/ibd_3_CMD_ConQuR_rel.csv"
# data_mat_conqur_rel, meta_data_conqur_rel = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(data_mat_conqur_rel, meta_data_conqur_rel, "study_name", output_dir_path + '/output_ibd_3_CMD_ConQuR_rel/ibd_3_CMD_ConQuR_rel', "disease", 30, ['gender', 'age'], 'Sam_id')

# # ### percentile_norm
# # # percentile_norm(overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/ibd_3_CMD_count_data.csv", overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/ibd_3_CMD_meta_data.csv", "disease", "IBD", "comma", overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/ibd_3_CMD/ibd_3_CMD", datatype = 'relab')
# # # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/ibd_3_CMD/ibd_3_CMD_percentile_norm.csv"
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/ibd_3_CMD/ibd_3_CMD_percentile_norm.csv"
# data_mat_percentile_norm, meta_data_percentile_norm = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(data_mat_percentile_norm, meta_data_percentile_norm, "study_name", output_dir_path + '/output_ibd_3_CMD_percentile_norm/ibd_3_CMD_percentile_norm', "disease", 30, ['gender', 'age'], 'Sam_id', datatype = 'relab')

# #### global evaluation
# input_frame_path = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/ibd_3_CMD_count_data.csv"
# bio_var = "disease"
# dataset_name = 'ibd_3_CMD'
# methods_list = ['nobc', 'harmony', 'limma', 'MMUPHin', 'ConQuR_rel', 'percentile_norm']
# global_eval_dataframe(input_frame_path, bio_var, dataset_name, methods_list, overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results", output_dir_path, datatype = 'relab', simulate = False)

# #### multi-method plot
# # df_l = [data_mat, res_h, data_mat_limma, data_mat_mmuphin, data_mat_conqur_rel, data_mat_conqur_tune_rel, data_mat_percentile_norm]
# # methods = ["nobc", "harmony", "limma", "MMUPhin", "ConQuR_rel", "Tune_ConQuR_rel", "percentile_norm"]
# # meta_data_l = [meta_data, meta_data_h, meta_data_limma, meta_data_mmuphin, meta_data_conqur_rel, meta_data_conqur_tune_rel, meta_data_percentile_norm]
# df_l = [data_mat, res_h, data_mat_limma, data_mat_mmuphin, data_mat_conqur_rel, data_mat_percentile_norm]
# methods = ["nobc", "harmony", "limma", "MMUPhin", "ConQuR_rel", "percentile_norm"]
# meta_data_l = [meta_data, meta_data_h, meta_data_limma, meta_data_mmuphin, meta_data_conqur_rel, meta_data_percentile_norm]
# plot_PCOA_multiple('ibd_3_CMD', df_l, methods, meta_data_l, used_var="study_name", output_root= output_dir_path + '/', datatype = "relab")


##############################################################################
# # CRC_8_CMD
# output_dir_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/outputs/CRC_8_CMD'
# address_directory = overall_path+'/mic_bc_benchmark/data/CRC_8_CMD'
# output_root = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/CRC_8_CMD"
# vars_use = ["study_name"]
# IDCol = 'Sam_id'

# # benchmarking other methods: 
# address_Y = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/CRC_8_CMD_meta_data.csv"

# ### nobc
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/CRC_8_CMD_count_data.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(data_mat, meta_data, "study_name", output_dir_path + '/output_CRC_8_CMD_nobc/CRC_8_CMD_nobc', "disease", 30, ['gender', 'age'], 'Sam_id', datatype = 'relab')

# ### harmony
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/CRC_8_CMD/CRC_8_CMD_harmony_adjusted_count.csv"
# res_h, meta_data_h = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(res_h, meta_data_h, "study_name", output_dir_path + '/output_CRC_8_CMD_harmony/CRC_8_CMD_harmony', "disease", 30, ['gender', 'age'], 'Sam_id', datatype = 'relab')

# ## limma
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/CRC_8_CMD/CRC_8_CMD_limma.csv"
# data_mat_limma, meta_data_limma = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(data_mat_limma, meta_data_limma,  "study_name", output_dir_path + '/output_CRC_8_CMD_limma/CRC_8_CMD_limma', "disease", 30, ['gender', 'age'], 'Sam_id', datatype = 'relab')

# ### MMUPHin
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/CRC_8_CMD/CRC_8_CMD_MMUPHin.csv"
# data_mat_mmuphin, meta_data_mmuphin = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(data_mat_mmuphin, meta_data_mmuphin,  "study_name", output_dir_path + '/output_CRC_8_CMD_MMUPHin/CRC_8_CMD_MMUPHin', "disease", 30, ['gender', 'age'], 'Sam_id', datatype = 'relab')

# ### ConQuR_rel
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/CRC_8_CMD/CRC_8_CMD_ConQuR_rel.csv"
# data_mat_conqur_rel, meta_data_conqur_rel = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(data_mat_conqur_rel, meta_data_conqur_rel, "study_name", output_dir_path + '/output_CRC_8_CMD_ConQuR_rel/CRC_8_CMD_ConQuR_rel', "disease", 30, ['gender', 'age'], 'Sam_id', datatype = 'relab')


# # ### percentile_norm
# # # print("doing percentile_norm")
# # # # percentile_norm(overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/CRC_8_CMD_count_data.csv", overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/CRC_8_CMD_meta_data.csv", "disease", "IBD", "comma", overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/CRC_8_CMD/CRC_8_CMD")
# # # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/CRC_8_CMD/CRC_8_CMD_percentile_norm.csv"
# # # data_mat_percentile_norm, meta_data_percentile_norm = load_results_from_benchmarked_methods(address_X, address_Y)
# address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/CRC_8_CMD/CRC_8_CMD_percentile_norm.csv"
# data_mat_percentile_norm, meta_data_percentile_norm = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(data_mat_percentile_norm, meta_data_percentile_norm, "study_name", output_dir_path + '/output_CRC_8_CMD_percentile_norm/CRC_8_CMD_percentile_norm', "disease", 30, ['gender', 'age'], 'Sam_id', datatype = 'relab')

# #### global evaluation
# input_frame_path = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_data/CRC_8_CMD_count_data.csv"
# bio_var = "disease"
# dataset_name = 'CRC_8_CMD'
# # methods_list = ['nobc', 'harmony', 'limma', 'MMUPHin', 'ConQuR_rel', 'Tune_ConQuR_rel', 'percentile_norm']
# methods_list = ['nobc', 'harmony', 'limma', 'MMUPHin', 'ConQuR_rel', 'percentile_norm']
# global_eval_dataframe(input_frame_path, bio_var, dataset_name, methods_list, overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results", output_dir_path, datatype = 'relab', simulate = False)

# # #### multi-method plot
# # # df_l = [data_mat, res_h, data_mat_limma, data_mat_mmuphin, data_mat_conqur_rel, data_mat_conqur_tune_rel, data_mat_percentile_norm]
# # # methods = ["nobc", "harmony", "limma", "MMUPhin", "ConQuR_rel", "Tune_ConQuR_rel", "percentile_norm"]
# # # meta_data_l = [meta_data, meta_data_h, meta_data_limma, meta_data_mmuphin, meta_data_conqur_rel, meta_data_conqur_tune_rel, meta_data_percentile_norm]
# df_l = [data_mat, res_h, data_mat_limma, data_mat_mmuphin, data_mat_conqur_rel, data_mat_percentile_norm]
# methods = ["nobc", "harmony", "limma", "MMUPhin", "ConQuR_rel", "percentile_norm"]
# meta_data_l = [meta_data, meta_data_h, meta_data_limma, meta_data_mmuphin, meta_data_conqur_rel, meta_data_percentile_norm]
# plot_PCOA_multiple('CRC_8_CMD', df_l, methods, meta_data_l, used_var="study_name", output_root= output_dir_path + '/', datatype = "relab")

##############################################################################
output_dir_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/outputs'
methods = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
datasets = ["autism_2_microbiomeHD", "cdi_3_microbiomeHD"]
output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
visualize_simulation_stats('count_rw', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, simulate = False)

output_dir_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/outputs'
methods = ["harmony", "limma", "MMUPHin", "ConQuR_rel", "percentile_norm"]
datasets = ["ibd_3_CMD", "CRC_8_CMD"]
output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
visualize_simulation_stats('relab_rw', output_dir_l, datasets, methods, highlighted_method = "ConQuR_rel", line = True, count_l = [False, False], simulate = False)

output_dir_path = '/athena/linglab/scratch/chf4012/simulation_data_eval_small_norelation_080723'
methods = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
datasets = ["out_1_0_0", "out_1_0_0.099", "out_1_0_0.299", "out_1_0_0.499", "out_1_0_0.699", "out_1_0_0.899"]
output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
visualize_simulation_stats('sim_1_0_bio', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, count_l = [True]*len(datasets), simulate = True)

output_dir_path = '/athena/linglab/scratch/chf4012/simulation_data_eval_small_norelation_080723'
methods = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
datasets = ["out_1_0.099_0", "out_1_0.099_0.099", "out_1_0.099_0.299", "out_1_0.099_0.499", "out_1_0.099_0.699", "out_1_0.099_0.899"]
output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
visualize_simulation_stats('sim_1_0.099_bio', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, count_l = [True]*len(datasets), simulate = True)

output_dir_path = '/athena/linglab/scratch/chf4012/simulation_data_eval_small_norelation_080723'
methods = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
datasets = ["out_1_0.299_0", "out_1_0.299_0.099", "out_1_0.299_0.299", "out_1_0.299_0.499", "out_1_0.299_0.699"]
output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
visualize_simulation_stats('sim_1_0.299_bio', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, count_l = [True]*len(datasets), simulate = True, sim_num_iters = 5)

output_dir_path = '/athena/linglab/scratch/chf4012/simulation_data_eval_small_norelation_080723'
methods = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
datasets = ["out_1_0.499_0", "out_1_0.499_0.099", "out_1_0.499_0.299", "out_1_0.499_0.499"]
output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
visualize_simulation_stats('sim_1_0.499_bio', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, count_l = [True]*len(datasets), simulate = True, sim_num_iters = 5)

output_dir_path = '/athena/linglab/scratch/chf4012/simulation_data_eval_small_norelation_080723'
methods = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
datasets = ["out_1_0.699_0", "out_1_0.699_0.099", "out_1_0.699_0.299"]
output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
visualize_simulation_stats('sim_1_0.699_bio', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, count_l = [True]*len(datasets), simulate = True, sim_num_iters = 5)

output_dir_path = '/athena/linglab/scratch/chf4012/simulation_data_eval_small_norelation_080723'
methods = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
datasets = ["out_1_0.899_0", "out_1_0.899_0.099"]
output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
visualize_simulation_stats('sim_1_0.899_bio', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, count_l = [True]*len(datasets), simulate = True, sim_num_iters = 5)


output_dir_path = '/athena/linglab/scratch/chf4012/simulation_data_eval_small_norelation_080723'
methods = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
datasets = ["out_1.25_0_0", "out_1.25_0_0.099", "out_1.25_0_0.299", "out_1.25_0_0.499", "out_1.25_0_0.699", "out_1.25_0_0.899"]
output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
visualize_simulation_stats('sim_1.25_0_bio', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, count_l = [True]*len(datasets), simulate = True)

output_dir_path = '/athena/linglab/scratch/chf4012/simulation_data_eval_small_norelation_080723'
methods = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
datasets = ["out_1.25_0.099_0", "out_1.25_0.099_0.099", "out_1.25_0.099_0.299", "out_1.25_0.099_0.499", "out_1.25_0.099_0.699", "out_1.25_0.099_0.899"]
output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
visualize_simulation_stats('sim_1.25_0.099_bio', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, count_l = [True]*len(datasets), simulate = True)

output_dir_path = '/athena/linglab/scratch/chf4012/simulation_data_eval_small_norelation_080723'
methods = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
datasets = ["out_1.25_0.299_0", "out_1.25_0.299_0.099", "out_1.25_0.299_0.299", "out_1.25_0.299_0.499", "out_1.25_0.299_0.699"]
output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
visualize_simulation_stats('sim_1.25_0.299_bio', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, count_l = [True]*len(datasets), simulate = True, sim_num_iters = 5)

output_dir_path = '/athena/linglab/scratch/chf4012/simulation_data_eval_small_norelation_080723'
methods = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
datasets = ["out_1.25_0.499_0", "out_1.25_0.499_0.099", "out_1.25_0.499_0.299", "out_1.25_0.499_0.499"]
output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
visualize_simulation_stats('sim_1.25_0.499_bio', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, count_l = [True]*len(datasets), simulate = True, sim_num_iters = 5)

output_dir_path = '/athena/linglab/scratch/chf4012/simulation_data_eval_small_norelation_080723'
methods = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
datasets = ["out_1.25_0.699_0", "out_1.25_0.699_0.099", "out_1.25_0.699_0.299"]
output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
visualize_simulation_stats('sim_1.25_0.699_bio', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, count_l = [True]*len(datasets), simulate = True, sim_num_iters = 5)

output_dir_path = '/athena/linglab/scratch/chf4012/simulation_data_eval_small_norelation_080723'
methods = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
datasets = ["out_1.25_0.899_0", "out_1.25_0.899_0.099"]
output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
visualize_simulation_stats('sim_1.5_0.899_bio', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, count_l = [True]*len(datasets), simulate = True, sim_num_iters = 5)

# do the same thing for out_1.5 now
output_dir_path = '/athena/linglab/scratch/chf4012/simulation_data_eval_small_norelation_080723'
methods = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
datasets = ["out_1.5_0_0", "out_1.5_0_0.099", "out_1.5_0_0.299", "out_1.5_0_0.499", "out_1.5_0_0.699", "out_1.5_0_0.899"]
output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
visualize_simulation_stats('sim_1.5_0_bio', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, count_l = [True]*len(datasets), simulate = True)

output_dir_path = '/athena/linglab/scratch/chf4012/simulation_data_eval_small_norelation_080723'
methods = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
datasets = ["out_1.5_0.099_0", "out_1.5_0.099_0.099", "out_1.5_0.099_0.299", "out_1.5_0.099_0.499", "out_1.5_0.099_0.699", "out_1.5_0.099_0.899"]
output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
visualize_simulation_stats('sim_1.5_0.099_bio', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, count_l = [True]*len(datasets), simulate = True)

output_dir_path = '/athena/linglab/scratch/chf4012/simulation_data_eval_small_norelation_080723'
methods = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
datasets = ["out_1.5_0.299_0", "out_1.5_0.299_0.099", "out_1.5_0.299_0.299", "out_1.5_0.299_0.499", "out_1.5_0.299_0.699"]
output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
visualize_simulation_stats('sim_1.5_0.299_bio', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, count_l = [True]*len(datasets), simulate = True, sim_num_iters = 5)

output_dir_path = '/athena/linglab/scratch/chf4012/simulation_data_eval_small_norelation_080723'
methods = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
datasets = ["out_1.5_0.499_0", "out_1.5_0.499_0.099", "out_1.5_0.499_0.299", "out_1.5_0.499_0.499"]
output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
visualize_simulation_stats('sim_1.5_0.499_bio', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, count_l = [True]*len(datasets), simulate = True, sim_num_iters = 5)

output_dir_path = '/athena/linglab/scratch/chf4012/simulation_data_eval_small_norelation_080723'
methods = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
datasets = ["out_1.5_0.699_0", "out_1.5_0.699_0.099", "out_1.5_0.699_0.299"]
output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
visualize_simulation_stats('sim_1.5_0.699_bio', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, count_l = [True]*len(datasets), simulate = True, sim_num_iters = 5)

output_dir_path = '/athena/linglab/scratch/chf4012/simulation_data_eval_small_norelation_080723'
methods = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
datasets = ["out_1.5_0.899_0", "out_1.5_0.899_0.099"]
output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
visualize_simulation_stats('sim_1.5_0.899_bio', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, count_l = [True]*len(datasets), simulate = True, sim_num_iters = 5)



##############################################################################