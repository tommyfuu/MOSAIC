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
import statsmodels.api as sm
import matplotlib as mpl
import sys, yaml

ARGPARSE_SWITCH = True

if ARGPARSE_SWITCH:
    import argparse
    parser = argparse.ArgumentParser(
                prog='evaluate_bc_methods',
                description='two options: (1) generate batch corrected results for harmony and percentile normalization along w runtime; (2) evaluate the batch correction results.',
                epilog='Text at the bottom of help')

    parser.add_argument("-o", "--option", type=int, default = 10, help='either 1 or 2, option 1: generate batch corrected results for harmony and percentile normalization along w runtime; option 2: evaluate the batch correction results; option 2: evaluate the batch correction results for multiple methods.')
    parser.add_argument("-i", "--iteration", type=int, default = 1, help='the iteration number')
    parser.add_argument("-d", "--datatype", default = 'count', help='either count or relab')
    parser.add_argument("-r", "--related", default = 'no', help='whether the batch effect is related to library size')
    parser.add_argument("-a", "--agent", default = 'cond_1', help='binarizing agent for differential abundance evaluation')
    parser.add_argument("-p", "--overallpath", default = '/athena/linglab/scratch/chf4012', help='overall path for data interaction and saving')

    args = vars(parser.parse_args())

    overall_path = args['overallpath']


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
    # check if current path has benchmark in it
    if '/benchmark' in os.getcwd():
        r.source('./PERMANOVA_supporting.R')
    else:
        r.source('./benchmark/PERMANOVA_supporting.R')
    # r_used_var = meta_data[used_var]
    # r_bio_var = meta_data[bio_var]

    # initial graphics
    if datatype == 'count':
        height = 12
    else:
        height = 6
    r.pdf(output_root+dataset_name+"_multi_PCOA_both_batch.pdf", len(batch_corrected_df_l)*6, height=height)
    if datatype == 'count':
        r.par(mfrow = robjects.IntVector([2, len(batch_corrected_df_l)]))
    else:
        r.par(mfrow = robjects.IntVector([1, len(batch_corrected_df_l)]))

    # plot the subplots in order
    if datatype == 'count':
        for idx, batch_corrected_df in enumerate(batch_corrected_df_l):
            data = np.array(batch_corrected_df)
            # use the first len(bio_var_l) rows
            data = data[:len(meta_data_l[idx])]
            # check if na/+-inf/huge value is in the data
            if np.isnan(data).any() or np.isinf(data).any() or np.max(data) > 100000:
                print("combat family explosion case, default #s and skipped")
                data = np.zeros(data.shape)
                r.Plot_single_PCoA(data, meta_data_l[idx][used_var], dissimilarity="Aitch", bc_method = methods[idx])
            else:
                data = np.where(data<np.percentile(data.flatten(), 0.01), 0, data)
                data = data+np.abs(np.min(data))
                r_used_var = meta_data_l[idx][used_var]
                r.Plot_single_PCoA(data, r_used_var, dissimilarity="Aitch", bc_method = methods[idx])
    
    for idx, batch_corrected_df in enumerate(batch_corrected_df_l):
        data = np.array(batch_corrected_df)
        if np.isnan(data).any() or np.isinf(data).any() or np.max(data) > 100000:
            print("combat family explosion case, default #s and skipped")
            data = np.zeros(data.shape)
            r.Plot_single_PCoA(data, meta_data_l[idx][used_var], dissimilarity="Aitch", bc_method = methods[idx])

        else:
            data = np.where(data<np.percentile(data.flatten(), 0.01), 0, data)
            data = data+np.abs(np.min(data))
            r_used_var = meta_data_l[idx][used_var]
            r.Plot_single_PCoA(data, r_used_var, dissimilarity="Bray", bc_method = methods[idx])
    grdevices.dev_off()
    return

class Evaluate(object):
    def __init__(
            self, batch_corrected_df, meta_data, batch_var, output_root, bio_var = False, n_pc=30, covar_l = [], IDCol = None, test_percent = 0.2, poslabel = '', pipeline = 'default', taxa_gt = None, datatype = "count", binarizing_agent_biovar = 'H', method = 'combat_seq'
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
        self.binarizing_agent_biovar = binarizing_agent_biovar
        self.method = method
        print("Current method", self.method)
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

        self.alpha_diversity_and_tests_for_batches_and_biovars()
        self.predict_difference_RF()
        self.calculate_R_sq()
        ## TODO: add FP stuff for simulation

        # convert summary dict to df
        self.summary_df = pd.DataFrame.from_dict(self.short_summary_dict, orient='index')
        self.summary_df.to_csv(self.output_root+"_summary.csv")

    def diff_abund_test(self):
        ### if two categories then wilcox test, three or more then kruskal-wallis test
        print("DIFF ABUND TEST")
        df = self.batch_corrected_df
        
        # if datatype is count, conduct relative abundance calculation
        if self.datatype == "count":
            # check if negative value exists, if exists, first add everything with the abs(most neg)
            if not self.method in ['limma', 'harmony']:
                df = df.div(df.sum(axis=1), axis=0)
        taxa_names = df.columns
        meta_data = self.meta_data

        # check if na/+-inf exists (combat weird explosion case), OR if some very large number exists (combat-seq explosion case), if so, skip the whole
        if df.isnull().values.any() or np.isinf(df.values).any() or np.max(df.values) > 100000:
            self.short_summary_dict["FDR"] = 0
            self.short_summary_dict["power"] = 0
            # fill p_val and FDR with 1
            df["p_value"] = 1
            df["FDR_p_value"] = 1
            df.to_csv(self.output_root+"_diff_abund_test.csv")
            return

        # use the metadata to get indices of samples in each bio_var category
        bio_var_l = list(meta_data[self.bio_var].values)
        # binarize bio_var_l with biovar_binarizing_agent
        bio_var_l = [1 if i == self.binarizing_agent_biovar else 0 for i in bio_var_l]

        taxa_names = df.columns
        # iterate over columns and conduct statistical test
        d = {}
        for idx, i in enumerate(taxa_names):
            # get the data for the current taxa
            data = df[i]
            # add covariates to df if exists
            if self.covar_l != []:
                for covar in self.covar_l:
                    # if not numerical, one-hot encode numerically and concat
                    if meta_data[covar].dtype == 'object':
                        # one-hot encode
                        covar_df = pd.get_dummies(meta_data[covar])
                        # convert True/False to 1/0
                        covar_df = covar_df.astype(int)
                        # reset index
                        covar_df.index = df.index
                        data = pd.concat([data, covar_df], axis=1)
            # add constant to df
            # print("data", list(data))
            data = sm.add_constant(data)
            # use the first len(bio_var_l) rows
            data = data.iloc[:len(bio_var_l)]
            # fit ols model
            model_ols = sm.OLS(bio_var_l, data, missing = 'drop').fit()
            # save p value to dict
            d[f'{i}'] = model_ols.pvalues[i]

        # print(d)
        # only keep columns that are taxa names (drop constant and covariates)
        df = self.batch_corrected_df[taxa_names]
        # transpose the batch corrected dataframe and save the p-values alongside
        df = df.T
        # save p value to df
        df["p_value"] = df.index.map(d)
        # replace nan p_value with 1
        df["p_value"] = df["p_value"].fillna(1)

        # add FDR corrected p-values
        df["FDR_p_value"] = multipletests(df['p_value'].values, method='fdr_bh')[1]
        ### if we know the ground truth, so we can calculate the stats such as FDR, TPR, FPR, sensitivity, etc
        if self.taxa_gt is not None:
            df["ground_truth"] = [1 if i in self.taxa_gt else 0 for i, taxa in enumerate(taxa_names)]
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
        # check if na/+-inf exists (combat weird explosion case), OR if some very large number exists (combat-seq explosion case), if so, skip the whole
        df = self.batch_corrected_df
        if df.isnull().values.any() or np.isinf(df.values).any() or np.max(df.values) > 100000:
            self.short_summary_dict["batch_aitch_r2"] = 0
            self.short_summary_dict["batch_bray_r2"] = 0
            self.short_summary_dict["batch_aitch_pval"] = 0
            self.short_summary_dict["batch_bray_pval"] = 0
            self.short_summary_dict["biovar_aitch_r2"] = 0
            self.short_summary_dict["biovar_bray_r2"] = 0
            self.short_summary_dict["biovar_aitch_pval"] = 0
            self.short_summary_dict["biovar_bray_pval"] = 0
            print("combat family explosion case, default #s and skipped")
            return
            
        data = np.array(self.batch_corrected_df)
        data = np.where(data<np.percentile(data.flatten(), 0.01), 0, data)
        data = data+np.abs(np.min(data))

        # use the first len(bio_var_l) rows
        data = data[:len(self.meta_data[self.IDCol])]

        # attempting rpy2
        #Must be activated
        pandas2ri.activate()
        numpy2ri.activate()

        # import r packages/functions
        r = robjects.r
        if '/benchmark' in os.getcwd():
            print("AAAA")
            print("AAAA")
            print("AAAA")
            print("AAAA")
            print("AAAA")
            print(os.getcwd())
            r.source('./PERMANOVA_supporting.R')
        else:
            r.source('./benchmark/PERMANOVA_supporting.R')
        # r.source('./PERMANOVA_supporting.R')
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

         # check if na/+-inf exists (combat weird explosion case), OR if some very large number exists (combat-seq explosion case), if so, skip the whole
        df = self.batch_corrected_df
        if df.isnull().values.any() or np.isinf(df.values).any() or np.max(df.values) > 100000:
            print("combat family explosion case, default #s and skipped")
            return 1, None
        data = np.array(self.batch_corrected_df)
        data = np.where(data<np.percentile(data.flatten(), 0.01), 0, data)
        data = data+np.abs(np.min(data))
        ids = list(self.meta_data[self.IDCol])
        # use the first len(bio_var_l) rows
        data = data[:len(ids)]
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
            # check if dfs are none
            if shannon_df_batch is not None:
                shannon_df_batch.to_csv(self.output_root+"_shannon_df_batch.csv")
                shannon_df_biovar.to_csv(self.output_root+"_shannon_df_biovar.csv")

        # save to summary dict
        self.short_summary_dict["batches_shannon_pval"] = shannon_global_pval_batch
        self.short_summary_dict["biovar_shannon_pval"] = shannon_global_pval_biovar
        return
    
    def predict_difference_RF(self):
        # check if na/+-inf exists (combat weird explosion case), OR if some very large number exists (combat-seq explosion case), if so, skip the whole
        df = self.batch_corrected_df
        if df.isnull().values.any() or np.isinf(df.values).any() or np.max(df.values) > 100000:
            self.short_summary_dict["rf_baseline_likelihood"] = 0
            self.short_summary_dict["rf_average_acc"] = 0
            self.short_summary_dict["rf_macro_precision"] = 0
            self.short_summary_dict["weighted_precision"] = 0
            self.short_summary_dict["macro_recall"] = 0
            self.short_summary_dict["weighted_recall"] = 0
            self.short_summary_dict["macro_f1"] = 0
            self.short_summary_dict["weighted_f1"] = 0
            self.short_summary_dict["auc"] = 0
            print("combat family explosion case, default #s and skipped")
            return

        # get the data
        used_x = self.batch_corrected_df.copy() # note that standard scaler is already conducted
        # use the first len(bio_var_l) rows
        used_x = used_x.iloc[:len(self.meta_data[self.IDCol])]
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
        current_method_dict.update({"batch_shannon_pval": shannon_pval_batch})
        current_method_dict.update({"biovar_shannon_pval": shannon_pval_biovar})
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

    current_runtime_kw_path = [result for result in benchmarked_results_dir_files if "runtime.txt" in result][0]
    with open(benchmarked_results_dir+'/'+current_runtime_kw_path) as f:
        lines = f.readlines()

    # print(method_dict)
    for line in lines:
        if "combat" in line and "combat" in method_dict.keys():
            method_dict["combat"]["runtime"] = float(line.split(" ")[-2])
        if "combat_seq" in line:
            method_dict["combat_seq"]["runtime"] = float(line.split(" ")[-2])
        if "limma" in line and "limma" in method_dict.keys():
            method_dict["limma"]["runtime"] = float(line.split(" ")[-2])
        if "MMUPHin" in line and "MMUPHin" in method_dict.keys():
            method_dict["MMUPHin"]["runtime"] = float(line.split(" ")[-2])
        if "ConQuR" in line and "Tune" not in line and 'libsize' not in line and "ConQuR" in method_dict.keys():
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

    if 'harmony' in method_dict.keys():
        benchmarked_data_harmony_dir = benchmarked_results_dir
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
                print("harmony lines")
                print(lines)
                method_dict["harmony"]["runtime"] = float(lines[0].split(" ")[-2])

    if 'percentile_norm' in method_dict.keys():
        benchmarked_data_harmony_dir = benchmarked_results_dir
        for time_file in current_runtime_kw_paths:
            time_file = benchmarked_data_harmony_dir + "/"+time_file
            if "percentile_norm_elapsed_time" in time_file:
                with open(time_file) as f:
                    lines = f.readlines()
                method_dict["percentile_norm"]["runtime"] = float(lines[0].split(" ")[0])
    if 'nobc' in method_dict:
        method_dict['nobc']['runtime'] = 'NA'
    # print into a pandas df where rows are datasets and columns and different stats
    results_df = pd.DataFrame.from_dict(method_dict, orient ='index') 
    results_df.T.to_csv(output_dir_path+"/global_benchmarking_stats_"+dataset_name+".csv")
    return pd.DataFrame.from_dict(method_dict, orient ='index') 
    
def visualize_simulation_stats(output_root, output_dir_l, datasets, methods, highlighted_method, simulate = False, sim_num_iters = 1000, dimensions = (5, 5), taxa_gt = None, line = True, count_l = [True, True, False, False], 
    marker_dict = {'nobc': 'o', "combat": 'v', "combat_seq": "D", "limma": '^', "MMUPHin": '<', "ConQuR": '>', "ConQuR_libsize": 's', "ConQuR_rel": 'p', "harmony":"+", "percentile_norm": "*"},
    method_colors_dict = {'nobc': 'black', "combat": 'red', "combat_seq": 'yellow', "limma": 'blue', "MMUPHin": 'green', "ConQuR": 'orange', "ConQuR_libsize": 'purple', "ConQuR_rel": 'pink', "harmony":"turquoise", "percentile_norm": "gray"},
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
    global_methods_batch_shannon_pval_rejection_proportions_l_dict = {}
    global_methods_biovar_shannon_pval_rejection_proportions_l_dict = {}

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

        # if 'count' in output_dir_path and 'out_1.5_0.25_0.75' in output_dir:
        #     # if combat_seq runtime is nan, print
        #     if np.isnan(methods_runtime['combat_seq']):
        #         print(iter, "combat_seq runtime is nan")
            # print(iter, "combat_seq and MMUPHin runtime", methods_runtime['combat_seq'], methods_runtime['MMUPHin'])

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
                    global_methods_batch_shannon_pval_rejection_proportions_l_dict[method] = []
                    global_methods_biovar_shannon_pval_rejection_proportions_l_dict[method] = []
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
                # for shannon pvals, calculate the rejection proportions for each method at each parameter combination instead
                global_methods_batch_shannon_pval_rejection_proportions_l_dict[method].append(np.mean([(cross_iter_batch_shannon_pval_dict[iter][method] < 0.05)*1 for iter in range(sim_num_iters)]))
                global_methods_biovar_shannon_pval_rejection_proportions_l_dict[method].append(np.mean([(cross_iter_biovar_shannon_pval_dict[iter][method] < 0.05)*1 for iter in range(sim_num_iters)]))
    print("plotting")
    def plot_stats(stats_summary_name, stats_name_l, stats_dict_1, stats_dict_2 = {}, postfix = '.png', ylim=[], pvallines = [False, False], line=True):
        mpl.rcParams['pdf.fonttype'] = 42 # ensure exported pdf has edited text
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
        for method in methods:
            color = method_colors_dict[method]
            if method != highlighted_method:
                alpha = 1
            else:
                alpha = 1
            # if count_l[0] or 'FDR' in stats_summary_name:
            ax1.tick_params(axis='both', which='major', labelsize=14)
            ax1.spines["top"].set_visible(False)
            ax1.spines["right"].set_visible(False)
            if line:
                ax1.plot(datasets, stats_dict_1[method], label=method, color=color, alpha=alpha, marker=marker_dict[method], linestyle=linestyle)
            else:
                dataset_to_x_dict = {datasets[i]: i for i in range(len(datasets))}
                # jitter the xs of each point using a normal distribution
                jitter = 0.03
                for dataset in datasets:
                    dataset_to_x_dict[dataset] += np.random.normal(0, jitter)
                ax1.plot(dataset_to_x_dict.values(), stats_dict_1[method], label=method, color=color, alpha=alpha, marker=marker_dict[method], linestyle=linestyle)
            ax1.set_title(stats_name_l[0])
            if len(datasets) == 2:
                ax1.set_xticks([0, 1], datasets)
            else:
                ax1.set_xticks(datasets)                    
            # set xticks to be verticle
            for tick in ax1.get_xticklabels():
                tick.set_rotation(90)

            if 'FDR' in stats_name_l:
                # plot log(1e-4, 0.2, 0.4, 0.6, 0.8, 1) on y axis while keeping the original values in the legend
                ax1.set_yticks([np.log2(1e-4), np.log2(0.2), np.log2(0.4), np.log2(0.6), np.log2(0.8), np.log2(1)], ["~0", "0.2", "0.4", "0.6", "0.8", "1"])

            if 'runtime' in stats_name_l:
                ax1.set_yticks([np.log2(1e-4), np.log2(10), np.log2(25), np.log2(50), np.log2(100), np.log2(200), np.log2(350), np.log2(500)], ["~0", "10", "25", "50", "100", "200", "350", "500"])
            
            if stats_dict_2 != {}:
                ax2.tick_params(axis='both', which='major', labelsize=14)
                ax2.spines["top"].set_visible(False)
                ax2.spines["right"].set_visible(False)
                if line:
                    ax2.plot(datasets, stats_dict_2[method], label=method, color=color, alpha=alpha, marker=marker_dict[method], linestyle=linestyle)
                else:
                    dataset_to_x_dict = {datasets[i]: i for i in range(len(datasets))}
                    # jitter the xs of each point using a normal distribution
                    jitter = 0.03
                    for dataset in datasets:
                        dataset_to_x_dict[dataset] += np.random.normal(0, jitter)
                    ax2.plot(dataset_to_x_dict.values(), stats_dict_2[method], label=method, color=color, alpha=alpha, marker=marker_dict[method], linestyle=linestyle)
                ax2.set_title(stats_name_l[1])
                if len(datasets) == 2:
                    ax2.set_xticks([0, 1], datasets)
                else:
                    ax2.set_xticks(datasets)
                # set xticks to be verticle
                for tick in ax2.get_xticklabels():
                    tick.set_rotation(90)

                if 'FDR' in stats_name_l:
                    # plot log(1e-4, 0.2, 0.4, 0.6, 0.8, 1) on y axis while keeping the original values in the legend
                    ax2.set_yticks([np.log2(1e-4), np.log2(0.2), np.log2(0.4), np.log2(0.6), np.log2(0.8), np.log2(1)], ["~0", "0.2", "0.4", "0.6", "0.8", "1"])

        plt.subplots_adjust(right=0.8)
        if stats_dict_2 != {}:
            ax2.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        else:
            ax1.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)

        # optionally set ylim
        if ylim != [] and line:
            ax1.set_ylim(ylim)
            if stats_dict_2 != {}:
                ax2.set_ylim(ylim)
        
        # optionally add 0.05 significance line
        print(stats_summary_name)
        if pvallines[0]:
            print("adding pval line for left")
            if 'FDR' in stats_summary_name:
                pvalline_val = np.log2(0.05)
            else:
                pvalline_val = 0.05
            ax1.axhline(y=pvalline_val, color='r', linestyle='--')
            if stats_dict_2 != {} and pvallines[1]:
                print("adding pval line for right")
                ax2.axhline(y=pvalline_val, color='r', linestyle='--')

        # find parent dir of output_dir_l[0]
        plt.savefig(output_root+"_"+stats_summary_name+postfix, bbox_inches="tight")
        plt.clf()
        plt.close()

    # plot        
    if taxa_gt is not None:
        plot_stats('FDR_sensitivity', ["FDR", "Sensitivity"], global_methods_FDR_r2_l_dict, global_methods_sensitivity_r2_l_dict, postfix=postfix, ylim=[np.log(1e-6), 0], line=line, pvallines=[True, False])
    
    plot_stats('runtime', ["runtime"], global_methods_runtime_l_dict, postfix=postfix)
    plot_stats('auc and weighted f1', ["auc", "weighted f1"], global_methods_rf_auc_l_dict, global_methods_rf_f1_l_dict, postfix=postfix, ylim=[0.4, 1], line=line)
    plot_stats('weighted precision and weighted recall', ["weighted precision", "weighted recall"], global_methods_rf_precision_l_dict, global_methods_rf_recall_l_dict, postfix=postfix, ylim=[0.4, 1], line=line)
    plot_stats('shannon_pval', ["PERMANOVA batch Shannon pval", "PERMANOVA biovar Shannon pval"], global_methods_batch_shannon_pval_l_dict, global_methods_biovar_shannon_pval_l_dict, postfix=postfix, ylim=[0, 1], pvallines=[True, True],line=line)
    plot_stats('shannon_pval_rejection_proportions', ["PERMANOVA batch Shannon pval rejection proportion", "PERMANOVA biovar Shannon pval rejection proportion"], global_methods_batch_shannon_pval_rejection_proportions_l_dict, global_methods_biovar_shannon_pval_rejection_proportions_l_dict, postfix=postfix, ylim=[0, 1], line=line, pvallines=[True, False])

    if not demonstrate:
        if count_l[0]:
            plot_stats('PERMANOVA_batch_R2', ["PERMANOVA batch R2 (Aitchinson)", "PERMANOVA batch R2 (Bray-Curtis)"], global_methods_batch_aitch_r2_l_dict, global_methods_batch_bray_r2_l_dict, postfix=postfix, ylim=[0, 0.3], line = True, pvallines=[True, True])
            plot_stats('PERMANOVA_biovar_R2', ["PERMANOVA biovar R2 (Aitchinson)", "PERMANOVA biovar R2 (Bray-Curtis)"], global_methods_biovar_aitch_r2_l_dict, global_methods_biovar_bray_r2_l_dict, postfix=postfix, ylim=[0, 0.3], line=line, pvallines=[True, True])
        else:
            plot_stats('PERMANOVA_batch_R2', ["PERMANOVA batch R2 (Bray-Curtis)"], global_methods_batch_bray_r2_l_dict, postfix=postfix, ylim=[0, 0.3], line=line, pvallines=[True, True])
            plot_stats('PERMANOVA_biovar_R2', ["PERMANOVA biovar R2 (Bray-Curtis)"], global_methods_biovar_bray_r2_l_dict, postfix=postfix, ylim=[0, 0.3], line=line, pvallines=[True, True])
    else:
        if count_l[0]:
            plot_stats('PERMANOVA_batch_R2', ["PERMANOVA batch R2 (Aitchinson)", "PERMANOVA batch R2 (Bray-Curtis)"], global_methods_batch_aitch_r2_l_dict, global_methods_batch_bray_r2_l_dict, postfix=postfix, line=line, pvallines=[True, True])
            plot_stats('PERMANOVA_biovar_R2', ["PERMANOVA biovar R2 (Aitchinson)", "PERMANOVA biovar R2 (Bray-Curtis)"], global_methods_biovar_aitch_r2_l_dict, global_methods_biovar_bray_r2_l_dict, postfix=postfix, line=line, pvallines=[True, True])
        else:
            plot_stats('PERMANOVA_batch_R2', ["PERMANOVA batch R2 (Bray-Curtis)"], global_methods_batch_bray_r2_l_dict, postfix=postfix, line=line, pvallines=[True, True])
            plot_stats('PERMANOVA_biovar_R2', ["PERMANOVA biovar R2 (Bray-Curtis)"], global_methods_biovar_bray_r2_l_dict, postfix=postfix, line=line, pvallines=[True, True])        
    return



## simulation evaluation - MIDAS
## null data
## STEP 1. GENERATE DATA FROM DATABASE


def iterative_methods_running_evaluate(run_or_evaluate, datatype, or_l, cond_effect_val_l, batch_effect_val_l, 
            iter = 1,
            address_XY_dir_path = overall_path+'/simulation_data_updated_MIDAS_yesrelation_090723', 
            output_dir_path = overall_path+"/simulation_data_updated_output_count_yesrelation_090723", 
            eval_dir_path = overall_path+"/simulation_data_updated_eval_count_yesrelation_090723",
            methods_list = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"],
            binarizing_agent_biovar = 'H'):

    IDCol = 'subjectid_text'
    for odds_ratio in or_l:
        for cond_effect_val in cond_effect_val_l:
            for batch_effect_val in batch_effect_val_l:
                if cond_effect_val + batch_effect_val <= 1:
                    print("odds_ratio", odds_ratio)
                    print("cond_effect_val", cond_effect_val)
                    print("batch_effect_val", batch_effect_val)
                    print("iter", iter)

                    address_Y = address_XY_dir_path + '/ibd_150_meta_' + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val) + '_iter_' + str(iter) + '.csv'

                    if run_or_evaluate == "run":
                        ## STEP 2. RUNNING METHODS (FOR THE TWO RUNNING IN PYTHON)
                        if datatype == "count":
                            address_X = address_XY_dir_path + '/ibd_150_count_' + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val) + '_iter_' + str(iter) + '.csv'
                        else:
                            address_X = address_XY_dir_path + '/ibd_150_relab_' + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val) + '_iter_' + str(iter) + '.csv'
                        # address_Y = address_XY_dir + '/ibd_150_meta_' + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val) + '_iter_' + str(iter) + '.csv'
                        data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
                        output_root = output_dir_path + "/out_"+str(odds_ratio)+"_"+ str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter) 
                        if not os.path.exists(output_root):
                            os.makedirs(output_root)
                        if not os.path.exists(output_root+"/ibd_"+str(odds_ratio)+"_"+ str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter) + "_harmony_adjusted_count.csv") or os.path.getsize(output_root+"/ibd_"+str(odds_ratio)+"_"+ str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter) + "_harmony_adjusted_count.csv") == 0:
                            res_h, meta_data_h = generate_harmony_results(data_mat, meta_data, IDCol, ["batchid"],  output_root+"/ibd_"+str(odds_ratio)+ "_" +str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter) +"_harmony")
                        if not os.path.exists(output_root+"/ibd_"+str(odds_ratio)+"_"+ str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter) + "_percentile_norm.csv") or os.path.getsize(output_root+"/ibd_"+str(odds_ratio)+"_"+ str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter) + "_percentile_norm.csv") == 0:
                            percentile_norm(address_X, address_Y, "cond", "cond_1", "comma", output_root+"/ibd_"+str(odds_ratio)+ "_"+str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter), simulate=True)
                
                        # res_h, meta_data_h = generate_harmony_results(data_mat, meta_data, IDCol, ["batchid"],  output_root+"/ibd_"+str(odds_ratio)+ "_" +str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter) +"_harmony")
                        # percentile_norm(address_X, address_Y, "cond", "cond_1", "comma", output_root+"/ibd_"+str(odds_ratio)+ "_"+str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter), simulate=True)
                    elif run_or_evaluate == "evaluate":    
                        # # global stats dataframe generation
                        # input_frame_path = address_XY_dir_path+"/ibd_150_count_"+ str(odds_ratio)+ "_"+ str(cond_effect_val)+ "_"+ str(batch_effect_val)+ '_iter_'+ str(iter)+ ".csv"
                        # bio_var = "cond"
                        # dataset_name = "batch_0"
                        # global_eval_dataframe(input_frame_path, bio_var, dataset_name, methods_list, output_dir_path+"/out_" + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val)+ '_iter_' + str(iter)+"/",
                        #     output_dir_path = eval_dir_path+'/out_' + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val)+ '_iter_' + str(iter), simulate = True, taxa_gt = True,
                        #     datatype = datatype)
                        # check if the whole iteration finished
                        if not os.path.isfile(eval_dir_path+'/out_' + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val)+ '_iter_' + str(iter)+"/ibd_150_"+str(odds_ratio) + "_" + str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter)+'_multi_PCOA_both_batch.pdf'):
                            if cond_effect_val == 0:
                                taxa_gt = []
                            else:
                                # read taxa_gt as a list from taxa_gt_path
                                taxa_gt_path = address_XY_dir_path+'/ibd_150_id_cond_' + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val) + '_iter_' + str(iter) + '.txt'
                                f = open(taxa_gt_path, "r")
                                taxa_gt = f.readline().split(' ')[:-1]
                                taxa_gt = [int(taxa)-1 for taxa in taxa_gt]
                
                            ### STEP 3. EVALUATE ALL THE METHODS - counts, yes relation
                            
                            if iter == 1:
                                pipeline = "default"
                            else:
                                pipeline = "scaled"

                            df_l = []
                            metadata_l = []

                            # check if already exists
                            output_dir = eval_dir_path+'/out_' + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val)+ '_iter_' + str(iter)
                            if not os.path.exists(output_dir):
                                os.makedirs(output_dir)

                            # evaluate individual methods
                            for method in methods_list:   
                                print("___________________________________________________________")
                                print("odds_ratio", odds_ratio)
                                print("cond_effect_val", cond_effect_val)
                                print("batch_effect_val", batch_effect_val)
                                print("iter", iter)
                                print("method", method)                     
                                # check if already exists
                                if not os.path.exists(output_dir+'/'+method):
                                    os.makedirs(output_dir+'/'+method)
                                
                                if method == 'harmony':
                                    address_X = output_dir_path + "/out_"+str(odds_ratio)+"_"+ str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter) + "/ibd_"+ str(odds_ratio)+ "_" + str(cond_effect_val)+ "_"+ str(batch_effect_val)+ '_iter_'+ str(iter)+ "_harmony_adjusted_count.csv"
                                elif method == 'percentile_norm':
                                    address_X = output_dir_path + "/out_"+str(odds_ratio)+"_"+ str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter) + "/ibd_"+ str(odds_ratio)+ "_" + str(cond_effect_val)+ "_"+ str(batch_effect_val)+ '_iter_'+ str(iter)+ "_percentile_norm.csv"
                                elif method == 'nobc':
                                    address_X = address_XY_dir_path+"/ibd_150_count_"+ str(odds_ratio)+ "_"+ str(cond_effect_val)+ "_"+ str(batch_effect_val)+ '_iter_'+ str(iter)+ ".csv"
                                else:
                                    address_X = output_dir_path + "/out_"+str(odds_ratio)+"_"+ str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter) +"/ibd_" + str(odds_ratio)+ "_"+ str(cond_effect_val)+ "_"+ str(batch_effect_val)+ '_iter_'+ str(iter)+ "_"+method+".csv"
                                data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
                                df_l.append(data_mat)
                                metadata_l.append(meta_data)

                                if not os.path.exists(output_dir+'/'+method+'/' + method +'_summary.csv'):
                                    Evaluate(data_mat, meta_data, 'batchid', output_dir+'/'+method+'/'+method,
                                            "cond", 30, IDCol=IDCol, pipeline = pipeline, taxa_gt = taxa_gt, datatype = datatype, binarizing_agent_biovar = binarizing_agent_biovar, method = method)
                                else:
                                    print("already exists", output_dir+'/'+method+'/out_' + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val) + '_iter_' + str(iter) + '_summary.csv')
                            
                            # global stats dataframe generation
                            input_frame_path = address_XY_dir_path+"/ibd_150_count_"+ str(odds_ratio)+ "_"+ str(cond_effect_val)+ "_"+ str(batch_effect_val)+ '_iter_'+ str(iter)+ ".csv"
                            bio_var = "cond"
                            dataset_name = "batch_0"
                            global_eval_dataframe(input_frame_path, bio_var, dataset_name, methods_list, output_dir_path+"/out_" + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val)+ '_iter_' + str(iter)+"/",
                                output_dir_path = eval_dir_path+'/out_' + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val)+ '_iter_' + str(iter), simulate = True, taxa_gt = True,
                                datatype = datatype)

                            # plot global PCOA vis
                            if not os.path.exists(eval_dir_path+'/out_' + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val)+ '_iter_' + str(iter)+"/ibd_150_"+str(odds_ratio) + "_" + str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter)+'_multi_PCOA_both_batch.pdf'):
                                plot_PCOA_multiple("ibd_150_"+str(odds_ratio) + "_" + str(cond_effect_val) + "_" + str(batch_effect_val) + "_iter_" + str(iter), df_l, methods_list, metadata_l, used_var="batchid", 
                                    output_root= eval_dir_path+'/out_' + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val)+ '_iter_' + str(iter)+"/",
                                    datatype = datatype)
                        else:
                            # global stats dataframe generation
                            input_frame_path = address_XY_dir_path+"/ibd_150_count_"+ str(odds_ratio)+ "_"+ str(cond_effect_val)+ "_"+ str(batch_effect_val)+ '_iter_'+ str(iter)+ ".csv"
                            bio_var = "cond"
                            dataset_name = "batch_0"
                            global_eval_dataframe(input_frame_path, bio_var, dataset_name, methods_list, output_dir_path+"/out_" + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val)+ '_iter_' + str(iter)+"/",
                                output_dir_path = eval_dir_path+'/out_' + str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val)+ '_iter_' + str(iter), simulate = True, taxa_gt = True,
                                datatype = datatype)
                            print("current iteration is already done: "+str(odds_ratio) + '_' + str(cond_effect_val) + '_' + str(batch_effect_val))

    return

### parameter combination set-up
or_l = [1, 1.25, 1.5]
cond_effect_val_l = [0, 0.25, 0.5, 0.75, 1]
batch_effect_val_l = [0, 0.25, 0.5, 0.75, 1]

# get current path
current_path = os.path.dirname(os.path.abspath(__file__))
print("current_path", current_path)
overall_path = current_path + "/../.."
with open(f'{current_path}/../config.yml') as file:
    config_data = yaml.load(file, Loader=yaml.FullLoader)

used_R_methods = list(config_data['used_R_methods'])
used_Python_methods = list(config_data['used_Python_methods'])
methods = used_R_methods
methods.extend(used_Python_methods)


if ARGPARSE_SWITCH:
    GLOBAL_DATATYPE = args['datatype']
    # methods_list_dict = {'count': ["nobc", "harmony", "combat_seq", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"],
    #                     'relab': ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR_rel", "percentile_norm"]}
    related = args['related']
    binarizing_agent_biovar = args['agent']
    if args['option'] == 1:
        # run two python methods on the simulation data 
        iterative_methods_running_evaluate(run_or_evaluate = 'run', datatype = GLOBAL_DATATYPE, iter = int(args['iteration']), or_l = or_l, cond_effect_val_l = cond_effect_val_l, batch_effect_val_l = batch_effect_val_l, 
                    address_XY_dir_path = overall_path+f'/simulation_outputs/simulation_data_MIDAS_1000_{related}relation_102023', 
                    output_dir_path = overall_path+f"/simulation_outputs/simulation_data_output_{GLOBAL_DATATYPE}_{related}relation_102023", 
                    eval_dir_path = overall_path+f"/simulation_outputs/simulation_data_eval_{GLOBAL_DATATYPE}_{related}relation_102023",
                    methods_list = methods, binarizing_agent_biovar = binarizing_agent_biovar)
        print("DONEEEEEE")
    elif args['option'] == 2:
        # run simulation evaluation
        iterative_methods_running_evaluate(run_or_evaluate = 'evaluate', datatype = GLOBAL_DATATYPE, iter = int(args['iteration']), or_l = or_l, cond_effect_val_l = cond_effect_val_l, batch_effect_val_l = batch_effect_val_l, 
                    address_XY_dir_path = overall_path+f'/simulation_outputs/simulation_data_MIDAS_1000_{related}relation_102023', 
                    output_dir_path = overall_path+f"/simulation_outputs/simulation_data_output_{GLOBAL_DATATYPE}_{related}relation_102023", 
                    eval_dir_path = overall_path+f"/simulation_outputs/simulation_data_eval_{GLOBAL_DATATYPE}_{related}relation_102023",
                    methods_list = methods, binarizing_agent_biovar = binarizing_agent_biovar)
    elif args['option'] == 3:
        # visualize simulation stats
        eval_dir_path = overall_path+f"/simulation_outputs/simulation_data_eval_{GLOBAL_DATATYPE}_{related}relation_102023"
        output_dir_path = overall_path+f"/simulation_outputs/simulation_data_eval_{GLOBAL_DATATYPE}_{related}relation_102023"

        if not os.path.exists(eval_dir_path+f'/line_plots_{GLOBAL_DATATYPE}_{related}'):
            os.makedirs(eval_dir_path+f'/line_plots_{GLOBAL_DATATYPE}_{related}')

        # methods = methods_list_dict[GLOBAL_DATATYPE]
        # odds ratio == 1
        datasets = ["out_1_0_0", "out_1_0.25_0", "out_1_0.5_0", "out_1_0.75_0", "out_1_1_0", "out_1_0_0.25", "out_1_0.25_0.25", "out_1_0.5_0.25", 
                    "out_1_0.75_0.25", "out_1_0_0.5", "out_1_0.25_0.5", "out_1_0.5_0.5", "out_1_0_0.75", "out_1_0.25_0.75", "out_1_0_1"]
        output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
        counts_l = [GLOBAL_DATATYPE=='count']*len(datasets)
        visualize_simulation_stats(eval_dir_path+f'/line_plots_{GLOBAL_DATATYPE}_{related}/sim_1_all_bio', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, count_l = counts_l, simulate = True, dimensions = (20, 10), taxa_gt = True, postfix = '.pdf', sim_num_iters=num_iters)

        datasets = ["out_1_0.25_0", "out_1_0.5_0", "out_1_0.75_0", "out_1_1_0", "out_1_0.25_0.25", "out_1_0.5_0.25", 
                    "out_1_0.75_0.25", "out_1_0.25_0.5", "out_1_0.5_0.5", "out_1_0.25_0.75"]
        output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
        counts_l = [GLOBAL_DATATYPE=='count']*len(datasets)
        visualize_simulation_stats(eval_dir_path+f'/line_plots_{GLOBAL_DATATYPE}_{related}/sim_1_all_bio_alwaysbio', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, count_l = counts_l, simulate = True, dimensions = (20, 10), taxa_gt = True, postfix = '.pdf', sim_num_iters=num_iters)


        # odds ratio == 1.25
        datasets = ["out_1.25_0_0", "out_1.25_0.25_0", "out_1.25_0.5_0", "out_1.25_0.75_0", "out_1.25_1_0", "out_1.25_0_0.25", "out_1.25_0.25_0.25", "out_1.25_0.5_0.25",
                    "out_1.25_0.75_0.25", "out_1.25_0_0.5", "out_1.25_0.25_0.5", "out_1.25_0.5_0.5", "out_1.25_0_0.75", "out_1.25_0.25_0.75", "out_1.25_0_1"]
        output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
        counts_l = [GLOBAL_DATATYPE=='count']*len(datasets)
        visualize_simulation_stats(eval_dir_path+f'/line_plots_{GLOBAL_DATATYPE}_{related}/sim_1.25_all_bio', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, count_l = counts_l, simulate = True, dimensions = (20, 10), taxa_gt = True, postfix = '.pdf', sim_num_iters=num_iters)

        datasets = ["out_1.25_0.25_0", "out_1.25_0.5_0", "out_1.25_0.75_0", "out_1.25_1_0", "out_1.25_0.25_0.25", "out_1.25_0.5_0.25",
                    "out_1.25_0.75_0.25", "out_1.25_0.25_0.5", "out_1.25_0.5_0.5", "out_1.25_0.25_0.75"]
        output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
        counts_l = [GLOBAL_DATATYPE=='count']*len(datasets)
        visualize_simulation_stats(eval_dir_path+f'/line_plots_{GLOBAL_DATATYPE}_{related}/sim_1.25_all_bio_alwaysbio', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, count_l = counts_l, simulate = True, dimensions = (20, 10), taxa_gt = True, postfix = '.pdf', sim_num_iters=num_iters)


        # odds ratio == 1.5
        datasets = ["out_1.5_0_0", "out_1.5_0.25_0", "out_1.5_0.5_0", "out_1.5_0.75_0", "out_1.5_1_0", "out_1.5_0_0.25", "out_1.5_0.25_0.25", "out_1.5_0.5_0.25",
                    "out_1.5_0.75_0.25", "out_1.5_0_0.5", "out_1.5_0.25_0.5", "out_1.5_0.5_0.5", "out_1.5_0_0.75", "out_1.5_0.25_0.75", "out_1.5_0_1"]
        output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
        counts_l = [GLOBAL_DATATYPE=='count']*len(datasets)
        visualize_simulation_stats(eval_dir_path+f'/line_plots_{GLOBAL_DATATYPE}_{related}/sim_1.5_all_bio', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, count_l = counts_l, simulate = True, dimensions = (20, 10), taxa_gt = True, postfix = '.pdf', sim_num_iters=num_iters)
    
        datasets = ["out_1.5_0.25_0", "out_1.5_0.5_0", "out_1.5_0.75_0", "out_1.5_1_0", "out_1.5_0.25_0.25", "out_1.5_0.5_0.25",
                    "out_1.5_0.75_0.25", "out_1.5_0.25_0.5", "out_1.5_0.5_0.5",  "out_1.5_0.25_0.75" ]
        output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
        counts_l = [GLOBAL_DATATYPE=='count']*len(datasets)
        visualize_simulation_stats(eval_dir_path+f'/line_plots_{GLOBAL_DATATYPE}_{related}/sim_1.5_all_bio_alwaysbio', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, count_l = counts_l, simulate = True, dimensions = (20, 10), taxa_gt = True, postfix = '.pdf', sim_num_iters=num_iters)

    elif args['option'] == 4:
        IDCol = 'Sam_id'
        dataset_name = config_data['dataset_name']
        if 'microbiomeHD' in dataset_name:
            vars_use = ["Dataset"]
            condition_var = "DiseaseState"
        elif 'CMD' in dataset_name:
            vars_use = ["study_name"]
            condition_var = "disease"
        
        address_X = f'{config_data["src"]}/cleaned_data/{dataset_name}/{dataset_name}_count_data.csv'
        address_Y = f'{config_data["src"]}/cleaned_data/{dataset_name}/{dataset_name}_meta_data.csv'
        data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
        res_h, meta_data_h = generate_harmony_results(data_mat, meta_data, IDCol, vars_use, f'{config_data["post_integration_outputs"]}/{dataset_name}/{dataset_name}_harmony')
        percentile_norm(address_X, address_Y, condition_var, config_data['condition_value'], 'comma', f'{config_data["post_integration_outputs"]}/{dataset_name}/{dataset_name}')

    elif args['option'] == 5:

        IDCol = 'Sam_id'
        dataset_name = config_data['dataset_name']
        if 'microbiomeHD' in dataset_name:
            vars_use = "Dataset"
            condition_var = "DiseaseState"
        elif 'CMD' in dataset_name:
            vars_use = "study_name"
            condition_var = "disease"

        address_Y = f'{config_data["src"]}/cleaned_data/{dataset_name}/{dataset_name}_meta_data.csv'

        methods.append('nobc')
        print("methods: ", methods)
        print("methods: ", methods)
        print("methods: ", methods)
        print("methods: ", methods)
        print("methods: ", methods)
        print("methods: ", methods)
        df_l = []
        meta_data_l = []
        for method in methods:
            if method == 'nobc':
                address_X = f'{config_data["src"]}/cleaned_data/{dataset_name}/{dataset_name}_count_data.csv'
            elif method == 'harmony':
                address_X = f'{config_data["post_integration_outputs"]}/{dataset_name}/{dataset_name}_{method}_adjusted_count.csv'
            else:
                address_X = f'{config_data["post_integration_outputs"]}/{dataset_name}/{dataset_name}_{method}.csv'

            data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
            df_l.append(data_mat)
            meta_data_l.append(meta_data)
            Evaluate(data_mat, meta_data, vars_use, f'{config_data["evaluation_outputs"]}/{dataset_name}/output_{dataset_name}_{method}/{dataset_name}_{method}', condition_var, 30, IDCol=IDCol, datatype = config_data['datatype'], binarizing_agent_biovar = config_data['binarizing_agent'], method = method)
        
        ### global evaluation
        input_frame_path = f'{config_data["src"]}/cleaned_data/{dataset_name}/{dataset_name}_count_data.csv'
        global_eval_dataframe(input_frame_path, condition_var, dataset_name, methods, f'{config_data["post_integration_outputs"]}', f'{config_data["evaluation_outputs"]}/{dataset_name}', simulate = False, datatype = config_data['datatype'])
        
        ### multi-method plot
        plot_PCOA_multiple(dataset_name, df_l, methods, meta_data_l, used_var=vars_use, output_root= f'{config_data["evaluation_outputs"]}/{dataset_name}/', datatype = config_data['datatype'])

                # ## run two left-over python methods for rw datasets
        # ## RUN HARMONY/PERCENTILE_NORM FOR RW
        # # autism 2 microbiomeHD
        # ################################################################################
        # vars_use = ["Dataset"]
        # IDCol = 'Sam_id'
        # address_X = overall_path+"/mic_bc_benchmark/data/cleaned_data/autism_2_microbiomeHD/autism_2_microbiomeHD_count_data.csv"
        # address_Y = overall_path+"/mic_bc_benchmark/data/cleaned_data/autism_2_microbiomeHD/autism_2_microbiomeHD_meta_data.csv"
        # data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
        # res_h, meta_data_h = generate_harmony_results(data_mat, meta_data, IDCol, vars_use, overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/autism_2_microbiomeHD/"+"autism_2_microbiomeHD_harmony")
        # percentile_norm(address_X, address_Y, "DiseaseState", "ASD", "comma", overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD")

        # # cdi 3 microbiomeHD
        # ################################################################################
        # vars_use = ["Dataset"]
        # IDCol = 'Sam_id'
        # address_X = overall_path+"/mic_bc_benchmark/data/cleaned_data/cdi_3_microbiomeHD/cdi_3_microbiomeHD_count_data.csv"
        # address_Y = overall_path+"/mic_bc_benchmark/data/cleaned_data/cdi_3_microbiomeHD/cdi_3_microbiomeHD_meta_data.csv"
        # data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
        # res_h, meta_data_h = generate_harmony_results(data_mat, meta_data, IDCol, vars_use, overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/cdi_3_microbiomeHD/"+"cdi_3_microbiomeHD_harmony")
        # percentile_norm(address_X, address_Y, "DiseaseState", "CDI", "comma", overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD")

        # # ibd_3_CMD
        # ################################################################################
        # vars_use = ["study_name"]
        # IDCol = 'Sam_id'
        # address_X = overall_path+"/mic_bc_benchmark/data/cleaned_data/ibd_3_CMD/ibd_3_CMD_count_data.csv"
        # address_Y = overall_path+"/mic_bc_benchmark/data/cleaned_data/ibd_3_CMD/ibd_3_CMD_meta_data.csv"
        # data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
        # res_h, meta_data_h = generate_harmony_results(data_mat, meta_data, IDCol, vars_use, overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/ibd_3_CMD/"+"ibd_3_CMD_harmony")
        # percentile_norm(address_X, address_Y, "disease", "IBD", "comma", overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/ibd_3_CMD/ibd_3_CMD")

        # # crc_8_CMD
        # ################################################################################
        # vars_use = ["study_name"]
        # IDCol = 'Sam_id'
        # address_X = overall_path+"/mic_bc_benchmark/data/cleaned_data/crc_8_CMD/crc_8_CMD_count_data.csv"
        # address_Y = overall_path+"/mic_bc_benchmark/data/cleaned_data/crc_8_CMD/crc_8_CMD_meta_data.csv"
        # data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
        # res_h, meta_data_h = generate_harmony_results(data_mat, meta_data, IDCol, vars_use, overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/crc_8_CMD/"+"crc_8_CMD_harmony")
        # percentile_norm(address_X, address_Y, "disease", "CRC", "comma", overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/crc_8_CMD/crc_8_CMD")


        # # evaluate the method results on rw datasets
        # ## EVALUATE METHODS ON REAL-WORLD DATASET
        # ################################################################################
        # # autism 2 microbiomeHD
        # output_dir_path = overall_path+'/outputs/autism_2_microbiomeHD'
        # address_directory = overall_path+'/mic_bc_benchmark/data/cleaned_data/autism_2_microbiomeHD'
        # vars_use = ["Dataset"]
        # IDCol = 'Sam_id'

        # address_Y = overall_path+"/mic_bc_benchmark/data/cleaned_data/autism_2_microbiomeHD/autism_2_microbiomeHD_meta_data.csv"

        # # nobc
        # address_X = overall_path+"/mic_bc_benchmark/data/cleaned_data/autism_2_microbiomeHD/autism_2_microbiomeHD_count_data.csv"
        # data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
        # ### nobc
        # Evaluate(data_mat, meta_data, 'Dataset', output_dir_path + '/output_autism_2_microbiomeHD_nobc/autism_2_microbiomeHD_nobc', "DiseaseState", 30, [], 'Sam_id', method = 'nobc')

        # ### harmony
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD_harmony_adjusted_count.csv"
        # res_h, meta_data_h = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(res_h, meta_data_h, 'Dataset', output_dir_path + '/output_autism_2_microbiomeHD_harmony/autism_2_microbiomeHD_harmony', "DiseaseState", 30, [], 'Sam_id', method = 'harmony')

        # # benchmarking other methods: 
        # ### combat (combat_seq)
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD_combat_seq.csv"
        # data_mat_combat, meta_data_combat = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(data_mat_combat, meta_data_combat, 'Dataset', output_dir_path + '/output_autism_2_microbiomeHD_combat_seq/autism_2_microbiomeHD_combat_seq', "DiseaseState", 30, [], 'Sam_id', method = 'combat_seq')

        # ### limma
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD_limma.csv"
        # data_mat_limma, meta_data_limma = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(data_mat_limma, meta_data_limma, 'Dataset', output_dir_path + '/output_autism_2_microbiomeHD_limma/autism_2_microbiomeHD_limma', "DiseaseState", 30, [], 'Sam_id', method = 'limma')

        # ### MMUPHin
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD_MMUPHin.csv"
        # data_mat_mmuphin, meta_data_mmuphin = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(data_mat_mmuphin, meta_data_mmuphin, 'Dataset', output_dir_path + '/output_autism_2_microbiomeHD_MMUPHin/autism_2_microbiomeHD_MMUPHin', "DiseaseState", 30, [], 'Sam_id', method = 'MMUPHin')

        # ### ConQuR
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD_ConQuR.csv"
        # data_mat_conqur, meta_data_conqur = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(data_mat_conqur, meta_data_conqur, 'Dataset', output_dir_path + '/output_autism_2_microbiomeHD_ConQuR/autism_2_microbiomeHD_ConQuR', "DiseaseState", 30, [], 'Sam_id', method = 'ConQuR')

        # ### ConQuR_libsize
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD_ConQuR_libsize.csv"
        # data_mat_conqur_libsize, meta_data_conqur_libsize = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(data_mat_conqur_libsize, meta_data_conqur_libsize, 'Dataset', output_dir_path + '/output_autism_2_microbiomeHD_ConQuR_libsize/autism_2_microbiomeHD_ConQuR_libsize', "DiseaseState", 30, [], 'Sam_id', method = 'ConQuR_libsize')

        # ### percentile_norm
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD_percentile_norm.csv"
        # data_mat_percentile_norm, meta_data_percentile_norm = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(data_mat_percentile_norm, meta_data_percentile_norm, 'Dataset', output_dir_path + '/output_autism_2_microbiomeHD_percentile_norm/autism_2_microbiomeHD_percentile_norm', "DiseaseState", 30, [], 'Sam_id', method = 'percentile_norm')

        # ### global evaluation
        # input_frame_path = overall_path+"/mic_bc_benchmark/data/cleaned_data/autism_2_microbiomeHD/autism_2_microbiomeHD_count_data.csv"
        # bio_var = "DiseaseState"
        # dataset_name = "autism_2_microbiomeHD"
        # methods_list = ["nobc", "harmony", "combat_seq", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
        # global_eval_dataframe(input_frame_path, bio_var, dataset_name, methods_list, overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results", output_dir_path, simulate = False)

        # ### multi-method plot
        # df_l = [data_mat, res_h, data_mat_combat, data_mat_limma, data_mat_mmuphin, data_mat_conqur, data_mat_conqur_libsize, data_mat_percentile_norm]
        # methods = ["nobc", "harmony", "combat_seq", "limma", "MMUPhin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
        # meta_data_l = [meta_data, meta_data_h, meta_data_combat, meta_data_limma, meta_data_mmuphin, meta_data_conqur, meta_data_conqur_libsize, meta_data_percentile_norm]
        # plot_PCOA_multiple('autism_2_microbiomeHD', df_l, methods, meta_data_l, used_var="Dataset", output_root= output_dir_path + '/')

        # ##############################################################################
        # # cdi 3 microbiomeHD
        # output_dir_path = overall_path+'/mic_bc_benchmark/outputs/cdi_3_microbiomeHD'
        # address_directory = overall_path+'/mic_bc_benchmark/data/cleaned_data/cdi_3_microbiomeHD'
        # vars_use = ["Dataset"]
        # IDCol = 'Sam_id'

        # address_Y = overall_path+"/mic_bc_benchmark/data/cleaned_data/cdi_3_microbiomeHD/cdi_3_microbiomeHD_meta_data.csv"

        # # nobc
        # address_X = overall_path+"/mic_bc_benchmark/data/cleaned_data/cdi_3_microbiomeHD/cdi_3_microbiomeHD_count_data.csv"
        # data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
        # ### nobc
        # Evaluate(data_mat, meta_data, 'Dataset', output_dir_path + '/output_cdi_3_microbiomeHD_nobc/cdi_3_microbiomeHD_nobc', "DiseaseState", 30, [], 'Sam_id', method = 'nobc')

        # ### harmony
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD_harmony_adjusted_count.csv"
        # res_h, meta_data_h = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(res_h, meta_data_h, 'Dataset', output_dir_path + '/output_cdi_3_microbiomeHD_harmony/cdi_3_microbiomeHD_harmony', "DiseaseState", 30, [], 'Sam_id', method = 'harmony')

        # # benchmarking other methods:
        # ### combat (combat_seq)
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD_combat_seq.csv"
        # data_mat_combat, meta_data_combat = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(data_mat_combat, meta_data_combat, 'Dataset', output_dir_path + '/output_cdi_3_microbiomeHD_combat_seq/cdi_3_microbiomeHD_combat_seq', "DiseaseState", 30, [], 'Sam_id', method = 'combat_seq')

        # ### limma
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD_limma.csv"
        # data_mat_limma, meta_data_limma = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(data_mat_limma, meta_data_limma, 'Dataset', output_dir_path + '/output_cdi_3_microbiomeHD_limma/cdi_3_microbiomeHD_limma', "DiseaseState", 30, [], 'Sam_id', method = 'limma')

        # ### MMUPHin
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD_MMUPHin.csv"
        # data_mat_mmuphin, meta_data_mmuphin = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(data_mat_mmuphin, meta_data_mmuphin, 'Dataset', output_dir_path + '/output_cdi_3_microbiomeHD_MMUPHin/cdi_3_microbiomeHD_MMUPHin', "DiseaseState", 30, [], 'Sam_id', method = 'MMUPHin')

        # ### ConQuR
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD_ConQuR.csv"
        # data_mat_conqur, meta_data_conqur = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(data_mat_conqur, meta_data_conqur, 'Dataset', output_dir_path + '/output_cdi_3_microbiomeHD_ConQuR/cdi_3_microbiomeHD_ConQuR', "DiseaseState", 30, [], 'Sam_id', method = 'ConQuR')

        # ### ConQuR_libsize
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD_ConQuR_libsize.csv"
        # data_mat_conqur_libsize, meta_data_conqur_libsize = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(data_mat_conqur_libsize, meta_data_conqur_libsize, 'Dataset', output_dir_path + '/output_cdi_3_microbiomeHD_ConQuR_libsize/cdi_3_microbiomeHD_ConQuR_libsize', "DiseaseState", 30, [], 'Sam_id', method = 'ConQuR_libsize')

        # ### percentile_norm
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD_percentile_norm.csv"
        # data_mat_percentile_norm, meta_data_percentile_norm = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(data_mat_percentile_norm, meta_data_percentile_norm, 'Dataset', output_dir_path + '/output_cdi_3_microbiomeHD_percentile_norm/cdi_3_microbiomeHD_percentile_norm', "DiseaseState", 30, [], 'Sam_id', method = 'percentile_norm')

        # ### global evaluation
        # input_frame_path = overall_path+"/mic_bc_benchmark/data/cleaned_data/cdi_3_microbiomeHD/cdi_3_microbiomeHD_count_data.csv"
        # bio_var = "DiseaseState"
        # dataset_name = "cdi_3_microbiomeHD"
        # methods_list = ["nobc", "harmony", "combat_seq", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
        # global_eval_dataframe(input_frame_path, bio_var, dataset_name, methods_list, overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results", output_dir_path, simulate = False)

        # ### multi-method plot
        # df_l = [data_mat, res_h, data_mat_combat, data_mat_limma, data_mat_mmuphin, data_mat_conqur, data_mat_conqur_libsize, data_mat_percentile_norm]
        # methods = ["nobc", "harmony", "combat_seq", "limma", "MMUPhin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
        # meta_data_l = [meta_data, meta_data_h, meta_data_combat, meta_data_limma, meta_data_mmuphin, meta_data_conqur, meta_data_conqur_libsize, meta_data_percentile_norm]
        # plot_PCOA_multiple('cdi_3_microbiomeHD', df_l, methods, meta_data_l, used_var="Dataset", output_root= output_dir_path + '/')

        # ##############################################################################
        # # ibd_3_CMD
        # output_dir_path = overall_path+'/mic_bc_benchmark/outputs/ibd_3_CMD'
        # address_directory = overall_path+'/mic_bc_benchmark/data/cleaned_data/ibd_3_CMD'
        # vars_use = ["study_name"]
        # IDCol = 'Sam_id'

        # address_Y = overall_path+"/mic_bc_benchmark/data/cleaned_data/ibd_3_CMD/ibd_3_CMD_meta_data.csv"

        # # nobc
        # address_X = overall_path+"/mic_bc_benchmark/data/cleaned_data/ibd_3_CMD/ibd_3_CMD_count_data.csv"
        # data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
        # ### nobc
        # Evaluate(data_mat, meta_data, "study_name", output_dir_path + '/output_ibd_3_CMD_nobc/ibd_3_CMD_nobc', "disease", 30, ['gender', 'age_category'], 'Sam_id', datatype = 'relab', method = 'nobc', binarizing_agent_biovar = 'IBD')

        # ### harmony
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/ibd_3_CMD/ibd_3_CMD_harmony_adjusted_count.csv"
        # res_h, meta_data_h = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(res_h, meta_data_h, "study_name", output_dir_path + '/output_ibd_3_CMD_harmony/ibd_3_CMD_harmony', "disease", 30, ['gender', 'age_category'], 'Sam_id', datatype = 'relab', method = 'harmony', binarizing_agent_biovar = 'IBD')

        # # benchmarking other methods:
        # ### combat
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/ibd_3_CMD/ibd_3_CMD_combat.csv"
        # data_mat_combat, meta_data_combat = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(data_mat_combat, meta_data_combat, "study_name", output_dir_path + '/output_ibd_3_CMD_combat/ibd_3_CMD_combat', "disease", 30, ['gender', 'age_category'], 'Sam_id', datatype = 'relab', method = 'combat', binarizing_agent_biovar = 'IBD')

        # ### limma
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/ibd_3_CMD/ibd_3_CMD_limma.csv"
        # data_mat_limma, meta_data_limma = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(data_mat_limma, meta_data_limma, "study_name", output_dir_path + '/output_ibd_3_CMD_limma/ibd_3_CMD_limma', "disease", 30, ['gender', 'age_category'], 'Sam_id', datatype = 'relab', method = 'limma', binarizing_agent_biovar = 'IBD')

        # ### MMUPHin
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/ibd_3_CMD/ibd_3_CMD_MMUPHin.csv"
        # data_mat_mmuphin, meta_data_mmuphin = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(data_mat_mmuphin, meta_data_mmuphin, "study_name", output_dir_path + '/output_ibd_3_CMD_MMUPHin/ibd_3_CMD_MMUPHin', "disease", 30, ['gender', 'age_category'], 'Sam_id', datatype = 'relab', method = 'MMUPHin', binarizing_agent_biovar = 'IBD')

        # ### ConQuR_rel
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/ibd_3_CMD/ibd_3_CMD_ConQuR_rel.csv"
        # data_mat_conqur_rel, meta_data_conqur_rel = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(data_mat_conqur_rel, meta_data_conqur_rel, "study_name", output_dir_path + '/output_ibd_3_CMD_ConQuR_rel/ibd_3_CMD_ConQuR_rel', "disease", 30, ['gender', 'age_category'], 'Sam_id', datatype = 'relab', method = 'ConQuR', binarizing_agent_biovar = 'IBD')

        # ### percentile_norm
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/ibd_3_CMD/ibd_3_CMD_percentile_norm.csv"
        # data_mat_percentile_norm, meta_data_percentile_norm = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(data_mat_percentile_norm, meta_data_percentile_norm, "study_name", output_dir_path + '/output_ibd_3_CMD_percentile_norm/ibd_3_CMD_percentile_norm', "disease", 30, ['gender', 'age_category'], 'Sam_id', datatype = 'relab', method = 'percentile_norm', binarizing_agent_biovar = 'IBD')

        # ### global evaluation
        # input_frame_path = overall_path+"/mic_bc_benchmark/data/cleaned_data/ibd_3_CMD/ibd_3_CMD_count_data.csv"
        # bio_var = "disease"
        # dataset_name = "ibd_3_CMD"
        # methods_list = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR_rel", "percentile_norm"]
        # global_eval_dataframe(input_frame_path, bio_var, dataset_name, methods_list, overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results", output_dir_path, simulate = False)

        # ### multi-method plot
        # df_l = [data_mat, res_h, data_mat_combat, data_mat_limma, data_mat_mmuphin, data_mat_conqur_rel, data_mat_percentile_norm]
        # methods = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR_rel", "percentile_norm"]
        # meta_data_l = [meta_data, meta_data_h, meta_data_combat, meta_data_limma, meta_data_mmuphin, meta_data_conqur_rel, meta_data_percentile_norm]
        # plot_PCOA_multiple('ibd_3_CMD', df_l, methods, meta_data_l, used_var="study_name", output_root= output_dir_path + '/', datatype = 'relab')

        # ###############################################################################
        # # crc 8 CMD
        # output_dir_path = overall_path+'/mic_bc_benchmark/outputs/crc_8_CMD'
        # address_directory = overall_path+'/mic_bc_benchmark/data/cleaned_data/crc_8_CMD'
        # vars_use = ["study_name"]
        # IDCol = 'Sam_id'

        # address_Y = overall_path+"/mic_bc_benchmark/data/cleaned_data/crc_8_CMD/crc_8_CMD_meta_data.csv"

        # # nobc
        # address_X = overall_path+"/mic_bc_benchmark/data/cleaned_data/crc_8_CMD/crc_8_CMD_count_data.csv"
        # data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
        # ### nobc
        # Evaluate(data_mat, meta_data, "study_name", output_dir_path + '/output_crc_8_CMD_nobc/crc_8_CMD_nobc', "disease", 30, ['gender', 'age'], 'Sam_id', datatype = 'relab', method = 'nobc', binarizing_agent_biovar = 'adenoma')

        # ### harmony
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/crc_8_CMD/crc_8_CMD_harmony_adjusted_count.csv"
        # res_h, meta_data_h = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(res_h, meta_data_h, "study_name", output_dir_path + '/output_crc_8_CMD_harmony/crc_8_CMD_harmony', "disease", 30, ['gender', 'age'], 'Sam_id', datatype = 'relab', method = 'harmony', binarizing_agent_biovar = 'adenoma')

        # # benchmarking other methods:
        # ### combat
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/crc_8_CMD/crc_8_CMD_combat.csv"
        # data_mat_combat, meta_data_combat = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(data_mat_combat, meta_data_combat, "study_name", output_dir_path + '/output_crc_8_CMD_combat/crc_8_CMD_combat', "disease", 30, ['gender', 'age'], 'Sam_id', datatype = 'relab', method = 'combat', binarizing_agent_biovar = 'adenoma')

        # ### limma
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/crc_8_CMD/crc_8_CMD_limma.csv"
        # data_mat_limma, meta_data_limma = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(data_mat_limma, meta_data_limma, "study_name", output_dir_path + '/output_crc_8_CMD_limma/crc_8_CMD_limma', "disease", 30, ['gender', 'age'], 'Sam_id', datatype = 'relab', method = 'limma', binarizing_agent_biovar = 'adenoma')

        # ### MMUPHin
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/crc_8_CMD/crc_8_CMD_MMUPHin.csv"
        # data_mat_mmuphin, meta_data_mmuphin = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(data_mat_mmuphin, meta_data_mmuphin, "study_name", output_dir_path + '/output_crc_8_CMD_MMUPHin/crc_8_CMD_MMUPHin', "disease", 30, ['gender', 'age'], 'Sam_id', datatype = 'relab', method = 'MMUPHin', binarizing_agent_biovar = 'adenoma')

        # ### ConQuR
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/crc_8_CMD/crc_8_CMD_ConQuR_rel.csv"
        # data_mat_conqur_rel, meta_data_conqur_rel = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(data_mat_conqur_rel, meta_data_conqur_rel, "study_name", output_dir_path + '/output_crc_8_CMD_ConQuR_rel/crc_8_CMD_ConQuR_rel', "disease", 30, ['gender', 'age'], 'Sam_id', datatype = 'relab', method = 'ConQuR', binarizing_agent_biovar = 'adenoma')

        # ### percentile_norm
        # address_X = overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results/crc_8_CMD/crc_8_CMD_percentile_norm.csv"
        # data_mat_percentile_norm, meta_data_percentile_norm = load_results_from_benchmarked_methods(address_X, address_Y)
        # Evaluate(data_mat_percentile_norm, meta_data_percentile_norm, "study_name", output_dir_path + '/output_crc_8_CMD_percentile_norm/crc_8_CMD_percentile_norm', "disease", 30, ['gender', 'age'], 'Sam_id', datatype = 'relab', method = 'percentile_norm', binarizing_agent_biovar = 'adenoma')

        # ### global evaluation
        # input_frame_path = overall_path+"/mic_bc_benchmark/data/cleaned_data/crc_8_CMD/crc_8_CMD_count_data.csv"
        # bio_var = "disease"
        # dataset_name = "crc_8_CMD"
        # methods_list = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR_rel", "percentile_norm"]
        # global_eval_dataframe(input_frame_path, bio_var, dataset_name, methods_list, overall_path+"/mic_bc_benchmark/benchmark/benchmarked_results", output_dir_path, simulate = False)

        # ### multi-method plot
        # df_l = [data_mat, res_h, data_mat_combat, data_mat_limma, data_mat_mmuphin, data_mat_conqur_rel, data_mat_percentile_norm]
        # methods = ["nobc", "harmony", "combat", "limma", "MMUPHin", "ConQuR_rel", "percentile_norm"]
        # meta_data_l = [meta_data, meta_data_h, meta_data_combat, meta_data_limma, meta_data_mmuphin, meta_data_conqur_rel, meta_data_percentile_norm]
        # plot_PCOA_multiple('crc_8_CMD', df_l, methods, meta_data_l, used_var="study_name", output_root= output_dir_path + '/', datatype = 'relab')

# # ## VISUALIZE LINE PLOTS FOR 2 COUNT-TYPE RW DATASETS and 2 RELAB-TYPE RW DATASETS
# ##############################################################################
# output_dir_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/outputs'
# methods = ["nobc", "harmony", "combat_seq", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "percentile_norm"]
# datasets = ["autism_2_microbiomeHD", "cdi_3_microbiomeHD"]
# output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
# visualize_simulation_stats('/athena/linglab/scratch/chf4012/mic_bc_benchmark/outputs/rw_data_plots/count_rw', output_dir_l, datasets, methods, highlighted_method = "ConQuR", line = True, simulate = False, count_l = [True, True], postfix = '.pdf')

# output_dir_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/outputs'
# methods = ["nobc", "combat", "harmony", "limma", "MMUPHin", "ConQuR_rel", "percentile_norm"]
# datasets = ["ibd_3_CMD", "crc_8_CMD"]
# output_dir_l = [output_dir_path+'/'+dataset for dataset in datasets]
# visualize_simulation_stats('/athena/linglab/scratch/chf4012/mic_bc_benchmark/outputs/rw_data_plots/relab_rw', output_dir_l, datasets, methods, highlighted_method = "ConQuR_rel", line = True, count_l = [False, False], simulate = False, postfix = '.pdf')