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

from preprocessing import *
from harmonicMic import run_harmonicMic
from harmony import run_harmony
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
from sklearn.decomposition import PCA
import seaborn as sns
from sklearn.preprocessing import StandardScaler
import itertools
from scipy import stats
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.stats.composition import clr
from skbio.stats.ordination import pcoa
from scipy.spatial import distance
# from sklearn.metrics.pairwise import euclidean_distances
from scipy.spatial.distance import pdist, squareform
import os
import time
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, KFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, f1_score, precision_score, recall_score
from skbio.stats.distance import DissimilarityMatrix, permanova
from sklearn.preprocessing import OneHotEncoder
import time

# data_mat, meta_data = load_data(address_X, address_Y, IDCol, index_col, PCA_first = False)
def generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root, option = "harmonicMic", PCA_first = False, diversity_weight = None):
    start_time = time.time()
    PCA_first_str = ""
    if option == "harmonicMic":
        ho = run_harmonicMic(data_mat, meta_data, vars_use, diversity_weight = diversity_weight)
    elif option == 'harmony':
        if PCA_first:
            PCA_first_str = "_PCA_first"
            pca_results  = PCA(n_components = 30, random_state=np.random.RandomState(88))
            pca_results.fit(data_mat)
            data_mat_1 = pca_results.transform(data_mat)
            data_mat = pd.DataFrame(data_mat_1, index = data_mat.index, columns= ['PC'+str(i) for i in range(30)])
        ho = run_harmony(data_mat, meta_data, vars_use)
    else:
        raise ValueError('Currently only support harmonicMic and harmony.')
    res = pd.DataFrame(ho.Z_corr)
    elapsed = time.time() - start_time
    print("time")
    print(elapsed)
    # give the index back
    res.index = data_mat.columns
    res.columns = list(meta_data[IDCol])
    
    if diversity_weight is not None:
        weight_str = "_weighted"
    else:
        weight_str = ""
    with open(output_root+"_"+option+PCA_first_str+weight_str+"_elapsed_time.txt", "w") as text_file:
        print(option, PCA_first_str, str(elapsed), "seconds", file=text_file)
        print("\n", file=text_file)

    res.T.to_csv(output_root+"_"+option+PCA_first_str+weight_str+"_adjusted_count.csv")
    return res.T, meta_data


def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    # source: https://matplotlib.org/stable/gallery/statistics/confidence_ellipse.html
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)


def run_eval(batch_corrected_df, meta_data, batch_var, output_root, bio_var = False, n_pc=30, covar = False, IDCol = None):
    a = Evaluate(batch_corrected_df, meta_data, batch_var, output_root, bio_var, n_pc, covar, IDCol = None)
    return

class Evaluate(object):
    def __init__(
            self, batch_corrected_df, meta_data, batch_var, output_root, bio_var = False, n_pc=30, covar = False, IDCol = None, test_percent = 0.2
    ):
        self.batch_corrected_df = batch_corrected_df
        self.meta_data = meta_data
        self.batch_var = batch_var
        self.output_root = output_root
        self.bio_var = bio_var # var retaining the biological info
        self.n_pc = n_pc
        self.covar = covar
        self.rng = np.random.RandomState(100)
        self.IDCol = IDCol
        self.test_percent = test_percent
        # make directory if not exists
        directory_path = "/".join(output_root.split("/")[:-1])
        if not os.path.exists(directory_path):
            os.makedirs(directory_path)
        # functions executed
        ## version 1, for batch evaluation
        self.alpha_beta_diversity_and_tests(self.batch_var) # has to happen before standard scaler

        self.calculate_R_sq(self.bio_var)
        self.calculate_R_sq(self.batch_var)

        self.standard_scaler()
        df = self.PCA_vis()
        if covar != False:
            self.covar == covar
            self.allPCs_covar_kw_test(df)
        self.PC01_kw_tests(df)
        self.predict_difference_RF()
        global_df = pd.DataFrame({"PC_0_p": [self.global_PC0_p], "PC_1_p":[self.global_PC1_p], "PC_0_bc_distance": [self.bc_global_pval_l[0]], "PC_1_bc_distance": [self.bc_global_pval_l[1]]})
        global_df.to_csv(self.output_root+"_global_batches_eval.csv")

    def standard_scaler(self):
        source_df = self.batch_corrected_df
        scaler = StandardScaler()
        data = source_df.T.to_numpy()
        scaler.fit(data)
        ss_scaled = scaler.transform(data).T
        self.batch_corrected_df = pd.DataFrame(ss_scaled, columns = source_df.columns, index=source_df.index)
        
    def PCA_vis(self):
        source_df = self.batch_corrected_df.copy()

        # zero mean for PCA
        source_df = source_df-source_df.mean()
        # fit PCA
        pca_results  = PCA(n_components = self.n_pc, random_state=self.rng)
        pca_results.fit(source_df)
        variances = pca_results.explained_variance_ratio_

        # ecdf plot
        ecdf_var = [variances[:i].sum() for i in range(len(variances))]
        fig, ax1 = plt.subplots(1,1,figsize = (20, 8))
        ax1.plot(range(len(ecdf_var)),ecdf_var)
        plt.xticks(range(len(ecdf_var)), range(len(ecdf_var)))
        plt.yticks([0.1*i for i in range(11)])
        plt.show()
        plt.savefig(self.output_root+"_ecdf.png")

        # pca visualization
        transformed = pca_results.transform(source_df)
        df = pd.DataFrame()
        for i in range(self.n_pc):
            df[f"PC{i}"] = transformed[:, i] 
        df['batches'] = list(self.meta_data[self.batch_var])
        df[self.bio_var] = list(self.meta_data[self.bio_var])
        ## pca visualization for bio_vars
        ### fetch info for bio_va
        all_colors_used = self.rng.uniform(0, 1, 3*len(np.unique(list(df[self.bio_var])))).tolist()
        colors = {var: all_colors_used[idx*3:idx*3+3] for idx,var in enumerate(np.unique(list(df[self.bio_var])))} # unlikely to get repeated colors
        fig, ax =  plt.subplots(1, 1, sharex=True, sharey=True, figsize=(20, 10))
        sns.scatterplot(df["PC0"] , df["PC1"], hue = df[self.bio_var], hue_order = np.unique(list(df[self.bio_var])), style = df["batches"], s = 100,ax = ax, cmap = "tab10", x_jitter = True, palette = colors)
        ax.legend(bbox_to_anchor=(1.04,1), loc="upper left")
        for index, current_biovar in enumerate(np.unique(list(df[self.bio_var]))):
            currentDF = df.loc[df[self.bio_var] == current_biovar]
            confidence_ellipse(currentDF['PC0'], currentDF['PC1'], ax, n_std=1.5, edgecolor=list(colors.values())[index])
        plt.savefig(self.output_root+"_PCA_"+self.bio_var+".png")

        ## pca visualization for batches
        fig, ax =  plt.subplots(1, 1, sharex=True, sharey=True, figsize=(20, 10))
        sns.scatterplot(df["PC0"] , df["PC1"], hue = df[self.bio_var], hue_order = np.unique(list(df[self.bio_var])), style = df["batches"], s = 100,ax = ax, cmap = "tab10", x_jitter = True, palette = colors)
        ax.legend(bbox_to_anchor=(1.04,1), loc="upper left")
        ### fetch info for bio_var
        for index, current_batch in enumerate(np.unique(list(df["batches"]))):
            currentDF = df.loc[df["batches"] == current_batch]
            confidence_ellipse(currentDF['PC0'], currentDF['PC1'], ax, n_std=1.5, edgecolor="grey")
        plt.savefig(self.output_root+"_PCA_batches.png")
        
        if self.covar != False:
            df[self.covar] = list(self.meta_data[self.covar])
        return df

    def PC01_kw_tests(self, df):
        # to investigate whether the bio_var info is retained
        bio_var_l = list(df[self.bio_var])
        bio_var_l = [x for x in bio_var_l if str(x) != 'nan']
        bio_var_l = list(np.unique(bio_var_l))
        
        # generate global metrics
        with open(self.output_root+"_"+self.bio_var+"_kw_pvals.txt", "w") as text_file:
            print("global batch kw p-vals \n", file=text_file)
            print("\n", file=text_file)

            data_PC0 = [df.loc[df[self.bio_var]==var]["PC0"].values for var in bio_var_l]
            data_PC1 = [df.loc[df[self.bio_var]==var]["PC1"].values for var in bio_var_l]
            
            self.global_PC0_p = stats.kruskal(*data_PC0)[1]
            self.global_PC1_p = stats.kruskal(*data_PC1)[1]
            print(". PC0", "across all biovar options, p-val = ", str(self.global_PC0_p), "\n",file=text_file)
            print(". PC1", "across all biovar options, p-val = ", str(self.global_PC1_p), "\n",file=text_file)
            
        dim = len(np.unique(bio_var_l))
        kw_heatmap_array = np.full((dim, dim), 1, dtype=float)
        for pair in itertools.combinations(bio_var_l, 2):
            data_for_kw_PC0 = [df.loc[df[self.bio_var]==var]["PC0"].values for var in pair]
            data_for_kw_PC1 = [df.loc[df[self.bio_var]==var]["PC1"].values for var in pair]
            PC0_p = stats.kruskal(*data_for_kw_PC0)[1]
            PC1_p = stats.kruskal(*data_for_kw_PC1)[1]
            print("PC0", pair[0], pair[1], PC0_p)
            print("PC1", pair[0], pair[1], PC1_p)
            print([bio_var_l.index(pair[0]), bio_var_l.index(pair[1])], kw_heatmap_array[bio_var_l.index(pair[0]), bio_var_l.index(pair[1])])
            kw_heatmap_array[bio_var_l.index(pair[0]), bio_var_l.index(pair[1])]=PC1_p
            kw_heatmap_array[bio_var_l.index(pair[1]), bio_var_l.index(pair[0])]=PC0_p
        fig, ax = plt.subplots()
        ax.plot([0, 1], [0, 1], transform=ax.transAxes, color="red")
        im = ax.imshow((-np.log2(kw_heatmap_array)).T, cmap="Oranges", rasterized=True,origin='lower')
        # Show all ticks and label them with the respective list entries
        ax.set_xticks(np.arange(len(bio_var_l)), labels=bio_var_l)
        ax.set_xlabel("PC0 (lower diagonal)")
        ax.set_yticks(np.arange(len(bio_var_l)), labels=bio_var_l)
        ax.set_ylabel("PC1 (upper diagonal)")
        # Rotate the tidata["PC0"].valuesck labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                rotation_mode="anchor")

        # Loop over data dimensionsif i!=j: and create text annotations.
        for i in range(len(bio_var_l)):
            for j in range(len(bio_var_l)):
                if i!=j:
                    text = ax.text(i, j, round(-np.log2(kw_heatmap_array)[i, j], 5),
                                ha="center", va="center", color="black", fontsize=8)

        ax.set_title("Kruskal-Wallis -log2(p-val) Heatmap")
        fig.tight_layout()
        plt.savefig(self.output_root+"_PC01_kw_tests.png")


        # to investigate whether batch effect is eliminated thouroughly
        batches_l = list(df["batches"])
        batches_l = [x for x in batches_l if str(x) != 'nan']
        batches_l = list(np.unique(batches_l))
        with open(self.output_root+"_batches_kw_pvals.txt", "w") as text_file:
            print("global batch kw p-vals \n", file=text_file)
            print("\n", file=text_file)

            data_PC0 = [df.loc[df["batches"]== batch]["PC0"].values for batch in batches_l]
            data_PC1 = [df.loc[df["batches"]== batch]["PC1"].values for batch in batches_l]
            
            self.global_PC0_p = stats.kruskal(*data_PC0)[1]
            self.global_PC1_p = stats.kruskal(*data_PC1)[1]
            print(". PC0", "across all batches, p-val = ", str(self.global_PC0_p), "\n",file=text_file)
            print(". PC1", "across all batches, p-val = ", str(self.global_PC1_p), "\n",file=text_file)
            # global_PC_0_bc_distance = 0
            # global_PC_1_bc_distance = 0

            print("global pair-wise batch kw p-vals \n", file=text_file)
        
            # num_of_pairs = 0
            for pair_batch in itertools.combinations(batches_l, 2):
                # num_of_pairs += 1
                print("current batch pair", str(pair_batch), '\n', file=text_file)
                current_df = df.loc[df["batches"].isin(pair_batch)]
                data_batch1 = current_df.loc[current_df["batches"]==pair_batch[0]]
                data_batch2 = current_df.loc[current_df["batches"]==pair_batch[1]]

                PC0_p = stats.kruskal(data_batch1["PC0"].values,data_batch2["PC0"].values)[1]
                PC1_p = stats.kruskal(data_batch1["PC1"].values,data_batch2["PC1"].values)[1]
                # global_PC_0_bc_distance += distance.braycurtis(data_batch1["PC0"].values,data_batch2["PC0"].values)
                # global_PC_1_bc_distance += distance.braycurtis(data_batch1["PC1"].values,data_batch2["PC1"].values)

                print(". PC0", str(pair_batch), PC0_p, "# samples in batch "+pair_batch[0]+":", len(data_batch1["PC0"].values), "# samples in batch "+pair_batch[1]+":", len(data_batch2["PC0"].values), file=text_file)
                print(". PC1", str(pair_batch), PC1_p, "# samples in batch "+pair_batch[0]+":", len(data_batch1["PC1"].values), "# samples in batch "+pair_batch[1]+":", len(data_batch2["PC1"].values), file=text_file)

            # global_PC_0_bc_distance = global_PC_0_bc_distance/num_of_pairs
            # global_PC_1_bc_distance = global_PC_1_bc_distance/num_of_pairs
            
            # global_df = pd.DataFrame({"PC_0_p": [global_PC0_p], "PC_1_p":[global_PC1_p], "PC_0_bc_distance": [global_PC_0_bc_distance], "PC_1_bc_distance": [global_PC_1_bc_distance]})
            # global_df.to_csv(self.output_root+"_global_batches_eval.csv")
            
            print("batch kw p-vals by bio_var \n", file=text_file)
            for var in bio_var_l:
                print(". bio_var == "+var, file=text_file)
                current_df_biovar = df.loc[df[self.bio_var] == var]
                for pair_batch in itertools.combinations(batches_l, 2):
                    current_df = current_df_biovar.loc[current_df_biovar["batches"].isin(pair_batch)]
                    data_batch1 = current_df.loc[current_df["batches"]==pair_batch[0]]
                    data_batch2 = current_df.loc[current_df["batches"]==pair_batch[1]]

                    PC0_p = stats.kruskal(data_batch1["PC0"].values,data_batch2["PC0"].values)[1]
                    PC1_p = stats.kruskal(data_batch1["PC1"].values,data_batch2["PC1"].values)[1]

                    print(".   PC0", str(pair_batch), PC0_p, "# new samples:", len(data_batch1["PC0"].values), "# old samples:", len(data_batch2["PC0"].values), file=text_file)
                    print(".   PC1", str(pair_batch), PC1_p, "# new samples:", len(data_batch1["PC1"].values), "# old samples:", len(data_batch2["PC1"].values), file=text_file)

    def allPCs_covar_kw_test(self, df):
        # to debug
        # lets do 30PCs for now
        fig, axes = plt.subplots(5, 6, sharex=True, sharey=True, figsize=(35, 35))
        current_ax_index = [0, 0]
        covar_unique = list(df[self.covar].values)
        covar_unique = [x for x in covar_unique if str(x) != 'nan']
        covar_unique = np.unique(covar_unique).tolist()
        for i in range(30):
            a = sns.scatterplot(df["PC0"],df["PC"+str(i)], hue = df[self.covar], style = df["batches"], s = 100,ax = axes[current_ax_index[0], current_ax_index[1]], cmap = "tab10", x_jitter = True)
            data_for_kw = [df.loc[df[self.covar]==var]["PC"+str(i)].values for var in covar_unique]
            PC1_p = stats.kruskal(*data_for_kw)[1]
            a.set_title("PC"+str(i)+" kw_pval="+str(round(PC1_p, 5)))
            if current_ax_index[1] == 5:
                current_ax_index = [current_ax_index[0]+1, 0]
            else:
                current_ax_index[1] +=1
            print("PC"+str(i), PC1_p)
        
        plt.savefig(self.output_root+"_allPCs_covar_kw_tests.png")

    def calculate_R_sq(self, R_sq_var):
        data = np.array(self.batch_corrected_df)
        data = np.where(data<np.percentile(data.flatten(), 0.01), 0, data)
        data = data+np.abs(np.min(data))
        ids = list(self.meta_data[self.IDCol])
        np.savetxt('data.out', data, delimiter=',')
        bc_div_bray = beta_diversity("braycurtis", data, ids)
        # compute aitchison distance
        data_clr = clr(data)
        data_clr = np.nan_to_num(data_clr)
        np.savetxt('data_clr.out', data_clr, delimiter=',')
        # bc_div_aitch = beta_diversity("euclidean", data_clr, ids)
        bc_div_aitch = pdist(data_clr)
        bc_div_aitch = squareform(bc_div_aitch)
        np.savetxt('bc_div_aitch.out', bc_div_aitch, delimiter=',')
        print(bc_div_aitch.shape)
        # bc_df_bray = pd.DataFrame(bc_div_bray.data, index=list(self.meta_data[self.IDCol]), columns=list(self.meta_data[self.IDCol]))
        
        # calculate R^2
        permanova_results_bray = permanova(bc_div_bray, self.meta_data, column=R_sq_var, permutations=999)
        # permanova_results_aitch = permanova(DissimilarityMatrix(bc_div_aitch), self.meta_data, column=R_sq_var, permutations=999)

        bray_r2 = 1 - 1 / (1 + permanova_results_bray[4] * permanova_results_bray[3] / (permanova_results_bray[2] - permanova_results_bray[3] - 1))
        # aitch_r2 = 1 - 1 / (1 + permanova_results_aitch[4] * permanova_results_aitch[3] / (permanova_results_aitch[2] - permanova_results_aitch[3] - 1))
        print(R_sq_var, "bray_r2", bray_r2)
        # print("aitch_r2", aitch_r2)
        return
    
    def alpha_beta_diversity_and_tests(self, test_var):
        data = np.array(self.batch_corrected_df)
        data = np.where(data<np.percentile(data.flatten(), 0.01), 0, data)
        data = data+np.abs(np.min(data))
        ids = list(self.meta_data[self.IDCol])
        print(data.shape)
        print("ids", len(ids))
        shannon_div = alpha_diversity('shannon', data, ids)
        bc_div = beta_diversity("braycurtis", data, ids)
        bc_df = pd.DataFrame(bc_div.data, index=list(self.meta_data[self.IDCol]), columns=list(self.meta_data[self.IDCol]))
        bc_pc = pcoa(bc_div)

        # visualize and save
        import matplotlib.pyplot as plt
        shannon_df = shannon_div.to_frame()
        shannon_df.columns = ["shannon"]
        shannon_df[test_var] = list(self.meta_data[test_var])
        shannon_fig = shannon_df.boxplot(column='shannon', by=test_var)
        shannon_fig.figure.savefig(self.output_root+"_"+test_var+"_alpha_shannon_pcoa.png")
        bc_metadata = self.meta_data
        bc_metadata.index = bc_metadata[self.IDCol]
        bc_metadata[test_var] = bc_metadata[test_var].fillna("unknown")
        bc_fig = bc_pc.plot(bc_metadata, test_var,
                 axis_labels=('PC 1', 'PC 2', 'PC 3'),
                 title='Samples colored by '+test_var, cmap='jet', s=50)
        bc_fig.figure.savefig(self.output_root+"_"+test_var+"_beta_bc_pcoa.png")

        # get list of unique test_var options
        test_var_col_l = list(self.meta_data[test_var])
        test_var_l = list(np.unique(self.meta_data[test_var].dropna()))

        # conduct statistical tests
        ## 1. alpha diversity
        alpha_kruskal_data = []
        alpha_pairwise_pval_dict = {}
        for var in test_var_l:
            removed_indices = [i for i, e in enumerate(test_var_col_l) if e != var]
            removed_indices = [var for i, var in enumerate(list(self.meta_data[self.IDCol])) if i in removed_indices]
            current_data_condition = list(shannon_div.drop(index = removed_indices))
            alpha_kruskal_data.append(current_data_condition)
        # calculate global kw p_val
        shannon_global_pval = stats.kruskal(*alpha_kruskal_data)[1]
        print("shannon_global_pval", shannon_global_pval)
        # calculate pairwise kw p_val
        kw_data_pair_l = list(itertools.combinations(alpha_kruskal_data, 2))
        for index, pair in enumerate(itertools.combinations(test_var_l, 2)):
            alpha_pairwise_pval_dict[pair] = stats.kruskal(*kw_data_pair_l[index])[1]
        print("shannon_pairwise_pval", alpha_pairwise_pval_dict)

        ## 2. beta diversity
        bc_global_pval_l = []
        bc_pairwise_pval_l = []
        for pc in ['PC1', 'PC2', 'PC3']:
            beta_kruskal_data = []
            bc_pairwise_pval_dict = {}
            current_data = bc_pc.samples[pc]
            for var in test_var_l:
                removed_indices = [i for i, e in enumerate(test_var_col_l) if e != var]
                removed_indices = [var for i, var in enumerate(list(self.meta_data[self.IDCol])) if i in removed_indices]
                current_data_condition = list(current_data.drop(index = removed_indices))
                beta_kruskal_data.append(current_data_condition)
            # calculate global kw p_val
            bc_global_pval_l.append(stats.kruskal(*beta_kruskal_data)[1])
            # calculate pairwise kw p_val
            kw_data_pair_l = list(itertools.combinations(beta_kruskal_data, 2))
            for index, pair in enumerate(itertools.combinations(test_var_l, 2)):
                bc_pairwise_pval_dict[pair] = stats.kruskal(*kw_data_pair_l[index])[1]
            bc_pairwise_pval_l.append(bc_pairwise_pval_dict)
        print("bray_curtis_global_pval, by PCs", bc_global_pval_l)
        print("bray_curtis_pairwise_pval, by PCs", bc_pairwise_pval_l)
        
        # return dataframes to csv for self checks
        ## alpha and beta diversity data themselves
        shannon_df.to_csv(self.output_root+"_shannon_df.csv")
        bc_df.to_csv(self.output_root+"_"+test_var+"_bray_curtis_dissimilaritymatrix.csv")
        ## kw significance testing to txt file
        with open(self.output_root+"_"+test_var+"_shannon_kw_pvals.txt", "w") as text_file:
            print("shannon_global_pval", shannon_global_pval, "\n", file=text_file)
            print("shannon_pairwise_pval \n", file=text_file)
            print(alpha_pairwise_pval_dict, file=text_file)
        with open(self.output_root+"_"+test_var+"_bray_curtis_kw_pvals.txt", "w") as text_file:
            print("bray_curtis_global_pvals, by PCs", bc_global_pval_l, "\n", file=text_file)
            print("bray_curtis_pairwise_pval \n", file=text_file)
            for index, pval_dict in enumerate(bc_pairwise_pval_l):
                print("PC"+str(index+1)+"\n", file=text_file)
                print(pval_dict, file=text_file)
                print("\n", file=text_file)

        if test_var == self.batch_var:
            self.bc_global_pval_l = bc_global_pval_l
        return 
    
    def predict_difference_RF(self):
        # get the data
        used_x = self.batch_corrected_df.copy() # note that standard scaler is already conducted
        used_y = list(self.meta_data[self.bio_var])

        # Creating the Training and Test set from data
        # X_train, X_test, y_train, y_test = train_test_split(used_x, used_y, test_size = self.test_percent, random_state = self.rng)
        kf = KFold(n_splits=5, random_state=self.rng, shuffle = True)# Define the split - into n_splits folds 
        kf.get_n_splits(used_x, used_y) # returns the number of splitting iterations in the cross-validator

        # train random forest classifier and evaluate
        used_x.reset_index(drop=True, inplace=True)
        eval_df = pd.DataFrame(columns = ["fold", "average_acc", "macro_prec", "weighted_prec", "macro_recall", "weighted_recall", "macro_f1", "weighted_f1", "auc",  "baseline_likelihood"])
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
            average_acc = sum(y_pred==y_test)/len(y_test)

            macro_precision = precision_score(y_test, y_pred, average = 'macro')
            weighted_precision = precision_score(y_test, y_pred, average = 'weighted')
            macro_recall = recall_score(y_test, y_pred, average = 'macro')
            weighted_recall = recall_score(y_test, y_pred, average = 'weighted')
            macro_f1 = f1_score(y_test, y_pred, average = 'macro')
            weighted_f1 = f1_score(y_test, y_pred, average = 'weighted')
            
            enc = OneHotEncoder(handle_unknown='ignore')
            enc.fit(np.array(y_test).reshape(-1, 1))
            y_test_oh = enc.transform(np.array(y_test).reshape(-1, 1)).toarray()
            y_pred_oh = enc.transform(np.array(y_pred).reshape(-1, 1)).toarray()
            auc = roc_auc_score(y_test_oh,  y_pred_oh)
            # find the most common element in y_test
            most_common_element = max(set(y_test), key = y_test.count)
            baseline_likelihood = y_test.count(most_common_element)/len(y_test)
            eval_df.loc[len(eval_df)] = [i+1, average_acc, macro_precision, weighted_precision, macro_recall, weighted_recall, macro_f1, weighted_f1, auc, baseline_likelihood]
        eval_df.to_csv(self.output_root+"_"+self.bio_var+"_rf_evaluate.csv")
        return
        

def global_eval_dataframe(input_frame_path, bio_var, dataset_name, methods_list, output_dir_path = "."):
    # fetch dimension (number of samples and number of features)
    data_mat = pd.read_csv(input_frame_path, index_col=0)
    num_of_samples = data_mat.shape[0]
    num_of_taxa = data_mat.shape[0]

    # rest of the stats
    method_dict = {}
    for method in methods_list:
        current_method_dict = {}
        current_dir_path = output_dir_path + "/output_" + dataset_name+"_"+method 
        current_dir_file_names = os.listdir(current_dir_path)

        # fetch between batch kw p-val (all methods + nbc)
        current_batch_kw_path = [result for result in current_dir_file_names if "batches_kw_pvals" in result][0]
        with open(current_dir_path+'/'+current_batch_kw_path) as f:
            lines = f.readlines()
        
        global_batch_p_val_PC0 = float([line for line in lines if "PC0 across all batches" in line][0].split(" ")[-2])
        global_batch_p_val_PC1 = float([line for line in lines if "PC1 across all batches" in line][0].split(" ")[-2])

        # fetch between biological variable kw pval (all methods + nbc)
        current_biovar_kw_path = [result for result in current_dir_file_names if bio_var+"_kw_pvals" in result][0]
        with open(current_dir_path+'/'+current_biovar_kw_path) as f:
            lines = f.readlines()
        
        global_biovar_p_val_PC0 = float([line for line in lines if "PC0 across all biovar options" in line][0].split(" ")[-2])
        global_biovar_p_val_PC1 = float([line for line in lines if "PC1 across all biovar options" in line][0].split(" ")[-2])

        current_method_dict.update({"global_batch_p_val_PC0": global_batch_p_val_PC0})
        current_method_dict.update({"global_batch_p_val_PC1": global_batch_p_val_PC1})
        current_method_dict.update({"global_biovar_p_val_PC0": global_biovar_p_val_PC0})
        current_method_dict.update({"global_biovar_p_val_PC1": global_biovar_p_val_PC1})
        method_dict[method] = current_method_dict

    if "combat" not in method_dict:
        method_dict["combat"] = {"global_batch_p_val_PC0": "NA", "global_batch_p_val_PC1": "NA", "global_biovar_p_val_PC0": "NA", "global_biovar_p_val_PC1": "NA", }
    # fetch time spent running
    benchmarked_results_dir = output_dir_path + "/benchmarked_results/" + dataset_name
    benchmarked_results_dir = os.listdir(benchmarked_results_dir)
    current_runtime_kw_path = [result for result in benchmarked_results_dir if "runtime copy.txt" in result][0]
    print(current_runtime_kw_path)
    with open(output_dir_path + "/benchmarked_results/" + dataset_name+'/'+current_runtime_kw_path) as f:
        lines = f.readlines()

    for line in lines:
        if "combat" in line:
            method_dict["combat"]["runtime"] = float(line.split(" ")[-2])
        if "limma" in line:
            method_dict["limma"]["runtime"] = float(line.split(" ")[-2])
        if "MMUPHin" in line:
            method_dict["MMUPHin"]["runtime"] = float(line.split(" ")[-2])
        if "ConquR" in line:
            method_dict["ConQuR"]["runtime"] = float(line.split(" ")[-2])
        if "ConquR_libsize" in line:
            method_dict["ConQuR_libsize"]["runtime"] = float(line.split(" ")[-2])
        if "percentile_normalization" in line:
            method_dict["Percentile_norm"]["runtime"] = float(line.split(" ")[-2])
    
    benchmarked_data_harmony_dir = output_dir_path + "/benchmarked_data"
    benchmarked_data_harmony_dir = os.listdir(benchmarked_data_harmony_dir)
    current_runtime_kw_paths = [result for result in benchmarked_data_harmony_dir if "elapsed_time.txt" in result]

    for time_file in current_runtime_kw_paths:
        time_file = output_dir_path + "/benchmarked_data/"+time_file
        if "harmonicMic_elapsed_time" in time_file:
            with open(time_file) as f:
                lines = f.readlines()
            method_dict["harmonicMic"]["runtime"] = float(lines[0].split(" ")[-2])
        if "harmonicMic_weighted_elapsed_time" in time_file:
            with open(time_file) as f:
                lines = f.readlines()
            method_dict["harmonicMic_weighted"]["runtime"] = float(lines[0].split(" ")[-2])
        if "harmony_elapsed_time" in time_file:
            with open(time_file) as f:
                lines = f.readlines()
            method_dict["harmony"]["runtime"] = float(lines[0].split(" ")[-2])
        if "harmony_PCA_first_elapsed_time" in time_file:
            with open(time_file) as f:
                lines = f.readlines()
            method_dict["harmony_PCs"]["runtime"] = float(lines[0].split(" ")[-2])
        
    # print into a pandas df where rows are datasets and columns and different stats
    results_df = pd.DataFrame.from_dict(method_dict, orient ='index') 
    results_df.T.to_csv(output_dir_path+"/global_benchmarking_stats_"+dataset_name+".csv")
    return pd.DataFrame.from_dict(method_dict, orient ='index') 

def PC01_kw_tests_perbatch(df, bio_var, batch, output_root):
    # to investigate whether the bio_var info is retained
    bio_var_l = list(df[bio_var])
    bio_var_l = [x for x in bio_var_l if str(x) != 'nan']
    bio_var_l = list(np.unique(bio_var_l))
    
    # generate global metrics
    with open(output_root+"_"+bio_var+"_kw_pvals.txt", "w") as text_file:
        print("global batch kw p-vals \n", file=text_file)
        print("\n", file=text_file)

        data_PC0 = [df.loc[df[bio_var]==var]["PC0"].values for var in bio_var_l]
        data_PC1 = [df.loc[df[bio_var]==var]["PC1"].values for var in bio_var_l]
        
        global_PC0_p = stats.kruskal(*data_PC0)[1]
        global_PC1_p = stats.kruskal(*data_PC1)[1]
        print(". PC0", "across all biovar options, p-val = ", str(global_PC0_p), "\n",file=text_file)
        print(". PC1", "across all biovar options, p-val = ", str(global_PC1_p), "\n",file=text_file)
        
    dim = len(np.unique(bio_var_l))
    kw_heatmap_array = np.full((dim, dim), 1, dtype=float)
    for pair in itertools.combinations(bio_var_l, 2):
        data_for_kw_PC0 = [df.loc[df[bio_var]==var]["PC0"].values for var in pair]
        data_for_kw_PC1 = [df.loc[df[bio_var]==var]["PC1"].values for var in pair]
        PC0_p = stats.kruskal(*data_for_kw_PC0)[1]
        PC1_p = stats.kruskal(*data_for_kw_PC1)[1]
        print("PC0", pair[0], pair[1], PC0_p)
        print("PC1", pair[0], pair[1], PC1_p)
        print([bio_var_l.index(pair[0]), bio_var_l.index(pair[1])], kw_heatmap_array[bio_var_l.index(pair[0]), bio_var_l.index(pair[1])])
        kw_heatmap_array[bio_var_l.index(pair[0]), bio_var_l.index(pair[1])]=PC1_p
        kw_heatmap_array[bio_var_l.index(pair[1]), bio_var_l.index(pair[0])]=PC0_p
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1], transform=ax.transAxes, color="red")
    im = ax.imshow((-np.log2(kw_heatmap_array)).T, cmap="Oranges", rasterized=True,origin='lower')
    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(bio_var_l)), labels=bio_var_l)
    ax.set_xlabel("PC0 (lower diagonal)")
    ax.set_yticks(np.arange(len(bio_var_l)), labels=bio_var_l)
    ax.set_ylabel("PC1 (upper diagonal)")
    # Rotate the tidata["PC0"].valuesck labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
            rotation_mode="anchor")

    # Loop over data dimensionsif i!=j: and create text annotations.
    for i in range(len(bio_var_l)):
        for j in range(len(bio_var_l)):
            if i!=j:
                text = ax.text(i, j, round(-np.log2(kw_heatmap_array)[i, j], 5),
                            ha="center", va="center", color="black", fontsize=8)

    ax.set_title("Kruskal-Wallis -log2(p-val) Heatmap")
    fig.tight_layout()
    plt.savefig(output_root+"_"+batch+"_PC01_kw_tests.png")

def PCA_vis_for_each_batch(source_df, meta_data, output_root, batch_var, bio_var, n_pc=30):
    scaler = StandardScaler()
    data = source_df.T.to_numpy()
    scaler.fit(data)
    ss_scaled = scaler.transform(data).T
    source_df = pd.DataFrame(ss_scaled, columns = source_df.columns, index=source_df.index)
    
    list_batches = list(np.unique(meta_data[batch_var]))
    rng = np.random.RandomState(88)
    for batch in list_batches:
        # get current_source_df and current_meta_data
        current_meta_data = meta_data[meta_data[batch_var]==batch]
        current_Sam_ids = list(current_meta_data['Sam_id'])
        current_source_df = source_df[source_df.index.isin(current_Sam_ids)]

        # zero mean for PCA
        current_source_df = current_source_df-current_source_df.mean()
        # fit PCA
        pca_results  = PCA(n_components = n_pc, random_state=rng)
        pca_results.fit(current_source_df)
        variances = pca_results.explained_variance_ratio_

        # ecdf plot
        ecdf_var = [variances[:i].sum() for i in range(len(variances))]
        fig, ax1 = plt.subplots(1,1,figsize = (20, 8))
        ax1.plot(range(len(ecdf_var)),ecdf_var)
        plt.xticks(range(len(ecdf_var)), range(len(ecdf_var)))
        plt.yticks([0.1*i for i in range(11)])
        plt.show()
        plt.savefig(output_root+"_ecdf.png")

        # pca visualization
        transformed = pca_results.transform(current_source_df)
        df = pd.DataFrame()
        for i in range(n_pc):
            df[f"PC{i}"] = transformed[:, i] 
        df[bio_var] = list(current_meta_data[bio_var])
        if len(np.unique(df[bio_var])) == 1:
            print(batch, "only has one bio_var option:", np.unique(df[bio_var])[0])
        else:
            ## pca visualization for bio_vars
            ### fetch info for bio_va
            all_colors_used = rng.uniform(0, 1, 3*len(np.unique(list(df[bio_var])))).tolist()
            colors = {var: all_colors_used[idx*3:idx*3+3] for idx,var in enumerate(np.unique(list(df[bio_var])))} # unlikely to get repeated colors
            fig, ax =  plt.subplots(1, 1, sharex=True, sharey=True, figsize=(20, 10))
            sns.scatterplot(df["PC0"] , df["PC1"], hue = df[bio_var], hue_order = np.unique(list(df[bio_var])), s = 100,ax = ax, cmap = "tab10", x_jitter = True, palette = colors)
            ax.legend(bbox_to_anchor=(1.04,1), loc="upper left")
            for index, current_biovar in enumerate(np.unique(list(df[bio_var]))):
                currentDF = df.loc[df[bio_var] == current_biovar]
                confidence_ellipse(currentDF['PC0'], currentDF['PC1'], ax, n_std=1.5, edgecolor=list(colors.values())[index])
            plt.savefig(output_root+"_PCA_"+bio_var+"_"+batch+".png")

            PC01_kw_tests_perbatch(df, bio_var, batch, output_root)


# Glickman dataset 
#################################################################################
# IDCol = 'Sam_id'
# index_col = "Unnamed: 0"
# vars_use = ["Dataset", "Sex"]
# output_root = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/Glickman"

# address_X = "/home/fuc/HRZE_TB/tom_organized_codes/batch_correction_PCA/1021_microbiome_batchcorrection/microbiome_merged_intersect_1023.csv"
# address_Y = "/home/fuc/HRZE_TB/tom_organized_codes/batch_correction_PCA/1021_microbiome_batchcorrection/intersect_metadata_1023.csv"
# data_mat, meta_data = load_data(address_X, address_Y, IDCol, index_col, output_root)
# # # PCA_vis_for_each_batch(data_mat, meta_data, output_root,"batch",  "Visit", n_pc=30)

# res, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmonicMic_theta0.5", option = "harmonicMic")
# Evaluate(res, meta_data, 'Dataset', './output_Glickman_harmonicMic_trial2/Glickman_harmonicMic_theta_trial2_1201', "Visit", 30, 'Sex', 'Sam_id')
# Evaluate(res, meta_data, 'Dataset', './output_Glickman_harmonicMic_theta0.5/Glickman_harmonicMic_theta0.5_1201', "Visit", 30, 'Sex', 'Sam_id')
# res, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmonicMic_weighted", option = "harmonicMic", diversity_weight=0.3)
# Evaluate(res, meta_data, 'Dataset', './output_Glickman_harmonicMic_weighted/Glickman_harmonicMic_weighted_1201', "Visit", 30, 'Sex', 'Sam_id')
# res_h, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmony", option = "harmony")
# Evaluate(res_h, meta_data, 'Dataset', './output_Glickman_harmony/Glickman_harmony_1201', "Visit", 30, 'Sex', 'Sam_id')
# res_h, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmony_PCs", option = "harmony", PCA_first=True)
# Evaluate(res_h, meta_data, 'Dataset', './output_Glickman_harmony_PCs/Glickman_harmony_PCs', "Visit", 30, 'Sex', 'Sam_id')
# Evaluate(data_mat, meta_data, 'Dataset', './output_Glickman_nobc/Glickman_nobc_1201', "Visit", 30, 'Sex', 'Sam_id')

# ## benchmarking other methods:
# address_Y = "/home/fuc/HRZE_TB/tom_organized_codes/batch_correction_PCA/1021_microbiome_batchcorrection/intersect_metadata_1023.csv"

# ### combat
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/Glickman/Glickman_combat.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, 'Dataset', './output_Glickman_combat/Glickman_combat_1201', "Visit", 30, 'Sex', 'Sam_id')

# ### ConQuR
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/Glickman/Glickman_ConQuR.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, 'Dataset', './output_Glickman_ConQuR/Glickman_ConQuR_1201', "Visit", 30, 'Sex', 'Sam_id')

# ### ConQuR_libsize
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/Glickman/Glickman_ConQuR_libsize.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, 'Dataset', './output_Glickman_ConQuR_libsize/Glickman_ConQuR_libsize_1201', "Visit", 30, 'Sex', 'Sam_id')

# ### limma
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/Glickman/Glickman_limma.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, 'Dataset', './output_Glickman_limma/Glickman_limma_1201', "Visit", 30, 'Sex', 'Sam_id')

# ### MMUPHin
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/Glickman/Glickman_MMUPHin.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, 'Dataset', './output_Glickman_MMUPHin/Glickman_MMUPHin_1201', "Visit", 30, 'Sex', 'Sam_id')

# ### Percentile_normalization
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/Glickman/Glickman_percentile_norm.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, 'Dataset', './output_Glickman_Percentile_norm/Glickman_Percentile_norm_1201', "Visit", 30, 'Sex', 'Sam_id')

# input_frame_path = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/Glickman_count_data.csv"
# bio_var = "Visit"
# dataset_name = "Glickman"
# methods_list = ["combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "Percentile_norm", "harmony", "harmonicMic", "harmony_PCs", "harmonicMic_weighted"]
# global_eval_dataframe(input_frame_path, bio_var, dataset_name, methods_list, output_dir_path = ".")

# #################################################################################



# # autism 2 microbiomeHD
# #################################################################################
# address_directory = '/home/fuc/harmonicMic/data/autism_2_microbiomeHD'
# output_root = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/autism_2_microbiomeHD"
# data_mat, meta_data = load_data_microbiomeHD(address_directory, output_root)
# vars_use = ["Dataset"]
# IDCol = 'Sam_id'
# res, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmonicMic", option = "harmonicMic")
# res_h, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmony", option = "harmony")
# Evaluate(res, meta_data, 'Dataset', './output_autism_2_microbiomeHD_harmonicMic/autism_2_microbiomeHD_harmonicMic_1201', "DiseaseState", 30, False, 'Sam_id')
# Evaluate(res_h, meta_data, 'Dataset', './output_autism_2_microbiomeHD_harmony/autism_2_microbiomeHD_harmony_1201', "DiseaseState", 30, False, 'Sam_id')
# res_h, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmony_PCs", option = "harmony", PCA_first=True)
# Evaluate(res_h, meta_data, 'Dataset', './output_autism_2_microbiomeHD_harmony_PCs/autism_2_microbiomeHD_harmony_PCs', "DiseaseState", 30, False, 'Sam_id')
# Evaluate(data_mat, meta_data, 'Dataset', './output_autism_2_microbiomeHD_nobc/autism_2_microbiomeHD_nobc_1201', "DiseaseState", 30, False, 'Sam_id')
# res, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmonicMic_weighted", option = "harmonicMic", diversity_weight=0.3)
# Evaluate(res, meta_data, 'Dataset', './output_autism_2_microbiomeHD_harmonicMic_weighted/autism_2_microbiomeHD_harmonicMic_weighted_1201', "DiseaseState", 30, False, 'Sam_id')

# # benchmarking other methods:DONE
# address_Y = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/autism_2_microbiomeHD_meta_data.csv"

# ### combat
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD_combat.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, 'Dataset', './output_autism_2_microbiomeHD_combat/autism_2_microbiomeHD_combat_1201', "DiseaseState", 30,  False, 'Sam_id')

# ### ConQuR
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD_ConQuR.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, 'Dataset', './output_autism_2_microbiomeHD_ConQuR/autism_2_microbiomeHD_ConQuR_1201', "DiseaseState", 30, False, 'Sam_id')

# ### ConQuR_libsize
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD_ConQuR_libsize.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, 'Dataset', './output_autism_2_microbiomeHD_ConQuR_libsize/autism_2_microbiomeHD_ConQuR_libsize_1201', "DiseaseState", 30,  False, 'Sam_id')

# ### limma
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD_limma.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, 'Dataset', './output_autism_2_microbiomeHD_limma/autism_2_microbiomeHD_limma_1201', "DiseaseState", 30,  False,'Sam_id')

# ### MMUPHin
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD_MMUPHin.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, 'Dataset', './output_autism_2_microbiomeHD_MMUPHin/autism_2_microbiomeHD_MMUPHin_1201', "DiseaseState" , 30, False,'Sam_id')

# ### Percentile_normalization
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD_percentile_norm.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, 'Dataset', './output_autism_2_microbiomeHD_Percentile_norm/autism_2_microbiomeHD_Percentile_norm_1201', "DiseaseState", 30, False, 'Sam_id')

# input_frame_path = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/autism_2_microbiomeHD_count_data.csv"
# bio_var = "DiseaseState"
# dataset_name = "autism_2_microbiomeHD"
# methods_list = ["combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "Percentile_norm", "harmony", "harmonicMic", "harmony_PCs", "harmonicMic_weighted"]
# global_eval_dataframe(input_frame_path, bio_var, dataset_name, methods_list, output_dir_path = ".")

# #################################################################################

# # cdi 3 microbiomeHD
# #################################################################################
# address_directory = '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/data/cdi_3_microbiomeHD'
# output_root = "/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/cdi_3_microbiomeHD"
# data_mat, meta_data = load_data_microbiomeHD(address_directory, output_root)
# PCA_vis_for_each_batch(data_mat, meta_data, output_root, "Dataset", "DiseaseState", n_pc=10)
# vars_use = ["Dataset"]
# IDCol = 'Sam_id'
# res, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmonicMic", option = "harmonicMic")
# res_h, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmony", option = "harmony")
# Evaluate(res, meta_data, 'Dataset', './output_cdi_3_microbiomeHD_harmonicMic/cdi_3_microbiomeHD_harmonicMic_1201', "DiseaseState", 30, False, 'Sam_id')
# Evaluate(res_h, meta_data, 'Dataset', './output_cdi_3_microbiomeHD_harmony/cdi_3_microbiomeHD_harmony_1201', "DiseaseState", 30, False, 'Sam_id')
# res_h, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmony_PCs", option = "harmony", PCA_first=True)
# Evaluate(res_h, meta_data, 'Dataset', './output_cdi_3_microbiomeHD_harmony_PCs/cdi_3_microbiomeHD_harmony_PCs', "DiseaseState", 30, False, 'Sam_id')
# Evaluate(data_mat, meta_data, 'Dataset', './output_cdi_3_microbiomeHD_nobc/cdi_3_microbiomeHD_nobc_1201', "DiseaseState", 30, False, 'Sam_id')
# res, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmonicMic_weighted", option = "harmonicMic", diversity_weight=0.3)
# Evaluate(res, meta_data, 'Dataset', './output_cdi_3_microbiomeHD_harmonicMic_weighted/cdi_3_microbiomeHD_harmonicMic_weighted_1201', "DiseaseState", 30, False, 'Sam_id')


# # benchmarking other methods: 
# address_Y = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/cdi_3_microbiomeHD_meta_data.csv"

# ### combat
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD_combat.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, 'Dataset', './output_cdi_3_microbiomeHD_combat/cdi_3_microbiomeHD_combat_1201', "DiseaseState", 30,  False, 'Sam_id')

# ### ConQuR
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD_ConQuR.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, 'Dataset', './output_cdi_3_microbiomeHD_ConQuR/cdi_3_microbiomeHD_ConQuR_1201', "DiseaseState", 30, False, 'Sam_id')

# ### ConQuR_libsize
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD_ConQuR_libsize.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, 'Dataset', './output_cdi_3_microbiomeHD_ConQuR_libsize/cdi_3_microbiomeHD_ConQuR_libsize_1201', "DiseaseState", 30,  False, 'Sam_id')

# ### limma
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD_limma.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, 'Dataset', './output_cdi_3_microbiomeHD_limma/cdi_3_microbiomeHD_limma_1201', "DiseaseState", 30,  False,'Sam_id')

# ### MMUPHin
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD_MMUPHin.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, 'Dataset', './output_cdi_3_microbiomeHD_MMUPHin/cdi_3_microbiomeHD_MMUPHin_1201', "DiseaseState" , 30, False,'Sam_id')

# ### Percentile_normalization
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD_percentile_norm.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, 'Dataset', './output_cdi_3_microbiomeHD_Percentile_norm/cdi_3_microbiomeHD_Percentile_norm_1201', "DiseaseState", 30, False, 'Sam_id')

# input_frame_path = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/cdi_3_microbiomeHD_count_data.csv"
# bio_var = "DiseaseState"
# dataset_name = "cdi_3_microbiomeHD"
# methods_list = ["combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "Percentile_norm", "harmony", "harmonicMic", "harmony_PCs", "harmonicMic_weighted"]
# global_eval_dataframe(input_frame_path, bio_var, dataset_name, methods_list, output_dir_path = ".")

# #################################################################################


# # ibd 3 CMD
# #################################################################################
# address_directory = '/home/fuc/harmonicMic/data/ibd_3_CMD'
# output_root = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/ibd_3_CMD"
# data_mat, meta_data = load_data_CMD(address_directory, output_root)
# vars_use = ["study_name"]
# IDCol = 'Sam_id'
# res, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmonicMic", option = "harmonicMic")
# res_h, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmony", option = "harmony")
# Evaluate(res, meta_data, "study_name", './output_ibd_3_CMD_harmonicMic/ibd_3_CMD_harmonicMic_1201', "disease", 30, False, 'Sam_id')
# Evaluate(res_h, meta_data, "study_name", './output_ibd_3_CMD_harmony/ibd_3_CMD_harmony_1201', "disease", 30, False, 'Sam_id')
# res_h, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmony_PCs", option = "harmony", PCA_first=True)
# Evaluate(res_h, meta_data, "study_name", './output_ibd_3_CMD_harmony_PCs/ibd_3_CMD_harmony_PCs', "disease", 30, False, 'Sam_id')
# Evaluate(data_mat, meta_data, "study_name", './output_ibd_3_CMD_nobc/ibd_3_CMD_nobc_1201', "disease", 30, False, 'Sam_id')
# res, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmonicMic_weighted", option = "harmonicMic", diversity_weight=0.3)
# Evaluate(res, meta_data, "study_name", './output_ibd_3_CMD_harmonicMic_weighted/ibd_3_CMD_harmonicMic_weighted_1201', "disease", 30, False, 'Sam_id')

# # benchmarking other methods: # TODO:
# address_Y = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/ibd_3_CMD_meta_data.csv"

# ### combat -> PROBLEMATIC
# # address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/ibd_3_CMD/ibd_3_CMD_combat.csv"
# # data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# # Evaluate(data_mat, meta_data, "study_name", './output_ibd_3_CMD_combat/ibd_3_CMD_combat_1201', "disease", 30,  False, 'Sam_id')

# ### ConQuR
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/ibd_3_CMD/ibd_3_CMD_ConQuR.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, "study_name", './output_ibd_3_CMD_ConQuR/ibd_3_CMD_ConQuR_1201', "disease", 30, False, 'Sam_id')

# ### ConQuR_libsize
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/ibd_3_CMD/ibd_3_CMD_ConQuR_libsize.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, "study_name", './output_ibd_3_CMD_ConQuR_libsize/ibd_3_CMD_ConQuR_libsize_1201', "disease", 30,  False, 'Sam_id')

# ### limma
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/ibd_3_CMD/ibd_3_CMD_limma.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, "study_name", './output_ibd_3_CMD_limma/ibd_3_CMD_limma_1201', "disease", 30,  False,'Sam_id')

# ### MMUPHin
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/ibd_3_CMD/ibd_3_CMD_MMUPHin.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, "study_name", './output_ibd_3_CMD_MMUPHin/ibd_3_CMD_MMUPHin_1201', "disease" , 30, False,'Sam_id')

# ### Percentile_normalization
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/ibd_3_CMD/ibd_3_CMD_percentile_norm.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, "study_name", './output_ibd_3_CMD_Percentile_norm/ibd_3_CMD_Percentile_norm_1201', "disease", 30, False, 'Sam_id')

# input_frame_path = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/ibd_3_CMD_count_data.csv"
# bio_var = "disease"
# dataset_name = "ibd_3_CMD"
# methods_list = ["limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "Percentile_norm", "harmony", "harmonicMic", "harmony_PCs", "harmonicMic_weighted"]
# global_eval_dataframe(input_frame_path, bio_var, dataset_name, methods_list, output_dir_path = ".")

# #################################################################################


# # adenoma 5 CMD
# #################################################################################
# address_directory = '/home/fuc/harmonicMic/data/adenoma_5_CMD'
# output_root = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/adenoma_5_CMD"
# data_mat, meta_data = load_data_CMD(address_directory, output_root)
# vars_use = ["study_name"]
# IDCol = 'Sam_id'

# res, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmonicMic", option = "harmonicMic")
# res_h, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmony", option = "harmony")
# Evaluate(res, meta_data, "study_name", './output_adenoma_5_CMD_harmonicMic/adenoma_5_CMD_harmonicMic_1201', "disease", 30, False, 'Sam_id')
# Evaluate(res_h, meta_data, "study_name", './output_adenoma_5_CMD_harmony/adenoma_5_CMD_harmony_1201', "disease", 30, False, 'Sam_id')
# res_h, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmony_PCs", option = "harmony", PCA_first=True)
# Evaluate(res_h, meta_data, "study_name", './output_adenoma_5_CMD_harmony_PCs/adenoma_5_CMD_harmony_PCs', "disease", 30, False, 'Sam_id')
# Evaluate(data_mat, meta_data, "study_name", './output_adenoma_5_CMD_nobc/adenoma_5_CMD_nobc_1201', "disease", 30, False, 'Sam_id')
# res, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmonicMic_weighted", option = "harmonicMic", diversity_weight=0.3)
# Evaluate(res, meta_data, "study_name", './output_adenoma_5_CMD_harmonicMic_weighted/adenoma_5_CMD_harmonicMic_weighted_1201', "disease", 30, False, 'Sam_id')

# # benchmarking other methods: # TODO:
# address_Y = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/adenoma_5_CMD_meta_data.csv"

# ### combat 
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/adenoma_5_CMD/adenoma_5_CMD_combat.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, "study_name", './output_adenoma_5_CMD_combat/adenoma_5_CMD_combat_1201', "disease", 30,  False, 'Sam_id')

# ### ConQuR
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/adenoma_5_CMD/adenoma_5_CMD_ConQuR.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, "study_name", './output_adenoma_5_CMD_ConQuR/adenoma_5_CMD_ConQuR_1201', "disease", 30, False, 'Sam_id')

# ### ConQuR_libsize
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/adenoma_5_CMD/adenoma_5_CMD_ConQuR_libsize.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, "study_name", './output_adenoma_5_CMD_ConQuR_libsize/adenoma_5_CMD_ConQuR_libsize_1201', "disease", 30,  False, 'Sam_id')

# ### limma
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/adenoma_5_CMD/adenoma_5_CMD_limma.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, "study_name", './output_adenoma_5_CMD_limma/adenoma_5_CMD_limma_1201', "disease", 30,  False,'Sam_id')

# ### MMUPHin
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/adenoma_5_CMD/adenoma_5_CMD_MMUPHin.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, "study_name", './output_adenoma_5_CMD_MMUPHin/adenoma_5_CMD_MMUPHin_1201', "disease" , 30, False,'Sam_id')

# ### Percentile_normalization
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/adenoma_5_CMD/adenoma_5_CMD_percentile_norm.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, "study_name", './output_adenoma_5_CMD_Percentile_norm/adenoma_5_CMD_Percentile_norm_1201', "disease", 30, False, 'Sam_id')

# input_frame_path = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/adenoma_5_CMD_count_data.csv"
# bio_var = "disease"
# dataset_name = "adenoma_5_CMD"
# methods_list = ["combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "Percentile_norm", "harmony", "harmonicMic", "harmony_PCs", "harmonicMic_weighted"]
# global_eval_dataframe(input_frame_path, bio_var, dataset_name, methods_list, output_dir_path = ".")

#################################################################################

# # CRC_8_CMD
# #################################################################################
# address_directory = '/home/fuc/harmonicMic/data/CRC_8_CMD'
# output_root = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/CRC_8_CMD"
# data_mat, meta_data = load_data_CMD(address_directory, output_root)
# vars_use = ["study_name"]
# IDCol = 'Sam_id'

# res, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmonicMic", option = "harmonicMic")
# res_h, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmony", option = "harmony")
# Evaluate(res, meta_data, "study_name", './output_CRC_8_CMD_harmonicMic/CRC_8_CMD_harmonicMic_1201', "disease", 30, False, 'Sam_id')
# Evaluate(res_h, meta_data, "study_name", './output_CRC_8_CMD_harmony/CRC_8_CMD_harmony_1201', "disease", 30, False, 'Sam_id')
# res_h, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmony_PCs", option = "harmony", PCA_first=True)
# Evaluate(res_h, meta_data, "study_name", './output_CRC_8_CMD_harmony_PCs/CRC_8_CMD_harmony_PCs', "disease", 30, False, 'Sam_id')
# Evaluate(data_mat, meta_data, "study_name", './output_CRC_8_CMD_nobc/CRC_8_CMD_nobc_1201', "disease", 30, False, 'Sam_id')
# res, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmonicMic_weighted", option = "harmonicMic", diversity_weight=0.3)
# Evaluate(res, meta_data, "study_name", './output_CRC_8_CMD_harmonicMic_weighted/CRC_8_CMD_harmonicMic_weighted_1201', "disease", 30, False, 'Sam_id')

# # benchmarking other methods: 
# address_Y = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/CRC_8_CMD_meta_data.csv"

# ### combat # TODO: to fix: contains NAN
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/CRC_8_CMD/CRC_8_CMD_combat.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, "study_name", './output_CRC_8_CMD_combat/CRC_8_CMD_combat_1201', "disease", 30,  False, 'Sam_id')

# ### ConQuR
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/CRC_8_CMD/CRC_8_CMD_ConQuR.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, "study_name", './output_CRC_8_CMD_ConQuR/CRC_8_CMD_ConQuR_1201', "disease", 30, False, 'Sam_id')

# ### ConQuR_libsize
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/CRC_8_CMD/CRC_8_CMD_ConQuR_libsize.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, "study_name", './output_CRC_8_CMD_ConQuR_libsize/CRC_8_CMD_ConQuR_libsize_1201', "disease", 30,  False, 'Sam_id')

# ### limma
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/CRC_8_CMD/CRC_8_CMD_limma.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, "study_name", './output_CRC_8_CMD_limma/CRC_8_CMD_limma_1201', "disease", 30,  False,'Sam_id')

# ### MMUPHin
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/CRC_8_CMD/CRC_8_CMD_MMUPHin.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, "study_name", './output_CRC_8_CMD_MMUPHin/CRC_8_CMD_MMUPHin_1201', "disease" , 30, False,'Sam_id')

# ### Percentile_normalization
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/CRC_8_CMD/CRC_8_CMD_percentile_norm.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, "study_name", './output_CRC_8_CMD_Percentile_norm/CRC_8_CMD_Percentile_norm_1201', "disease", 30, False, 'Sam_id')

# input_frame_path = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/CRC_8_CMD_count_data.csv"
# bio_var = "disease"
# dataset_name = "CRC_8_CMD"
# methods_list = ["combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "Percentile_norm", "harmony", "harmonicMic", "harmony_PCs", "harmonicMic_weighted"]
# global_eval_dataframe(input_frame_path, bio_var, dataset_name, methods_list, output_dir_path = ".")


#################################################################################

# # T2D 10 CMD 
# #################################################################################
# address_directory = '/home/fuc/harmonicMic/data/T2D_10_CMD'
# output_root = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/T2D_10_CMD"
# data_mat, meta_data = load_data_CMD(address_directory, output_root)
# PCA_vis_for_each_batch(data_mat, meta_data, output_root, "study_name", "disease", n_pc=10)
# vars_use = ["study_name"]
# IDCol = 'Sam_id'


# res, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmonicMic", option = "harmonicMic", diversity_weight=0.1)
# res_h, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmony", option = "harmony")
# Evaluate(res, meta_data, "study_name", './output_T2D_10_CMD_harmonicMic/T2D_10_CMD_harmonicMic_1201', "disease", 30, False, 'Sam_id')
# Evaluate(res, meta_data, "study_name", './output_T2D_10_CMD_harmonicMic_trial/T2D_10_CMD_harmonicMic_trial_1201', "disease", 30, False, 'Sam_id')
# Evaluate(res, meta_data, "study_name", './output_T2D_10_CMD_harmonicMic_trial_2/T2D_10_CMD_harmonicMic_trial_2_1201', "disease", 30, False, 'Sam_id')
# Evaluate(res_h, meta_data, "study_name", './output_T2D_10_CMD_harmony/T2D_10_CMD_harmony_1201', "disease", 30, False, 'Sam_id')
# res_h, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmony_PCs", option = "harmony", PCA_first=True)
# Evaluate(res_h, meta_data, "study_name", './output_T2D_10_CMD_harmony_PCs/T2D_10_CMD_harmony_PCs', "disease", 30, False, 'Sam_id')
# Evaluate(data_mat, meta_data, "study_name", './output_T2D_10_CMD_nobc/T2D_10_CMD_nobc_1201', "disease", 30, False, 'Sam_id')
# res, meta_data = generate_harmonicMic_results(data_mat, meta_data, IDCol, vars_use, output_root+"harmonic_weighted", option = "harmonicMic", diversity_weight=0.3)
# Evaluate(res, meta_data, "study_name", './output_T2D_10_CMD_harmonicMic_weighted/T2D_10_CMD_harmonicMic_weighted_1201', "disease", 30, False, 'Sam_id')

# # benchmarking other methods: # TODO:
# address_Y = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/T2D_10_CMD_meta_data.csv"

# ### combat 
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/T2D_10_CMD/T2D_10_CMD_combat.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, "study_name", './output_T2D_10_CMD_combat/T2D_10_CMD_combat_1201', "disease", 30,  False, 'Sam_id')

# ### ConQuR
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/T2D_10_CMD/T2D_10_CMD_ConQuR.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, "study_name", './output_T2D_10_CMD_ConQuR/T2D_10_CMD_ConQuR_1201', "disease", 30, False, 'Sam_id')

# ### ConQuR_libsize
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/T2D_10_CMD/T2D_10_CMD_ConQuR_libsize.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, "study_name", './output_T2D_10_CMD_ConQuR_libsize/T2D_10_CMD_ConQuR_libsize_1201', "disease", 30,  False, 'Sam_id')

# ### limma
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/T2D_10_CMD/T2D_10_CMD_limma.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, "study_name", './output_T2D_10_CMD_limma/T2D_10_CMD_limma_1201', "disease", 30,  False,'Sam_id')

# ### MMUPHin
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/T2D_10_CMD/T2D_10_CMD_MMUPHin.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, "study_name", './output_T2D_10_CMD_MMUPHin/T2D_10_CMD_MMUPHin_1201', "disease" , 30, False,'Sam_id')

# ### Percentile_normalization
# address_X = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/T2D_10_CMD/T2D_10_CMD_percentile_norm.csv"
# data_mat, meta_data = load_results_from_benchmarked_methods(address_X, address_Y)
# Evaluate(data_mat, meta_data, "study_name", './output_T2D_10_CMD_Percentile_norm/T2D_10_CMD_Percentile_norm_1201', "disease", 30, False, 'Sam_id')

# input_frame_path = "/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/T2D_10_CMD_count_data.csv"
# bio_var = "disease"
# dataset_name = "T2D_10_CMD"
# methods_list = ["combat", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize", "Percentile_norm", "harmony", "harmonicMic", "harmony_PCs", "harmonicMic_weighted"]
# global_eval_dataframe(input_frame_path, bio_var, dataset_name, methods_list, output_dir_path = ".")
# #################################################################################


address_directory = '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/data/hanninganGD'
output_root = "/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/hanninganGD"
data_mat, meta_data = load_data_CMD(address_directory, output_root, id = 'patient_visit_id')
Evaluate(data_mat, meta_data, 'location', './output_hanninganGD_nobc/hanninganGD_nobc_050223', "disease", 30,  False, 'patient_visit_id')
