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
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
from sklearn.decomposition import PCA
import seaborn as sns
import random
from sklearn.preprocessing import StandardScaler
import itertools
from scipy import stats

def generate_harmonicMic_results(address_X, address_Y, IDCol, vars_use, index_col = False):
    data_mat, meta_data = load_data(address_X, address_Y, IDCol, index_col)
    print("aa")
    print(data_mat.shape)
    print(meta_data.shape)
    ho = run_harmonicMic(data_mat, meta_data, vars_use)
    res = pd.DataFrame(ho.Z_corr)

    # give the index back
    res.index = data_mat.columns
    return res, meta_data

def generate_harmony_results(address_X, address_Y, IDCol, vars_use, index_col = False):
    data_mat, meta_data = load_data(address_X, address_Y, IDCol, index_col)
    print("aa")
    print(data_mat.shape)
    print(meta_data.shape)
    ho = run_harmony(data_mat, meta_data, vars_use)
    res = pd.DataFrame(ho.Z_corr)

    # give the index back
    res.index = data_mat.columns
    return res, meta_data
# address_X = "/home/fuc/HRZE_TB/tom_organized_codes/batch_correction_PCA/1021_microbiome_batchcorrection/microbiome_merged_intersect_1023.csv"
# address_Y = "/home/fuc/HRZE_TB/tom_organized_codes/batch_correction_PCA/1021_microbiome_batchcorrection/intersect_metadata_1023.csv"
# IDCol = 'Sam_id'
# index_col = "Unnamed: 0"
# vars_use = ["Dataset", "Sex"]
# res, meta_data = generate_harmonicMic_results(address_X, address_Y, IDCol, vars_use, index_col)

# res_h, meta_data = generate_harmony_results(address_X, address_Y, IDCol, vars_use, index_col)
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


def run_eval(batch_corrected_df, meta_data, batch_var, output_root, bio_var = False, n_pc=30, covar = False):
    a = Evaluate(batch_corrected_df, meta_data, batch_var, output_root, bio_var, n_pc, covar)
    return

# Evaluate(res, meta_data, 'Dataset', 'Glickman_harmonicMic', "Visit", 30, 'Sex')
# Evaluate(res_h, meta_data, 'Dataset', 'Glickman_harmony', "Visit", 30, 'Sex')
class Evaluate(object):
    def __init__(
            self, batch_corrected_df, meta_data, batch_var, output_root, bio_var = False, n_pc=30, covar = False
    ):
        self.batch_corrected_df = batch_corrected_df
        self.meta_data = meta_data
        self.batch_var = batch_var
        self.output_root = output_root
        self.bio_var = bio_var # var retaining the biological info
        print("bio_var", bio_var)
        self.n_pc = n_pc
        self.covar = covar
        self.rng = np.random.RandomState(88)
        # functions executed
        self.standard_scaler()
        df = self.PCA_vis()
        self.PC01_kw_tests(df)
        self.covar == covar
        if covar != False:
            self.covar == covar
            self.allPCs_covar_kw_test(df)

    def standard_scaler(self):
        source_df = self.batch_corrected_df
        scaler = StandardScaler()
        data = source_df.T.to_numpy()
        scaler.fit(data)
        ss_scaled = scaler.transform(data).T
        self.batch_corrected_df = pd.DataFrame(ss_scaled, columns = source_df.columns, index=source_df.index)
    
    def PCA_vis(self):
        source_df = self.batch_corrected_df

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
        df['batches'] = self.meta_data[self.batch_var]
        df[self.bio_var] = self.meta_data[self.bio_var]

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
            df[self.covar] = self.meta_data[self.covar]
        return df

    def PC01_kw_tests(self, df):
        # TO DEBUG
        bio_var_l = list(df[self.bio_var])
        bio_var_l = [x for x in bio_var_l if str(x) != 'nan']
        bio_var_l = list(np.unique(bio_var_l))
        dim = len(np.unique(bio_var_l))
        kw_heatmap_array = np.full((dim, dim), 1, dtype=float)
        for pair in itertools.combinations(bio_var_l, 2):
            data_for_kw_PC0 = [df.loc[df[self.bio_var]==var]["PC0"].values for var in pair]
            data_for_kw_PC1 = [df.loc[df[self.bio_var]==var]["PC1"].values for var in pair]
            # data1 = df.loc[df['Day'] == pair[0]]
            # data2 = df.loc[df['Day'] == pair[1]]
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
        plt.savefig(self.output_root+"_OC01_kw_tests.png")

    def allPCs_covar_kw_test(self, df):
        # to debug
        # lets do 30PCs for now
        fig, axes = plt.subplots(5, 6, sharex=True, sharey=True, figsize=(35, 35))
        current_ax_index = [0, 0]
        covar_unique = list(df[self.covar].values)
        covar_unique = [x for x in covar_unique if str(x) != 'nan']
        covar_unique = np.unique(covar_unique).tolist()
        for i in range(30):
            print(current_ax_index)
            a = sns.scatterplot(df["PC0"],df["PC"+str(i)], hue = df[self.covar], style = df["batches"], s = 100,ax = axes[current_ax_index[0], current_ax_index[1]], cmap = "tab10", x_jitter = True)
            a.set_title("PC"+str(i))
            data_for_kw = [df.loc[df[self.covar]==var]["PC"+str(i)].values for var in covar_unique]
            # data_m = df.loc[df[self.covar]=="Male"]
            # data_f = df.loc[df['Sex']=="Female"]
            PC1_p = stats.kruskal(*data_for_kw)[1]
            if current_ax_index[1] == 5:
                current_ax_index = [current_ax_index[0]+1, 0]
            else:
                current_ax_index[1] +=1
            print("PC"+str(i), PC1_p)
        
        plt.savefig(self.output_root+"_allPCs_covar_kw_tests.png")
