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
# kmeans does not always return k centroids, but kmeans2 does
from scipy.cluster.vq import kmeans2
import logging
import itertools
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa
from scipy import stats

# create logger
logger = logging.getLogger('harmonypy')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

# from IPython.core.debugger import set_trace

def run_harmonicMic(
    data_mat: np.ndarray,
    meta_data: pd.DataFrame,
    vars_use,
    theta = None,
    lamb = None,
    sigma = 0.1, 
    nclust = None,
    tau = 0,
    block_size = 0.05, 
    max_iter_harmony = 10,
    max_iter_kmeans = 20,
    epsilon_cluster = 1e-5,
    epsilon_harmony = 1e-4, 
    plot_convergence = False,
    verbose = True,
    reference_values = None,
    cluster_prior = None,
    random_state = 0
):
    """Run Harmony.
    """

    # theta = None
    # lamb = None
    # sigma = 0.1
    # nclust = None
    # tau = 0
    # block_size = 0.05
    # epsilon_cluster = 1e-5
    # epsilon_harmony = 1e-4
    # plot_convergence = False
    # verbose = True
    # reference_values = None
    # cluster_prior = None
    # random_state = 0

    N = meta_data.shape[0]
    if data_mat.shape[1] != N:
        data_mat = data_mat.T

    assert data_mat.shape[1] == N, \
       "data_mat and meta_data do not have the same number of cells" 

    if nclust is None:
        nclust = np.min([np.round(N / 30.0), 100]).astype(int)

    if type(sigma) is float and nclust > 1:
        sigma = np.repeat(sigma, nclust)

    if isinstance(vars_use, str):
        vars_use = [vars_use]

    phi = pd.get_dummies(meta_data[vars_use]).to_numpy().T
    phi_n = meta_data[vars_use].describe().loc['unique'].to_numpy().astype(int)

    ### NOTE: ADDED phi_dict
    phi_dict = {}
    for var in vars_use:
        print(var, np.unique(list(meta_data[var])))
        print(pd.get_dummies(meta_data[var]).to_numpy().T)
        phi_dict.update({var: pd.get_dummies(meta_data[var]).to_numpy().T})
    # print("phi_n", phi_n)
    if theta is None:
        theta = np.repeat([1] * len(phi_n), phi_n)
    elif isinstance(theta, float) or isinstance(theta, int):
        theta = np.repeat([theta] * len(phi_n), phi_n)
    elif len(theta) == len(phi_n):
        theta = np.repeat([theta], phi_n)

    assert len(theta) == np.sum(phi_n), \
        "each batch variable must have a theta"

    if lamb is None:
        lamb = np.repeat([1] * len(phi_n), phi_n)
    elif isinstance(lamb, float) or isinstance(lamb, int):
        lamb = np.repeat([lamb] * len(phi_n), phi_n)
    elif len(lamb) == len(phi_n):
        lamb = np.repeat([lamb], phi_n)

    assert len(lamb) == np.sum(phi_n), \
        "each batch variable must have a lambda"

    # Number of items in each category.
    N_b = phi.sum(axis = 1)
    # Proportion of items in each category.
    Pr_b = N_b / N

    if tau > 0:
        theta = theta * (1 - np.exp(-(N_b / (nclust * tau)) ** 2))

    lamb_mat = np.diag(np.insert(lamb, 0, 0))

    phi_moe = np.vstack((np.repeat(1, N), phi))

    np.random.seed(random_state)
    print(data_mat.shape)
    ho = HarmonicMic(
        data_mat, phi, phi_moe, Pr_b, sigma, theta, max_iter_harmony, max_iter_kmeans,
        epsilon_cluster, epsilon_harmony, nclust, block_size, lamb_mat, verbose, vars_use, phi_dict
    )

    return ho

class HarmonicMic(object):
    def __init__(
            self, Z, Phi, Phi_moe, Pr_b, sigma,
            theta, max_iter_harmony, max_iter_kmeans, 
            epsilon_kmeans, epsilon_harmony, K, block_size,
            lamb, verbose, vars_use, phi_dict
    ):
        self.Z_corr = np.array(Z)
        self.Z_orig = np.array(Z)

        self.Z_cos = self.Z_orig / self.Z_orig.max(axis=0)
        self.Z_cos = self.Z_cos / np.linalg.norm(self.Z_cos, ord=2, axis=0)

        self.Phi             = Phi
        self.Phi_moe         = Phi_moe

        # NOTE: added vars_use and phi_dict
        self.vars_use        = vars_use
        self.phi_dict        = phi_dict
        self.N               = self.Z_corr.shape[1]
        self.Pr_b            = Pr_b
        self.B               = self.Phi.shape[0] # number of batch variables
        self.d               = self.Z_corr.shape[0]
        self.window_size     = 3
        self.epsilon_kmeans  = epsilon_kmeans
        self.epsilon_harmony = epsilon_harmony

        self.lamb            = lamb
        self.sigma           = sigma
        self.sigma_prior     = sigma
        self.block_size      = block_size
        self.K               = K                # number of clusters
        self.max_iter_harmony = max_iter_harmony
        self.max_iter_kmeans = max_iter_kmeans
        self.verbose         = verbose
        self.theta           = theta
        # NOTE: added variable - ratio, has the shape of #ofclusters*#ofsamples
        self.ratio           = None

        self.objective_harmony        = []
        self.objective_kmeans         = []
        self.objective_kmeans_dist    = []
        self.objective_kmeans_entropy = []
        self.objective_kmeans_cross   = []
        self.kmeans_rounds  = []

        self.allocate_buffers()
        self.init_cluster()
        self.harmonize(self.max_iter_harmony, self.verbose)

    def result(self):
        return self.Z_corr

    def allocate_buffers(self):
        self._scale_dist = np.zeros((self.K, self.N))
        self.dist_mat    = np.zeros((self.K, self.N))
        self.O           = np.zeros((self.K, self.B))
        self.E           = np.zeros((self.K, self.B))
        self.W           = np.zeros((self.B + 1, self.d))
        self.Phi_Rk      = np.zeros((self.B + 1, self.N))
        # NOTE: added ratio variable
        self.ratio       = np.zeros((self.K, self.Z_corr.shape[1]))

    def init_cluster(self):
        # Start with cluster centroids
        km_centroids, km_labels = kmeans2(self.Z_cos.T, self.K, minit='++')
        print("self.Z_cos.shape", self.Z_cos.shape)
        self.Y = km_centroids.T
        # (1) Normalize
        self.Y = self.Y / np.linalg.norm(self.Y, ord=2, axis=0)
        # (2) Assign cluster probabilities
        self.dist_mat = 2 * (1 - np.dot(self.Y.T, self.Z_cos))
        self.R = -self.dist_mat
        self.R = self.R / self.sigma[:,None]
        self.R -= np.max(self.R, axis = 0)
        self.R = np.exp(self.R)
        self.R = self.R / np.sum(self.R, axis = 0)
        # (3) Batch diversity statistics
        self.E = np.outer(np.sum(self.R, axis=1), self.Pr_b)
        self.O = np.inner(self.R , self.Phi)
        self.compute_objective()
        # Save results
        self.objective_harmony.append(self.objective_kmeans[-1])

    def compute_objective(self):
        kmeans_error = np.sum(np.multiply(self.R, self.dist_mat))
        # Entropy
        _entropy = np.sum(safe_entropy(self.R) * self.sigma[:,np.newaxis])
        # Cross Entropy
        x = (self.R * self.sigma[:,np.newaxis])
        y = np.tile(self.theta[:,np.newaxis], self.K).T
        z = np.log((self.O + 1) / (self.E + 1))
        ### NOTE: added ratio, remember that we are harmonizing for one PC/one feature at a time, so all we worry about is samples
        cluster_by_covar_prob = np.dot(self.R, self.Phi.T)
        cluster_sum = np.sum(cluster_by_covar_prob, axis=1)
        self.ratio = cluster_by_covar_prob/cluster_sum[:, np.newaxis]
        # self.ratio[self.ratio < 0.001*self.K*self.B] = 0.0 # note that e^0=1 so this would make the ratio be 1 and have no impact on the z
        # z = 1.01**(self.ratio) * z # element-wise operation
        z = np.exp(self.ratio) * z # element-wise operation
        # in the toy example, 6 is the cluster size, 175 is the sample size, 5 is the number of one-hot encoded covariate options
        w = np.dot(y * z, self.Phi)
        # the variable blocks
        _cross_entropy = np.sum(x * w) # element-wise multiplicationa and sum
        print(kmeans_error)
        print("x, y, z, w, O, E, K, _cross_entropy shape", x.shape, y.shape, z.shape, w.shape, self.O.shape, self.E.shape, self.K, _cross_entropy)
        print("self.theta", self.theta.shape)

        ### NOTE: added bray-curtis between sample similarity to objective function to minimize
        bray_curtis_sum = 0
        ## compute beta diversity using skbio
        # data = self.Z_corr.T+abs(np.min(self.Z_corr.T)) # to avoid negative values
        # 11/1 change: more scientific way of avoiding negative values
        data = np.where(self.Z_corr.T<np.percentile(self.Z_corr.flatten(), 0.01), 0, self.Z_corr.T)
        data = data+np.abs(np.min(data))
        ids = list(range(self.Z_corr.T.shape[0]))
        bc_dm = beta_diversity("braycurtis", data, ids)
        bc_pc = pcoa(bc_dm)
        ## use the first 3 PCs individually to calculate kruskal-wallis values

        for pc in ['PC1', 'PC2', 'PC3']:
            for var in self.phi_dict.keys():
                current_phi = self.phi_dict[var]
                current_data = bc_pc.samples[pc]
                kruskal_data = []
                for phi_one_condition in list(current_phi):
                    removed_indices = [i for i, e in enumerate(phi_one_condition) if e == 0]
                    current_data_condition = list(current_data.drop(index = removed_indices))
                    kruskal_data.append(current_data_condition)
                bray_curtis_sum += 1-stats.kruskal(*kruskal_data)[1] # the higher this value, the more influence beta diversity has on objective value
                print(bray_curtis_sum)

        # normalize
        bray_curtis_sum = bray_curtis_sum/len(self.phi_dict.keys())

        ### NOTE 1125 change, added shannon index within-sample diversity to objective function to minimize
        shannon_sum = 0
        for var in self.phi_dict.keys():
            current_phi = self.phi_dict[var]
            current_data = bc_pc.samples[pc]
            kruskal_data = []
            for phi_one_condition in list(current_phi):
                removed_indices = [i for i, e in enumerate(phi_one_condition) if e == 0]
                current_data_condition = list(current_data.drop(index = removed_indices))
                kruskal_data.append(current_data_condition)
            shannon_sum += 1-stats.kruskal(*kruskal_data)[1]
         # normalize
        shannon_sum = shannon_sum/len(self.phi_dict.keys())

        # Save results
        # self.objective_kmeans.append(kmeans_error + _entropy + _cross_entropy)
        self.objective_kmeans.append(kmeans_error + _entropy + _cross_entropy + bray_curtis_sum + shannon_sum)
        self.objective_kmeans_dist.append(kmeans_error)
        self.objective_kmeans_entropy.append(_entropy)
        self.objective_kmeans_cross.append(_cross_entropy)
    
    def harmonize(self, iter_harmony=10, verbose=True):
        converged = False
        for i in range(1, iter_harmony + 1):
            if verbose:
                logger.info("Iteration {} of {}".format(i, iter_harmony))
            # STEP 1: Clustering
            self.cluster()
            # STEP 2: Regress out covariates
            # self.moe_correct_ridge()
            self.Z_cos, self.Z_corr, self.W, self.Phi_Rk = moe_correct_ridge(
                self.Z_orig, self.Z_cos, self.Z_corr, self.R, self.W, self.K,
                self.Phi_Rk, self.Phi_moe, self.lamb
            )
            # STEP 3: Check for convergence
            converged = self.check_convergence(1)
            if converged:
                if verbose:
                    logger.info(
                        "Converged after {} iteration{}"
                        .format(i, 's' if i > 1 else '')
                    )
                break
        if verbose and not converged:
            logger.info("Stopped before convergence")
        return 0

    def cluster(self):
        # Z_cos has changed
        # R is assumed to not have changed
        # Update Y to match new integrated data
        self.dist_mat = 2 * (1 - np.dot(self.Y.T, self.Z_cos))
        print("self.dist_mat.shape", self.dist_mat.shape)
        for i in range(self.max_iter_kmeans):
            # print("kmeans {}".format(i))
            # STEP 1: Update Y
            ### NOTE: ANOTHER CHANGE: weight the Rs by feature count weights
            #### 1.1 compute normalized sample sums
            norm_sample_sums = np.sum(self.Z_cos, axis=0)
            #### 1.2 set 75% as the reference
            ref_sample_sum = np.percentile(norm_sample_sums, 75)
            ### 1.3 set R weights as the ratio between feature sum and the reference level
            R_weights = norm_sample_sums/ref_sample_sum
            print("R_weights.shape", R_weights.shape)

            self.Y = np.dot(self.Z_cos, (self.R*R_weights).T)
            self.Y = self.Y / np.linalg.norm(self.Y, ord=2, axis=0)
            # STEP 2: Update dist_mat -> cosine distance
            self.dist_mat = 2 * (1 - np.dot(self.Y.T, self.Z_cos))
            print("self.dist_mat.shape", self.dist_mat.shape)
            # STEP 3: Update R
            self.update_R()
            # maybe do stuff here instead
            # STEP 4: Check for convergence
            self.compute_objective()
            if i > self.window_size:
                converged = self.check_convergence(0)
                if converged:
                    break
        self.kmeans_rounds.append(i)
        self.objective_harmony.append(self.objective_kmeans[-1])
        return 0

    def update_R(self):
        self._scale_dist = -self.dist_mat
        self._scale_dist = self._scale_dist / self.sigma[:,None]
        self._scale_dist -= np.max(self._scale_dist, axis=0)
        self._scale_dist = np.exp(self._scale_dist)
        # Update cells in blocks
        update_order = np.arange(self.N)
        np.random.shuffle(update_order)
        n_blocks = np.ceil(1 / self.block_size).astype(int)
        blocks = np.array_split(update_order, n_blocks)
        print("blocks", len(blocks))
        print("R.shape", self.R.shape)
        print("Z_corr.shape", self.Z_corr.shape)
        for b in blocks:
            # STEP 1: Remove cells
            self.E -= np.outer(np.sum(self.R[:,b], axis=1), self.Pr_b)
            self.O -= np.dot(self.R[:,b], self.Phi[:,b].T)
            # STEP 2: Recompute R for removed cells
            self.R[:,b] = self._scale_dist[:,b]
            self.R[:,b] = np.multiply(
                self.R[:,b],
                np.dot(
                    np.power((self.E + 1) / (self.O + 1), self.theta),
                    self.Phi[:,b]
                )
            )
            self.R[:,b] = self.R[:,b] / np.linalg.norm(self.R[:,b], ord=1, axis=0)
            # STEP 3: Put cells back
            self.E += np.outer(np.sum(self.R[:,b], axis=1), self.Pr_b)
            self.O += np.dot(self.R[:,b], self.Phi[:,b].T)
        return 0

    def check_convergence(self, i_type):
        obj_old = 0.0
        obj_new = 0.0
        # Clustering, compute new window mean
        if i_type == 0:
            okl = len(self.objective_kmeans)
            for i in range(self.window_size):
                obj_old += self.objective_kmeans[okl - 2 - i]
                obj_new += self.objective_kmeans[okl - 1 - i]
            if abs(obj_old - obj_new) / abs(obj_old) < self.epsilon_kmeans:
                return True
            return False
        # Harmony
        if i_type == 1:
            obj_old = self.objective_harmony[-2]
            obj_new = self.objective_harmony[-1]
            if (obj_old - obj_new) / abs(obj_old) < self.epsilon_harmony:
                return True
            return False
        return True


def safe_entropy(x: np.array):
    y = np.multiply(x, np.log(x))
    y[~np.isfinite(y)] = 0.0
    return y

def moe_correct_ridge(Z_orig, Z_cos, Z_corr, R, W, K, Phi_Rk, Phi_moe, lamb):
    Z_corr = Z_orig.copy()
    for i in range(K):
        Phi_Rk = np.multiply(Phi_moe, R[i,:])
        x = np.dot(Phi_Rk, Phi_moe.T) + lamb
        W = np.dot(np.dot(np.linalg.inv(x), Phi_Rk), Z_orig.T)
        W[0,:] = 0 # do not remove the intercept
        W = W.astype("float64")
        Phi_Rk = Phi_Rk.astype("float64")
        Z_corr = Z_corr.astype("float64")
        Z_corr -= np.dot(W.T, Phi_Rk)
    Z_cos = Z_corr / np.linalg.norm(Z_corr, ord=2, axis=0)
    return Z_cos, Z_corr, W, Phi_Rk

