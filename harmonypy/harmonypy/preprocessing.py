# harmonicMic - A data alignment algorithm dedicated to microbiome data.
# Copyright (C) 2022  Chenlian (Tom) Fu <chf4012@med.cornell.edu; tfu@g.hmc.edu>
#               2019  Kamil Slowikowski <kslowikowski@gmail.com>
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

def load_data(address_X, address_Y, IDCol, index_col = False):
    if index_col != False:
        data_mat = pd.read_csv(address_X, index_col=index_col)
    else:
        data_mat = pd.read_csv(address_X)
    meta_data = pd.read_csv(address_Y)
    data_mat = preprocess(data_mat, meta_data, IDCol)
    data_mat = np.array(data_mat)

    return data_mat

def preprocess(data_mat, meta_data, IDCol):
    print(data_mat.shape)
    print(data_mat)
    # remove samples with all zeros
    data_mat = data_mat.loc[~(data_mat==0).all(axis=1)]
    kept_samples = data_mat.index
    meta_data = meta_data[meta_data[IDCol].isin(kept_samples)]

    print(data_mat.shape)
    # remove features where the sum of counts are below 0.01% compared to the total sum of all counts (include those features with all zeroes)
    # Arumugam, Manimozhiyan, Jeroen Raes, Eric Pelletier, Denis Le Paslier, Takuji Yamada, Daniel R Mende, Gabriel R Fernandes, et al. 2011. “Enterotypes of the Human Gut Microbiome.” Nature 473 (7346). Nature Publishing Group: 174.
    col_names = list(data_mat.columns)
    col_sums = data_mat.sum(axis = 0)
    total = data_mat.to_numpy().sum()
    col_sums_ratios = [col_sum/total for col_sum in col_sums]
    removable_feature_names = [col_names[index] for index, ratio in enumerate(col_sums_ratios) if ratio<0.0001]
    print(removable_feature_names)
    data_mat.drop(removable_feature_names, axis=1, inplace=True)
    print(data_mat.shape)
    return data_mat, meta_data


