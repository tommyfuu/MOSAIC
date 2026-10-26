"""Real-data ranges of the parameters varied in the simulation ablations (sample size, number of taxa, zero fraction,
condition balance, batch balance), per dataset and per study (batch), to anchor the ablation settings.

Writes realworld_parameter_ranges_by_dataset.csv and realworld_parameter_ranges_by_study.csv to
trial_rw/outputs/ablation_anchoring/.
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import REPO_ROOT

DATASETS = {  # name: (batch column, condition column, case label)
    'autism_2_microbiomeHD': ('Dataset', 'DiseaseState', 'ASD'),
    'cdi_3_microbiomeHD_clean': ('Dataset', 'DiseaseState', 'CDI'),
    'ibd_3_CMD': ('study_name', 'disease', 'IBD'),
    'crc_11_CMD': ('study_name', 'disease', 'CRC'),
}
OUT = os.path.join(REPO_ROOT, 'trial_rw', 'outputs', 'ablation_anchoring')
os.makedirs(OUT, exist_ok=True)

ds_rows, st_rows = [], []
for ds, (bcol, ccol, case) in DATASETS.items():
    d = os.path.join(REPO_ROOT, 'data', 'cleaned_data', ds)
    X = pd.read_csv(f'{d}/{ds}_count_data.csv', index_col=0)
    M = pd.read_csv(f'{d}/{ds}_meta_data.csv')
    assert len(M) == len(X)
    Z = X.values == 0
    batches = M[bcol].astype(str).values
    case_frac = (M[ccol].astype(str) == case).mean()
    bshare = pd.Series(batches).value_counts(normalize=True)
    lib = X.sum(axis=1)
    ds_rows.append(dict(dataset=ds, n_samples=len(X), n_taxa=X.shape[1], n_studies=len(bshare),
                        zero_fraction=Z.mean(), case_fraction=case_frac,
                        smallest_study_share=bshare.min(), largest_study_share=bshare.max(),
                        median_library_size=np.median(lib)))
    for b in bshare.index:
        idx = batches == b
        Xb = X.values[idx]
        st_rows.append(dict(dataset=ds, study=b, n_samples=int(idx.sum()), study_share=bshare[b],
                            n_taxa_present=int((Xb.sum(axis=0) > 0).sum()), zero_fraction=(Xb == 0).mean(),
                            case_fraction=(M.loc[idx, ccol].astype(str) == case).mean()))
ds_t = pd.DataFrame(ds_rows); st_t = pd.DataFrame(st_rows)
ds_t.round(4).to_csv(f'{OUT}/realworld_parameter_ranges_by_dataset.csv', index=False)
st_t.round(4).to_csv(f'{OUT}/realworld_parameter_ranges_by_study.csv', index=False)
pd.set_option('display.width', 220)
print(ds_t.round(3).to_string(index=False)); print()
print(st_t.round(3).to_string(index=False)); print()
print('per-study ranges: n', st_t.n_samples.min(), '-', st_t.n_samples.max(), '| taxa present', st_t.n_taxa_present.min(), '-', st_t.n_taxa_present.max(),
      '| zero frac', round(st_t.zero_fraction.min(), 3), '-', round(st_t.zero_fraction.max(), 3),
      '| case frac', round(st_t.case_fraction.min(), 3), '-', round(st_t.case_fraction.max(), 3))
