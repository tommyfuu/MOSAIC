"""Build the strongly-relevant-taxa tables (preprint Supplemental Tables 3-4) from the key-term searches.

A taxon is strongly relevant if it has >=1 publication for all five key terms, or >=10 publications
for at least three of them.

Writes to trial_rw/outputs/pubtator_keyterms/:
    strongly_relevant_counts_by_method.csv   (Table 3: significant / strongly relevant taxa per method)
    {dataset}_taxa_relevance.csv             (Table 4: per-taxon term counts, methods detecting it)
"""
import glob
import os
import sys

import pandas as pd

sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import REPO_ROOT

IN = os.path.join(REPO_ROOT, 'trial_rw', 'outputs')
OUT = os.path.join(IN, 'pubtator_keyterms')
DATASETS = ['autism_2_microbiomeHD', 'cdi_3_microbiomeHD_clean', 'ibd_3_CMD', 'crc_11_CMD']

counts_rows = []
for ds in DATASETS:
    union = pd.read_csv(f'{IN}/{ds}_taxa_union_wilcoxon.csv', index_col=0)
    methods = [c for c in union.columns if c != 'cleaned_taxa']
    terms = {}
    for f in sorted(glob.glob(f'{OUT}/{ds}__*.csv')):
        t = pd.read_csv(f)
        terms[t['condition'].iloc[0] if len(t) else os.path.basename(f)] = dict(zip(t['search_term'], t['pubtator_count']))
    if len(terms) != 5:
        print(f'{ds}: expected 5 key-term files, found {len(terms)}; skipping')
        continue
    tab = union[['cleaned_taxa']].copy()
    for term, d in terms.items():
        tab[term] = tab['cleaned_taxa'].map(d)
    tc = tab[list(terms)]
    tab['n_terms_ge1'] = (tc >= 1).sum(axis=1)
    tab['n_terms_ge10'] = (tc >= 10).sum(axis=1)
    tab['strongly_relevant'] = (tab['n_terms_ge1'] == 5) | (tab['n_terms_ge10'] >= 3)
    tab['methods_detecting'] = union[methods].apply(lambda r: ';'.join(m for m in methods if bool(r[m])), axis=1)
    tab.to_csv(f'{OUT}/{ds}_taxa_relevance.csv')
    for m in methods:
        sig = union[m].astype(bool)
        counts_rows.append(dict(dataset=ds, method=m, n_significant=int(sig.sum()),
                                n_strongly_relevant=int((sig & tab['strongly_relevant']).sum())))
    print(f'{ds}: {len(tab)} taxa, {int(tab["strongly_relevant"].sum())} strongly relevant')

counts = pd.DataFrame(counts_rows)
counts.to_csv(f'{OUT}/strongly_relevant_counts_by_method.csv', index=False)
print(counts.pivot(index='method', columns='dataset', values='n_strongly_relevant').to_string())
