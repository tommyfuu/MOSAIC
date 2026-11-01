"""Build the strongly-relevant-taxa tables (preprint Supplemental Tables 3-4) from the key-term searches.

Key terms come from each dataset's config YAML (key_terms.py). With five terms, a taxon is strongly relevant if it
has >=1 publication for all five, or >=10 publications for at least three of them; with more than five, only the
second criterion is used.

Writes to trial_rw/outputs/pubtator_keyterms/:
    strongly_relevant_counts_by_method.csv   (Table 3: significant / strongly relevant taxa per method)
    {dataset}_taxa_relevance.csv             (Table 4: per-taxon term counts, methods detecting it)
"""
import os
import sys

import pandas as pd

sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import REPO_ROOT
from key_terms import key_terms, slug, strongly_relevant

IN = os.path.join(REPO_ROOT, 'trial_rw', 'outputs')
OUT = os.path.join(IN, 'pubtator_keyterms')
# datasets to process: given on the command line, or the four real-world datasets of the manuscript
DATASETS = sys.argv[1:] or ['autism_2_microbiomeHD', 'cdi_3_microbiomeHD_clean', 'ibd_3_CMD', 'crc_11_CMD']

counts_rows = []
for ds in DATASETS:
    union = pd.read_csv(f'{IN}/{ds}_taxa_union_wilcoxon.csv', index_col=0)
    methods = [c for c in union.columns if c != 'cleaned_taxa']
    terms = {}
    for term in sorted(key_terms(ds), key=slug):  # same column order as the earlier sorted file listing
        f = f'{OUT}/{ds}__{slug(term)}.csv'
        if not os.path.exists(f):
            raise SystemExit(f'{ds}: missing {os.path.basename(f)}; run run_pubtator_keyterms.sh first')
        t = pd.read_csv(f)
        terms[term] = dict(zip(t['search_term'], t['pubtator_count']))
    tab = union[['cleaned_taxa']].copy()
    for term, d in terms.items():
        tab[term] = tab['cleaned_taxa'].map(d)
    tc = tab[list(terms)]
    tab['n_terms_ge1'] = (tc >= 1).sum(axis=1)
    tab['n_terms_ge10'] = (tc >= 10).sum(axis=1)
    tab['strongly_relevant'] = strongly_relevant(tab['n_terms_ge1'], tab['n_terms_ge10'], len(terms))
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
