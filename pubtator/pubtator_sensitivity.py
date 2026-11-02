"""Sensitivity analysis of the PubTator "strongly relevant taxa" evaluation (Supplemental Tables 3, 5 and 6).

Two kinds of variation of the default rule (with five key terms, a taxon is strongly relevant if it has >=1
publication for all five, or >=10 publications for at least three of them; with more than five, only the second
criterion is used; see key_terms.py):
  1. thresholds: stricter/looser count cutoffs and number of terms;
  2. popularity normalization: each key-term count divided by the taxon's total PubTator publication count
     (searched alone, {dataset}_taxon_total.csv), i.e. the share of the taxon's literature that mentions the term,
     so heavily studied taxa do not qualify on publication volume alone.
For every rule, dataset and method: number of significant taxa, number strongly relevant, and the proportion
(precision). Method rankings under each rule are compared to the default rule by Spearman correlation.
Also reports the popularity (median taxon publication count) of each method's significant taxa.

Writes to trial_rw/outputs/pubtator_keyterms/sensitivity/:
    sensitivity_long.csv          rule x dataset x method counts and proportions
    sensitivity_rank_agreement.csv Spearman rho vs. the default rule, per dataset and rule
    method_popularity.csv         median / mean log10 publication count of each method's significant taxa
Run after run_pubtator_keyterms.sh (needs a {dataset}__{term}.csv file per key term in the dataset's config YAML,
and {dataset}_taxon_total.csv).
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import REPO_ROOT
from key_terms import key_terms, slug, strongly_relevant

IN = os.path.join(REPO_ROOT, 'trial_rw', 'outputs')
KT = os.path.join(IN, 'pubtator_keyterms')
OUT = os.path.join(KT, 'sensitivity')
os.makedirs(OUT, exist_ok=True)
# datasets to process: given on the command line, or the four real-world datasets of the manuscript
DATASETS = sys.argv[1:] or ['autism_2_microbiomeHD', 'cdi_3_microbiomeHD_clean', 'ibd_3_CMD', 'crc_11_CMD']
# strategies shown in the main figures only (same set as SUCCINCT_DROP_METHODS in benchmark/evaluate.py);
# taxa found only by the dropped strategies are left out of the union
DROP_METHODS = {'bat', 'rbe', 'bmc', 'ruv3', 'svd'}


def n_terms(x, cut):
    return (x >= cut).sum(axis=1)


def make_rules(k):
    """Rule name -> function(counts table C [taxa x k terms], fraction table F) -> boolean Series. The first rule
    is the default; with k = 5 the names and rules are those of the manuscript."""
    default = 'default (all 5 >=1, or >=3 terms >=10)' if k == 5 else 'default (>=3 terms >=10)'
    return {
    default: lambda C, F: strongly_relevant(n_terms(C, 1), n_terms(C, 10), k),
    f'all {k} terms >=1': lambda C, F: n_terms(C, 1) == k,
    f'all {k} terms >=5': lambda C, F: n_terms(C, 5) == k,
    '>=2 terms >=10': lambda C, F: n_terms(C, 10) >= 2,
    '>=3 terms >=10': lambda C, F: n_terms(C, 10) >= 3,
    '>=4 terms >=10': lambda C, F: n_terms(C, 10) >= 4,
    '>=3 terms >=20': lambda C, F: n_terms(C, 20) >= 3,
    '>=3 terms >=50': lambda C, F: n_terms(C, 50) >= 3,
    'popularity-normalized: >=3 terms >=1% of taxon pubs': lambda C, F: n_terms(F, 0.01) >= 3,
    'popularity-normalized: >=3 terms >=5% of taxon pubs': lambda C, F: n_terms(F, 0.05) >= 3,
    'popularity-normalized: >=3 terms >=10% of taxon pubs': lambda C, F: n_terms(F, 0.10) >= 3,
    'popularity-normalized: >=1 term >=10% of taxon pubs': lambda C, F: n_terms(F, 0.10) >= 1,
    }

long_rows, pop_rows, rank_rows = [], [], []
for ds in DATASETS:
    union = pd.read_csv(f'{IN}/{ds}_taxa_union_wilcoxon.csv', index_col=0)
    methods = [c for c in union.columns if c != 'cleaned_taxa' and c not in DROP_METHODS]
    union = union[union[methods].astype(bool).any(axis=1)][['cleaned_taxa'] + methods]
    terms = sorted(key_terms(ds), key=slug)  # same column order as the earlier sorted file listing
    RULES = make_rules(len(terms))
    DEFAULT = next(iter(RULES))
    C = pd.DataFrame(index=union.index)
    for term in terms:
        f = f'{KT}/{ds}__{slug(term)}.csv'
        if not os.path.exists(f):
            raise SystemExit(f'{ds}: missing {os.path.basename(f)}; run run_pubtator_keyterms.sh first')
        t = pd.read_csv(f)
        C[term] = union['cleaned_taxa'].map(dict(zip(t['search_term'], t['pubtator_count'])))
    C = C.astype(float)
    tot_f = f'{KT}/{ds}_taxon_total.csv'
    if os.path.exists(tot_f):
        t = pd.read_csv(tot_f)
        total = union['cleaned_taxa'].map(dict(zip(t['search_term'], t['pubtator_count']))).astype(float)
    else:
        print(f'{ds}: no {os.path.basename(tot_f)}; popularity-normalized rules skipped')
        total = pd.Series(np.nan, index=union.index)
    F = C.div(total.where(total > 0), axis=0)  # share of the taxon's publications mentioning each term

    rel = {}
    for name, rule in RULES.items():
        if name.startswith('popularity') and total.isna().all():
            continue
        r = rule(C.fillna(0), F.fillna(0))
        rel[name] = r
        for m in methods:
            sig = union[m].astype(bool)
            n_sig, n_rel = int(sig.sum()), int((sig & r).sum())
            long_rows.append(dict(dataset=ds, rule=name, method=m, n_significant=n_sig, n_strongly_relevant=n_rel,
                                  proportion_relevant=n_rel / n_sig if n_sig else np.nan))
    for m in methods:
        sig = union[m].astype(bool)
        lp = np.log10(total[sig] + 1)
        pop_rows.append(dict(dataset=ds, method=m, n_significant=int(sig.sum()),
                             median_taxon_pubs=float(total[sig].median()) if sig.any() else np.nan,
                             mean_log10_taxon_pubs=float(lp.mean()) if sig.any() else np.nan))
    pop_rows.append(dict(dataset=ds, method='(all tested taxa in union)', n_significant=len(union),
                         median_taxon_pubs=float(total.median()), mean_log10_taxon_pubs=float(np.log10(total + 1).mean())))

    # agreement of method rankings with the default rule (counts and proportions)
    d = pd.DataFrame([r for r in long_rows if r['dataset'] == ds])
    base = d[d.rule == DEFAULT].set_index('method')
    for name in rel:
        cur = d[d.rule == name].set_index('method').reindex(base.index)
        row = dict(dataset=ds, rule=name, n_taxa_relevant=int(rel[name].sum()), n_taxa=len(union))
        for col in ('n_strongly_relevant', 'proportion_relevant'):
            ok = base[col].notna() & cur[col].notna()
            row[f'spearman_{col}'] = (spearmanr(base.loc[ok, col], cur.loc[ok, col])[0]
                                      if ok.sum() > 2 and base.loc[ok, col].nunique() > 1 and cur.loc[ok, col].nunique() > 1
                                      else np.nan)
        rank_rows.append(row)
    print(f'{ds}: {len(union)} taxa; strongly relevant under default = {int(rel[DEFAULT].sum())}')

long = pd.DataFrame(long_rows)
long.to_csv(f'{OUT}/sensitivity_long.csv', index=False)
ranks = pd.DataFrame(rank_rows)
ranks.to_csv(f'{OUT}/sensitivity_rank_agreement.csv', index=False)
pd.DataFrame(pop_rows).to_csv(f'{OUT}/method_popularity.csv', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30); pd.set_option('display.max_colwidth', 55)
print(ranks.round(2).to_string(index=False))
for ds in long.dataset.unique():
    print(f'\n== {ds}: proportion of significant taxa that are strongly relevant')
    print(long[long.dataset == ds].pivot(index='method', columns='rule', values='proportion_relevant').round(2).T.to_string())
print(OUT)
