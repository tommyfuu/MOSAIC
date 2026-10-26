"""Excel version of the PubTator sensitivity analysis (output of pubtator_sensitivity.py).

Sheets: README, rank_agreement, one counts and one proportions sheet per dataset (rules x methods), method_popularity.
Writes trial_rw/outputs/pubtator_keyterms/sensitivity/SuppTable_pubtator_sensitivity.xlsx. Needs openpyxl.
Usage: python3 make_sensitivity_xlsx.py
"""
import os
import sys

import pandas as pd

sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import REPO_ROOT

D = os.path.join(REPO_ROOT, 'trial_rw', 'outputs', 'pubtator_keyterms', 'sensitivity')
long = pd.read_csv(f'{D}/sensitivity_long.csv')
ranks = pd.read_csv(f'{D}/sensitivity_rank_agreement.csv')
pop = pd.read_csv(f'{D}/method_popularity.csv')
SHORT = {'autism_2_microbiomeHD': 'autism', 'cdi_3_microbiomeHD_clean': 'cdi', 'ibd_3_CMD': 'ibd', 'crc_11_CMD': 'crc'}
readme = pd.DataFrame({'note': [
    'Sensitivity of the "strongly relevant taxa" evaluation (Supplemental Tables 3-4) to the relevance rule.',
    'Default rule: >=1 publication for all five key terms, or >=10 publications for at least three of them.',
    'Strategies: those shown in the main figures (original data, ComBat-seq/ComBat, limma, MMUPHin, ConQuR, ConQuR-libsize/ConQuR-rel, Harmony, Percentile Normalization, DEBIAS-M); taxa found only by other strategies are excluded.',
    'Threshold rules vary the publication cutoffs and the number of key terms required.',
    'Popularity-normalized rules divide each key-term count by the taxon\'s total PubTator publication count, '
    'so a taxon qualifies by the share of its literature that mentions the term rather than by volume.',
    'rank_agreement: Spearman correlation of the method ranking under each rule with the default rule '
    '(by number and by proportion of significant taxa that are strongly relevant).',
    'counts_* / proportions_*: per dataset, rules (rows) x methods (columns).',
    'method_popularity: median and mean log10 total publication count of each method\'s significant taxa.',
    'IBD key terms: inflammatory bowel disease, Crohn\'s disease, ulcerative colitis, abdominal pain, intestinal inflammation.',
    'CRC key terms: colon cancer, colorectal cancer, tumorigenesis, rectal bleeding, MSI subtype.']})
with pd.ExcelWriter(f'{D}/SuppTable_pubtator_sensitivity.xlsx') as xw:
    readme.to_excel(xw, sheet_name='README', index=False)
    ranks.round(3).to_excel(xw, sheet_name='rank_agreement', index=False)
    for ds, s in SHORT.items():
        sub = long[long.dataset == ds]
        if sub.empty:
            continue
        sub.pivot(index='rule', columns='method', values='n_strongly_relevant').to_excel(xw, sheet_name=f'counts_{s}')
        sub.pivot(index='rule', columns='method', values='proportion_relevant').round(3).to_excel(xw, sheet_name=f'proportions_{s}')
    pop.round(3).to_excel(xw, sheet_name='method_popularity', index=False)
print(f'{D}/SuppTable_pubtator_sensitivity.xlsx')
