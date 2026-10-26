# PubTator literature check (Supplemental Tables 3-4)

Run from `pubtator/`, after the real-world evaluations (`trial_rw/outputs`) are up to date:

1. `python3 rw_data_pubtator_prep.py`: for each real-world dataset, takes the union of taxa with
   FDR < 0.05 in any method's `*_diff_abund_test.csv`, cleans the taxon names (genus level for
   MicrobiomeHD, species level for CMD) and writes `trial_rw/outputs/{dataset}_taxa_union_wilcoxon.csv`.
2. `sbatch run_pubtator_keyterms.sh`: searches every taxon together with each of the five clinical key
   terms per disease (as in the preprint Methods) using `pubtator_api_search.py`, and writes
   `trial_rw/outputs/pubtator_keyterms/{dataset}__{term}.csv`. Sequential, ~1 h for ~1,000 taxa
   (PubTator3 allows ~3 requests/sec); rerunning resumes from existing files.
3. `python3 pubtator_relevance_tables.py`: strongly relevant taxa (>=1 publication for all 5 terms, or
   >=10 publications for at least 3) per method (`strongly_relevant_counts_by_method.csv`, Supp Table 3)
   and per taxon (`{dataset}_taxa_relevance.csv`).
4. `python3 make_supp_table4.py [--succinct]`: Supp Table 4 as an Excel workbook with live formulas
   (`SuppTable4_pubtator_keyterms[_succinct].xlsx`); `--succinct` drops bat/ruv3/svd/rbe/bmc. Needs `openpyxl`.
5. `python3 pubtator_sensitivity.py`: sensitivity of the strongly-relevant rule to the publication thresholds and to
   normalization by taxon popularity (total publications per taxon, `{dataset}_taxon_total.csv`); writes
   `pubtator_keyterms/sensitivity/`. `python3 make_sensitivity_xlsx.py` turns it into an Excel table (needs `openpyxl`).

IBD key terms were revised in October 2026 (previously "Irritable Bowel Syndrome", "FODMAP diet", "inflammation");
the old results are kept in `pubtator_keyterms/archive_ibd_terms_v1_IBS/`. The CRC term "diarrhea" was replaced by
"tumorigenesis" at the same time (old results in `pubtator_keyterms/archive_crc_terms_v1_diarrhea/`).

`run_pubtator_search.sh` runs a single condition term per dataset (quick check, not used for the tables).
`archive/` holds the older Selenium-based scraper and earlier versions of the prep script.
