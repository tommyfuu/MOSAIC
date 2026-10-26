#!/usr/bin/env python3
"""Build the union of significant taxa (FDR < 0.05) across batch-correction methods for
each real-world dataset, as input for the PubTator search.
"""
import os
import glob
import pandas as pd
import sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import REPO_ROOT, SIM_OUTPUT_ROOT, BENCH_DIR, DATA_DIR

REPO = REPO_ROOT

DATASETS = {
    "autism_2_microbiomeHD": f"{REPO}/trial_rw/outputs",
    "cdi_3_microbiomeHD_clean": f"{REPO}/trial_rw/outputs",
    "ibd_3_CMD": f"{REPO}/trial_rw/outputs",
    "crc_11_CMD": f"{REPO}/trial_rw/outputs",
}

MICROBIOMEHD_DATASETS = {"autism_2_microbiomeHD", "cdi_3_microbiomeHD_clean"}

for dataset_name, dir_path in DATASETS.items():
    prefix = f"output_{dataset_name}_"
    method_dirs = sorted(glob.glob(f"{dir_path}/{dataset_name}/{prefix}*"))
    methods = [os.path.basename(d)[len(prefix):] for d in method_dirs if os.path.isdir(d)]
    print(f"{dataset_name}: methods found = {methods}")

    current_dataset_method_to_taxa_dict = {}
    for method in methods:
        test_path = f"{dir_path}/{dataset_name}/output_{dataset_name}_{method}/{dataset_name}_{method}_diff_abund_test.csv"
        if not os.path.exists(test_path):
            print(f"  {method}: no diff_abund_test.csv found, skipping")
            continue
        runtime_df = pd.read_csv(test_path, index_col=0)
        runtime_df_sig = runtime_df[runtime_df["FDR_p_value"] < 0.05]
        print(f"  {method}: {runtime_df_sig.shape[0]} significant taxa (FDR<0.05)")
        taxa_l = list(runtime_df_sig.index)

        if dataset_name in MICROBIOMEHD_DATASETS:
            cleaned_l = [" ".join(taxa.split("_")) for taxa in taxa_l]
        else:
            cleaned_l = [" ".join(taxa.split("|")[-1].split("g__")[-1].split("__")[1].split("_")) for taxa in taxa_l]
        current_dataset_method_to_taxa_dict[method] = dict(zip(taxa_l, cleaned_l))

    original_taxa_to_genus_dict = {}
    for method in current_dataset_method_to_taxa_dict:
        original_taxa_to_genus_dict.update(current_dataset_method_to_taxa_dict[method])

    taxa_union = set()
    for method in current_dataset_method_to_taxa_dict:
        taxa_union.update(current_dataset_method_to_taxa_dict[method].keys())
    taxa_union = sorted(taxa_union)

    taxa_union_df = pd.DataFrame(index=taxa_union)
    taxa_union_df["cleaned_taxa"] = [original_taxa_to_genus_dict[t] for t in taxa_union_df.index]
    for method in current_dataset_method_to_taxa_dict:
        taxa_union_df[method] = taxa_union_df.index.isin(current_dataset_method_to_taxa_dict[method].keys())

    out_path = f"{dir_path}/{dataset_name}_taxa_union_wilcoxon.csv"
    taxa_union_df.to_csv(out_path)
    print(f"  -> wrote {out_path} ({len(taxa_union)} unique significant taxa across {len(current_dataset_method_to_taxa_dict)} methods)")

print("DONE")
