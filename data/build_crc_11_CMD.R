library(dplyr)
library(mia)
library(phyloseq)
library(curatedMetagenomicData)
library(yaml)

src <- '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data'
dataset_name <- 'crc_11_CMD'
output_dir <- paste0(src, '/pruned_', dataset_name)
dir.create(output_dir, showWarnings = FALSE)

studies <- c("FengQ_2015", "GuptaA_2019", "HanniganGD_2017",
             "ThomasAM_2018a", "ThomasAM_2018b", "ThomasAM_2019_c",
             "VogtmannE_2016", "WirbelJ_2018", "YachidaS_2019",
             "YuJ_2015", "ZellerG_2014")
conditions <- c("CRC", "healthy")

cat("Fetching data from curatedMetagenomicData...\n")
phylo_obj <-
    filter(sampleMetadata, disease %in% conditions & study_name %in% studies) |>
    select(where(~ !all(is.na(.x)))) |>
    returnSamples("relative_abundance") |>
    makePhyloseqFromTreeSummarizedExperiment(assay.type = "relative_abundance")

cat("Phyloseq object before pruning:\n")
print(phylo_obj)

# Sample-wise cleaning: remove samples below 5th percentile libsize
libsize_threshold_num <- quantile(sample_sums(phylo_obj), 0.05)
phylo_obj <- prune_samples(sample_sums(phylo_obj) > libsize_threshold_num, phylo_obj)

# Taxon-wise cleaning: remove taxa with very low total abundance
phylo_obj <- prune_taxa(taxa_sums(phylo_obj) > 0.05, phylo_obj)

cat("Phyloseq object after pruning:\n")
print(phylo_obj)

# Save otu table, sample metadata, taxonomy table
study_l_str <- paste(studies, collapse = "|")
write.csv(otu_table(phylo_obj),
          paste0(output_dir, "/otu_table_", study_l_str, ".csv"),
          row.names = TRUE)
write.csv(as.matrix(sample_data(phylo_obj)),
          paste0(output_dir, "/sample_table_", study_l_str, ".csv"),
          row.names = TRUE)
write.csv(as.matrix(tax_table(phylo_obj)),
          paste0(output_dir, "/tax_table_", study_l_str, ".csv"),
          row.names = TRUE)

cat("\nSummary of conditions:\n")
print(table(sample_data(phylo_obj)$disease))
cat("\nSummary of studies:\n")
print(table(sample_data(phylo_obj)$study_name))
cat("\nDone. Pruned data saved to:", output_dir, "\n")
