library(dplyr)
library(readr)
library(mia)
library(phyloseq)
library(curatedMetagenomicData)

# do argparse for R that takes in overall_path, study_name, libsize_threshold, and relab_threshold
args <- commandArgs(trailingOnly = TRUE)
print(args)
if (length(args) == 2) {
    dataset_type <- args[1] # either 'microbiomeHD' or 'CMD'
    study_name <- args[2]
} else {
    print("Usage: Rscript step1_clean_and_prune.R overall_path study_name")
    quit()
}




load_phyloseq_from_merged_microbiomeHD <- function (otu_mat_path, meta_data_path){
    # data loading
    otu_mat<- t(read.csv(otu_mat_path, row.names=1))
    
    # split otu_mat's rownames(otu) into a dataframe by splitting with '.' and get the existing taxa classes/levels
    taxa_l = lapply(rownames(otu_mat), function(x) strsplit(x, "[.]")[[1]])
    taxa_classes <- lapply(taxa_l[[1]], function(x) strsplit(x, '__')[[1]][1])
    taxa_l = lapply(taxa_l, function(x) lapply(x, function(y) strsplit(y, '__')[[1]][2][1]))
    taxa_mat = as.data.frame(t(as.data.frame(do.call(cbind, taxa_l))))
    colnames(taxa_mat) <- taxa_classes
    
    # assign otu_mat's rownames as the first column of taxa_mat
    rownames(otu_mat) = rownames(taxa_mat)

    # clean taxa_mat and otu_mat for taxa that do not have all levels
    cleaned_taxa_mat <- taxa_mat[!is.na(taxa_mat$'g'), ]
    dim(cleaned_taxa_mat)
    # get otu_mat for cleaned taxa_mat
    cleaned_otu_mat <- otu_mat[rownames(cleaned_taxa_mat), ]

    # reasigns the rownames of otu_mat as sp+number
    rownames(cleaned_otu_mat) <- paste0('sp', 1:nrow(cleaned_otu_mat))

    # load metadata
    meta_data = read.csv(meta_data_path)
    meta_data <- meta_data %>% 
        tibble::column_to_rownames("Sam_id") 

    # turn data into phyloseq object
    phyloseq_dataset <- phyloseq(otu_table(cleaned_otu_mat, taxa_are_rows = TRUE), sample_data(meta_data), tax_table(cleaned_taxa_mat))

    return(phyloseq_dataset)
}

load_phyloseq_from_merged_CMD <- function (list_of_studies, list_of_conditions){
    study_l_str <- paste(list_of_studies,collapse="|")
    study_l_abund_str <- paste("(", study_l_str, ").", "relative_abundance", sep="")

    current_phylo_obj <-
        filter(sampleMetadata, disease %in% list_of_conditions & study_name %in% list_of_studies) |>
        select(where(~ !all(is.na(.x)))) |>
        returnSamples("relative_abundance") |>
        makePhyloseqFromTreeSummarizedExperiment(assay.type = "relative_abundance")

    return (list(current_phylo_obj, study_l_str))
}

clean_prune_save_phyloseq <- function (phyloseq_dataset, out_str, libsize_threshold, relab_threshold, save = FALSE, save_to = NULL){
    # sample-wise cleaning: remove samples with less than 0.05% libsize_threshold percentile
    print(phyloseq_dataset)
    libsize_threshold_num = quantile(sample_sums(phyloseq_dataset), libsize_threshold)
    phyloseq_dataset <- prune_samples(sample_sums(phyloseq_dataset) > libsize_threshold_num, phyloseq_dataset)
    print(phyloseq_dataset)
    # taxon-wise cleaning: remove taxa with less than relab_threshold relative abundance
    phyloseq_dataset <- prune_taxa(taxa_sums(phyloseq_dataset) > relab_threshold, phyloseq_dataset)
    print(phyloseq_dataset)
    if(save){
        # currently can only save otu table, metadata, and taxonomy info
        write.csv(otu_table(phyloseq_dataset), paste(save_to, "/otu_table_", out_str, '.csv', sep=""), row.names = TRUE)
        write.csv(as.matrix(sample_data(phyloseq_dataset)), paste(save_to, "/sample_table_", out_str, '.csv', sep=""), row.names = TRUE)
        write.csv(tax_table(phyloseq_dataset), paste(save_to, "/tax_table_", out_str, '.csv', sep=""), row.names = TRUE)
    }
    return(phyloseq_dataset)
}


# overall_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/'
# autism_phyloseq_obj <- load_phyloseq_from_merged_microbiomeHD(paste0(overall_path, 'intermediate_autism_2_microbiomeHD/autism_2_microbiomeHD_count_data.csv'), paste0(overall_path, 'intermediate_autism_2_microbiomeHD/autism_2_microbiomeHD_meta_data.csv'))
# clean_prune_save_phyloseq(autism_phyloseq_obj, 'autism_2_microbiomeHD', 0.05, 0.05, save = TRUE, save_to = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/pruned_autism_2_microbiomeHD')

# cdi_phyloseq_obj <- load_phyloseq_from_merged_microbiomeHD(paste0(overall_path, 'intermediate_cdi_3_microbiomeHD/cdi_3_microbiomeHD_count_data.csv'), paste0(overall_path, 'intermediate_cdi_3_microbiomeHD/cdi_3_microbiomeHD_meta_data.csv'))
# clean_prune_save_phyloseq(cdi_phyloseq_obj, 'cdi_2_microbiomeHD', 0.05, 0.05, save = TRUE, save_to = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/pruned_cdi_3_microbiomeHD')

# ibd_phyloseq_obj <- load_phyloseq_from_merged_CMD(c("HMP_2019_ibdmdb", "LiJ_2014", "NielsenHB_2014"), c("IBD", 'healthy'))
# clean_prune_save_phyloseq(ibd_phyloseq_obj[[1]], ibd_phyloseq_obj[[2]], 0.05, 0.05, save = TRUE, save_to = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/pruned_ibd_3_CMD')

# crc_phyloseq_obj <- load_phyloseq_from_merged_CMD(c("FengQ_2015", "HanniganGD_2017", "ThomasAM_2018a", "YachidaS_2019", "ZellerG_2014"), c("adenoma", "healthy"))
# clean_prune_save_phyloseq(crc_phyloseq_obj[[1]], crc_phyloseq_obj[[2]], 0.05, 0.05, save = TRUE, save_to = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/pruned_crc_8_CMD')

# find the path to this Rscript

## set working directory to the location of this Rscript


library(yaml)
current_path = getwd()
if (grepl("/data", current_path)){
    config_object <- yaml.load_file('../config.yml')
} else{
    config_object <- yaml.load_file('./config.yml')
}
# config_object <- yaml.load_file('./config.yml')
src_path <- config_object$src
libsize_threshold <- config_object$RW_LIBSIZE_THRESHOLD
relab_threshold <- config_object$RW_RELAB_THRESHOLD

# make sure they are numeric
libsize_threshold <- as.numeric(libsize_threshold)
relab_threshold <- as.numeric(relab_threshold)
print(libsize_threshold)
print(relab_threshold)

# mkdir
output_dir <- paste0(src_path, '/pruned_', study_name)
dir.create(output_dir, showWarnings = FALSE)

if(dataset_type == 'microbiomeHD'){
    phyloseq_obj <- load_phyloseq_from_merged_microbiomeHD(paste0(src_path, '/intermediate_', study_name, '/', study_name, '_count_data.csv'), paste0(src_path, '/intermediate_', study_name, '/', study_name, '_meta_data.csv'))
    clean_prune_save_phyloseq(phyloseq_obj, study_name, libsize_threshold, relab_threshold, save = TRUE, save_to = output_dir)
} else if(dataset_type == 'CMD'){
    study_list <- config_object$CMD_datasets
    condition_list <- config_object$CMD_conditions
    phyloseq_obj <- load_phyloseq_from_merged_CMD(list_of_studies = study_list, list_of_conditions = condition_list)
    clean_prune_save_phyloseq(phyloseq_obj[[1]], phyloseq_obj[[2]], libsize_threshold, relab_threshold, save = TRUE, save_to = output_dir)

} else {
    print("dataset_type must be either 'microbiomeHD' or 'CMD'")
    quit()
}
