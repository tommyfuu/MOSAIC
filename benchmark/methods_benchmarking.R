library(xtable) # table
# library(mixOmics)
library(sva) # ComBat
library(limma) # removeBatchEffect (LIMMA)
library(vegan) # RDA
library(MMUPHin)
library(FDboost)
# source("./adjust_batch.R")
# print("HIIII")
# source("./continuous_discover.R")
# source("./discrete_discover.R")
# source("./controls.R")
# source("./helpers_adjust_batch.R")
# source("./helpers_continuous_discover.R")
# source("./helpers_discrete_discover.R")
# source("./helpers_lm_meta.R")
# source("./helpers.R")
library(ConQuR)
library(doParallel) 
library(dplyr)
library(readr)
library(tibble)
library(mixOmics)

run_methods <- function(data_mat_path, meta_data_path, output_root, batch_ref, dataset = "Dataset", covar = NULL, controlled = FALSE, Sam_id = 'Sam_id', transpose = FALSE, count = FALSE) {
    sink_file_name = paste(output_root, "_runtime.txt", sep="")
    
    # r split and join everything except for last element
    output_dir = strsplit(output_root, "/")
    output_dir = paste(output_dir[[1]][1:(length(output_dir[[1]])-1)], collapse = "/")
    dir.create(file.path(output_dir))

    # data loading <- load data_mat and meta_data, output of the preprocessing.py file
    if(transpose == TRUE) {
        count_data = read.table(data_mat_path, sep=",",header=T,check.names = F)
        # count_data = t(count_data)
    }
    else {
        count_data = read.table(data_mat_path, sep=",",header=T,row.names=1,check.names = F)
    }
    # save count_data
    # write.csv(count_data, paste(output_root, "_count_data.csv", sep=""), row.names = TRUE)
    # print(count_data)
    metadata = as.data.frame(read_csv(meta_data_path))
    # sink(paste(output_root, "_runtime.txt", sep=""))
    
    ## TODO:  potentially need preprocessing such as log and +1
    count_data.clr <- logratio.transfo(count_data+1, logratio = 'CLR')
    cat("runtime documenting...", file=sink_file_name, append=FALSE)

    ## combat
    start_time <- Sys.time()
    batch_info <- as.factor(setNames(as.character(metadata[, dataset]), metadata[[Sam_id]]))
    if(is.null(covar)) {
        count_data.combat <- t(ComBat(t(count_data.clr), batch = batch_info, par.prior=FALSE, mod=NULL))
    }
    else {
        if(length(covar)>1){
            covar = covar[1]
        }
        covar_df = factor(metadata[, covar])
        count_data.combat <- t(ComBat(t(count_data.clr), batch = batch_info, mod = covar_df, par.prior=FALSE, mean.only=TRUE))
    }
    # if(count == FALSE) {
    #     count_data.combat = clr(count_data.combat, inverse = TRUE)
    # }
    write.csv(count_data.combat, paste(output_root, "_combat.csv", sep=""), row.names = TRUE)
    end_time <- Sys.time()
    cat(c("combat runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
    cat('\n', file=sink_file_name, append=TRUE)

    ## combat_seq
    start_time <- Sys.time()
    batch_info <- as.factor(setNames(as.character(metadata[, dataset]), metadata[[Sam_id]]))
    if(is.null(covar)) {
        count_data.combat <- t(ComBat_seq(t(count_data), batch = batch_info))
    }
    else {
        if(length(covar)>1){
            covar = covar[1]
        }
        covar_df = factor(metadata[, covar])
        count_data.combat <- t(ComBat_seq(t(count_data), batch = batch_info, covar_mod = covar_df))
    }
    write.csv(count_data.combat, paste(output_root, "_combat_seq.csv", sep=""), row.names = TRUE)
    end_time <- Sys.time()
    cat(c("combat_seq runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
    cat('\n', file=sink_file_name, append=TRUE)

    ## limma
    ## note that the 'covariates' argument here cannot be used because limma requires continuous covariates
    start_time <- Sys.time()
    count_data.limma <- t(removeBatchEffect(t(count_data.clr), batch = batch_info))
    if(count == TRUE) {
        count_data.limma = clr(count_data.limma, inverse = TRUE)
    }
    write.csv(count_data.limma,paste(output_root, "_limma.csv", sep=""), row.names = TRUE)
    end_time <- Sys.time()
    cat(c("limma runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
    cat('\n', file=sink_file_name, append=TRUE)

    ## MMUPHin
    start_time <- Sys.time()
    if (count == FALSE) {
        count_data_t_relab = t(t(count_data)/rowSums(t(count_data)))
        if(any(is.na(t(count_data_t_relab)))) {
            count_data_t_relab = t(t(count_data+0.001)/rowSums(t(count_data+0.001)))
        }
    }
    else{
        count_data_t_relab = count_data
    }
    metadata_mupphin <- metadata
    row.names(metadata_mupphin) <- metadata[[Sam_id]]
    feature_abd = t(count_data_t_relab)
    # print(t(count_data_t_relab))
    colnames(feature_abd) = row.names(metadata_mupphin)
    if(is.null(covar)) {
        fit_adjust_batch <- adjust_batch(feature_abd = feature_abd,
                                    batch = dataset,
                                    data = metadata_mupphin,
                                    control = list(verbose = FALSE))
    }
    else {
        fit_adjust_batch <- adjust_batch(feature_abd = feature_abd,
                                    batch = dataset,
                                    covariates = covar,
                                    data = metadata_mupphin,
                                    control = list(verbose = FALSE))
    }
    count_data.MMUPHin <- fit_adjust_batch$feature_abd_adj
    write.csv(t(count_data.MMUPHin),paste(output_root, "_MMUPHin.csv", sep=""), row.names = TRUE)
    end_time <- Sys.time()
    cat(c("MMUPHin runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
    cat('\n', file=sink_file_name, append=TRUE)

    ## ConquR
    ## need to put batchid and covar into countdata dataframe
    # count_data.conqur_pre <- count_data
    start_time <- Sys.time()
    batchid <- factor(metadata[, dataset])
    if (is.null(covar)){
        count_data.conqur = ConQuR(tax_tab=count_data, batchid=batchid, covariates=NULL, simple_match = T, batch_ref = batch_ref,
                         logistic_lasso=T, quantile_type="lasso", interplt=F)
    }
    else {
        covar_df = factor(metadata[, covar])
        count_data.conqur = ConQuR(tax_tab=count_data, batchid=batchid, covariates=covar_df, batch_ref = batch_ref,
                         logistic_lasso=T, quantile_type="lasso", interplt=F)
    }
    write.csv(count_data.conqur,paste(output_root, "_ConQuR.csv", sep=""), row.names = TRUE)
    end_time <- Sys.time()
    cat(c("ConquR runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
    cat('\n', file=sink_file_name, append=TRUE)

    ## ConQuR_libsize
    ## need to put batchid and covar into countdata dataframe
    start_time <- Sys.time()
    batchid <- factor(metadata[, dataset])
    if (is.null(covar)){
        count_data.conqur_libsize = ConQuR_libsize(tax_tab=count_data, batchid=batchid, covariates=NULL, simple_match = T, batch_ref = batch_ref,
                         logistic_lasso=T, quantile_type="lasso", interplt=F)
    }
    else {
        covar_df = factor(metadata[, covar])
        count_data.conqur_libsize = ConQuR_libsize(tax_tab=count_data, batchid=batchid, covariates=covar_df, batch_ref = batch_ref,
                         logistic_lasso=T, quantile_type="lasso", interplt=F, num_core=5)
    }
    write.csv(count_data.conqur_libsize,paste(output_root, "_ConQuR_libsize.csv", sep=""), row.names = TRUE)
    end_time <- Sys.time()
    cat(c("ConquR_libsize runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)

    # ## attempt Tune_ConQuR
    # start_time <- Sys.time()
    # batchid <- factor(metadata[, dataset])
    # # here for batch_ref_pool attempt all the potential batch_refs
    # # unique values in the batchid column
    # batch_ref_pool <- unique(metadata[, dataset])

    # print("batch ref pool")
    # print(batch_ref_pool)

    # if (is.null(covar)){
    #     print("Running ConQuR_finetuned version now")
    #     count_data.conqur_finetuned = Tune_ConQuR(tax_tab=count_data, batchid=batchid, covariates=NULL,
    #                        batch_ref_pool=batch_ref_pool,
    #                        logistic_lasso_pool=T, 
    #                        quantile_type_pool=c("standard", "lasso"),
    #                        simple_match_pool=T,
    #                        lambda_quantile_pool=c(NA, "2p/n"),
    #                        interplt_pool=F,
    #                        frequencyL=0,
    #                        frequencyU=1)
    # }
    # else {
    #     covar_df = factor(metadata[, covar])
    #     count_data.conqur_finetuned = Tune_ConQuR(tax_tab=count_data, batchid=batchid, covariates=covar_df,
    #                        batch_ref_pool=batch_ref_pool,
    #                        logistic_lasso_pool=c(T, F), 
    #                        quantile_type_pool=c("standard", "lasso"),
    #                        simple_match_pool=T,
    #                        lambda_quantile_pool=c(NA, "2p/n"),
    #                        interplt_pool=F,
    #                        frequencyL=0,
    #                        frequencyU=1)
    # }
    # write.csv(count_data.conqur_finetuned,paste(output_root, "_ConQuR_finetuned.csv", sep=""), row.names = TRUE)
    # end_time <- Sys.time()
    # cat(c("ConquR_finetuned runtime", toString(end_time - start_time), "seconds"))

   
    
}

# glickman
# run_methods('/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/Glickman_count_data.csv',
# '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/Glickman_meta_data.csv',
# '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_results/Glickman/Glickman',
# dataset = "Dataset",
# batch_ref = 'Old',
# covar='Sex'
# )


# autism 2 microbiomeHD
run_methods('/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/autism_2_microbiomeHD_count_data.csv',
'/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/autism_2_microbiomeHD_meta_data.csv',
'/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD',
dataset = "Dataset",
batch_ref = 'asd_son',
covar = '')

# cdi 3 microbiomeHD
run_methods('/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/cdi_3_microbiomeHD_count_data.csv',
'/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/cdi_3_microbiomeHD_meta_data.csv',
'/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD',
dataset = "Dataset",
batch_ref = 'cdi_schubert')

# # ibd 3 CMD
# run_methods('/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/ibd_3_CMD_count_data.csv',
# '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/ibd_3_CMD_meta_data.csv',
# '/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/ibd_3_CMD/ibd_3_CMD',
# dataset = "study_name",
# batch_ref = 'HMP_2019_ibdmdb')

# adenoma 5 CMD
# run_methods('/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/adenoma_5_CMD_count_data.csv',
# '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/adenoma_5_CMD_meta_data.csv',
# '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_results/adenoma_5_CMD/adenoma_5_CMD',
# dataset = "study_name",
# batch_ref = 'FengQ_2015',
# covar = c("gender")
# )

# # CRC_8_CMD
# run_methods('/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/CRC_8_CMD_count_data.csv',
# '/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/CRC_8_CMD_meta_data.csv',
# '/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/CRC_8_CMD/CRC_8_CMD',
# dataset = "study_name",
# batch_ref = 'FengQ_2015',
# # covar = c("gender")
# )

# # T2D 10 CMD
# run_methods('/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/T2D_10_CMD_count_data.csv',
# '/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/T2D_10_CMD_meta_data.csv',
# '/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/T2D_10_CMD/T2D_10_CMD',
# dataset = "study_name",
# batch_ref = 'Castro-NallarE_2015',
# # covar = c("gender")
# )

# IBD_MDB study
# for (i in c('0.0', '1.0', '2.0', '3.0', '4.0')){
#     run_methods(paste('/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/ibdmdb_interval_', i, '_count_data.csv', sep=""),
#         paste('/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/ibdmdb_interval_', i, '_meta_data.csv', sep=""),
#         paste('/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_results/ibdmdb_interval_', i, '/ibdmdb_interval_', i, sep=""),
#         #     '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/ibdmdb_interval_0.0_count_data.csv',
#         # '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/ibdmdb_interval_0.0_meta_data.csv',
#         # '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_results/ibdmdb_interval_0.0/ibdmdb_interval_0.0',
#         dataset = 'location',
#         batch_ref = 'Los_Angeles',
#         Sam_id = 'patient_visit_id'
#     )
# }

# for (i in c('0.0')){
#     run_methods(paste('/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/ibdmdb_interval_', i, '_count_data.csv', sep=""),
#         paste('/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/ibdmdb_interval_', i, '_meta_data.csv', sep=""),
#         paste('/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_results/ibdmdb_interval_', i, '/ibdmdb_interval_', i, sep=""),
#         #     '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/ibdmdb_interval_0.0_count_data.csv',
#         # '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/ibdmdb_interval_0.0_meta_data.csv',
#         # '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_results/ibdmdb_interval_0.0/ibdmdb_interval_0.0',
#         dataset = 'location',
#         batch_ref = 'Los_Angeles',
#         Sam_id = 'patient_visit_id'
#     )
# }

# # hanningan study
# run_methods('/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/hanninganGD_noBoston_count_data.csv',
# '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/hanninganGD_noBoston_meta_data.csv',
# '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_results/hanninganGD_noBoston/hanninganGD_noBoston',
# dataset = 'location',
# batch_ref = 'Toronto',
# Sam_id = 'patient_visit_id'
# )

# # NOTICED THAT the finetuning version requires at least one covariate
# run_methods('/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/hanninganGD_noBoston_count_data.csv',
# '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/hanninganGD_noBoston_meta_data.csv',
# '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_results/hanninganGD_noBoston_gender/hanninganGD_noBoston_gender',
# dataset = 'location',
# batch_ref = 'Toronto',
# Sam_id = 'patient_visit_id',
# covar = c("gender")
# )


# bin_corr_val_l = c(0, 0.1, 0.3, 0.5, 0.7, 0.9)
# cond_effect_val_l = c(0, 0.099, 0.299, 0.499, 0.699, 0.899)
# batch_effect_val_l = c(0, 0.099, 0.299, 0.499, 0.699, 0.899)
bin_corr_val_l = c(0.3)
cond_effect_val_l = c(0, 0.099, 0.899)
batch_effect_val_l = c(0, 0.099, 0.899)
scaled_midas_methods_bencharking <- function(bin_corr_val_l, cond_effect_val_l, batch_effect_val_l, num_iter){   
  for (bin_corr_val in bin_corr_val_l) {
    for (cond_effect_val in cond_effect_val_l) {
      for (batch_effect_val in batch_effect_val_l) {
        if (cond_effect_val + batch_effect_val <= 1) {
          for (iter in seq(1, num_iter)){
            # print(bin_corr_val)
            # print(cond_effect_val)
            # print(batch_effect_val)
            
            output_file_path_count = paste0("/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/data/simulation_data_midas/ibd_150_count_", bin_corr_val, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter, ".csv")
            # run the methods on this
            run_methods(output_file_path_count,
                        paste0("/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/data/simulation_data_midas/ibd_150_meta_", bin_corr_val, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter, ".csv"),
                        paste0("/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_results/simulation_MIDAS/ibd_150_", bin_corr_val, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter),
                        dataset = 'batchid',
                        batch_ref = "batch_0",
                        Sam_id = 'subjectid_text',
                        transpose = TRUE,
                        count = TRUE,
            )
          }
        }
      }
    }
  }
}

scaled_midas_methods_bencharking(bin_corr_val_l, cond_effect_val_l, batch_effect_val_l, 1)
