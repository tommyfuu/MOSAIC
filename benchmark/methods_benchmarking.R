library(xtable) # table
# library(mixOmics)
library(sva) # ComBat
library(limma) # removeBatchEffect (LIMMA)
library(vegan) # RDA
library(MMUPHin)
library(ConQuR)
library(doParallel) 
library(dplyr)
library(readr)
library(tibble)
library(mixOmics)

run_methods <- function(data_mat_path, meta_data_path, output_root, batch_ref, dataset = "Dataset", covar = NULL, controlled = FALSE, Sam_id = 'Sam_id') {
    # r split and join everything except for last element
    output_dir = strsplit(output_root, "/")
    output_dir = paste(output_dir[[1]][1:(length(output_dir[[1]])-1)], collapse = "/")
    dir.create(file.path(output_dir))

    # data loading <- load data_mat and meta_data, output of the preprocessing.py file
    count_data = read.table(data_mat_path, sep=",",header=T,row.names=1,check.names = F)

    # print(count_data)
    metadata = as.data.frame(read_csv(meta_data_path))
    sink(paste(output_root, "_runtime.txt", sep=""))

    ## TODO:  potentially need preprocessing such as log and +1
    count_data.clr <- logratio.transfo(count_data+1, logratio = 'CLR')

    
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
    write.csv(count_data.combat, paste(output_root, "_combat.csv", sep=""), row.names = TRUE)
    end_time <- Sys.time()
    cat(c("combat runtime", toString(end_time - start_time), "seconds"))

    ## limma
    ## note that the 'covariates' argument here cannot be used because limma requires continuous covariates
    start_time <- Sys.time()
    count_data.limma <- t(removeBatchEffect(t(count_data.clr), batch = batch_info))
    write.csv(count_data.limma,paste(output_root, "_limma.csv", sep=""), row.names = TRUE)
    end_time <- Sys.time()
    cat(c("limma runtime", toString(end_time - start_time), "seconds"))
    cat('\n')

    ## MMUPHin
    start_time <- Sys.time()
    count_data_t_relab = t(t(count_data)/rowSums(t(count_data)))
    write.csv(count_data_t_relab, "WTF.csv")
    if(any(is.na(t(count_data_t_relab)))) {
        count_data_t_relab = t(t(count_data+0.001)/rowSums(t(count_data+0.001)))
    }
    metadata_mupphin <- metadata
    row.names(metadata_mupphin) <- metadata[[Sam_id]]
    if(is.null(covar)) {
        fit_adjust_batch <- adjust_batch(feature_abd = t(count_data_t_relab),
                                    batch = dataset,
                                    data = metadata_mupphin,
                                    control = list(verbose = FALSE))
    }
    else {
        fit_adjust_batch <- adjust_batch(feature_abd = t(count_data_t_relab),
                                    batch = dataset,
                                    covariates = covar,
                                    data = metadata_mupphin,
                                    control = list(verbose = FALSE))
    }
    count_data.MMUPHin <- fit_adjust_batch$feature_abd_adj
    write.csv(t(count_data.MMUPHin),paste(output_root, "_MMUPHin.csv", sep=""), row.names = TRUE)
    end_time <- Sys.time()
    cat(c("MMUPHin runtime", toString(end_time - start_time), "seconds"))
    cat('\n')

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
    cat(c("ConquR runtime", toString(end_time - start_time), "seconds"))
    cat('\n')

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
    cat(c("ConquR_libsize runtime", toString(end_time - start_time), "seconds"))

    ## attempt Tune_ConQuR
    start_time <- Sys.time()
    batchid <- factor(metadata[, dataset])
    # here for batch_ref_pool attempt all the potential batch_refs
    # unique values in the batchid column
    batch_ref_pool <- unique(metadata[, dataset])

    print("batch ref pool")
    print(batch_ref_pool)

    if (is.null(covar)){
        print("Running ConQuR_finetuned version now")
        count_data.conqur_finetuned = Tune_ConQuR(tax_tab=count_data, batchid=batchid, covariates=NULL,
                           batch_ref_pool=batch_ref_pool,
                           logistic_lasso_pool=T, 
                           quantile_type_pool=c("standard", "lasso"),
                           simple_match_pool=T,
                           lambda_quantile_pool=c(NA, "2p/n"),
                           interplt_pool=F,
                           frequencyL=0,
                           frequencyU=1)
    }
    else {
        covar_df = factor(metadata[, covar])
        count_data.conqur_finetuned = Tune_ConQuR(tax_tab=count_data, batchid=batchid, covariates=covar_df,
                           batch_ref_pool=batch_ref_pool,
                           logistic_lasso_pool=c(T, F), 
                           quantile_type_pool=c("standard", "lasso"),
                           simple_match_pool=c(T, F),
                           lambda_quantile_pool=c(NA, "2p/n"),
                           interplt_pool=F,
                           frequencyL=0,
                           frequencyU=1)
    }
    write.csv(count_data.conqur_finetuned,paste(output_root, "_ConQuR_finetuned.csv", sep=""), row.names = TRUE)
    end_time <- Sys.time()
    cat(c("ConquR_finetuned runtime", toString(end_time - start_time), "seconds"))

   
    
}

# glickman
# run_methods('/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/Glickman_count_data.csv',
# '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/Glickman_meta_data.csv',
# '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_results/Glickman/Glickman',
# dataset = "Dataset",
# batch_ref = 'Old',
# covar='Sex'
# )


# # autism 2 microbiomeHD
# run_methods('/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/autism_2_microbiomeHD_count_data.csv',
# '/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/autism_2_microbiomeHD_meta_data.csv',
# '/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD',
# dataset = "Dataset",
# batch_ref = 'asd_son')

# # cdi 3 microbiomeHD
# run_methods('/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/cdi_3_microbiomeHD_count_data.csv',
# '/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/cdi_3_microbiomeHD_meta_data.csv',
# '/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD',
# dataset = "Dataset",
# batch_ref = 'cdi_schubert')

# # ibd 3 CMD
# run_methods('/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/ibd_3_CMD_count_data.csv',
# '/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/ibd_3_CMD_meta_data.csv',
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

## NOTICED THAT the finetuning version requires at least one covariate
run_methods('/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/hanninganGD_noBoston_count_data.csv',
'/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/hanninganGD_noBoston_meta_data.csv',
'/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_results/hanninganGD_noBoston_gender/hanninganGD_noBoston_gender',
dataset = 'location',
batch_ref = 'Toronto',
Sam_id = 'patient_visit_id',
covar = c("gender")
)