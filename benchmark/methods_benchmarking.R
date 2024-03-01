#!/usr/bin/env Rscript
library(xtable) # table
# library(mixOmics)
library(sva) # ComBat
library(limma) # removeBatchEffect (LIMMA)
library(vegan) # RDA
library(MMUPHin)
library(FDboost)  
library(doParallel) 
library(dplyr)
library(readr)
library(tibble)
library(stringr)
library(mixOmics)

args = commandArgs(trailingOnly=TRUE)

print(args)
if (length(args)==0 || length(args)>1 ) {
#   stop("There has to be exactly one argument supplied with this script for simulation runs", call.=FALSE)
    option = 5
    print("not running on real-world dataset")
} else if (length(args)==1) {
  # default output file
  option= args[1]
  print(option)
}

## load ConQuR
# current_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR'
current_path = getwd()
conqur_path = str_replace(current_path, "mic_bc_benchmark/benchmark", "mic_bc_benchmark/ConQuR")
source(paste0(conqur_path, "/ConQuR_help_functions.R"))
source(paste0(conqur_path, "/ConQuR_main_tune.R"))
source(paste0(conqur_path, "/ConQuR_help_functions_libsize_old.R"))
source(paste0(conqur_path, "/ConQuR_main_tune_libsize_old.R"))
source(paste0(conqur_path, "/ConQuR_help_functions_rel.R"))
source(paste0(conqur_path, "/ConQuR_main_tune_rel.R"))
source(paste0(conqur_path, "/supporting_functions.R"))


run_methods <- function(data_mat_path, meta_data_path, output_root, batch_ref, dataset = "Dataset", covar = NULL, controlled = FALSE, Sam_id = 'Sam_id', transpose = FALSE, count = FALSE, 
                        used_methods =  c("combat_seq", "limma", "MMUPHin", 'ConQuR', 'ConQuR_libsize')) {
    

    
    print(output_root)
    # set up file save address
    sink_file_name = paste(output_root, "_runtime.txt", sep="")
    print(data_mat_path)
    output_dir = strsplit(output_root, "/")
    output_dir = paste(output_dir[[1]][1:(length(output_dir[[1]])-1)], collapse = "/")
    if(!dir.create(file.path(output_dir))){
        dir.create(file.path(output_dir))
    }

    # data loading <- load data_mat and meta_data, output of the preprocessing.py file
    if(transpose == TRUE) {
        count_data = read.table(data_mat_path, sep=",",header=T,check.names = F)
    }
    else {
        count_data = read.table(data_mat_path, sep=",",header=T,row.names=1,check.names = F)
    }
    print("check point1")

    metadata = as.data.frame(read_csv(meta_data_path))

    # count the number of samples and taxa to determine whether to use 10 cores or 2
    print(nrow(count_data))
    print(ncol(count_data))
    if(nrow(count_data)>300){
        num_core = 10
    }
    else {
        num_core = 2
    }

    print("check point2")



    # for ConQuR
    batchid <- factor(metadata[, dataset])
    
    
    ## CASE 1. count data (combat will not work)
    if(count == TRUE) {
        count_data.clr <- logratio.transfo(count_data+1, logratio = 'CLR')
        # if ConQuR_libsize output does not exist and its runtime documented File.readlines(paste(output_root, "_runtime.txt", sep="")).any?{ |l| l['ConQuR_libsize'] }
        if(!file.exists(paste(output_root, "_ConQuR_libsize.csv", sep=""))){
            run_count = TRUE
        }
        else if (file.exists(paste(output_root, "_runtime.txt", sep="")) & length(grep('ConQuR_libsize', readLines(paste(output_root, "_runtime.txt", sep=""))))==0){
            run_count = TRUE
        } 
        else {
            run_count = FALSE
        }

        if (run_count == TRUE){  
            cat("runtime documenting...\n", file=sink_file_name, append=FALSE)
            ### 1.1 combat_seq (in place of combat)
            if ('combat_seq' %in% used_methods) {
                print("in count mode, can only use combat seq")
                batch_info <- as.factor(setNames(as.character(metadata[, dataset]), metadata[[Sam_id]]))
                start_time <- Sys.time()
                ## check for covars
                if(is.null(covar)) {
                    count_data.combat_seq <- t(ComBat_seq(t(count_data), batch = batch_info))
                }
                else if(length(covar)==1){
                    covar_df = as.numeric(factor(metadata[, covar]))
                    count_data.combat_seq <- t(ComBat_seq(t(count_data), batch = batch_info, group = covar_df))
                }
                else{
                    covar_df = apply(metadata[, covar], 2, function(x) as.numeric(factor(x)))
                    count_data.combat_seq <- t(ComBat_seq(t(count_data), batch = batch_info, covar_mod = covar_df))
                }
                end_time <- Sys.time()
                cat(c("combat_seq runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.combat_seq, paste(output_root, "_combat_seq.csv", sep=""), row.names = TRUE)
            }

            ### 1.2 limma
            if ('limma' %in% used_methods) {
                batch_info <- as.factor(setNames(as.character(metadata[, dataset]), metadata[[Sam_id]]))
                start_time <- Sys.time()
                ## check for covars
                count_data.limma <-removeBatchEffect(t(count_data.clr), batch = batch_info)
                count_data.limma <- t(count_data.limma)
                if(is.null(covar)) {
                    count_data.limma <- t(removeBatchEffect(t(count_data.clr), batch = batch_info))
                }
                else if(length(covar)==1){
                    covar_df = matrix(as.numeric(factor(metadata[, covar])))
                    # design0 <- model.matrix(factor(metadata[, covar]))
                    count_data.limma <- t(removeBatchEffect(t(count_data.clr), batch = batch_info, design = covar_df))
                }
                else{
                    covar_df = apply(metadata[, covar], 2, function(x) as.numeric(factor(x)))
                    # design0 <- model.matrix(factor(metadata[, covar]))
                    count_data.limma <- t(removeBatchEffect(t(count_data.clr), batch = batch_info, design = covar_df))
                }
                ## the following line is count specific
                count_data.limma = clr(count_data.limma, inverse = TRUE)
                ## the above line is count specific
                end_time <- Sys.time()
                cat(c("limma runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.limma, paste(output_root, "_limma.csv", sep=""), row.names = TRUE)
            }

            ### 1.3 MMUPHin - count/relab agnostic
            if ('MMUPHin' %in% used_methods) {
                metadata_mupphin <- metadata
                feature_abd = count_data
                row.names(metadata_mupphin) <- metadata[[Sam_id]]
                feature_abd = as.data.frame(t(feature_abd))
                colnames(feature_abd) <- rownames(metadata_mupphin)
                start_time <- Sys.time()
                fit_adjust_batch <- adjust_batch(feature_abd = feature_abd,
                                        batch = dataset,
                                        covariates = covar,
                                        data = metadata_mupphin,
                                        control = list(verbose = FALSE))
                count_data.MMUPHin <- t(fit_adjust_batch$feature_abd_adj)                        
                end_time <- Sys.time()
                cat(c("MMUPHin runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.MMUPHin, paste(output_root, "_MMUPHin.csv", sep=""), row.names = TRUE)
            }
            
            ### 1.4 ConQuR (count version) - simple
            if ('ConQuR' %in% used_methods) {
                print("running ConQuR")
                start_time <- Sys.time()
                count_data.ConQuR = ConQuR(tax_tab=count_data, batchid=batchid, covariates=covar_df, batch_ref=batch_ref, num_core = num_core)                       
                end_time <- Sys.time()
                cat(c("ConQuR runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.ConQuR, paste(output_root, "_ConQuR.csv", sep=""), row.names = TRUE)
            }

            ### 1.5 ConQuR_libsize (only runnable for counts bc relab does not have libsize) - simple
            if ('ConQuR_libsize' %in% used_methods) {
                print("running ConQuR_libsize")
                start_time <- Sys.time()
                count_data.ConQuR_libsize = ConQuR_libsize(tax_tab=count_data, batchid=batchid, covariates=covar_df, batch_ref=batch_ref, num_core = num_core)                       
                end_time <- Sys.time()
                cat(c("ConQuR_libsize runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.ConQuR_libsize, paste(output_root, "_ConQuR_libsize.csv", sep=""), row.names = TRUE)
            }

            ### 1.6 Tune_ConQuR (count version)
            if ('Tune_ConQuR' %in% used_methods) {
                start_time <- Sys.time()
                count_data.Tune_ConQuR = Tune_ConQuR(tax_tab=count_data, batchid=batchid, covariates=covar_df,
                            batch_ref_pool=batch_ref,
                            logistic_lasso_pool=c(T, F),
                            quantile_type_pool=c("standard", "lasso", "composite"),
                            simple_match_pool=c(T, F),
                            lambda_quantile_pool=c(NA, "2p/n", "2p/logn"),
                            interplt_pool=c(T, F),
                            frequencyL=0,
                            frequencyU=1,
                            cutoff=0.25,
                            num_core = num_core)
                end_time <- Sys.time()
                cat(c("Tune_ConQuR runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.Tune_ConQuR$tax_final, paste(output_root, "_Tune_ConQuR.csv", sep=""), row.names = TRUE)
            }

            ### 1.7 Tune_ConQuR_libsize (only runnable for counts bc relab does not have libsize) - simple
            if ('Tune_ConQuR_libsize' %in% used_methods) {
                start_time <- Sys.time()
                count_data.Tune_ConQuR_libsize = Tune_ConQuR_libsize(tax_tab=count_data, batchid=batchid, covariates=covar_df,
                            batch_ref_pool=batch_ref,
                            logistic_lasso_pool=c(T, F),
                            quantile_type_pool=c("standard", "lasso", "composite"),
                            simple_match_pool=c(T, F),
                            lambda_quantile_pool=c(NA, "2p/n", "2p/logn"),
                            interplt_pool=c(T, F),
                            frequencyL=0,
                            frequencyU=1,
                            cutoff=0.25,
                            num_core = num_core)
                end_time <- Sys.time()
                cat(c("Tune_ConQuR_libsize runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.Tune_ConQuR_libsize$tax_final, paste(output_root, "_Tune_ConQuR_libsize.csv", sep=""), row.names = TRUE)
            }
        }
        else {
            print("already done")
        }
    }

    ## CASE 2. relative abundance data
    else {
        count_data.clr <- logratio.transfo(count_data+min(count_data[count_data>0]), logratio = 'CLR')
        print("check whether there are NAs in count_data.clr")
        print(sum(is.na(count_data.clr)))
        print("check whether there are Inf in count_data.clr")
        print(sum(is.infinite(count_data.clr)))
        if(!file.exists(paste(output_root, "_ConQuR_rel.csv", sep=""))){
            run_relab = TRUE
        }
        else if (file.exists(paste(output_root, "_runtime.txt", sep="")) & length(grep('ConQuR_rel', readLines(paste(output_root, "_runtime.txt", sep=""))))==0){
            run_relab = TRUE
        } 
        else {
            run_relab = FALSE
        }

        if (run_relab == TRUE){  
            cat("runtime documenting...\n", file=sink_file_name, append=FALSE)
            ### 1.1 combat
            if ('combat' %in% used_methods) {
                print("in relative abundance mode, can only use combat")
                batch_info <- as.factor(setNames(as.character(metadata[, dataset]), metadata[[Sam_id]]))
                start_time <- Sys.time()
                ## check for covars
                if(is.null(covar)) {
                    # count_data.combat <- t(ComBat(t(count_data*100), batch = batch_info, par.prior=FALSE, mod=NULL))
                    count_data.combat <- t(ComBat(t(count_data.clr), batch = batch_info, par.prior=FALSE, mod=NULL))
                }
                else if(length(covar)==1){
                    covar_df = as.numeric(factor(metadata[, covar]))
                    # count_data.combat <- t(ComBat(t(count_data*100), batch = batch_info, par.prior=FALSE, mod=covar_df))
                    count_data.combat <- t(ComBat(t(count_data.clr), batch = batch_info, par.prior=FALSE, mod=covar_df))
                }
                else{
                    print(colnames(metadata))
                    covar_df = apply(metadata[, covar], 2, function(x) as.numeric(factor(x)))
                    # count_data.combat <- t(ComBat(t(count_data*100), batch = batch_info, par.prior=FALSE, mod=covar_df))
                    count_data.combat <- t(ComBat(t(count_data.clr), batch = batch_info, par.prior=FALSE, mod=covar_df))
                }
                # now normalize it back to relative abundance
                # count_data.combat <- count_data.combat/rowSums(count_data.combat)
                end_time <- Sys.time()
                cat(c("combat runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.combat, paste(output_root, "_combat.csv", sep=""), row.names = TRUE)
            }

            ### 1.2 limma
            if ('limma' %in% used_methods) {
                batch_info <- as.factor(setNames(as.character(metadata[, dataset]), metadata[[Sam_id]]))
                start_time <- Sys.time()
                ## check for covars
                # print("AAAAAAA")
                # print("BBBBBB")
                if(is.null(covar)) {
                    count_data.limma <- t(removeBatchEffect(t(count_data.clr), batch = batch_info))
                }
                else if(length(covar)==1){
                    # covar_df = matrix(as.numeric(factor(metadata[, covar])))
                    # design0 <- model.matrix(factor(metadata[, covar]))
                    count_data.limma <- t(removeBatchEffect(t(count_data.clr), batch = batch_info))
                }
                else{
                    covar_df = apply(metadata[, covar], 2, function(x) as.numeric(factor(x)))
                    # design0 <- model.matrix(factor(metadata[, covar]))
                    count_data.limma <- t(removeBatchEffect(t(count_data.clr), batch = batch_info, design = covar_df))
                }
                end_time <- Sys.time()
                cat(c("limma runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.limma, paste(output_root, "_limma.csv", sep=""), row.names = TRUE)
            }

            ### 1.3 MMUPHin - count/relab agnostic
            if ('MMUPHin' %in% used_methods) {
                metadata_mupphin <- metadata
                feature_abd = count_data
                row.names(metadata_mupphin) <- metadata[[Sam_id]]
                feature_abd = as.data.frame(t(feature_abd))
                colnames(feature_abd) <- rownames(metadata_mupphin)
                start_time <- Sys.time()
                fit_adjust_batch <- adjust_batch(feature_abd = feature_abd,
                                        batch = dataset,
                                        covariates = covar,
                                        data = metadata_mupphin,
                                        control = list(verbose = FALSE))
                count_data.MMUPHin <- t(fit_adjust_batch$feature_abd_adj)                        
                end_time <- Sys.time()
                cat(c("MMUPHin runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.MMUPHin, paste(output_root, "_MMUPHin.csv", sep=""), row.names = TRUE)
            }
            
            ### 1.4 ConQuR_rel (rel_ab version) - simple
            if ('ConQuR_rel' %in% used_methods) {
                start_time <- Sys.time()
                count_data.ConQuR_rel = ConQuR_rel(tax_tab=count_data, batchid=batchid, covariates=covar_df, batch_ref=batch_ref, num_core = num_core)                       
                end_time <- Sys.time()
                cat(c("ConQuR_rel runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.ConQuR_rel, paste(output_root, "_ConQuR_rel.csv", sep=""), row.names = TRUE)
            }

            ### 1.5 Tune_ConQuR_rel (rel_ab version)
            if ('Tune_ConQuR_rel' %in% used_methods) {
                start_time <- Sys.time()
                count_data.Tune_ConQuR_rel = Tune_ConQuR_rel(tax_tab=count_data, batchid=batchid, covariates=covar_df,
                            batch_ref_pool=batch_ref,
                            logistic_lasso_pool=c(T, F),
                            quantile_type_pool=c("standard", "lasso", "composite"),
                            simple_match_pool=c(T, F),
                            lambda_quantile_pool=c(NA, "2p/n", "2p/logn"),
                            interplt_pool=c(T, F),
                            frequencyL=0,
                            frequencyU=1,
                            cutoff=0.25, num_core = num_core)
                end_time <- Sys.time()
                cat(c("Tune_ConQuR_rel runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.Tune_ConQuR_rel$tax_final, paste(output_root, "_Tune_ConQuR_rel.csv", sep=""), row.names = TRUE)
            }
        
        }
        else {
            print("already done")
        }
    }
    
}



if(option == 1){
    # autism 2 microbiomeHD
    current_root = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/cleaned_data/autism_2_microbiomeHD/autism_2_microbiomeHD'
    output_root = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD'
    run_methods(paste0(current_root, "_count_data.csv"),
        paste0(current_root, "_meta_data.csv"),
        output_root,
        dataset = "Dataset",
        batch_ref = 'asd_son',
        covar = c("DiseaseState"),
        count = TRUE,
        used_methods =  c("combat_seq", "limma", "MMUPHin", 'ConQuR', 'ConQuR_libsize')
    )
} else if(option == 2){
    print("benchmarking cdi")
    # cdi 3 microbiomeHD
    current_root = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/cleaned_data/cdi_3_microbiomeHD/cdi_3_microbiomeHD'
    output_root = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD'
    run_methods(paste0(current_root, "_count_data.csv"),
        paste0(current_root, "_meta_data.csv"),
        output_root,
        dataset = "Dataset",
        batch_ref = 'cdi_youngster',
        covar = c("DiseaseState"),
        count = TRUE,
        used_methods = c("combat_seq", "limma", "MMUPHin", 'ConQuR', 'ConQuR_libsize')
    )
} else if(option ==3){
    print("benchmarking ibd")
    # ibd 3 CMD
    current_root = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/cleaned_data/ibd_3_CMD/ibd_3_CMD'
    output_root = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark/benchmarked_results/ibd_3_CMD/ibd_3_CMD'
    run_methods(paste0(current_root, "_count_data.csv"),
        paste0(current_root, "_meta_data.csv"),
        output_root,
        dataset = "study_name",
        batch_ref = 'HMP_2019_ibdmdb',
        covar = c('disease', 'gender', 'age_category'),
        used_methods = c("combat", "limma", "MMUPHin", 'ConQuR_rel')
    )
} else {
    # crc_8_CMD
    current_root = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/cleaned_data/crc_8_CMD/crc_8_CMD'
    output_root = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark/benchmarked_results/crc_8_CMD/crc_8_CMD'
    run_methods(paste0(current_root, "_count_data.csv"),
        paste0(current_root, "_meta_data.csv"),
        output_root,
        dataset = "study_name",
        batch_ref = 'FengQ_2015',
        covar = c("disease", "gender", "age"),
        used_methods = c("combat", "limma", "MMUPHin", 'ConQuR_rel')
    )
}



# # ibd 3 CMD
# current_root = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/cleaned_data/ibd_3_CMD/ibd_3_CMD'
# output_root = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark/benchmarked_results/ibd_3_CMD/ibd_3_CMD'
# run_methods(paste0(current_root, "_count_data.csv"),
#     paste0(current_root, "_meta_data.csv"),
#     output_root,
#     dataset = "study_name",
#     batch_ref = 'HMP_2019_ibdmdb',
#     covar = c("disease", "gender", "age"),
#     used_methods = c("combat", "limma", "MMUPHin", 'ConQuR_rel')
# )



# # crc_8_CMD
# current_root = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/cleaned_data/crc_8_CMD/crc_8_CMD'
# output_root = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark/benchmarked_results/crc_8_CMD/crc_8_CMD'
# run_methods(paste0(current_root, "_count_data.csv"),
#     paste0(current_root, "_meta_data.csv"),
#     output_root,
#     dataset = "study_name",
#     batch_ref = 'FengQ_2015',
#     covar = c("disease", "gender", "age"),
#     used_methods = c("combat", "limma", "MMUPHin", 'ConQuR_rel')
# )

