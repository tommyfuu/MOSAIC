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
library(mixOmics)

args = commandArgs(trailingOnly=TRUE)

print(args)
if (length(args)==0 || length(args)>1 ) {
    GLOBAL_ITER = 1
    print("default running iteration 1")
} else if (length(args)==1) {
    # default output file
    GLOBAL_ITER = args[1]
}

## load ConQuR
conqur_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR'
source(paste0(conqur_path, "/ConQuR_help_functions.R"))
source(paste0(conqur_path, "/ConQuR_main_tune.R"))
source(paste0(conqur_path, "/ConQuR_help_functions_libsize_old.R"))
source(paste0(conqur_path, "/ConQuR_main_tune_libsize_old.R"))
source(paste0(conqur_path, "/ConQuR_help_functions_rel.R"))
source(paste0(conqur_path, "/ConQuR_main_tune_rel.R"))
source(paste0(conqur_path, "/supporting_functions.R"))


run_methods <- function(data_mat_path, meta_data_path, output_root, batch_ref, dataset = "Dataset", covar = NULL, controlled = FALSE, Sam_id = 'Sam_id', transpose = FALSE, count = FALSE, 
                        used_methods = c("combat_seq", "limma", "MMUPHin", 'conqur_libsize', "conqur")) {
    

    
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
        print(">>>>>")
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
                start_time <- Sys.time()
                count_data.ConQuR = ConQuR(tax_tab=count_data, batchid=batchid, covariates=covar_df, batch_ref=batch_ref, num_core = num_core)                       
                end_time <- Sys.time()
                cat(c("ConQuR runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.ConQuR, paste(output_root, "_ConQuR.csv", sep=""), row.names = TRUE)
            }

            ### 1.5 ConQuR_libsize (only runnable for counts bc relab does not have libsize) - simple
            if ('ConQuR' %in% used_methods) {
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


# just to speed up
# An mc-version of the sapply function.
library(parallel)
mcsapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
  FUN <- match.fun(FUN)
  answer <- parallel::mclapply(X = X, FUN = FUN, ...)
  if (USE.NAMES && is.character(X) && is.null(names(answer))) 
    names(answer) <- X
  if (!isFALSE(simplify) && length(answer)) 
    simplify2array(answer, higher = (simplify == "array"))
  else answer
}

run_methods_per_iter <- function(iter, overall_path, output_dir, or, cond_effect_val, batch_effect_val, used_methods, batch_ref = "batch_0", dataset = 'batchid', covar = c("cond"), Sam_id = 'Sam_id', transpose = TRUE, count = TRUE) {
                            if(count == TRUE){
                                output_file_path_count = paste0(overall_path, "/ibd_150_count_", or, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter, ".csv")
                            }
                            else{
                                output_file_path_count = paste0(overall_path, "/ibd_150_relab_", or, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter, ".csv")
                            }
                            # output_file_path_count = paste0(overall_path, "/ibd_150_count_", or, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter, ".csv")
                            run_methods(output_file_path_count,
                                paste0(overall_path, "/ibd_150_meta_", or, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter, ".csv"),
                                paste0(output_dir, "/out_", or, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter, "/ibd_", or, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter),
                                dataset = dataset,
                                batch_ref = batch_ref,
                                covar = covar,
                                Sam_id = Sam_id,
                                transpose = transpose,
                                count = count,
                                used_methods = used_methods,
                            )
                        }
    

scaled_midas_methods_bencharking <- function(output_dir, overall_path, method_l, or_l, cond_effect_val_l, batch_effect_val_l, num_iter, count = TRUE){   
  for (or in or_l) {
    for (cond_effect_val in cond_effect_val_l) {
      for (batch_effect_val in batch_effect_val_l) {
        if (cond_effect_val + batch_effect_val <= 1) {
            print(output_dir)
            print(or)
            print(cond_effect_val)
            print(batch_effect_val)
            print(iter)
            mcsapply(seq(1, num_iter), function(iter) run_methods_per_iter(iter, overall_path, output_dir, or, cond_effect_val, batch_effect_val, used_methods = method_l, count = count),
                mc.cores=5)
            print("WHAT'S HAPPENING??")
          
        }
      }
    }
  }
}


scaled_slurm_methods_bencharking <- function(output_dir, overall_path, method_l, or_l, cond_effect_val_l, batch_effect_val_l, iter, count = TRUE){   
  for (or in or_l) {
    for (cond_effect_val in cond_effect_val_l) {
      for (batch_effect_val in batch_effect_val_l) {
        if (cond_effect_val + batch_effect_val <= 1) {
            print(output_dir)
            print(or)
            print(cond_effect_val)
            print(batch_effect_val)
            print(iter)
            if(count == TRUE){
                output_file_path_count = paste0(overall_path, "/ibd_150_count_", or, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter, ".csv")
            }
            else{
                output_file_path_count = paste0(overall_path, "/ibd_150_relab_", or, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter, ".csv")
            }
            # output_file_path_count = paste0(overall_path, "/ibd_150_count_", or, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter, ".csv")
            run_methods(output_file_path_count,
                paste0(overall_path, "/ibd_150_meta_", or, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter, ".csv"),
                paste0(output_dir, "/out_", or, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter, "/ibd_", or, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter),
                dataset = "batchid",
                batch_ref = "batch_0", 
                covar = c("cond"),
                Sam_id = 'Sam_id',
                transpose = TRUE,
                count = count,
                used_methods = method_l,
            )
            print("WHAT'S HAPPENING??")
          
        }
      }
    }
  }
}

# # uncomment these codes for scaled slurm runs
# or_l = c(1, 1.25, 1.5)
# cond_effect_val_l = c(0, 0.25, 0.5, 0.75, 1)
# batch_effect_val_l = c(0, 0.25, 0.5, 0.75, 1)

# print("count no relation")
# overall_path = '/athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_MIDAS_1000_norelation_102023'
# output_dir = '/athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_output_count_norelation_102023'
# method_l = c("combat_seq", "limma", "MMUPHin", 'ConQuR', 'ConQuR_libsize')
# scaled_slurm_methods_bencharking(output_dir, overall_path, method_l, or_l, cond_effect_val_l, batch_effect_val_l, GLOBAL_ITER, count = TRUE)

# print("relab no relation")
# overall_path = '/athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_MIDAS_1000_norelation_102023'
# output_dir = '/athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_output_relab_norelation_102023'
# method_l = c("combat", "limma", "MMUPHin", 'ConQuR_rel')
# scaled_slurm_methods_bencharking(output_dir, overall_path, method_l, or_l, cond_effect_val_l, batch_effect_val_l, GLOBAL_ITER, count = FALSE)

# print("count yes relation")
# overall_path = '/athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_MIDAS_1000_yesrelation_102023'
# output_dir = '/athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_output_count_yesrelation_102023'
# method_l = c("combat_seq", "limma", "MMUPHin", 'ConQuR', 'ConQuR_libsize')
# scaled_slurm_methods_bencharking(output_dir, overall_path, method_l, or_l, cond_effect_val_l, batch_effect_val_l, GLOBAL_ITER, count = TRUE)

# print("relab yes relation")
# overall_path = '/athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_MIDAS_1000_yesrelation_102023'
# output_dir = '/athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_output_relab_yesrelation_102023'
# method_l = c("combat", "limma", "MMUPHin", 'ConQuR_rel')
# scaled_slurm_methods_bencharking(output_dir, overall_path, method_l, or_l, cond_effect_val_l, batch_effect_val_l, GLOBAL_ITER, count = FALSE)

