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

## load ConQuR
conqur_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR'
source(paste0(conqur_path, "/ConQuR_help_functions.R"))
source(paste0(conqur_path, "/ConQuR_main_tune.R"))
# source(paste0(conqur_path, "/ConQuR_help_functions_libsize.R"))
# source(paste0(conqur_path, "/ConQuR_main_tune_libsize.R"))
source(paste0(conqur_path, "/ConQuR_help_functions_libsize_old.R"))
source(paste0(conqur_path, "/ConQuR_main_tune_libsize_old.R"))
source(paste0(conqur_path, "/ConQuR_help_functions_rel.R"))
source(paste0(conqur_path, "/ConQuR_main_tune_rel.R"))
source(paste0(conqur_path, "/supporting_functions.R"))

run_methods <- function(data_mat_path, meta_data_path, output_root, batch_ref, dataset = "Dataset", covar = NULL, controlled = FALSE, Sam_id = 'Sam_id', transpose = FALSE, count = FALSE, 
                        used_methods = c("combat_seq", "limma", "MMUPHin", 'conqur_libsize', "conqur")) {
    
    # set up file save address
    sink_file_name = paste(output_root, "_runtime.txt", sep="")
    print(data_mat_path)
    output_dir = strsplit(output_root, "/")
    output_dir = paste(output_dir[[1]][1:(length(output_dir[[1]])-1)], collapse = "/")
    dir.create(file.path(output_dir))

    # data loading <- load data_mat and meta_data, output of the preprocessing.py file
    if(transpose == TRUE) {
        count_data = read.table(data_mat_path, sep=",",header=T,check.names = F)
    }
    else {
        count_data = read.table(data_mat_path, sep=",",header=T,row.names=1,check.names = F)
    }
    print("check point1")

    metadata = as.data.frame(read_csv(meta_data_path))

    print("check point2")
    
    ## TODO:  potentially need preprocessing such as log and +1
    count_data.clr <- logratio.transfo(count_data+1, logratio = 'CLR')
    # count_data.clr <- mutate_all(count_data.clr, function(x) as.numeric(as.character(x)))
    cat("runtime documenting...\n", file=sink_file_name, append=FALSE)

    
    ## CASE 1. count data (combat will not work)
    if(count == TRUE) {

        ### 1.1 combat_seq (in place of combat)
        if ('combat' %in% used_methods) {
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
            cat(c("combat runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
            cat('\n', file=sink_file_name, append=TRUE)
            write.csv(count_data.combat_seq, paste(output_root, "_combat.csv", sep=""), row.names = TRUE)
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
            cat(c("limma runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
            cat('\n', file=sink_file_name, append=TRUE)
            write.csv(count_data.limma, paste(output_root, "_limma.csv", sep=""), row.names = TRUE)
        }

        ### 1.3 MMUPHin - count/relab agnostic
        if ('MMUPHin' %in% used_methods) {
            metadata_mupphin <- metadata
            row.names(metadata_mupphin) <- metadata[[Sam_id]]
            feature_abd = t(count_data)
            start_time <- Sys.time()
            fit_adjust_batch <- adjust_batch(feature_abd = feature_abd,
                                    batch = dataset,
                                    covariates = covar,
                                    data = metadata_mupphin,
                                    control = list(verbose = FALSE))
            count_data.MMUPHin <- fit_adjust_batch$feature_abd_adj                        
            end_time <- Sys.time()
            cat(c("MMUPHin runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
            cat('\n', file=sink_file_name, append=TRUE)
            write.csv(count_data.MMUPHin, paste(output_root, "_MMUPHin.csv", sep=""), row.names = TRUE)
        }
        
        ### 1.4 ConQuR (count version) - simple
        if ('ConQuR' %in% used_methods) {
            start_time <- Sys.time()
            batchid <- factor(metadata[, dataset])
            count_data.ConQuR = ConQuR(tax_tab=count_data, batchid=batchid, covariates=covar_df, batch_ref=batch_ref)                       
            end_time <- Sys.time()
            cat(c("ConQuR runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
            cat('\n', file=sink_file_name, append=TRUE)
            write.csv(count_data.ConQuR, paste(output_root, "_ConQuR.csv", sep=""), row.names = TRUE)
        }

        ### 1.5 ConQuR_libsize (only runnable for counts bc relab does not have libsize) - simple
        if ('ConQuR' %in% used_methods) {
            start_time <- Sys.time()
            count_data.ConQuR_libsize = ConQuR_libsize(tax_tab=count_data, batchid=batchid, covariates=covar_df, batch_ref=batch_ref)                       
            end_time <- Sys.time()
            cat(c("ConQuR_libsize runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
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
                         cutoff=0.25)
            end_time <- Sys.time()
            cat(c("Tune_ConQuR runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
            cat('\n', file=sink_file_name, append=TRUE)
            write.csv(count_data.Tune_ConQuR, paste(output_root, "_Tune_ConQuR.csv", sep=""), row.names = TRUE)
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
                         cutoff=0.25)
            end_time <- Sys.time()
            cat(c("Tune_ConQuR_libsize runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
            cat('\n', file=sink_file_name, append=TRUE)
            write.csv(count_data.Tune_ConQuR_libsize, paste(output_root, "_Tune_ConQuR_libsize.csv", sep=""), row.names = TRUE)
        }

    }

    ## CASE 2. relative abundance data
    else {

        ### 1.1 combat
        if ('combat' %in% used_methods) {
            print("in relative abundance mode, can only use combat")
            batch_info <- as.factor(setNames(as.character(metadata[, dataset]), metadata[[Sam_id]]))
            start_time <- Sys.time()
            ## check for covars
            if(is.null(covar)) {
                count_data.combat <- t(ComBat(t(count_data*100), batch = batch_info, par.prior=FALSE, mod=NULL))
            }
            else if(length(covar)==1){
                covar_df = as.numeric(factor(metadata[, covar]))
                count_data.combat <- t(ComBat(t(count_data*100), batch = batch_info, par.prior=FALSE, mod=covar_df))
            }
            else{
                covar_df = apply(metadata[, covar], 2, function(x) as.numeric(factor(x)))
                count_data.combat <- t(ComBat(t(count_data*100), batch = batch_info, par.prior=FALSE, mod=covar_df))
            }
            # now normalize it back to relative abundance
            count_data.combat <- count_data.combat/rowSums(count_data.combat)
            end_time <- Sys.time()
            cat(c("combat runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
            cat('\n', file=sink_file_name, append=TRUE)
            write.csv(count_data.combat, paste(output_root, "_combat.csv", sep=""), row.names = TRUE)
        }

        ### 1.2 limma
        if ('limma' %in% used_methods) {
            batch_info <- as.factor(setNames(as.character(metadata[, dataset]), metadata[[Sam_id]]))
            start_time <- Sys.time()
            ## check for covars
            # print("AAAAAAA")
            print(removeBatchEffect(t(count_data.clr), batch = batch_info))
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
            cat(c("limma runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
            cat('\n', file=sink_file_name, append=TRUE)
            write.csv(count_data.limma, paste(output_root, "_limma.csv", sep=""), row.names = TRUE)
        }

        ### 1.3 MMUPHin - count/relab agnostic
        if ('MMUPHin' %in% used_methods) {
            metadata_mupphin <- metadata
            row.names(metadata_mupphin) <- metadata[[Sam_id]]
            feature_abd = t(count_data)
            start_time <- Sys.time()
            fit_adjust_batch <- adjust_batch(feature_abd = feature_abd,
                                    batch = dataset,
                                    covariates = covar,
                                    data = metadata_mupphin,
                                    control = list(verbose = FALSE))
            count_data.MMUPHin <- t(fit_adjust_batch$feature_abd_adj)                        
            end_time <- Sys.time()
            cat(c("MMUPHin runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
            cat('\n', file=sink_file_name, append=TRUE)
            write.csv(count_data.MMUPHin, paste(output_root, "_MMUPHin.csv", sep=""), row.names = TRUE)
        }
        
        ### 1.4 ConQuR_rel (rel_ab version) - simple
        if ('ConQuR_rel' %in% used_methods) {
            start_time <- Sys.time()
            count_data.ConQuR_rel = ConQuR_rel(tax_tab=count_data, batchid=batchid, covariates=covar_df, batch_ref=batch_ref)                       
            end_time <- Sys.time()
            cat(c("ConQuR_simple runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
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
                         cutoff=0.25)
            end_time <- Sys.time()
            cat(c("Tune_ConQuR_rel runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
            cat('\n', file=sink_file_name, append=TRUE)
            write.csv(count_data.Tune_ConQuR_rel, paste(output_root, "_Tune_ConQuR_rel.csv", sep=""), row.names = TRUE)
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
                         cutoff=0.25)
            end_time <- Sys.time()
            cat(c("Tune_ConQuR_libsize runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
            cat('\n', file=sink_file_name, append=TRUE)
            write.csv(count_data.Tune_ConQuR_libsize, paste(output_root, "_Tune_ConQuR_libsize.csv", sep=""), row.names = TRUE)
        }

    }








    # ## combat/combat_seq
    # if ('combat' %in% used_methods) {
    # start_time <- Sys.time()
    # print("check point3")
    # batch_info <- as.factor(setNames(as.character(metadata[, dataset]), metadata[[Sam_id]]))
    # print("check point4")
    # if(is.null(covar)) {
    #     if(count == TRUE) {
    #         count_data.combat <- t(ComBat(t(count_data.clr), batch = batch_info, par.prior=FALSE, mod=NULL))
    #         }
    #     else {
    #         count_data.combat <- t(ComBat(t(count_data*100), batch = batch_info, par.prior=FALSE, mod=NULL))
    #         # now normalize it back to relative abundance
    #         count_data.combat <- count_data.combat/rowSums(count_data.combat)
    #         }
    #     }
    # else {
    #     if(length(covar)==1){
    #         covar_df = as.numeric(factor(metadata[, covar]))
    #     }
    #     else{
    #         # turn each column into a factor
    #         print("doing more than 1 covariate")
    #         covar_df = apply(metadata[, covar], 2, function(x) as.numeric(factor(x)))
    #     }
    #     # covar_df = as.numeric(factor(metadata[, covar]))
    #     # print(design)
    #     if(count == TRUE) {
    #         count_data.combat <- t(ComBat(t(count_data.clr), batch = batch_info, mod = covar_df, par.prior=FALSE, mean.only=TRUE))
    #         }
    #     else {
    #         count_data.combat <- t(ComBat(t(count_data*100), batch = batch_info, mod = covar_df, par.prior=FALSE, mean.only=TRUE))
    #         # now normalize it back to relative abundance
    #         count_data.combat <- count_data.combat/rowSums(count_data.combat)
    #         }
    #     }
    # write.csv(count_data.combat, paste(output_root, "_combat.csv", sep=""), row.names = TRUE)
    # end_time <- Sys.time()
    # cat(c("combat runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
    # cat('\n', file=sink_file_name, append=TRUE)
    # }

    # ## combat_seq
    # ### CAN ONLY WORK IF count == TRUE
    # if ('combat_seq' %in% used_methods) {
    # start_time <- Sys.time()
    # batch_info <- as.factor(setNames(as.character(metadata[, dataset]), metadata[[Sam_id]]))
    # if(is.null(covar)) {
    #     if(count == TRUE) {
    #         count_data.combat_seq <- t(ComBat_seq(t(count_data), batch = batch_info))
            
    #         write.csv(count_data.combat_seq, paste(output_root, "_combat_seq.csv", sep=""), row.names = TRUE)
    #         end_time <- Sys.time()
    #         cat(c("combat_seq runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
    #         cat('\n', file=sink_file_name, append=TRUE)
    #     }
    #     else {
    #         print("count must be TRUE for combat_seq")
    #     }
    # }
    # else {
    #     if(length(covar)==1){
    #         covar_df = as.numeric(factor(metadata[, covar]))
    #     }
    #     else{
    #         print("doing more than 1 covariate")
    #         # turn each column into a factor
    #         covar_df = apply(metadata[, covar], 2, function(x) as.numeric(factor(x)))
    #     }
    #     if(count == TRUE) {
    #         count_data.combat_seq <- t(ComBat_seq(t(count_data), batch = batch_info, covar_mod = covar_df))
    #         write.csv(count_data.combat_seq, paste(output_root, "_combat_seq.csv", sep=""), row.names = TRUE)
    #         end_time <- Sys.time()
    #         cat(c("combat_seq runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
    #         cat('\n', file=sink_file_name, append=TRUE)
    #         }
    #     else {
    #      print("count must be TRUE for combat_seq")
    #         }
    #     }
    # }
    # ## limma
    # ## note that the 'covariates' argument here cannot be used because limma requires continuous covariates
    # if ('limma' %in% used_methods) {
    # start_time <- Sys.time()
    # count_data.limma <- t(removeBatchEffect(t(count_data.clr), batch = batch_info))
    # if(count == TRUE) {
    #     count_data.limma = clr(count_data.limma, inverse = TRUE)
    # }
    # write.csv(count_data.limma,paste(output_root, "_limma.csv", sep=""), row.names = TRUE)
    # end_time <- Sys.time()
    # cat(c("limma runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
    # cat('\n', file=sink_file_name, append=TRUE)
    # }

    # ## MMUPHin
    # if ('mmuphin' %in% used_methods) {
    # start_time <- Sys.time()
    # metadata_mupphin <- metadata
    # row.names(metadata_mupphin) <- metadata[[Sam_id]]
    # feature_abd = t(count_data)
    # # print(t(count_data_t_relab))
    # colnames(feature_abd) = row.names(metadata_mupphin)
    # if(is.null(covar)) {
    #     fit_adjust_batch <- adjust_batch(feature_abd = feature_abd,
    #                                 batch = dataset,
    #                                 data = metadata_mupphin,
    #                                 control = list(verbose = FALSE))
    # }
    # else {
    #     fit_adjust_batch <- adjust_batch(feature_abd = feature_abd,
    #                                 batch = dataset,
    #                                 covariates = covar,
    #                                 data = metadata_mupphin,
    #                                 control = list(verbose = FALSE))
    # }
    # count_data.MMUPHin <- fit_adjust_batch$feature_abd_adj
    # write.csv(t(count_data.MMUPHin),paste(output_root, "_MMUPHin.csv", sep=""), row.names = TRUE)
    # end_time <- Sys.time()
    # cat(c("MMUPHin runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
    # cat('\n', file=sink_file_name, append=TRUE)
    # }

    # ## ConquR
    # if ('conqur' %in% used_methods) {
    # start_time <- Sys.time()
    # batchid <- factor(metadata[, dataset])
    # if (is.null(covar)){
    #     count_data.conqur = ConQuR(tax_tab=count_data, batchid=batchid, covariates=NULL, simple_match = T, batch_ref = batch_ref,
    #                      logistic_lasso=T, quantile_type="lasso", interplt=F)
    #     # count_data.conqur = Tune_ConQuR(tax_tab=count_data, batchid=batchid, covariates=NULL,
    #     #                    batch_ref_pool=batch_ref,  
    #     #                    logistic_lasso_pool=c(T, F),
    #     #                    quantile_type_pool=c("standard", "lasso", "composite"),
    #     #                    simple_match_pool=c(T, F),
    #     #                    lambda_quantile_pool=c(NA, "2p/n", "2p/logn"),
    #     #                    interplt_pool=c(T, F),
    #     #                    frequencyL=0,
    #     #                    frequencyU=1,
    #     #                    cutoff=0.25)
    # }
    # else {
        
    #     # covar_df = data.frame(factor(metadata[, covar]))
    #     if(length(covar)==1){
    #         # covar = covar[1]
    #         covar_df = as.numeric(factor(metadata[, covar]))
    #     }
    #     else{
    #         print("ConQuR: doing more than 1 covariate")
    #         # turn each column into a factor
    #         covar_df = apply(metadata[, covar], 2, function(x) as.numeric(factor(x)))
    #     }
    #     count_data.conqur = ConQuR(tax_tab=count_data, batchid=batchid, covariates=covar_df, batch_ref = batch_ref,
    #                      logistic_lasso=T, quantile_type="lasso", interplt=F)
    #     # count_data.conqur = Tune_ConQuR(tax_tab=count_data, batchid=batchid, covariates=covar_df,
    #     #                    batch_ref_pool=batch_ref,  
    #     #                    logistic_lasso_pool=c(T, F),
    #     #                    quantile_type_pool=c("standard", "lasso", "composite"),
    #     #                    simple_match_pool=c(T, F),
    #     #                    lambda_quantile_pool=c(NA, "2p/n", "2p/logn"),
    #     #                    interplt_pool=c(T, F),
    #     #                    frequencyL=0,
    #     #                    frequencyU=1,
    #     #                    cutoff=0.25)
    # }
    # write.csv(count_data.conqur,paste(output_root, "_ConQuR.csv", sep=""), row.names = TRUE)
    # end_time <- Sys.time()
    # cat(c("ConquR runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
    # cat('\n', file=sink_file_name, append=TRUE)
    # }

    # ## ConQuR_libsize
    # ## need to put batchid and covar into countdata dataframe
    # if ('conqur_libsize' %in% used_methods) {
    # start_time <- Sys.time()
    # batchid <- factor(metadata[, dataset])
    # if (is.null(covar)){
    #     if(count == TRUE){
    #         count_data.conqur_libsize = ConQuR_libsize(tax_tab=count_data, batchid=batchid, covariates=NULL, simple_match = T, batch_ref = batch_ref,
    #                      logistic_lasso=T, quantile_type="lasso", interplt=F)
    #         # count_data.conqur = Tune_ConQuR_libsize(tax_tab=count_data, batchid=batchid, covariates=NULL,
    #         #                batch_ref_pool=batch_ref,
    #         #                logistic_lasso_pool=c(T, F),
    #         #                quantile_type_pool=c("standard", "lasso", "composite"),
    #         #                simple_match_pool=c(T, F),
    #         #                lambda_quantile_pool=c(NA, "2p/n", "2p/logn"),
    #         #                interplt_pool=c(T, F),
    #         #                frequencyL=0,
    #         #                frequencyU=1,
    #         #                cutoff=0.25)           
    #         write.csv(count_data.conqur_libsize,paste(output_root, "_ConQuR_libsize.csv", sep=""), row.names = TRUE)
    #         end_time <- Sys.time()
    #         cat(c("ConquR_libsize runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
    #         }
    #     else{
    #         print("Due to unknown libsize, conqur_libsize cannot be run")
    #         } 
    #     }
    # else {
    #     if(length(covar)==1){
    #         covar_df = as.numeric(factor(metadata[, covar]))
    #     }
    #     else{
    #         print("ConQuR_libsize: doing more than 1 covariate")
    #         # turn each column into a factor
    #         covar_df = apply(metadata[, covar], 2, function(x) as.numeric(factor(x)))
    #     }
    #     # covar_df = data.frame(factor(metadata[, covar]))
    #     if(count == TRUE){
    #         count_data.conqur_libsize = ConQuR_libsize(tax_tab=count_data, batchid=batchid, covariates=covar_df, batch_ref = batch_ref,
    #                         logistic_lasso=T, quantile_type="lasso", interplt=F, num_core=5)
    #         # count_data.conqur = Tune_ConQuR_libsize(tax_tab=count_data, batchid=batchid, covariates=covar_df,
    #         #                batch_ref_pool=batch_ref,  
    #         #                logistic_lasso_pool=c(T, F),
    #         #                quantile_type_pool=c("standard", "lasso", "composite"),
    #         #                simple_match_pool=c(T, F),
    #         #                lambda_quantile_pool=c(NA, "2p/n", "2p/logn"),
    #         #                interplt_pool=c(T, F),
    #         #                frequencyL=0,
    #         #                frequencyU=1,
    #         #                cutoff=0.25) 
    #         write.csv(count_data.conqur_libsize,paste(output_root, "_ConQuR_libsize.csv", sep=""), row.names = TRUE)
    #         end_time <- Sys.time()
    #         cat(c("ConquR_libsize runtime", toString(end_time - start_time), "seconds"), file=sink_file_name, append=TRUE)
    #     }
    #     else{
    #         print("Due to unknown libsize, conqur_libsize cannot be run")
    #     }
    # }
    # }

    
}

overall_path = '/athena/linglab/scratch/chf4012/simulation_data_MIDAS_small_072623'

# autism 2 microbiomeHD
run_methods('/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark/benchmarked_data/autism_2_microbiomeHD_count_data.csv',
'/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark/benchmarked_data/autism_2_microbiomeHD_meta_data.csv',
'/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD',
dataset = "Dataset",
covar = c("DiseaseState"),
count = TRUE,
batch_ref = 'asd_son',
used_methods = c("combat", "limma", "MMUPHin", 'ConQuR', 'ConQuR_libsize', 'Tune_ConQuR', 'Tune_ConQuR_libsize')
)

# cdi 3 microbiomeHD
run_methods('/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark/benchmarked_data/cdi_3_microbiomeHD_count_data.csv',
'/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark/benchmarked_data/cdi_3_microbiomeHD_meta_data.csv',
'/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD',
dataset = "Dataset",
covar = c("DiseaseState"),
count = TRUE,
batch_ref = 'cdi_youngster',
used_methods = c("combat", "limma", "MMUPHin", 'ConQuR', 'ConQuR_libsize', 'Tune_ConQuR', 'Tune_ConQuR_libsize')
)

# # ibd 3 CMD
# run_methods('/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/ibd_3_CMD_count_data.csv',
# '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/ibd_3_CMD_meta_data.csv',
# '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_results/ibd_3_CMD/ibd_3_CMD',
# dataset = "study_name",
# batch_ref = 'HMP_2019_ibdmdb',
# covar = c("disease", "gender", "age"),
# used_methods = c("combat", "limma", "MMUPHin", 'ConQuR', 'ConQuR_libsize', 'Tune_ConQuR', 'Tune_ConQuR_libsize')
# )


# # CRC_8_CMD
# run_methods('/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/CRC_8_CMD_count_data.csv',
# '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/CRC_8_CMD_meta_data.csv',
# '/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_results/CRC_8_CMD/CRC_8_CMD',
# dataset = "study_name",
# batch_ref = 'FengQ_2015',
# covar = c("disease", "gender", "age"),
# used_methods = c("combat", "limma", "MMUPHin", 'ConQuR', 'ConQuR_libsize', 'Tune_ConQuR', 'Tune_ConQuR_libsize')
# )



or_l = c(1, 1.25, 1.5)
cond_effect_val_l = c(0, 0.099, 0.299, 0.499, 0.699, 0.899)
batch_effect_val_l = c(0, 0.099, 0.299, 0.499, 0.699, 0.899)
output_dir = '/athena/linglab/scratch/chf4012/simulation_data_output_small_072623'
scaled_midas_methods_bencharking <- function(method_l, or_l, cond_effect_val_l, batch_effect_val_l, num_iter){   
  for (or in or_l) {
    for (cond_effect_val in cond_effect_val_l) {
      for (batch_effect_val in batch_effect_val_l) {
        if (cond_effect_val + batch_effect_val <= 1) {
          for (iter in seq(1, num_iter)){
            print(or)
            print(cond_effect_val)
            print(batch_effect_val)
            
            output_file_path_count = paste0(overall_path, "/ibd_150_count_", or, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter, ".csv")
            # run the methods on this
            run_methods(output_file_path_count,
                        paste0(overall_path, "/ibd_150_meta_", or, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter, ".csv"),
                        paste0(output_dir, "/ibd_150_", or, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter),
                        dataset = 'batchid',
                        batch_ref = "batch_0",
                        covar = c("cond"),
                        Sam_id = 'subjectid_text',
                        transpose = TRUE,
                        count = TRUE,
                        used_methods = method_l,
            )
          }
        }
      }
    }
  }
}
method_l = c("combat", "limma", "MMUPHin", 'ConQuR', 'ConQuR_libsize', 'Tune_ConQuR', 'Tune_ConQuR_libsize')
scaled_midas_methods_bencharking(method_l, or_l, cond_effect_val_l, batch_effect_val_l, 5)
