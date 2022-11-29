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

run_methods <- function(data_mat_path, meta_data_path, output_root, batch_ref, dataset = "Dataset", covar = NULL) {
    # data loading <- load data_mat and meta_data, output of the preprocessing.py file
    count_data = read.table(data_mat_path, sep=",",header=T,row.names=1,check.names = F)
    metadata = as.data.frame(read_csv(meta_data_path))
    sink(paste(output_root, "_runtime.txt", sep=""))


    ## TODO:  potentially need preprocessing such as log and +1
    # count_data <- count_data + 1
    count_data.clr <- logratio.transfo(count_data+1, logratio = 'CLR')

    
    ## combat
    start_time <- Sys.time()
    batch_info <- as.factor(setNames(as.character(metadata[, dataset]), metadata$Sam_id))
    if(is.null(covar)) {
        count_data.combat <- t(ComBat(t(count_data.clr), batch = batch_info, par.prior=FALSE, mean.only=TRUE))
    }
    else {
        covar_df = factor(metadata[, covar])
        count_data.combat <- t(ComBat(t(count_data.clr), batch = batch_info, mod = covar_df, par.prior=FALSE, mean.only=TRUE))
    }
    write.csv(count_data.combat, paste(output_root, "_combat.csv", sep=""), row.names = TRUE)
    end_time <- Sys.time()
    cat(c("combat runtime", toString(end_time - start_time), "seconds"))
    cat('\n')

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
    metadata_mupphin <- metadata
    row.names(metadata_mupphin) <- metadata$Sam_id
    if(is.null(covar)) {
        fit_adjust_batch <- adjust_batch(feature_abd = t(count_data_t_relab ),
                                    batch = dataset,
                                    data = metadata_mupphin,
                                    control = list(verbose = FALSE))
    }
    else {
        fit_adjust_batch <- adjust_batch(feature_abd = t(count_data_t_relab ),
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
        # a <- rep(list("1"), dim(count_data)[1])
        # a <- seq(1, dim(count_data)[1])
        # covar_df <- factor(a)
        covar_df <- batchid 
    }
    else {
        covar_df = factor(metadata[, covar])
    }
    # print("ASACAVVC")
    # print(covar_df)
    count_data.conqur = ConQuR(tax_tab=count_data, batchid=batchid, covariates=covar_df, batch_ref = batch_ref,
                         logistic_lasso=T, quantile_type="lasso", interplt=T)
    write.csv(count_data.conqur,paste(output_root, "_ConQuR.csv", sep=""), row.names = TRUE)
    end_time <- Sys.time()
    cat(c("ConquR runtime", toString(end_time - start_time), "seconds"))
    cat('\n')
    
}

# glickman
# run_methods('/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/Glickman_count_data.csv',
# '/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/Glickman_meta_data.csv',
# '/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/Glickman/Glickman',
# dataset = "Dataset",
# batch_ref = 'Old',
# covar='Sex'
# )


# autism 2 microbiomeHD
run_methods('/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/autism_2_microbiomeHD_count_data.csv',
'/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/autism_2_microbiomeHD_meta_data.csv',
'/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/autism_2_microbiomeHD/autism_2_microbiomeHD',
dataset = "Dataset",
batch_ref = 'asd_son')

# cdi 3 microbiomeHD
run_methods('/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/cdi_3_microbiomeHD_count_data.csv',
'/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/cdi_3_microbiomeHD_meta_data.csv',
'/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/cdi_3_microbiomeHD/cdi_3_microbiomeHD',
dataset = "Dataset",
batch_ref = 'cdi_schubert')

# ibd 3 CMD
run_methods('/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/ibd_3_CMD_count_data.csv',
'/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/ibd_3_CMD_meta_data.csv',
'/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/ibd_3_CMD/ibd_3_CMD',
dataset = "study_name",
batch_ref = 'HMP_2019_ibdmdb')