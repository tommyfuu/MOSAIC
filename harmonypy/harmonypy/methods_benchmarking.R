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

run_methods <- function(data_mat_path, meta_data_path, output_root, covar = NULL) {
    # data loading <- load data_mat and meta_data, output of the preprocessing.py file
    count_data = read.table(data_mat_path, sep=",",header=T,row.names=1,check.names = F)
    metadata = read_csv(meta_data_path)
    count_data <- count_data %>%
        column_to_rownames('...1')


    ## TODO:  potentially need preprocessing such as log and +1
    # count_data <- count_data + 1
    count_data.clr <- logratio.transfo(count_data+1, logratio = 'clr')

    ## combat
    batch_info <- as.factor(setNames(as.character(metadata$batch), metadata$Sam_id))
    count_data.combat <- t(ComBat(t(count_data.clr), batch = batch_info, par.prior=FALSE, mean.only=TRUE))
    write.csv(count_data.combat, paste(output_root, "_combat.csv", sep=""), row.names = TRUE)
    
    ## limma
    count_data.limma <- t(removeBatchEffect(t(count_data.clr), batch = batch_info))
    write.csv(count_data.limma,paste(output_root, "_limma.csv", sep=""), row.names = TRUE)

    ## ConquR
    ## need to put batchid and covar into countdata dataframe
    count_data.conqur_pre <- count_data
    batchid <- meta_data$Dataset
    if (covar == NULL) {
        a <- rep(list("1"), dim(taxa)[1])
        covar_df <- factor(a)
    }
    else {
        covar_df = Sample_Data[, covar]
    }
    count_data.conqur = ConQuR(tax_tab=taxa, batchid=batchid, covariates=covar_df, batch_ref="0",
                         logistic_lasso=T, quantile_type="lasso", interplt=T)
    write.csv(count_data.conqur,paste(output_root, "_ConQuR.csv", sep=""), row.names = TRUE)
    
    ## MMUPHin
    count_data_t_relab = t(t(count_data)/rowSums(t(count_data)))
    metadata_mupphin <- metadata %>% column_to_rownames('Sam_id')
    if(covar == NULL) {
        fit_adjust_batch <- adjust_batch(feature_abd = t(count_data_t_relab ),
                                    batch = "Dataset",
                                    data = metadata_mupphin,
                                    control = list(verbose = FALSE))
    }
    else {
        fit_adjust_batch <- adjust_batch(feature_abd = t(count_data_t_relab ),
                                    batch = "Dataset",
                                    covariates = "Sex",
                                    data = metadata_mupphin,
                                    control = list(verbose = FALSE))
    }
    count_data.MMUPHin <- fit_adjust_batch$feature_abd_adj
    write.csv(count_data.MMUPHin,paste(output_root, "_MMUPHin.csv", sep=""), row.names = TRUE)
}

run_methods('/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/Glickman_count_data.csv',
'/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_data/Glickman_meta_data.csv',
'/home/fuc/harmonicMic/harmonypy/harmonypy/benchmarked_results/Glickman/Glickman',
'Sex'
)
## combat

## limma

## batch mean centering

## ConquR

## MMUPHin