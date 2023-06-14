library(xtable) # table
library(ConQuR)
library(doParallel) 
library(dplyr)
library(readr)
library(tibble)


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
                           simple_match_pool=T,
                           lambda_quantile_pool=c(NA, "2p/n"),
                           interplt_pool=F,
                           frequencyL=0,
                           frequencyU=1)
    }
    write.csv(count_data.conqur_finetuned,paste(output_root, "_ConQuR_finetuned.csv", sep=""), row.names = TRUE)
    end_time <- Sys.time()
    cat(c("ConquR_finetuned runtime", toString(end_time - start_time), "seconds"))
    
}



run_methods('/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/hanninganGD_noBoston_count_data.csv',
'/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_data/hanninganGD_noBoston_meta_data.csv',
'/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/benchmarked_results/hanninganGD_noBoston_gender/hanninganGD_noBoston_gender',
dataset = 'location',
batch_ref = 'Toronto',
Sam_id = 'patient_visit_id',
covar = c("gender")
)