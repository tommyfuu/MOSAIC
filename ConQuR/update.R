
overall_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR'
source(paste0(overall_path, "/ConQuR_help_functions.R"))
source(paste0(overall_path, /ConQuR_main_tune.R"))
source(paste0(overall_path, /ConQuR_help_functions_libsize.R"))
source(paste0(overall_path, /ConQuR_main_tune_libsize.R"))
# source(paste0(overall_path, /ConQuR_help_functions_libsize_old.R"))
# source(paste0(overall_path, /ConQuR_main_tune_libsize_old.R")
source(paste0(overall_path, /ConQuR_help_functions_rel.R"))
source(paste0(overall_path, /ConQuR_main_tune_rel.R")



###### relative abundance update ######

ibd_taxa = read.csv("/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark/benchmarked_data/ibd_3_CMD_count_data.csv")
ibd_meta = read.csv("/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark/benchmarked_data/ibd_3_CMD_meta_data.csv")

# length(unique(ibd_meta$subject_id))
# identical(ibd_meta$Sam_id, ibd_taxa$X)

rownames(ibd_taxa) = ibd_taxa$X
ibd_taxa = ibd_taxa[, -1]

apply(ibd_taxa, 1, sum)
summary( apply(ibd_taxa, 2, function(z){mean(z == 0)}) ) 

table(ibd_meta$study_name)

### the key point is to make sure the data consists of complete cases only ###
### in the future, may update the code to do sanity check and manipulation ###

id_to_delete = which(ibd_meta$gender == "" | is.na(ibd_meta$age))
ibd_meta = ibd_meta[-id_to_delete, ]
ibd_taxa = ibd_taxa[-id_to_delete, ]

batchid = factor(ibd_meta$study_name)
covar = data.frame(factor(ibd_meta$disease), factor(ibd_meta$gender), ibd_meta$age)

# ConQuR_rel
fit1 = ConQuR_rel(tax_tab=ibd_taxa, batchid=batchid, covariates=covar,
                  batch_ref="HMP_2019_ibdmdb")
fit2 = ConQuR_rel(tax_tab=ibd_taxa, batchid=batchid, covariates=covar,
                  batch_ref="HMP_2019_ibdmdb",
                  logistic_lasso=T, quantile_type="lasso")

output_root = 'idb_3_CMD'
write.csv(fit1, paste(output_root, "_relab_basic.csv", sep=""), row.names = TRUE)
write.csv(fit2, paste(output_root, "_relab_lasso.csv", sep=""), row.names = TRUE)

# ConQuR_rel tuned
## Tom: I didn't run the following fine-tuned version as it is very time-consuming. Please help run it on SCU after your account is well set-up ###
taxa_tuned_rel = Tune_ConQuR_rel(tax_tab=ibd_taxa, batchid=batchid, covariates=covar,
                                 batch_ref_pool="HMP_2019_ibdmdb",
                                 logistic_lasso_pool=c(T, F),
                                 quantile_type_pool=c("standard", "lasso", "composite"),
                                 simple_match_pool=c(T, F),
                                 lambda_quantile_pool=c(NA, "2p/n", "2p/logn"),
                                 interplt_pool=c(T, F),
                                 frequencyL=0,
                                 frequencyU=1,
                                 cutoff=0.25) 
# save(taxa_tuned_rel, file="relative_tuned_results.Rdata")
write.csv(taxa_tuned_rel, paste(output_root, "_relab_tuned.csv", sep=""), row.names = TRUE)

load(file="relative_results.Rdata")
# load(file="relative_tuned_results.Rdata")

pdf("ibd_relative.pdf", width=12, height=8)
par(mfrow=c(2, 3))
Plot_PCoA(TAX=ibd_taxa, factor=batchid)
# Plot_PCoA(TAX=taxa_tuned_rel$tax_final, factor=batchid)
Plot_PCoA(TAX=fit1, factor=batchid)
Plot_PCoA(TAX=fit2, factor=batchid)

Plot_PCoA(TAX=ibd_taxa, factor=batchid, dissimilarity="Aitch")
# Plot_PCoA(TAX=taxa_tuned_rel$tax_final, factor=batchid, dissimilarity="Aitch")
Plot_PCoA(TAX=fit1, factor=batchid, dissimilarity="Aitch")
Plot_PCoA(TAX=fit2, factor=batchid, dissimilarity="Aitch")
dev.off()



###### Tune update ######

autism_taxa = read.csv("/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark/benchmarked_data/autism_2_microbiomeHD_count_data.csv")
autism_meta = read.csv("/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark/benchmarked_data/autism_2_microbiomeHD_meta_data.csv")

rownames(autism_taxa) = autism_taxa$X
autism_taxa = autism_taxa[, -1]

# length(unique(autism_meta$Sam_id))
# identical(autism_meta$Sam_id, autism_taxa$X)

apply(autism_taxa, 1, sum)
summary( apply(autism_taxa, 2, function(z){mean(z == 0)}) ) 

table(autism_meta[, 1], autism_meta[, 2])

batchid = factor(autism_meta$Dataset)
covar = factor(autism_meta$DiseaseState)
  
# ConQuR
taxa_tuned = Tune_ConQuR(tax_tab=autism_taxa, batchid=batchid, covariates=covar,
                         batch_ref_pool="asd_son",
                         logistic_lasso_pool=c(T, F),
                         quantile_type_pool=c("standard", "lasso", "composite"),
                         simple_match_pool=c(T, F),
                         lambda_quantile_pool=c(NA, "2p/n", "2p/logn"),
                         interplt_pool=c(T, F),
                         frequencyL=0,
                         frequencyU=1,
                         cutoff=0.25)

fit_autism = ConQuR(tax_tab=autism_taxa, batchid=batchid, covariates=covar, batch_ref="asd_son")
fit_autism_lasso = ConQuR(tax_tab=autism_taxa, batchid=batchid, covariates=covar, batch_ref="asd_son",
                    logistic_lasso=T, quantile_type="lasso")
fit_autism_simple = ConQuR(tax_tab=autism_taxa, batchid=batchid, covariates=covar, batch_ref="asd_son", simple_match=T)

# save(taxa_tuned, fit_autism, fit_autism_lasso, fit_autism_simple, file="autism_results.Rdata")
output_root = 'autism_2_microbiomeHD'
write.csv(taxa_tuned, paste(output_root, "_count_tuned.csv", sep=""), row.names = TRUE)
write.csv(fit_autism, paste(output_root, "_count_basic.csv", sep=""), row.names = TRUE)
write.csv(fit_autism_simple, paste(output_root, "_count_simple_basic.csv", sep=""), row.names = TRUE)
write.csv(fit_autism_lasso, paste(output_root, "_count_lasso.csv", sep=""), row.names = TRUE)


# ConQuR_libsize 
taxa_tuned_libsize = Tune_ConQuR_libsize(tax_tab=autism_taxa, batchid=batchid, covariates=covar,
                                         batch_ref_pool="asd_son",
                                         logistic_lasso_pool=c(T, F),
                                         quantile_type_pool=c("standard", "lasso", "composite"),
                                         simple_match_pool=c(T, F),
                                         lambda_quantile_pool=c(NA, "2p/n", "2p/logn"),
                                         interplt_pool=c(T, F),
                                         frequencyL=0,
                                         frequencyU=1,
                                         cutoff=0.25) 

fit_autism_libsize = ConQuR_libsize(tax_tab=autism_taxa, batchid=batchid, covariates=covar, batch_ref="asd_son")
fit_autism_lasso_libsize = ConQuR_libsize(tax_tab=autism_taxa, batchid=batchid, covariates=covar, batch_ref="asd_son",
                                          logistic_lasso=T, quantile_type="lasso")
fit_autism_simple_libsize = ConQuR_libsize(tax_tab=autism_taxa, batchid=batchid, covariates=covar, batch_ref="asd_son", simple_match=T)
  
# save(taxa_tuned_libsize, fit_autism_libsize, fit_autism_lasso_libsize, fit_autism_simple_libsize, file="autism_results_libsize.Rdata")

write.csv(taxa_tuned_libsize, paste(output_root, "_count_libsize_tuned.csv", sep=""), row.names = TRUE)
write.csv(fit_autism_libsize, paste(output_root, "_count_libsize_basic.csv", sep=""), row.names = TRUE)
write.csv(fit_autism_simple_libsize, paste(output_root, "_count_libsize_basic_simple.csv", sep=""), row.names = TRUE)
write.csv(fit_autism_lasso_libsize, paste(output_root, "_count_libsize_lasso.csv", sep=""), row.names = TRUE)

# ConQuR_libsize old version
taxa_tuned_libsize = Tune_ConQuR_libsize(tax_tab=autism_taxa, batchid=batchid, covariates=covar,
                                         batch_ref_pool="asd_son",
                                         logistic_lasso_pool=c(T, F),
                                         quantile_type_pool=c("standard", "lasso", "composite"),
                                         simple_match_pool=c(T, F),
                                         lambda_quantile_pool=c(NA, "2p/n", "2p/logn"),
                                         interplt_pool=c(T, F),
                                         frequencyL=0,
                                         frequencyU=1,
                                         cutoff=0.25) 

fit_autism_libsize = ConQuR_libsize(tax_tab=autism_taxa, batchid=batchid, covariates=covar, batch_ref="asd_son")
fit_autism_lasso_libsize = ConQuR_libsize(tax_tab=autism_taxa, batchid=batchid, covariates=covar, batch_ref="asd_son",
                                          logistic_lasso=T, quantile_type="lasso")
fit_autism_simple_libsize = ConQuR_libsize(tax_tab=autism_taxa, batchid=batchid, covariates=covar, batch_ref="asd_son", simple_match=T)

# save(taxa_tuned_libsize, fit_autism_libsize, fit_autism_lasso_libsize, fit_autism_simple_libsize, file="autism_results_libsize_old.Rdata")
write.csv(taxa_tuned_libsize, paste(output_root, "_count_libsize_old_tuned.csv", sep=""), row.names = TRUE)
write.csv(fit_autism_libsize, paste(output_root, "_count_libsize_old_basic.csv", sep=""), row.names = TRUE)
write.csv(fit_autism_simple_libsize, paste(output_root, "_count_libsize_old_basic_simple.csv", sep=""), row.names = TRUE)
write.csv(fit_autism_lasso_libsize, paste(output_root, "_count_libsize_old_lasso.csv", sep=""), row.names = TRUE)



pdf("autism_count.pdf", width=20, height=8)
par(mfrow=c(2, 5))

# ConQuR, batchid corrected

load(file="autism_results.Rdata")
Plot_PCoA(TAX=autism_taxa, factor=batchid)
Plot_PCoA(TAX=taxa_tuned$tax_final, factor=batchid)
Plot_PCoA(TAX=fit_autism, factor=batchid)
Plot_PCoA(TAX=fit_autism_lasso, factor=batchid)
Plot_PCoA(TAX=fit_autism_simple, factor=batchid)

Plot_PCoA(TAX=autism_taxa, factor=batchid, dissimilarity="Aitch")
Plot_PCoA(TAX=taxa_tuned$tax_final, factor=batchid, dissimilarity="Aitch")
Plot_PCoA(TAX=fit_autism, factor=batchid, dissimilarity="Aitch")
Plot_PCoA(TAX=fit_autism_lasso, factor=batchid, dissimilarity="Aitch")
Plot_PCoA(TAX=fit_autism_simple, factor=batchid, dissimilarity="Aitch")

# ConQuR_libsize, batchid corrected, jitter before divide libsize and recover in proposed.method.fast.libsize

load(file="autism_results_libsize.Rdata")
Plot_PCoA(TAX=autism_taxa, factor=batchid)
Plot_PCoA(TAX=taxa_tuned_libsize$tax_final, factor=batchid)
Plot_PCoA(TAX=fit_autism_libsize, factor=batchid)
Plot_PCoA(TAX=fit_autism_lasso_libsize, factor=batchid)
Plot_PCoA(TAX=fit_autism_simple_libsize, factor=batchid)

Plot_PCoA(TAX=autism_taxa, factor=batchid, dissimilarity="Aitch")
Plot_PCoA(TAX=taxa_tuned_libsize$tax_final, factor=batchid, dissimilarity="Aitch")
Plot_PCoA(TAX=fit_autism_libsize, factor=batchid, dissimilarity="Aitch")
Plot_PCoA(TAX=fit_autism_lasso_libsize, factor=batchid, dissimilarity="Aitch")
Plot_PCoA(TAX=fit_autism_simple_libsize, factor=batchid, dissimilarity="Aitch")

# ConQuR_libsize, batchid corrected

load(file="autism_results_libsize.Rdata")
Plot_PCoA(TAX=autism_taxa, factor=batchid)
Plot_PCoA(TAX=taxa_tuned_libsize$tax_final, factor=batchid)
Plot_PCoA(TAX=fit_autism_libsize, factor=batchid)
Plot_PCoA(TAX=fit_autism_lasso_libsize, factor=batchid)
Plot_PCoA(TAX=fit_autism_simple_libsize, factor=batchid)

Plot_PCoA(TAX=autism_taxa, factor=batchid, dissimilarity="Aitch")
Plot_PCoA(TAX=taxa_tuned_libsize$tax_final, factor=batchid, dissimilarity="Aitch")
Plot_PCoA(TAX=fit_autism_libsize, factor=batchid, dissimilarity="Aitch")
Plot_PCoA(TAX=fit_autism_lasso_libsize, factor=batchid, dissimilarity="Aitch")
Plot_PCoA(TAX=fit_autism_simple_libsize, factor=batchid, dissimilarity="Aitch")

dev.off()


