
overall_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR'
source(paste0(overall_path, "/ConQuR_help_functions.R"))
source(paste0(overall_path, "/ConQuR_main_tune.R"))
source(paste0(overall_path, "/ConQuR_help_functions_libsize.R"))
source(paste0(overall_path, "/ConQuR_main_tune_libsize.R"))
# source(paste0(overall_path, "/ConQuR_help_functions_libsize_old.R"))
# source(paste0(overall_path, "/ConQuR_main_tune_libsize_old.R"))
source(paste0(overall_path, "/ConQuR_help_functions_rel.R"))
source(paste0(overall_path, "/ConQuR_main_tune_rel.R"))
source(paste0(overall_path, "/supporting_functions.R"))

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
output_root = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR/autism_2_microbiomeHD'
save(taxa_tuned, fit_autism, fit_autism_lasso, fit_autism_simple, file="/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR/autism_results.Rdata")


# ConQuR_libsize - new, batchid corrected, jitter before divide libsize and recover in proposed.method.fast.libsize
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
  
save(taxa_tuned_libsize, fit_autism_libsize, fit_autism_lasso_libsize, fit_autism_simple_libsize, file="./autism_results_libsize.Rdata")



# plot
load(file="/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR/autism_results.Rdata")

pdf("/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR/autism_count.pdf", width=20, height=8)
par(mfrow=c(2, 5))
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

dev.off()

# ConQuR_libsize, batchid corrected, jitter before divide libsize and recover in proposed.method.fast.libsize

# load(file="/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR/autism_results_libsize.Rdata")
# pdf("/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR/autism_count_new.pdf", width=20, height=8)
# par(mfrow=c(2, 5))
# Plot_PCoA(TAX=autism_taxa, factor=batchid)
# Plot_PCoA(TAX=taxa_tuned_libsize$tax_final, factor=batchid)
# Plot_PCoA(TAX=fit_autism_libsize, factor=batchid)
# Plot_PCoA(TAX=fit_autism_lasso_libsize, factor=batchid)
# Plot_PCoA(TAX=fit_autism_simple_libsize, factor=batchid)

# Plot_PCoA(TAX=autism_taxa, factor=batchid, dissimilarity="Aitch")
# Plot_PCoA(TAX=taxa_tuned_libsize$tax_final, factor=batchid, dissimilarity="Aitch")
# Plot_PCoA(TAX=fit_autism_libsize, factor=batchid, dissimilarity="Aitch")
# Plot_PCoA(TAX=fit_autism_lasso_libsize, factor=batchid, dissimilarity="Aitch")
# Plot_PCoA(TAX=fit_autism_simple_libsize, factor=batchid, dissimilarity="Aitch")

# dev.off()
# write.csv(taxa_tuned, paste(output_root, "_count_tuned.csv", sep=""), row.names = TRUE)
# write.csv(fit_autism, paste(output_root, "_count_basic.csv", sep=""), row.names = TRUE)
# write.csv(fit_autism_simple, paste(output_root, "_count_simple_basic.csv", sep=""), row.names = TRUE)
# write.csv(fit_autism_lasso, paste(output_root, "_count_lasso.csv", sep=""), row.names = TRUE)

# write.csv(taxa_tuned_libsize, paste(output_root, "_count_libsize_old_tuned.csv", sep=""), row.names = TRUE)
# write.csv(fit_autism_libsize, paste(output_root, "_count_libsize_old_basic.csv", sep=""), row.names = TRUE)
# write.csv(fit_autism_simple_libsize, paste(output_root, "_count_libsize_old_basic_simple.csv", sep=""), row.names = TRUE)
# write.csv(fit_autism_lasso_libsize, paste(output_root, "_count_libsize_old_lasso.csv", sep=""), row.names = TRUE)


