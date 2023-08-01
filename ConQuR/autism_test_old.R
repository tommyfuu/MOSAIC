
overall_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR'
source(paste0(overall_path, "/ConQuR_help_functions.R"))
source(paste0(overall_path, "/ConQuR_main_tune.R"))
source(paste0(overall_path, "/ConQuR_help_functions_libsize_old.R"))
source(paste0(overall_path, "/ConQuR_main_tune_libsize_old.R"))
source(paste0(overall_path, "/ConQuR_help_functions_rel.R"))
source(paste0(overall_path, "/ConQuR_main_tune_rel.R"))
source(paste0(overall_path, "/supporting_functions.R"))

###### Tune update ######

autism_taxa = read.csv("/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark/benchmarked_data/autism_2_microbiomeHD_count_data.csv")
autism_meta = read.csv("/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark/benchmarked_data/autism_2_microbiomeHD_meta_data.csv")

rownames(autism_taxa) = autism_taxa$X
autism_taxa = autism_taxa[, -1]


apply(autism_taxa, 1, sum)
summary( apply(autism_taxa, 2, function(z){mean(z == 0)}) ) 

table(autism_meta[, 1], autism_meta[, 2])

batchid = factor(autism_meta$Dataset)
covar = factor(autism_meta$DiseaseState)
  
# ConQuR
# taxa_tuned = Tune_ConQuR(tax_tab=autism_taxa, batchid=batchid, covariates=covar,
#                          batch_ref_pool="asd_son",
#                          logistic_lasso_pool=c(T, F),
#                          quantile_type_pool=c("standard", "lasso", "composite"),
#                          simple_match_pool=c(T, F),
#                          lambda_quantile_pool=c(NA, "2p/n", "2p/logn"),
#                          interplt_pool=c(T, F),
#                          frequencyL=0,
#                          frequencyU=1,
#                          cutoff=0.25)

fit_autism = ConQuR(tax_tab=autism_taxa, batchid=batchid, covariates=covar, batch_ref="asd_son")
fit_autism_lasso = ConQuR(tax_tab=autism_taxa, batchid=batchid, covariates=covar, batch_ref="asd_son",
                    logistic_lasso=T, quantile_type="lasso")
fit_autism_simple = ConQuR(tax_tab=autism_taxa, batchid=batchid, covariates=covar, batch_ref="asd_son", simple_match=T)

# save(taxa_tuned, fit_autism, fit_autism_lasso, fit_autism_simple, file="autism_results.Rdata")
output_root = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR/autism_2_microbiomeHD'
# save(taxa_tuned, fit_autism, fit_autism_lasso, fit_autism_simple, file="/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR/autism_results.Rdata")
save(fit_autism, fit_autism_lasso, fit_autism_simple, file="/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR/autism_results_1.Rdata")


# # ConQuR_libsize - old
# taxa_tuned_libsize_old = Tune_ConQuR_libsize(tax_tab=autism_taxa, batchid=batchid, covariates=covar,
#                                          batch_ref_pool="asd_son",
#                                          logistic_lasso_pool=c(T, F),
#                                          quantile_type_pool=c("standard", "lasso", "composite"),
#                                          simple_match_pool=c(T, F),
#                                          lambda_quantile_pool=c(NA, "2p/n", "2p/logn"),
#                                          interplt_pool=c(T, F),
#                                          frequencyL=0,
#                                          frequencyU=1,
#                                          cutoff=0.25) 

fit_autism_libsize_old = ConQuR_libsize(tax_tab=autism_taxa, batchid=batchid, covariates=covar, batch_ref="asd_son")
fit_autism_lasso_libsize_old = ConQuR_libsize(tax_tab=autism_taxa, batchid=batchid, covariates=covar, batch_ref="asd_son",
                                          logistic_lasso=T, quantile_type="lasso")
fit_autism_simple_libsize_old = ConQuR_libsize(tax_tab=autism_taxa, batchid=batchid, covariates=covar, batch_ref="asd_son", simple_match=T)
  
# save(taxa_tuned_libsize_old, fit_autism_libsize_old, fit_autism_lasso_libsize_old, fit_autism_simple_libsize_old, file="./autism_results_libsize_old.Rdata")
save(fit_autism_libsize_old, fit_autism_lasso_libsize_old, fit_autism_simple_libsize_old, file="./autism_results_libsize_old_1.Rdata")


# plot
# load(file="/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR/autism_results_libsize_old.Rdata")
load(file="/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR/autism_results_libsize_old_1.Rdata")
load(file="/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR/autism_results_1.Rdata")

pdf("/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR/autism_count_old_1.pdf", width=20, height=8)
# par(mfrow=c(2, 5))
par(mfrow=c(2, 5))
Plot_PCoA(TAX=autism_taxa, factor=batchid)
# Plot_PCoA(TAX=taxa_tuned_libsize_old$tax_final, factor=batchid)
Plot_PCoA(TAX=fit_autism, factor=batchid)
Plot_PCoA(TAX=fit_autism_libsize_old, factor=batchid)
Plot_PCoA(TAX=fit_autism_lasso_libsize_old, factor=batchid)
Plot_PCoA(TAX=fit_autism_simple_libsize_old, factor=batchid)

Plot_PCoA(TAX=autism_taxa, factor=batchid, dissimilarity="Aitch")
# Plot_PCoA(TAX=taxa_tuned_libsize_old$tax_final, factor=batchid, dissimilarity="Aitch")
Plot_PCoA(TAX=fit_autism, factor=batchid, dissimilarity="Aitch")
Plot_PCoA(TAX=fit_autism_libsize_old, factor=batchid, dissimilarity="Aitch")
Plot_PCoA(TAX=fit_autism_lasso_libsize_old, factor=batchid, dissimilarity="Aitch")
Plot_PCoA(TAX=fit_autism_simple_libsize_old, factor=batchid, dissimilarity="Aitch")

dev.off()
write.csv(taxa_tuned, paste(output_root, "_count_tuned.csv", sep=""), row.names = TRUE)
write.csv(fit_autism, paste(output_root, "_count_basic.csv", sep=""), row.names = TRUE)
write.csv(fit_autism_simple, paste(output_root, "_count_simple_basic.csv", sep=""), row.names = TRUE)
write.csv(fit_autism_lasso, paste(output_root, "_count_lasso.csv", sep=""), row.names = TRUE)

write.csv(taxa_tuned_libsize_old, paste(output_root, "_count_libsize_old_tuned.csv", sep=""), row.names = TRUE)
write.csv(fit_autism_libsize_old, paste(output_root, "_count_libsize_old_basic.csv", sep=""), row.names = TRUE)
write.csv(fit_autism_simple_libsize_old, paste(output_root, "_count_libsize_old_basic_simple.csv", sep=""), row.names = TRUE)
write.csv(fit_autism_lasso_libsize_old, paste(output_root, "_count_libsize_old_lasso.csv", sep=""), row.names = TRUE)


