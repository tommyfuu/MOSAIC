
overall_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR'
source(paste0(overall_path, "/ConQuR_help_functions.R"))
source(paste0(overall_path, "/ConQuR_main_tune.R"))
source(paste0(overall_path, "/ConQuR_help_functions_libsize.R"))
source(paste0(overall_path, "/ConQuR_main_tune_libsize.R"))
# source(paste0(overall_path, /ConQuR_help_functions_libsize_old.R"))
# source(paste0(overall_path, /ConQuR_main_tune_libsize_old.R")
source(paste0(overall_path, "/ConQuR_help_functions_rel.R"))
source(paste0(overall_path, "/ConQuR_main_tune_rel.R"))
source(paste0(overall_path, "/supporting_functions.R"))


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

print("checkpoint 1")
# ### the key point is to make sure the data consists of complete cases only ###
# ### in the future, may update the code to do sanity check and manipulation ###

# print("ibd_meta_0")
# print(ibd_meta)

id_to_delete = which(ibd_meta$gender == "" | is.na(ibd_meta$age))
print('id_to_delete')
print(id_to_delete)
ibd_meta = ibd_meta[-id_to_delete, ]
ibd_taxa = ibd_taxa[-id_to_delete, ]

batchid = factor(ibd_meta$study_name)
covar = data.frame(factor(ibd_meta$disease), factor(ibd_meta$gender), ibd_meta$age)
# print("checkpoint 2")
# # ConQuR_rel
# fit1 = ConQuR_rel(tax_tab=ibd_taxa, batchid=batchid, covariates=covar,
#                   batch_ref="HMP_2019_ibdmdb", num_core = 10)

# print("checkpoint 3")
# fit2 = ConQuR_rel(tax_tab=ibd_taxa, batchid=batchid, covariates=covar,
#                   batch_ref="HMP_2019_ibdmdb",
#                   logistic_lasso=T, quantile_type="lasso", num_core = 10)

# output_root = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR/idb_3_CMD'
# # write.csv(fit1, paste(output_root, "_relab_basic.csv", sep=""), row.names = TRUE)
# # write.csv(fit2, paste(output_root, "_relab_lasso.csv", sep=""), row.names = TRUE)
# print("checkpoint 4")
# # ConQuR_rel tuned
# ## Tom: I didn't run the following fine-tuned version as it is very time-consuming. Please help run it on SCU after your account is well set-up ###
# taxa_tuned_rel = Tune_ConQuR_rel(tax_tab=ibd_taxa, batchid=batchid, covariates=covar,
#                                  batch_ref_pool="HMP_2019_ibdmdb",
#                                  logistic_lasso_pool=c(T, F),
#                                  quantile_type_pool=c("standard", "lasso", "composite"),
#                                  simple_match_pool=c(T, F),
#                                  lambda_quantile_pool=c(NA, "2p/n", "2p/logn"),
#                                  interplt_pool=c(T, F),
#                                  frequencyL=0,
#                                  frequencyU=1,
#                                  cutoff=0.25, num_core = 10) 
# print("checkpoint 5")
# save(fit1, fit2, file="/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR/relative_results.Rdata")

# save(taxa_tuned_rel, file="/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR/relative_tuned_results.Rdata")
# # write.csv(taxa_tuned_rel, paste(output_root, "_relab_tuned.csv", sep=""), row.names = TRUE)

load(file="/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR/relative_results.Rdata")
load(file="relative_tuned_results.Rdata")
print("checkpoint 6")
pdf("/athena/linglab/scratch/chf4012/mic_bc_benchmark/ConQuR/ibd_relative.pdf", width=12, height=8)
par(mfrow=c(2, 3))
Plot_PCoA(TAX=ibd_taxa, factor=batchid)
# Plot_PCoA(TAX=taxa_tuned_rel$tax_final, factor=batchid)
Plot_PCoA(TAX=fit1, factor=batchid)
Plot_PCoA(TAX=fit2, factor=batchid)
Plot_PCoA(TAX=taxa_tuned_rel, factor=batchid)

Plot_PCoA(TAX=ibd_taxa, factor=batchid, dissimilarity="Aitch")
# Plot_PCoA(TAX=taxa_tuned_rel$tax_final, factor=batchid, dissimilarity="Aitch")
Plot_PCoA(TAX=fit1, factor=batchid, dissimilarity="Aitch")
Plot_PCoA(TAX=fit2, factor=batchid, dissimilarity="Aitch")
Plot_PCoA(TAX=taxa_tuned_rel, factor=batchid, dissimilarity="Aitch")
dev.off()

# print("checkpoint 7")
# write.csv(fit1, paste(output_root, "_relab_basic.csv", sep=""), row.names = TRUE)
# write.csv(fit2, paste(output_root, "_relab_lasso.csv", sep=""), row.names = TRUE)
# write.csv(taxa_tuned_rel, paste(output_root, "_relab_tuned.csv", sep=""), row.names = TRUE)


