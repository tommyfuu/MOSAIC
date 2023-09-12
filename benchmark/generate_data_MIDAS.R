# load data
load("/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark/ibd_150.Rdata")
library("bindata")
library("MIDAS")
library(tibble)

args = commandArgs(trailingOnly=TRUE)
print(args)
if (length(args)==0 || length(args)>1 ) {
  stop("There has to be exactly one argument supplied with this script for simulation runs", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  GLOBAL_ITER = args[1]
}

# instantiate variables
otu_original = t(otu)
p = ncol(otu_original)

# making sure we are generating from the right data
if (p==301) {
  print("p=301")
  n = nrow(otu_original) * 3 # 450 samples in this case
  # get library sizes of the original dataset
  print("library sizes of the original dataset")
  libsize_l = rowSums(otu_original)
  # print(libsize_l)
  print(length(libsize_l))
  sampled_libsize_l = sample(libsize_l, n, replace = TRUE)
} 

# function that turns from odds ratio to binary correlation (from ConQuR paper)
bincorr <- function(OR, p1, p2) {    
  if (OR==1) p11=p2-p2+p1*p2 else {
    p11_1=p2-(1/2/(1-OR)*(1-p1+OR*p1+p2-OR*p2-
                            sqrt((-1+p1-OR*p1-p2+OR*p2)^2-4*(1-OR)*(p2-p1*p2))))
    p11_2=p2-(1/2/(1-OR)*(1-p1+OR*p1+p2-OR*p2-
                            sqrt((-1+p1-OR*p1-p2+OR*p2)^2)-4*(1-OR)*(p2-p1*p2)))
    if (p11_1>0 && p11_1<p1 && p11_1<p2) p11=p11_1 else p11=p11_2
  }
  bincorr=(p11-p1*p2)/sqrt(p1*(1-p1)*p2*(1-p2))
  return(bincorr)
}

# ### get the 50 differentially abundant taxa for conditional effect and 250 be impacted by batch ###
print("checkpoint 1")
id_cond = NULL
id_batch = 1:301

zp = colMeans(otu_original==0)

for (tau in 0:24/25) {
  print(tau)
  id_cond = c(id_cond,order(abs(zp-tau))[1:2])
}


# function to generate from midas for each iteration
midas_generate_per_iter <- function(otu_original, output_root, or, cond_effect_val, batch_effect_val, iter, libsize_l, batch_libsize_related){
  output_file_path_count = paste0(output_root, "/ibd_150_count_", or, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter, ".csv")
  print(output_file_path_count)
  output_file_path_relab = paste0(output_root, "/ibd_150_relab_", or, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter, ".csv")
  output_file_path_meta = paste0(output_root, "/ibd_150_meta_", or, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter, ".csv")
  
  if (file.exists(output_file_path_relab)){
    print(output_file_path_relab)
    print("file exists")
    print("___")
  }
  else{
    print("cond_effect_val")
    print(cond_effect_val)
    midas_simulate(otu_original, n, or, cond_effect_val, batch_effect_val, output_file_path_count, output_file_path_relab, output_file_path_meta, libsize_l=libsize_l, batch_libsize_related=batch_libsize_related)
    print("___")
  }
}


# function to generate simulated data with MIDAs - linear combination
midas_simulate <- function(otu_original, n, or, cond_effect, batch_effect, out_count, out_relab, out_meta, p_cond = 0.5, p_batch = 0.5, libsize_l = NULL, batch_libsize_related = FALSE){
  if(is.null(libsize_l)){
    libsize_l = rep(10000,n)
  } 
  
  print(length(zp))
  # check if p_batch and libsize are related
  if(batch_libsize_related == TRUE){
    print("TRUE")
    p_batch = 1 / (1 + exp(-scale(libsize_l)))
    bin_corr = sapply(p_batch, function(x) bincorr(or, p_cond, x))
    
    cond_batchid_vec = matrix(ncol=2, nrow=n)
    for (i in 1:n){
      cond_batchid_vec[i, ] = rmvbin(1, c(p_cond, p_batch[i]), bincorr=(1-bin_corr[i])*diag(2)+bin_corr[i])
    }
  }
  else{
    print("FALSE")
    bin_corr = bincorr(or, p_cond, p_batch)
    print('OR, p_cond, p_batch, bin_corr')
    # print(p_cond)
    print(c(or, p_cond, p_batch))
    print(length(bin_corr))
    # print(c(or, p_cond, p_batch, bin_corr))
    cond_batchid_vec = rmvbin(n, c(p_cond, p_batch), bincorr=(1-bin_corr)*diag(2)+bin_corr)
  }
  # # turn this into bincorr with the bincorr function and document bincorr value as well
  # bin_corr = bincorr(or, p_cond, p_batch)
  # print('OR, p_cond, p_batch, bin_corr')
  # print(c(or, p_cond, p_batch, bin_corr))
  # # generate batch and condition id
  # cond_batchid_vec = rmvbin(n, c(p_cond, p_batch), bincorr=(1-bin_corr)*diag(2)+bin_corr)
  cond = cond_batchid_vec[, 1] # find the ids of samples that are affected by condition
  batchid = cond_batchid_vec[, 2] # find the ids of samples that are affected by batch

  # save metadata as a dataframe with an additional column being the subject id
  subjectid_text = sprintf("Subject_%d", 1:n) 
  current_metadata <- data.frame(subjectid_text, batchid, cond)

  current_metadata$batchid <- replace(current_metadata$batchid, current_metadata$batchid == 0, "batch_0")
  current_metadata$batchid <- replace(current_metadata$batchid, current_metadata$batchid == 1, "batch_1")
  
  current_metadata$cond <- replace(current_metadata$cond, current_metadata$cond == 0, "cond_0")
  current_metadata$cond <- replace(current_metadata$cond, current_metadata$cond == 1, "cond_1")
  write.csv(current_metadata, file = out_meta, row.names = FALSE)
  print("doing its thing")
  # print(cond)
  # simulate data
  fitted = Midas.setup(otu_original, fit.beta=FALSE)

  perm_batch = sample(id_batch, length(id_batch))
  perm_cond = sample(id_cond, length(id_cond))

  rel0 = fitted$mean.rel.abund
  rel_batch = fitted$mean.rel.abund
  rel_cond = fitted$mean.rel.abund
  rel_batch[id_batch] = rel_batch[perm_batch]
  rel_cond[id_cond] = rel_cond[perm_cond]

  rel10 = unlist(lapply(fitted$rel.abund.1,mean))
  rel1_batch = unlist(lapply(fitted$rel.abund.1,mean))
  rel1_cond = unlist(lapply(fitted$rel.abund.1,mean))
  rel1_batch[id_batch] = rel1_batch[perm_batch]
  rel1_cond[id_cond] = rel1_cond[perm_cond]


  fitted_c0b0 = Midas.modify(fitted,
                            lib.size = libsize_l)
                            # taxa.1.prop = 'same')
  fitted_c1b0 = Midas.modify(fitted, 
                            lib.size = libsize_l,
                            mean.rel.abund.1 = (1-cond_effect)*rel10 + cond_effect*rel1_cond,
                            mean.rel.abund = (1-cond_effect)*rel0 + cond_effect*rel_cond)
                            # taxa.1.prop = 'same')
  fitted_c0b1 = Midas.modify(fitted, 
                            lib.size = libsize_l,
                            mean.rel.abund.1 = (1-batch_effect)*rel10 + batch_effect*rel1_batch,
                            mean.rel.abund = (1-batch_effect)*rel0 + batch_effect*rel_batch)
                            # taxa.1.prop = 'same')
  fitted_c1b1 = Midas.modify(fitted, 
                            lib.size = libsize_l,
                            mean.rel.abund.1 = (1-cond_effect-batch_effect)*rel10 + cond_effect*rel1_cond + batch_effect*rel1_batch,
                            mean.rel.abund = (1-cond_effect-batch_effect)*rel0 + cond_effect*rel_cond + batch_effect*rel_batch)
                            # taxa.1.prop = 'same')


  otu_c0b0 = Midas.sim(fitted_c0b0)$sim_count
  otu_c1b0 = Midas.sim(fitted_c1b0)$sim_count
  otu_c0b1 = Midas.sim(fitted_c0b1)$sim_count
  otu_c1b1 = Midas.sim(fitted_c1b1)$sim_count

  otu = matrix(NA, nrow = n, ncol = p)
  # print(cond==0&batchid==0)
  otu[cond==0&batchid==0,] = otu_c0b0[cond==0&batchid==0,] # unadjusted ones
  otu[cond==1&batchid==0,] = otu_c1b0[cond==1&batchid==0,] # ones adjuetd only for condition
  otu[cond==0&batchid==1,] = otu_c0b1[cond==0&batchid==1,] # ones adjuetd only for batch
  otu[cond==1&batchid==1,] = otu_c1b1[cond==1&batchid==1,] # ones adjuetd for both condition and batch (correlated by the coefficient)

  # save otu
  write.csv(otu, file = out_count, row.names = FALSE)
  # relative abundance
  rela = otu / rowSums(otu)
  write.csv(rela, file = out_relab, row.names = FALSE)
}
print("checkpoint 2")


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


# generate simulated data using both set ups 08072023
scaled_midas_data_generation <- function(output_root, otu_original, n, or_l, cond_effect_val_l, batch_effect_val_l, num_iter, libsize_l, batch_libsize_related = FALSE){   
  for (or in or_l) {
    for (cond_effect_val in cond_effect_val_l) {
      for (batch_effect_val in batch_effect_val_l) {
        if (cond_effect_val + batch_effect_val <= 1) {
          mcsapply(seq(1, num_iter), function(iter) midas_generate_per_iter(otu_original, output_root, or, cond_effect_val, batch_effect_val, iter, libsize_l, batch_libsize_related), mc.cores = 10)
          # sapply(seq(1, num_iter), function(iter) midas_generate_per_iter(output_root, or, cond_effect_val, batch_effect_val, iter, libsize_l, batch_libsize_related))
        }
      }
    }
  }
}

# generate simulated data using both set ups 08072023
scaled_slurm_midas_data_generation <- function(output_root, otu_original, n, or_l, cond_effect_val_l, batch_effect_val_l, iter, libsize_l, batch_libsize_related = FALSE){   
  for (or in or_l) {
    for (cond_effect_val in cond_effect_val_l) {
      for (batch_effect_val in batch_effect_val_l) {
        if (cond_effect_val + batch_effect_val <= 1) {
          sink(file = paste0(output_root, "/ibd_150_id_cond_", or, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter, ".txt"))
          cat(id_cond, "\n")
          sink()
          sink(file = paste0(output_root, "/ibd_150_id_batch_", or, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter, ".txt"))
          cat(id_batch, "\n")
          sink()

          midas_generate_per_iter(otu_original, output_root, or, cond_effect_val, batch_effect_val, iter, libsize_l, batch_libsize_related)
        }
      }
    }
  }
}



# or_l = c(1, 1.25, 1.5)
# cond_effect_val_l = c(0, 0.25, 0.5, 0.75, 1)
# batch_effect_val_l = c(0, 0.25, 0.5, 0.75, 1)
# # output_root = '/athena/linglab/scratch/chf4012/simulation_data_updated_MIDAS_norelation_090723'
# # # scaled_midas_data_generation(output_root, otu_original, n, or_l, cond_effect_val_l, batch_effect_val_l, num_iter=1000, libsize_l=sampled_libsize_l, batch_libsize_related = TRUE)
# # scaled_slurm_midas_data_generation(output_root, otu_original, n, or_l, cond_effect_val_l, batch_effect_val_l, iter=GLOBAL_ITER, libsize_l=sampled_libsize_l, batch_libsize_related = TRUE)
# output_root = '/athena/linglab/scratch/chf4012/simulation_data_updated_MIDAS_yesrelation_090723'
# # scaled_midas_data_generation(output_root, otu_original, n, or_l, cond_effect_val_l, batch_effect_val_l, num_iter=1000, libsize_l=sampled_libsize_l, batch_libsize_related = FALSE)
# scaled_slurm_midas_data_generation(output_root, otu_original, n, or_l, cond_effect_val_l, batch_effect_val_l, iter=GLOBAL_ITER, libsize_l=sampled_libsize_l, batch_libsize_related = FALSE)


or_l = c(1, 1.25, 1.5)
cond_effect_val_l = c(0, 0.25, 0.5, 0.75, 1)
batch_effect_val_l = c(0, 0.25, 0.5, 0.75, 1)
# output_root = '/athena/linglab/scratch/chf4012/simulation_data_updated_MIDAS_norelation_090723'
# scaled_slurm_midas_data_generation(output_root, otu_original, n, or_l, cond_effect_val_l, batch_effect_val_l, iter=GLOBAL_ITER, libsize_l=sampled_libsize_l, batch_libsize_related = FALSE)
output_root = '/athena/linglab/scratch/chf4012/simulation_data_updated_MIDAS_yesrelation_090723'
scaled_slurm_midas_data_generation(output_root, otu_original, n, or_l, cond_effect_val_l, batch_effect_val_l, iter=GLOBAL_ITER, libsize_l=sampled_libsize_l, batch_libsize_related = TRUE)
