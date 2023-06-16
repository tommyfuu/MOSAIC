# load data
load("/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/ibd_150.Rdata") 
# otu = read.csv("/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/data/QinJ_2012_T2D/otu_table_QinJ_2012.csv", header = TRUE, row.names = 1)
library("bindata")
library("MIDAS")

# instantiate variables
otu_original = t(otu)
libsize = rowSums(otu_original)
bin_corr = 0.3
cond_effect = 0.3
batch_effect = 0.3
p = ncol(otu_original)


if (p==301) {
  print("p=301")
  n = nrow(otu_original) * 5
  id_batch = 101:250
  id_cond = 1:150
} 
if (p==233) {
  print("p=233")
  n = nrow(otu_original) 
  id_batch = c(21:40,51:70,81:160)
  id_cond = c(1:20,41:50,71:80,121:200)
}

print("checkpoint 1")
# function to generate simulated data with MIDAs - linear combination
midas_simulate <- function(otu_original, n, bin_corr, cond_effect, batch_effect, out_count, out_relab){

  # generate batch and condition id
  cond_batchid_vec = rmvbin(n, c(0.5, 0.5), bincorr=(1-bin_corr)*diag(2)+bin_corr)
  cond = cond_batchid_vec[, 1] # find the ids of taxa that are affected by condition
  batchid = cond_batchid_vec[, 2] # find the ids of taxa that are affected by batch

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
                            lib.size = rep(10000,n))
  fitted_c1b0 = Midas.modify(fitted, 
                            lib.size = rep(10000,n),
                            mean.rel.abund.1 = (1-cond_effect)*rel10 + cond_effect*rel1_cond,
                            mean.rel.abund = (1-cond_effect)*rel0 + cond_effect*rel_cond)
  fitted_c0b1 = Midas.modify(fitted, 
                            lib.size = rep(10000,n),
                            mean.rel.abund.1 = (1-batch_effect)*rel10 + batch_effect*rel1_batch,
                            mean.rel.abund = (1-batch_effect)*rel0 + batch_effect*rel_batch)
  fitted_c1b1 = Midas.modify(fitted, 
                            lib.size = rep(10000,n),
                            mean.rel.abund.1 = (1-cond_effect-batch_effect)*rel10 + cond_effect*rel1_cond + batch_effect*rel1_batch,
                            mean.rel.abund = (1-cond_effect-batch_effect)*rel0 + cond_effect*rel_cond + batch_effect*rel_batch)


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
# function to generate simulated data with MIDAs - linear combination
midas_FC_simulate <- function(otu_original, n, cond_effect_FC, batch_effect_FC, out_count, out_relab){
  # simulate data
  fitted = Midas.setup(otu_original, fit.beta=FALSE)

  ## generate raw simulated data without set library size
  fitted_null = Midas.modify(fitted)
  otu_null = Midas.sim(fitted_null)$sim_count

  ## generate r vector where an element is 1 if id_cond is 1 and 0 otherwise
  cond = rep(0,n)
  cond[id_cond] = 1
  ## repeat the same manipulation for batch
  batchid = rep(0,n)
  batchid[id_batch] = 1

  ## manipulate taxa by fold changes
  # otu_FC = matrix(NA, nrow = n, ncol = p)
  otu_FC = otu_null
  otu_FC = round(otu_FC,0)
  # print(cond==0&batchid==0)
  otu_FC[cond==0&batchid==0] = otu_null[cond==0&batchid==0] # unadjusted ones
  otu_FC[cond==1&batchid==0] = otu_null[cond==1&batchid==0]*cond_effect_FC # ones adjuetd only for condition
  otu_FC[cond==0&batchid==1] = otu_null[cond==0&batchid==1]*batch_effect_FC # ones adjuetd only for batch
  otu_FC[cond==1&batchid==1] = otu_null[cond==1&batchid==1]*cond_effect_FC*batch_effect_FC # ones adjuetd for both condition and batch (correlated by the coefficient)

  # save otu
  write.csv(otu_FC, file = out_count, row.names = FALSE)
  # relative abundance
  rela_FC = otu_FC / rowSums(otu_FC)
  write.csv(rela_FC, file = out_relab, row.names = FALSE)
}
print("checkpoint 3")
# generate simulated data using Jiuyao's set up + FC set up  
scaled_midas_data_generation <- function(otu_original, n, bin_corr_val_l, cond_effect_val_l, batch_effect_val_l, num_iter){   
  for (bin_corr_val in bin_corr_val_l) {
    for (cond_effect_val in cond_effect_val_l) {
      for (batch_effect_val in batch_effect_val_l) {
        if (cond_effect_val + batch_effect_val <= 1) {
          for (iter in seq(1, num_iter)){
            print(bin_corr_val)
            print(cond_effect_val)
            print(batch_effect_val)
            
            output_file_path_count = paste0("/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/data/simulation_data_midas_1/ibd_150_count_", bin_corr_val, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter, ".csv")
            print(output_file_path_count)
            output_file_path_relab = paste0("/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/data/simulation_data_midas_1/ibd_150_relab_", bin_corr_val, "_", cond_effect_val, "_", batch_effect_val, '_iter_', iter, ".csv")
            
            if (file.exists(output_file_path_count)){
              print(output_file_path_count)
              print("file exists")
              print("___")
            }
            else{
              midas_simulate(otu_original, n, bin_corr_val, cond_effect_val, batch_effect_val, output_file_path_count, output_file_path_relab)
              print("___")
            }
          }
        }
      }
    }
  }
}
print("checkpoint 4")
# repeat but for the FC version
scaled_midas_FC_data_generation <- function(otu_original, n, cond_effect_val_l, batch_effect_val_l, num_iter){  
  for (cond_effect_val in cond_effect_val_l) {
    for (batch_effect_val in batch_effect_val_l) {
      for (iter in seq(1, num_iter)){
        print(cond_effect_val)
        print(batch_effect_val)
        print(iter)
        output_file_path_count = paste0("/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/data/simulation_data_midas_FC/ibd_150_count_", cond_effect_val, "_", batch_effect_val, '_iter_', iter, ".csv")
        output_file_path_relab = paste0("/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/data/simulation_data_midas_FC/ibd_150_relab_", cond_effect_val, "_", batch_effect_val, '_iter_', iter, ".csv")
        
        if (file.exists(output_file_path_count)){
          print(output_file_path_count)
          print("file exists")
          print("___")
        }
        else{
          midas_FC_simulate(otu_original, n, cond_effect_val, batch_effect_val, output_file_path_count, output_file_path_relab)
          print("___")
        }
      }
    }
  }
}



# bin_corr_val_l = c(0, 0.1, 0.3, 0.5, 0.7, 0.9)
# cond_effect_val_l = c(0, 0.099, 0.299, 0.499, 0.699, 0.899)
# batch_effect_val_l = c(0, 0.099, 0.299, 0.499, 0.699, 0.899)
# scaled_midas_data_generation(otu_original, n, bin_corr_val_l, cond_effect_val_l, batch_effect_val_l, num_iter=10)

# test for effect
cond_effect_val_l = c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)
batch_effect_val_l = c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)
scaled_midas_FC_data_generation(otu_original, n, cond_effect_val_l, batch_effect_val_l, num_iter=10)


# midas_bc_biovar(otu_original, n, bin_corr, cond_effect, batch_effect, "./ibd_150_relab.csv")

# midas_bc_biovar(otu_original, n, bin_corr, cond_effect, batch_effect, "./ibd_150_relab.csv")
# midas_bc_biovar(otu_original, n, bin_corr, cond_effect, batch_effect, "./QinJ_2012_relab.csv")
# write relative abundance to csv
# write.csv(rela, file = "./ibd_150_relab.csv", row.names = FALSE, out = )



# 1. modify count:
## - 1.1 use midas to generate data
# fitted_c0b0 = Midas.modify(fitted,
#                             lib.size = rep(10000,n)) - can even remove the library size so we can use the original data's library size
## - 1.2 multiply certain taxa by the fold change number
### certain taxa can be randomly selected or selected using id_batch = 101:250
#### the taxa that are already really big in the original data should not be multiplied by a fold change that's too big
#### the taxa that are already really small in the original data would not be impacted by the fold change too much
##### or we can just normalize back so the rel abund is still summing to 1

# 2. when simulating with midas, 
## dirichlet multi-normal distribution, multiply mean factor by fc
