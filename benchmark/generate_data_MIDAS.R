load("/Users/chenlianfu/Documents/GitHub/mic_bc_benchmark/benchmark/ibd_150.Rdata") 
library("bindata")
otu_original = otu
libsize = rowSums(otu_original)

p = ncol(otu_original)


if (p==301) {
  n = nrow(otu_original) * 5
  id_batch = 101:250
  id_cond = 1:150
} 
if (p==233) {
  n = nrow(otu_original) 
  id_batch = c(21:40,51:70,81:160)
  id_cond = c(1:20,41:50,71:80,121:200)
}

cond_batchid_vec = rmvbin(n, c(0.5, 0.5), bincorr=(1-bin_corr)*diag(2)+bin_corr)
cond = cond_batchid_vec[, 1]
batchid = cond_batchid_vec[, 2]

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
otu[cond==0&batchid==0,] = otu_c0b0[cond==0&batchid==0,]
otu[cond==1&batchid==0,] = otu_c1b0[cond==1&batchid==0,]
otu[cond==0&batchid==1,] = otu_c0b1[cond==0&batchid==1,]
otu[cond==1&batchid==1,] = otu_c1b1[cond==1&batchid==1,]

rela = otu / rowSums(otu)