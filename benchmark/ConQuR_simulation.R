
###### starting data ######

library(phyloseq)
library(HMP2Data)

mom = momspi16S()
mom = prune_samples(momspi16S_samp$sample_body_site=="vagina" &
                      momspi16S_samp$visit_number==4,
                    mom)
mom = prune_samples(sample_sums(mom)>=4000, mom)
zp = rowMeans(otu_table(mom)==0)
mom = prune_taxa(zp <= 0.95, mom)

otu = matrix(otu_table(mom), ncol=nsamples(mom))
rownames(otu) = taxa_names(mom)
colnames(otu) = sample_names(mom)

shape = dirmult::dirmult(t(otu))$gamma
libsize = as.vector(colSums(otu))

id = NULL
zp = rowMeans(otu==0)
for (tau in 0:9/10) {
  id = c(id,order(abs(zp-tau))[1:2])
}

id1 = 1:5
id2 = 6:20

save(mom, otu, id, 
     shape, libsize,
     id1, id2,
     file="mom_270.Rdata")



###### simulation ######

CASE = 1 # 1 - 6, 6 scenarios
set.seed(CASE)

library(HMP)
library(dirmult)
library(vegan)
library(ade4)
library(compositions)
library(bindata)


load("mom_270.Rdata") 
p = nrow(otu)
n = ncol(otu)
batch_decrease_id = order(apply(otu, 1, function(z){ mean(z == 0)}))[2*(1:(p/2+1))-1]
batch_increase_id = order(apply(otu, 1, function(z){ mean(z == 0)}))[2*(1:(p/2))]
id_d = id[id1]
id_i = id[id2]

MC = 500
OR_cond_batchid = 1.25 
p_cond = 0.5

bincorr <- function(OR, p1, p2) {    #from odds ratio to binary correlation
  if (OR==1) p11=p2-p2+p1*p2 else {
    p11_1=p2-(1/2/(1-OR)*(1-p1+OR*p1+p2-OR*p2-
                            sqrt((-1+p1-OR*p1-p2+OR*p2)^2-4*(1-OR)*(p2-p1*p2))))
    p11_2=p2-(1/2/(1-OR)*(1-p1+OR*p1+p2-OR*p2-
                            sqrt((-1+p1-OR*p1-p2+OR*p2)^2)-4*(1-OR)*(p2-p1*p2)))
    if (p11_1>0 && p11_1<=p1 && p11_1<p2) p11=p11_1 else p11=p11_2
  }
  bincorr=(p11-p1*p2)/sqrt(p1*(1-p1)*p2*(1-p2))
  return(bincorr)
}


# scenario
effect_size_candidate = 0:2
libsize_influence_batchid_candidate = 0:1
eg = expand.grid(effect_size_candidate=effect_size_candidate, libsize_influence_batchid_candidate=libsize_influence_batchid_candidate)

e = eg[CASE,]

if (e[1] == 1){ # cond effect size > batch effect size
  cond_effect = 64
  batch_effect = 4
} else if (e[1] == 2){ # cond effect size < batch effect size
  cond_effect = 4
  batch_effect = 64
} else{
  cond_effect = 16
  batch_effect = 1 
}


otu0 = array(dim=c(MC,n,p))
cond0 = matrix(nrow=MC, ncol=n)
batchid0 = matrix(nrow=MC, ncol=n)

for (ii in 1:2){
  
  if (e[2] == 0) {
    bin_corr = bincorr(OR_cond_batchid, p_cond, 0.5)
    cond_batchid_vec = rmvbin(n, c(p_cond, 0.5), bincorr=(1-bin_corr)*diag(2)+bin_corr)
  } else {
    p_batchid = 1 / (1 + exp(-scale(libsize)))
    bin_corr = sapply(p_batchid, function(x) bincorr(OR_cond_batchid, p_cond, x))
    
    cond_batchid_vec = matrix(ncol=2, nrow=n)
    for (i in 1:n){
      cond_batchid_vec[i, ] = rmvbin(1, c(p_cond, p_batchid[i]), bincorr=(1-bin_corr[i])*diag(2)+bin_corr[i])
    }
  }
  cond = cond_batchid_vec[, 1]
  batchid = cond_batchid_vec[, 2]
  
  otu_tmp = matrix(nrow = p, ncol = n)
  for (i in 1:n) {
    otu_tmp[, i] = bayesm::rdirichlet(otu[, i] + 0.5) * sum(otu[, i])
  }
  
  for (i in which(cond==1)) {
    # add cond effect
    s1 = sum(otu_tmp[id_d, i])
    s2 = sum(otu_tmp[id_i, i])  
    
    if (s2 != 0){
      for (i1 in id_d) {
        otu_tmp[i1, i] = round(otu_tmp[i1, i] / cond_effect)
      }
      
      cond_effect2 = 1 + s1/s2*(cond_effect - 1)/cond_effect
      
      for (i1 in id_i) {
        otu_tmp[i1, i] = round(otu_tmp[i1, i] * cond_effect2)
      }
    }
  }
  
  for (i in which(batchid==1)){
    # add batch effect 
    s1 = sum(otu_tmp[batch_decrease_id, i])
    s2 = sum(otu_tmp[batch_increase_id, i])  
    
    if (s2 != 0){
      for (i1 in batch_decrease_id) {
        otu_tmp[i1, i] = round(otu_tmp[i1, i] / batch_effect)
      }
      
      batch_effect2 = 1 + s1/s2*(batch_effect - 1)/batch_effect
      
      for (i1 in batch_increase_id) {
        otu_tmp[i1, i] = round(otu_tmp[i1, i] * batch_effect2)
      }
    }
  }
  
  otu_tmp = round(otu_tmp)
  
  otu0[ii,,] = t(otu_tmp) 
  cond0[ii,] = cond
  batchid0[ii,] = batchid
  
}

save(otu0, cond0, batchid0, 
     file = paste0("simulated/", CASE,".Rdata"))



###### analysis ######

CASE = 1 # 1 - 3000, 6 scenarios x 500 MC

load(paste0("simulated/", ceiling(CASE / 500), ".Rdata"))
ii = CASE %% 500
if (ii == 0) ii = 500

taxa = otu0[ii,,] 

cond = factor( cond0[ii,] ) 
covar = data.frame(cond)

n = nrow(taxa)
batchid = factor( batchid0[ii,] )


# ConQuR
library(doParallel)
library(ConQuR)
taxa_tuned = Tune_ConQuR(tax_tab=taxa, batchid=batchid, covariates=covar,
                         batch_ref_pool="0",
                         logistic_lasso_pool=c(T, F),
                         quantile_type_pool=c("standard", "lasso", "composite"),
                         simple_match_pool=c(T, F),
                         lambda_quantile_pool=c(NA, "2p/n", "2p/logn"),
                         interplt_pool=c(T, F),
                         frequencyL=0,
                         frequencyU=1,
                         cutoff=0.25)

# ComBat-seq
library(sva)
adjusted_counts <- t( ComBat_seq(t(taxa), batch=batchid, group=NULL, covar_mod=covar, full_mod=FALSE) )

# MMUPHin
library(MMUPHin)
meta = data.frame(cond, batchid)
rownames(taxa) = 1:nrow(taxa)
fit_adjust_batch <- adjust_batch(feature_abd = t(taxa),
                                 batch = "batchid",
                                 covariates = "cond",
                                 data = meta,
                                 control = list(verbose = FALSE, diagnostic_plot = NULL))
MMUPHin_adj <- t(fit_adjust_batch$feature_abd_adj)

# Others...


