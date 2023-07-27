
library(quantreg)
library(cqrReg)
library(glmnet)
library(dplyr)
library(doParallel)
library(gplots)
library(vegan)
library(ade4)
library(compositions)
# library(randomForest)
# library(ROCR)
# library(ape)
# library(GUniFrac)
library(fastDummies)



### ConQuR_libsize main function 

ConQuR_libsize <- function(tax_tab, batchid, covariates,
                           libsize_tune=NULL,
                           batch_ref,
                           logistic_lasso=F, quantile_type="standard", simple_match=F,
                           lambda_quantile="2p/n", interplt=F,
                           delta=0.4999, taus=seq(0.005, 0.995, by=0.005), num_core=2){

  # relevel batch id
  batchid = relevel(batchid, ref=batch_ref)

  ###### libsize specific ######
  libsize = libsize_tune
  if (is.null(libsize)) libsize = apply(tax_tab, 1, sum)

  tax_tab_temp = tax_tab
  tax_tab_temp[tax_tab_temp > 0] = tax_tab_temp[tax_tab_temp > 0] + runif(length(tax_tab_temp[tax_tab_temp > 0]), min=0, max=1)
  for (ii in 1:nrow(tax_tab)){
    tax_tab_temp[ii, ] = tax_tab_temp[ii, ] / libsize[ii]
  }
  tax_tab = tax_tab_temp # transform count into relative abundance


  registerDoParallel(num_core)

  if (simple_match == T){

    #### by simple quantile-quantile matching is chosen ####

    ### correct each of the taxa ###
    tax_new = foreach (ll=1:ncol(tax_tab), .combine=cbind) %do%{
      y = as.numeric( tax_tab[, ll] )
      simple_QQ_libsize(y=y, batchid=batchid, batch_ref=batch_ref, taus=taus)
    }

  } else{

    #### otherwise, correct data via regression ####

    ### process data ###
    X = data.frame(covariates, batchid)
    X_span = model.matrix(~., X)[,-1]

    X_correct = X
    X_correct$batchid = batch_ref
    X_correct$batchid = factor(X_correct$batchid)

    X_span_correct = X_span
    X_span_correct[, grep( "batchid", colnames(X_span_correct) )] = 0

    ### correct each of the taxa ###
    tax_new = foreach (ll=1:ncol(tax_tab), .combine=cbind) %do%{
      y = as.numeric( tax_tab[, ll] )

      ###### libsize specific ######
      ConQuR_each_libsize(y=y, X=X, X_span=X_span, X_correct=X_correct, X_span_correct=X_span_correct, batchid=batchid, batch_ref=batch_ref,
                          libsize = libsize,
                          delta=delta, taus=taus, logistic_lasso=logistic_lasso, quantile_type=quantile_type, lambda_quantile=lambda_quantile, interplt=interplt)
    }

  }


  ###### libsize specific ######
  if (ncol(tax_tab) == 1) tax_new = matrix(tax_new, nrow=nrow(tax_tab))

  tax_count_new = matrix(ncol=ncol(tax_new), nrow=nrow(tax_new))
  for (ll in 1:nrow(tax_new)){
    tax_count_new[ll, ] = round( tax_new[ll, ] * libsize[ll], digits=0 ) # transform back to count by multiplying relative abundance by library size
  }

  rownames(tax_count_new) = rownames(tax_tab)
  colnames(tax_count_new) = colnames(tax_tab)
  return(tax_count_new)

}





### Tune ConQuR_libsize over variations

Tune_ConQuR_libsize <- function(tax_tab, batchid, covariates,
                                batch_ref_pool,
                                logistic_lasso_pool,
                                quantile_type_pool,
                                simple_match_pool,
                                lambda_quantile_pool,
                                interplt_pool,
                                frequencyL,
                                frequencyU,
                                cutoff=0.1, delta=0.4999, taus=seq(0.005, 0.995, by=0.005), num_core=2){

  # compute library size
  libsize_tune = apply(tax_tab, 1, sum)

  # prepare the method list, taxa pool corrresponding to the grid of frequency, and result table
  method1 = expand.grid(logistic=logistic_lasso_pool, quantile=quantile_type_pool, lambda=lambda_quantile_pool, interplt=interplt_pool)
  method1 = method1[-c( which(method1$quantile=="standard" & !is.na(method1$lambda)), which(method1$quantile!="standard" & is.na(method1$lambda)) ), ]
  method2 = any(simple_match_pool == T)

  freq_grid = seq(from=frequencyL, to=frequencyU, by=cutoff)
  prate = apply(tax_tab, 2, function(z){length( which(z > 0) ) / nrow(tax_tab) })

  cut_list = NULL
  for (ii in 1:(length(freq_grid)-1)){
    start = freq_grid[ii]
    end = freq_grid[ii+1]

    cut_list[[ii]] = which(prate > start & prate <= end)
  }

  R2_initial = data.frame(method=c("original", apply(method1, 1, paste, collapse="_")), batch_R2=NA)
  if (method2 == T) R2_initial = rbind(R2_initial, c("simple", NA))

  # greedy experiment on reference batch and on each cutoff of frequency
  tax_new_list = vector(mode="list", length=length(batch_ref_pool))
  method_chosen_list = vector(mode="list", length=length(batch_ref_pool))

  names(tax_new_list) = names(method_chosen_list) = batch_ref_pool

  # search on reference batch
  for (current_ref in 1:length(batch_ref_pool)){

    tax_new = matrix(ncol=ncol(tax_tab), nrow=nrow(tax_tab))
    method_chosen = NULL

    # search on each cutoff of frequency
    for (current_cutoff in 1:length(cut_list)){

      R2 = R2_initial
      tab_list = NULL

      if (length( cut_list[[current_cutoff]] ) == 0) next

      # benchmark -- on original data
      tab = tax_tab[, cut_list[[current_cutoff]], drop=F]
      tab_list[[1]] = tab

      index = which( apply(tab, 1, sum) > 0 )
      R2[1, 2] = adonis(formula = tab[index, ] ~ batchid[index])$aov.tab[1, 5]

      # do correction for all combination, record results and compute PERMANOVA R2
      for (current_method in 1:nrow(method1)){
        tab_list[[1+current_method]] = ConQuR_libsize(tax_tab=tab, batchid=batchid, covariates=covariates,
                                                      libsize_tune=libsize_tune,
                                                      batch_ref=as.character(batch_ref_pool[current_ref]),
                                                      logistic_lasso=method1[current_method, 'logistic'],
                                                      quantile_type=as.character(method1[current_method, 'quantile']),
                                                      simple_match=F,
                                                      lambda_quantile=as.character(method1[current_method, 'lambda']),
                                                      interplt=method1[current_method, 'interplt'],
                                                      delta=delta, taus=taus, num_core=num_core)

        index = which( apply(tab_list[[1+current_method]], 1, sum) > 0 )
        if (length(index) != 0  & length( unique(batchid[index]) ) != 1){
          R2[1+current_method, 2] = adonis(formula = tab_list[[1+current_method]][index, ] ~ batchid[index])$aov.tab[1, 5]
        }
      }

      if (method2 == T){
        tab_list[[1+nrow(method1)+1]] = ConQuR_libsize(tax_tab=tab, batchid=batchid, covariates=covariates,
                                                       libsize_tune=libsize_tune,
                                                       batch_ref=as.character(batch_ref_pool[current_ref]),
                                                       logistic_lasso=F,
                                                       quantile_type="standard",
                                                       simple_match=T,
                                                       lambda_quantile="2p/n",
                                                       interplt=F,
                                                       delta=delta, taus=taus, num_core=num_core)

        index = which( apply(tab_list[[1+nrow(method1)+1]], 1, sum) > 0 )
        if (length(index) != 0 & length( unique(batchid[index]) ) != 1){
          R2[1+nrow(method1)+1, 2] = adonis(formula = tab_list[[1+nrow(method1)+1]][index, ] ~ batchid[index])$aov.tab[1, 5]
        }
      }

      # find the optimal choice
      index_optimal = which.min(as.numeric( R2$batch_R2 ))

      for (col in 1:length( cut_list[[current_cutoff]] )){
        tax_new[, cut_list[[current_cutoff]][col]] = tab_list[[index_optimal]][, col]
      }

      method_chosen[current_cutoff] = R2$method[index_optimal]

    }

    # correct the remaining taxa by the default method -- standard logistic + standard quantile + no interpolation
    cut_list_remaining = which(prate <= frequencyL | prate > frequencyU)

    if (!length( cut_list_remaining ) == 0){

      tab = tax_tab[, cut_list_remaining, drop=F]

      tax_new_remaining = ConQuR_libsize(tax_tab=tab, batchid=batchid, covariates=covariates,
                                         libsize_tune=libsize_tune,
                                         batch_ref=as.character(batch_ref_pool[current_ref]),
                                         logistic_lasso=F,
                                         quantile_type="standard",
                                         simple_match=F,
                                         lambda_quantile="2p/n",
                                         interplt=F,
                                         delta=delta, taus=taus, num_core=num_core)

      for (col in 1:length( cut_list_remaining )){
        tax_new[, cut_list_remaining[col]] = tax_new_remaining[, col]
      }

    }

    tax_new_list[[current_ref]] = tax_new
    method_chosen_list[[current_ref]] = method_chosen

  }

  # determine the reference batch
  R2_ref = NULL
  for (current_ref in 1:length(batch_ref_pool)){
    index = which( apply(tax_new_list[[current_ref]], 1, sum) > 0 )
    R2_ref[current_ref] = adonis(formula = tax_new_list[[current_ref]][index, ] ~ batchid[index])$aov.tab[1, 5]
  }

  index_final = which.min(R2_ref)

  tax_final = tax_new_list[[index_final]]
  rownames(tax_final) = rownames(tax_tab)
  colnames(tax_final) = colnames(tax_tab)

  method_final = matrix(ncol=length(cut_list), nrow=7)
  rownames(method_final) = c("batch_ref", "original", "simple_match", "logistic_lasso", "quantile_type", "lambda", "interplt")
  colnames(method_final) = paste0(freq_grid[-length(freq_grid)], "-", freq_grid[-1])

  method_temp = method_chosen_list[[index_final]]
  for (ii in 1:length(cut_list)){
    if (is.na(method_temp[ii])){
      method_final[, ii] = rep(NA, 7)
    } else if (method_temp[ii] == "original"){
      method_final[, ii] = c(NA, T, rep(NA, 5))
    } else if (method_temp[ii] == "simple"){
      method_final[, ii] = c(names(method_chosen_list)[index_final], F, T, rep(NA, 4))
    } else{
      method_final[, ii] = c(names(method_chosen_list)[index_final], F, F, unlist( strsplit(method_temp[ii], '_') ) )
    }
  }

  return(list(tax_final=tax_final, method_final=method_final))

}




