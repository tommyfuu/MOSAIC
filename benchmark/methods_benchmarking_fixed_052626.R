#!/usr/bin/env Rscript
library(xtable) # table
# library(mixOmics)
library(sva) # ComBat
library(limma) # removeBatchEffect (LIMMA)
library(vegan) # RDA
library(MMUPHin)
library(FDboost)
library(doParallel)
library(dplyr)
library(readr)
library(tibble)
library(stringr)
library(mixOmics)

args = commandArgs(trailingOnly=TRUE)

## load ConQuR
current_path = getwd()
print("currnet path")
print(current_path)
if(grepl("mic_bc_benchmark/benchmark", current_path)){
    conqur_path = str_replace(current_path, "mic_bc_benchmark/benchmark", "mic_bc_benchmark/ConQuR")
} else {
    conqur_path = str_replace(current_path, "mic_bc_benchmark", "mic_bc_benchmark/ConQuR")
}
source(paste0(conqur_path, "/ConQuR_help_functions.R"))
source(paste0(conqur_path, "/ConQuR_main_tune.R"))
source(paste0(conqur_path, "/ConQuR_help_functions_libsize_old.R"))
source(paste0(conqur_path, "/ConQuR_main_tune_libsize_old.R"))
source(paste0(conqur_path, "/ConQuR_help_functions_rel.R"))
source(paste0(conqur_path, "/ConQuR_main_tune_rel.R"))
source(paste0(conqur_path, "/supporting_functions.R"))

check_metadata <- function(data, variables, no_missing = TRUE) {
  if(is.null(variables)) return(NULL)
  variables_absent <- setdiff(variables, colnames(data))
  if(length(variables_absent) > 0) {
    stop("Following variable(s) not present in data:\n",
         paste(variables_absent, collapse = ","))
  }
  variables_missing <- vapply(variables,
                              function(variable) {
                                any(is.na(data[[variable]]))
                              },
                              TRUE)
  if(any(variables_missing) & no_missing) {
    stop("Following variable(s) in data have missing values:\n",
         paste(variables[variables_missing], collapse = ","))
  }
  return(data[, variables, drop = FALSE])
}

check_batch <- function(x, min_n_batch = 2) {
  if(!is.factor(x)) {
    warning("Batch variable is not a factor as provided and will be converted ",
            "to one.")
    x <- as.factor(x)
  }
  if(nlevels(x) < min_n_batch)
    stop("Must have at least ", min_n_batch, " batches!")
  return(x)
}

construct_design <- function(data, with_intercept = TRUE) {
  if(is.null(data)) return(NULL)
  if(with_intercept)
    model.matrix(~ ., data = data)
  else
    model.matrix(~ . - 1, data = data)
}


run_methods <- function(data_mat_path, meta_data_path, output_root, batch_ref, dataset = "Dataset", covar = NULL, controlled = FALSE, Sam_id = 'Sam_id', transpose = FALSE, count = FALSE,
                        used_methods = c("combat_seq", "limma", "MMUPHin", 'ConQuR', 'ConQuR_libsize',
                                         'ruv3', 'bmc', 'bat', 'rbe', 'svd')) {

    print(output_root)
    sink_file_name = paste(output_root, "_runtime.txt", sep="")
    print(data_mat_path)
    output_dir = strsplit(output_root, "/")
    output_dir = paste(output_dir[[1]][1:(length(output_dir[[1]])-1)], collapse = "/")
    if(!dir.create(file.path(output_dir))){
        dir.create(file.path(output_dir))
    }

    if(transpose == TRUE) {
        count_data = read.table(data_mat_path, sep=",",header=T,check.names = F)
    }
    else {
        count_data = read.table(data_mat_path, sep=",",header=T,row.names=1,check.names = F)
    }
    print("check point1")

    metadata = as.data.frame(read_csv(meta_data_path))

    print(nrow(count_data))
    print(ncol(count_data))
    if(nrow(count_data)>300){
        num_core = 10
    }
    else {
        num_core = 2
    }

    print("check point2")

    # for ConQuR
    batchid <- factor(metadata[, dataset])

    # -------------------------------------------------------------------------
    # Helper: build covar_df from metadata (used by multiple methods)
    # -------------------------------------------------------------------------
    make_covar_df <- function() {
        if (is.null(covar)) return(NULL)
        if (length(covar) == 1) {
            matrix(as.numeric(factor(metadata[, covar])))
        } else {
            apply(metadata[, covar], 2, function(x) as.numeric(factor(x)))
        }
    }
    # define globally so ConQuR/MMUPHin etc. always have it (NULL when covar=NULL)
    covar_df <- make_covar_df()

    # -------------------------------------------------------------------------
    # Helper: determine whether a given block still needs to run.
    # Returns TRUE if any method in `methods_to_check` is absent from runtime.txt.
    # -------------------------------------------------------------------------
    needs_run <- function(methods_to_check) {
        if (length(methods_to_check) == 0) return(FALSE)
        rtf <- paste0(output_root, "_runtime.txt")
        if (!file.exists(rtf)) return(TRUE)
        runtime_lines <- readLines(rtf)
        all_done <- all(sapply(methods_to_check, function(m) any(grepl(m, runtime_lines))))
        return(!all_done)
    }

    ## CASE 1. count data (combat will not work)
    if(count == TRUE) {
        count_data.clr <- logratio.transfo(count_data+1, logratio = 'CLR')

        count_methods_all <- c("combat_seq", "limma", "MMUPHin",
                                "ConQuR", "ConQuR_libsize",
                                "Tune_ConQuR", "Tune_ConQuR_libsize",
                                "ruv3", "bmc", "bat", "rbe", "svd")
        count_methods_to_check <- intersect(used_methods, count_methods_all)
        run_count <- needs_run(count_methods_to_check)

        if (run_count == TRUE) {
            already_has_header <- file.exists(sink_file_name) && any(grepl("runtime documenting", readLines(sink_file_name)))
            if (!already_has_header) cat("runtime documenting...\n", file=sink_file_name, append=TRUE)

            ### 1.1 combat_seq (in place of combat)
            if ('combat_seq' %in% used_methods && !file.exists(paste0(output_root, "_combat_seq.csv"))) {
                print("in count mode, can only use combat seq")
                batch_info <- as.factor(setNames(as.character(metadata[, dataset]), metadata[[Sam_id]]))
                start_time <- Sys.time()
                if(is.null(covar)) {
                    count_data.combat_seq <- t(ComBat_seq(t(count_data), batch = batch_info))
                }
                else if(length(covar)==1){
                    covar_df = as.numeric(factor(metadata[, covar]))
                    count_data.combat_seq <- t(ComBat_seq(t(count_data), batch = batch_info, group = covar_df))
                }
                else{
                    covar_df = apply(metadata[, covar], 2, function(x) as.numeric(factor(x)))
                    count_data.combat_seq <- t(ComBat_seq(t(count_data), batch = batch_info, covar_mod = covar_df))
                }
                end_time <- Sys.time()
                cat(c("combat_seq runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.combat_seq, paste(output_root, "_combat_seq.csv", sep=""), row.names = TRUE)
            }

            ### 1.2 limma
            if ('limma' %in% used_methods && !file.exists(paste0(output_root, "_limma.csv"))) {
                batch_info <- as.factor(setNames(as.character(metadata[, dataset]), metadata[[Sam_id]]))
                start_time <- Sys.time()
                count_data.limma <-removeBatchEffect(t(count_data.clr), batch = batch_info)
                count_data.limma <- t(count_data.limma)
                if(is.null(covar)) {
                    count_data.limma <- t(removeBatchEffect(t(count_data.clr), batch = batch_info))
                }
                else if(length(covar)==1){
                    covar_df = matrix(as.numeric(factor(metadata[, covar])))
                    count_data.limma <- t(removeBatchEffect(t(count_data.clr), batch = batch_info, design = covar_df))
                }
                else{
                    covar_df = apply(metadata[, covar], 2, function(x) as.numeric(factor(x)))
                    count_data.limma <- t(removeBatchEffect(t(count_data.clr), batch = batch_info, design = covar_df))
                }
                ## the following line is count specific
                count_data.limma = clr(count_data.limma, inverse = TRUE)
                end_time <- Sys.time()
                cat(c("limma runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.limma, paste(output_root, "_limma.csv", sep=""), row.names = TRUE)
            }

            ### 1.3 MMUPHin - count/relab agnostic
            if ('MMUPHin' %in% used_methods && !file.exists(paste0(output_root, "_MMUPHin.csv"))) {
                metadata_mupphin <- metadata
                feature_abd = count_data
                row.names(metadata_mupphin) <- metadata[[Sam_id]]
                feature_abd = as.data.frame(t(feature_abd))
                colnames(feature_abd) <- rownames(metadata_mupphin)
                df_batch = check_metadata(metadata, dataset)
                batchmod <- construct_design(data = df_batch, with_intercept = FALSE)
                print(batchmod)
                var_batch <- check_batch(df_batch[[dataset]], min_n_batch = 2)
                n_batch <- nlevels(x = var_batch)
                print(n_batch)
                print(ncol(batchmod))
                start_time <- Sys.time()
                fit_adjust_batch <- adjust_batch(feature_abd = feature_abd,
                                        batch = dataset,
                                        covariates = covar,
                                        data = metadata_mupphin,
                                        control = list(verbose = FALSE))
                count_data.MMUPHin <- t(fit_adjust_batch$feature_abd_adj)
                end_time <- Sys.time()
                cat(c("MMUPHin runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.MMUPHin, paste(output_root, "_MMUPHin.csv", sep=""), row.names = TRUE)
            }

            ### 1.4 ConQuR (count version) - simple
            if ('ConQuR' %in% used_methods && !file.exists(paste0(output_root, "_ConQuR.csv"))) {
                print("running ConQuR")
                start_time <- Sys.time()
                count_data.ConQuR = ConQuR(tax_tab=count_data, batchid=batchid, covariates=covar_df, batch_ref=batch_ref, simple_match=is.null(covar_df), num_core = num_core)
                end_time <- Sys.time()
                cat(c("ConQuR runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.ConQuR, paste(output_root, "_ConQuR.csv", sep=""), row.names = TRUE)
            }

            ### 1.5 ConQuR_libsize (only runnable for counts) - simple
            if ('ConQuR_libsize' %in% used_methods && !file.exists(paste0(output_root, "_ConQuR_libsize.csv"))) {
                print("running ConQuR_libsize")
                start_time <- Sys.time()
                count_data.ConQuR_libsize = ConQuR_libsize(tax_tab=count_data, batchid=batchid, covariates=covar_df, batch_ref=batch_ref, simple_match=is.null(covar_df), num_core = num_core)
                end_time <- Sys.time()
                cat(c("ConQuR_libsize runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.ConQuR_libsize, paste(output_root, "_ConQuR_libsize.csv", sep=""), row.names = TRUE)
            }

            ### 1.6 Tune_ConQuR (count version)
            if ('Tune_ConQuR' %in% used_methods && !file.exists(paste0(output_root, "_Tune_ConQuR.csv"))) {
                start_time <- Sys.time()
                count_data.Tune_ConQuR = Tune_ConQuR(tax_tab=count_data, batchid=batchid, covariates=covar_df,
                            batch_ref_pool=batch_ref,
                            logistic_lasso_pool=c(T, F),
                            quantile_type_pool=c("standard", "lasso", "composite"),
                            simple_match_pool=c(T, F),
                            lambda_quantile_pool=c(NA, "2p/n", "2p/logn"),
                            interplt_pool=c(T, F),
                            frequencyL=0,
                            frequencyU=1,
                            cutoff=0.25,
                            num_core = num_core)
                end_time <- Sys.time()
                cat(c("Tune_ConQuR runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.Tune_ConQuR$tax_final, paste(output_root, "_Tune_ConQuR.csv", sep=""), row.names = TRUE)
            }

            ### 1.7 Tune_ConQuR_libsize (only runnable for counts) - simple
            if ('Tune_ConQuR_libsize' %in% used_methods && !file.exists(paste0(output_root, "_Tune_ConQuR_libsize.csv"))) {
                start_time <- Sys.time()
                count_data.Tune_ConQuR_libsize = Tune_ConQuR_libsize(tax_tab=count_data, batchid=batchid, covariates=covar_df,
                            batch_ref_pool=batch_ref,
                            logistic_lasso_pool=c(T, F),
                            quantile_type_pool=c("standard", "lasso", "composite"),
                            simple_match_pool=c(T, F),
                            lambda_quantile_pool=c(NA, "2p/n", "2p/logn"),
                            interplt_pool=c(T, F),
                            frequencyL=0,
                            frequencyU=1,
                            cutoff=0.25,
                            num_core = num_core)
                end_time <- Sys.time()
                cat(c("Tune_ConQuR_libsize runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.Tune_ConQuR_libsize$tax_final, paste(output_root, "_Tune_ConQuR_libsize.csv", sep=""), row.names = TRUE)
            }

            # ----------------------------------------------------------------
            # NEW METHODS (count case)  -- implementations follow MBECS v1.3.1
            # ----------------------------------------------------------------

            ### 1.8 ruv3 (RUV-III via ruv::RUVIII; mirrors mbecRUV3 in MBECS)
            # Negative controls: features non-significant in LM(feature ~ covariates)
            # at p > 0.05 (pseudo-negatives when no known controls are supplied).
            # Replicate structure: uses metadata$replicate when technical replicates
            # exist; falls back to biological-condition pseudo-replicates otherwise.
            if ('ruv3' %in% used_methods && !file.exists(paste0(output_root, "_ruv3.csv"))) {
                if (!requireNamespace("ruv", quietly = TRUE)) {
                    warning("ruv package not installed; skipping ruv3. Install with: install.packages('ruv')")
                    cat(c("ruv3 SKIPPED (ruv package not installed)"), file=sink_file_name, append=TRUE)
                    cat('\n', file=sink_file_name, append=TRUE)
                } else {
                    start_time <- Sys.time()
                    Y_ruv <- as.matrix(count_data.clr)
                    # Negative controls: features with no significant covariate association
                    ctl_ruv <- tryCatch({
                        if (!is.null(covar) && length(covar) > 0) {
                            pv <- sapply(seq_len(ncol(Y_ruv)), function(j) {
                                df_lm <- data.frame(y = Y_ruv[, j],
                                                    metadata[, covar, drop = FALSE])
                                fit <- lm(y ~ ., data = df_lm)
                                f   <- summary(fit)$fstatistic
                                if (is.null(f)) return(1.0)
                                pf(f[1], f[2], f[3], lower.tail = FALSE)
                            })
                            pv > 0.05
                        } else {
                            rep(TRUE, ncol(Y_ruv))
                        }
                    }, error = function(e) rep(TRUE, ncol(Y_ruv)))
                    if (sum(ctl_ruv) < 2) ctl_ruv <- rep(TRUE, ncol(Y_ruv))
                    # Replicate structure
                    if ('replicate' %in% colnames(metadata) &&
                        length(unique(metadata[['replicate']])) < nrow(metadata)) {
                        M_ruv <- ruv::replicate.matrix(metadata[['replicate']])
                    } else {
                        message("No technical replicates found; using biological-condition pseudo-replicates for RUV-III.")
                        if (!is.null(covar) && length(covar) > 0) {
                            cond_var <- factor(metadata[, covar[1]])
                            M_ruv <- model.matrix(~ cond_var - 1)
                        } else {
                            M_ruv <- matrix(1, nrow(Y_ruv), 1)
                        }
                    }
                    k_ruv <- max(1L, nlevels(batchid) - 1L)
                    max_k <- max(1L, nrow(Y_ruv) - ncol(M_ruv) - 1L)
                    k_ruv <- min(k_ruv, max_k)
                    Y_ruv[!is.finite(Y_ruv)] <- 0
                    tryCatch({
                        ruv3_out <- ruv::RUVIII(Y_ruv, M_ruv, ctl_ruv, k = k_ruv)
                        ruv3_corrected <- if (is.list(ruv3_out)) ruv3_out$newY else ruv3_out
                        count_data.ruv3 <- clr(ruv3_corrected, inverse = TRUE)
                        rownames(count_data.ruv3) <- rownames(count_data)
                        end_time <- Sys.time()
                        cat(c("ruv3 runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                        cat('\n', file=sink_file_name, append=TRUE)
                        write.csv(count_data.ruv3, paste0(output_root, "_ruv3.csv"), row.names = TRUE)
                    }, error = function(e) {
                        end_time <- Sys.time()
                        warning(paste("ruv3 SVD failed, writing uncorrected data:", conditionMessage(e)))
                        cat(c("ruv3 runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                        cat('\n', file=sink_file_name, append=TRUE)
                        cat("ruv3 note: SVD failed; output is uncorrected (fallback to original data)\n", file=sink_file_name, append=TRUE)
                        write.csv(count_data, paste0(output_root, "_ruv3.csv"), row.names = TRUE)
                    })
                }
            }

            ### 1.9 bmc (Batch Mean Centering; mirrors mbecBMC in MBECS)
            # Centers each batch to zero mean per feature using scale(center=TRUE,
            # scale=FALSE). No global mean is added back (MBECS convention).
            if ('bmc' %in% used_methods && !file.exists(paste0(output_root, "_bmc.csv"))) {
                start_time <- Sys.time()
                count_data.bmc <- as.matrix(count_data.clr)
                for (b in levels(batchid)) {
                    idx <- which(as.character(batchid) == b)
                    if (length(idx) == 0) next
                    count_data.bmc[idx, ] <- scale(
                        count_data.bmc[idx, , drop = FALSE],
                        center = TRUE, scale = FALSE
                    )
                }
                count_data.bmc <- clr(count_data.bmc, inverse = TRUE)
                rownames(count_data.bmc) <- rownames(count_data)
                end_time <- Sys.time()
                cat(c("bmc runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.bmc, paste0(output_root, "_bmc.csv"), row.names = TRUE)
            }

            ### 1.10 bat (ComBat on CLR; mirrors mbecBat in MBECS)
            # Applies sva::ComBat (non-parametric) to CLR-transformed counts.
            # Distinct from combat_seq which uses ComBat_seq on raw counts.
            if ('bat' %in% used_methods && !file.exists(paste0(output_root, "_bat.csv"))) {
                start_time <- Sys.time()
                if (is.null(covar)) {
                    bat_corrected <- sva::ComBat(
                        dat         = t(count_data.clr),
                        batch       = batchid,
                        mod         = matrix(1, nrow(count_data), 1),
                        par.prior   = FALSE,
                        prior.plots = FALSE
                    )
                } else {
                    bat_mod <- stats::model.matrix(~ as.factor(metadata[, covar[1]]))
                    bat_corrected <- sva::ComBat(
                        dat         = t(count_data.clr),
                        batch       = batchid,
                        mod         = bat_mod,
                        par.prior   = FALSE,
                        prior.plots = FALSE
                    )
                }
                count_data.bat <- clr(t(bat_corrected), inverse = TRUE)
                rownames(count_data.bat) <- rownames(count_data)
                end_time <- Sys.time()
                cat(c("bat runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.bat, paste0(output_root, "_bat.csv"), row.names = TRUE)
            }

            ### 1.11 rbe (removeBatchEffect on CLR; mirrors mbecRBE in MBECS)
            # Uses limma::removeBatchEffect on CLR data with an explicit intercept
            # design when no covariates are present, and model.matrix(~covar) otherwise.
            # Differs from the 'limma' block above in design matrix construction.
            if ('rbe' %in% used_methods && !file.exists(paste0(output_root, "_rbe.csv"))) {
                batch_info_rbe <- as.factor(setNames(as.character(metadata[, dataset]), metadata[[Sam_id]]))
                start_time <- Sys.time()
                if (is.null(covar)) {
                    rbe_design    <- matrix(1, nrow(count_data), 1)
                    rbe_corrected <- t(removeBatchEffect(
                        t(count_data.clr), batch = batch_info_rbe, design = rbe_design))
                } else {
                    rbe_design    <- stats::model.matrix(~ as.factor(metadata[, covar[1]]))
                    rbe_corrected <- t(removeBatchEffect(
                        t(count_data.clr), batch = batch_info_rbe, design = rbe_design))
                }
                count_data.rbe <- clr(rbe_corrected, inverse = TRUE)
                rownames(count_data.rbe) <- rownames(count_data)
                end_time <- Sys.time()
                cat(c("rbe runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.rbe, paste0(output_root, "_rbe.csv"), row.names = TRUE)
            }

            ### 1.12 svd (deflate first PC of CLR cross-product; mirrors mbecSVD in MBECS)
            # Centers and scales features, computes the p×p cross-product matrix,
            # runs SVD, deflates the dominant singular component (assumed = batch),
            # then rescales back to original mean and SD.
            if ('svd' %in% used_methods && !file.exists(paste0(output_root, "_svd.csv"))) {
                start_time <- Sys.time()
                Y_svd      <- as.matrix(count_data.clr)
                tmp.sd     <- apply(Y_svd, 2, sd)
                tmp.sd[tmp.sd == 0] <- 1
                tmp.mean   <- colMeans(Y_svd)
                Y_scaled   <- scale(Y_svd, center = TRUE, scale = TRUE)
                Y_cross    <- crossprod(Y_scaled)           # p x p
                svd.res    <- svd(Y_cross)
                svd.u1     <- svd.res$u[, 1]
                tmp.t1     <- Y_scaled %*% svd.u1 / drop(sqrt(crossprod(svd.u1)))
                tmp.c1     <- crossprod(Y_scaled, tmp.t1) / drop(crossprod(tmp.t1))
                Y_defl     <- Y_scaled - tmp.t1 %*% t(tmp.c1)
                corrected_svd <- sweep(sweep(Y_defl, 2, tmp.sd, `*`), 2, tmp.mean, `+`)
                count_data.svd <- clr(corrected_svd, inverse = TRUE)
                rownames(count_data.svd) <- rownames(count_data)
                end_time <- Sys.time()
                cat(c("svd runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.svd, paste0(output_root, "_svd.csv"), row.names = TRUE)
            }
        }
        else {
            print("already done")
        }
    }

    ## CASE 2. relative abundance data
    else {
        count_data.clr <- logratio.transfo(count_data+min(count_data[count_data>0]), logratio = 'CLR')
        print("check whether there are NAs in count_data.clr")
        print(sum(is.na(count_data.clr)))
        print("check whether there are Inf in count_data.clr")
        print(sum(is.infinite(count_data.clr)))

        relab_methods_all <- c("combat", "limma", "MMUPHin",
                                "ConQuR_rel", "Tune_ConQuR_rel",
                                "ruv3", "bmc", "bat", "rbe", "svd")
        relab_methods_to_check <- intersect(used_methods, relab_methods_all)
        run_relab <- needs_run(relab_methods_to_check)

        print("run_relab")
        print(run_relab)
        if (run_relab == TRUE) {
            already_has_header <- file.exists(sink_file_name) && any(grepl("runtime documenting", readLines(sink_file_name)))
            if (!already_has_header) cat("runtime documenting...\n", file=sink_file_name, append=TRUE)

            ### 2.1 combat
            if ('combat' %in% used_methods && !file.exists(paste0(output_root, "_combat.csv"))) {
                print("in relative abundance mode, can only use combat")
                batch_info <- as.factor(setNames(as.character(metadata[, dataset]), metadata[[Sam_id]]))
                start_time <- Sys.time()
                if(is.null(covar)) {
                    count_data.combat <- t(ComBat(t(count_data.clr), batch = batch_info, par.prior=FALSE, mod=NULL))
                }
                else if(length(covar)==1){
                    covar_df = as.numeric(factor(metadata[, covar]))
                    count_data.combat <- t(ComBat(t(count_data.clr), batch = batch_info, par.prior=FALSE, mod=covar_df))
                }
                else{
                    print(colnames(metadata))
                    covar_df = apply(metadata[, covar], 2, function(x) as.numeric(factor(x)))
                    print(covar_df)
                    count_data.combat <- t(ComBat(t(count_data.clr), batch = batch_info, par.prior=FALSE, mod=covar_df))
                }
                end_time <- Sys.time()
                cat(c("combat runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.combat, paste(output_root, "_combat.csv", sep=""), row.names = TRUE)
            }

            ### 2.2 limma
            if ('limma' %in% used_methods && !file.exists(paste0(output_root, "_limma.csv"))) {
                batch_info <- as.factor(setNames(as.character(metadata[, dataset]), metadata[[Sam_id]]))
                start_time <- Sys.time()
                if(is.null(covar)) {
                    count_data.limma <- t(removeBatchEffect(t(count_data.clr), batch = batch_info))
                }
                else if(length(covar)==1){
                    count_data.limma <- t(removeBatchEffect(t(count_data.clr), batch = batch_info))
                }
                else{
                    covar_df = apply(metadata[, covar], 2, function(x) as.numeric(factor(x)))
                    count_data.limma <- t(removeBatchEffect(t(count_data.clr), batch = batch_info, design = covar_df))
                }
                end_time <- Sys.time()
                cat(c("limma runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.limma, paste(output_root, "_limma.csv", sep=""), row.names = TRUE)
            }

            ### 2.3 MMUPHin - count/relab agnostic
            if ('MMUPHin' %in% used_methods && !file.exists(paste0(output_root, "_MMUPHin.csv"))) {
                metadata_mupphin <- metadata
                feature_abd = count_data
                row.names(metadata_mupphin) <- metadata[[Sam_id]]
                feature_abd = as.data.frame(t(feature_abd))
                colnames(feature_abd) <- rownames(metadata_mupphin)
                start_time <- Sys.time()
                fit_adjust_batch <- adjust_batch(feature_abd = feature_abd,
                                        batch = dataset,
                                        covariates = covar,
                                        data = metadata_mupphin,
                                        control = list(verbose = FALSE))
                count_data.MMUPHin <- t(fit_adjust_batch$feature_abd_adj)
                end_time <- Sys.time()
                cat(c("MMUPHin runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.MMUPHin, paste(output_root, "_MMUPHin.csv", sep=""), row.names = TRUE)
            }

            ### 2.4 ConQuR_rel (rel_ab version) - simple
            if ('ConQuR_rel' %in% used_methods && !file.exists(paste0(output_root, "_ConQuR_rel.csv"))) {
                start_time <- Sys.time()
                count_data.ConQuR_rel = ConQuR_rel(tax_tab=count_data, batchid=batchid, covariates=covar_df, batch_ref=batch_ref, simple_match=is.null(covar_df), num_core = num_core)
                end_time <- Sys.time()
                cat(c("ConQuR_rel runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.ConQuR_rel, paste(output_root, "_ConQuR_rel.csv", sep=""), row.names = TRUE)
            }

            ### 2.5 Tune_ConQuR_rel (rel_ab version)
            if ('Tune_ConQuR_rel' %in% used_methods && !file.exists(paste0(output_root, "_Tune_ConQuR_rel.csv"))) {
                start_time <- Sys.time()
                count_data.Tune_ConQuR_rel = Tune_ConQuR_rel(tax_tab=count_data, batchid=batchid, covariates=covar_df,
                            batch_ref_pool=batch_ref,
                            logistic_lasso_pool=c(T, F),
                            quantile_type_pool=c("standard", "lasso", "composite"),
                            simple_match_pool=c(T, F),
                            lambda_quantile_pool=c(NA, "2p/n", "2p/logn"),
                            interplt_pool=c(T, F),
                            frequencyL=0,
                            frequencyU=1,
                            cutoff=0.25, num_core = num_core)
                end_time <- Sys.time()
                cat(c("Tune_ConQuR_rel runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.Tune_ConQuR_rel$tax_final, paste(output_root, "_Tune_ConQuR_rel.csv", sep=""), row.names = TRUE)
            }

            # ----------------------------------------------------------------
            # NEW METHODS (relab case)  -- implementations follow MBECS v1.3.1
            # All CLR-based methods return CLR-space output (same as combat / limma).
            # ----------------------------------------------------------------

            ### 2.6 ruv3 (RUV-III on CLR; returns CLR-space output)
            if ('ruv3' %in% used_methods && !file.exists(paste0(output_root, "_ruv3.csv"))) {
                if (!requireNamespace("ruv", quietly = TRUE)) {
                    warning("ruv package not installed; skipping ruv3. Install with: install.packages('ruv')")
                    cat(c("ruv3 SKIPPED (ruv package not installed)"), file=sink_file_name, append=TRUE)
                    cat('\n', file=sink_file_name, append=TRUE)
                } else {
                    start_time <- Sys.time()
                    Y_ruv <- as.matrix(count_data.clr)
                    ctl_ruv <- tryCatch({
                        if (!is.null(covar) && length(covar) > 0) {
                            pv <- sapply(seq_len(ncol(Y_ruv)), function(j) {
                                df_lm <- data.frame(y = Y_ruv[, j],
                                                    metadata[, covar, drop = FALSE])
                                fit <- lm(y ~ ., data = df_lm)
                                f   <- summary(fit)$fstatistic
                                if (is.null(f)) return(1.0)
                                pf(f[1], f[2], f[3], lower.tail = FALSE)
                            })
                            pv > 0.05
                        } else {
                            rep(TRUE, ncol(Y_ruv))
                        }
                    }, error = function(e) rep(TRUE, ncol(Y_ruv)))
                    if (sum(ctl_ruv) < 2) ctl_ruv <- rep(TRUE, ncol(Y_ruv))
                    if ('replicate' %in% colnames(metadata) &&
                        length(unique(metadata[['replicate']])) < nrow(metadata)) {
                        M_ruv <- ruv::replicate.matrix(metadata[['replicate']])
                    } else {
                        message("No technical replicates found; using biological-condition pseudo-replicates for RUV-III.")
                        if (!is.null(covar) && length(covar) > 0) {
                            cond_var <- factor(metadata[, covar[1]])
                            M_ruv <- model.matrix(~ cond_var - 1)
                        } else {
                            M_ruv <- matrix(1, nrow(Y_ruv), 1)
                        }
                    }
                    k_ruv <- max(1L, nlevels(batchid) - 1L)
                    max_k <- max(1L, nrow(Y_ruv) - ncol(M_ruv) - 1L)
                    k_ruv <- min(k_ruv, max_k)
                    Y_ruv[!is.finite(Y_ruv)] <- 0
                    tryCatch({
                        ruv3_out <- ruv::RUVIII(Y_ruv, M_ruv, ctl_ruv, k = k_ruv)
                        ruv3_corrected <- if (is.list(ruv3_out)) ruv3_out$newY else ruv3_out
                        count_data.ruv3 <- ruv3_corrected  # stay in CLR space for relab
                        rownames(count_data.ruv3) <- rownames(count_data)
                        end_time <- Sys.time()
                        cat(c("ruv3 runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                        cat('\n', file=sink_file_name, append=TRUE)
                        write.csv(count_data.ruv3, paste0(output_root, "_ruv3.csv"), row.names = TRUE)
                    }, error = function(e) {
                        end_time <- Sys.time()
                        warning(paste("ruv3 SVD failed, writing uncorrected data:", conditionMessage(e)))
                        cat(c("ruv3 runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                        cat('\n', file=sink_file_name, append=TRUE)
                        cat("ruv3 note: SVD failed; output is uncorrected (fallback to original data)\n", file=sink_file_name, append=TRUE)
                        write.csv(count_data.clr, paste0(output_root, "_ruv3.csv"), row.names = TRUE)
                    })
                }
            }

            ### 2.7 bmc (Batch Mean Centering on CLR; returns CLR-space output)
            if ('bmc' %in% used_methods && !file.exists(paste0(output_root, "_bmc.csv"))) {
                start_time <- Sys.time()
                count_data.bmc <- as.matrix(count_data.clr)
                for (b in levels(batchid)) {
                    idx <- which(as.character(batchid) == b)
                    if (length(idx) == 0) next
                    count_data.bmc[idx, ] <- scale(
                        count_data.bmc[idx, , drop = FALSE],
                        center = TRUE, scale = FALSE
                    )
                }
                rownames(count_data.bmc) <- rownames(count_data)
                end_time <- Sys.time()
                cat(c("bmc runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.bmc, paste0(output_root, "_bmc.csv"), row.names = TRUE)
            }

            ### 2.8 bat (ComBat on CLR; returns CLR-space output)
            if ('bat' %in% used_methods && !file.exists(paste0(output_root, "_bat.csv"))) {
                start_time <- Sys.time()
                if (is.null(covar)) {
                    bat_corrected <- sva::ComBat(
                        dat         = t(count_data.clr),
                        batch       = batchid,
                        mod         = matrix(1, nrow(count_data), 1),
                        par.prior   = FALSE,
                        prior.plots = FALSE
                    )
                } else {
                    bat_mod <- stats::model.matrix(~ as.factor(metadata[, covar[1]]))
                    bat_corrected <- sva::ComBat(
                        dat         = t(count_data.clr),
                        batch       = batchid,
                        mod         = bat_mod,
                        par.prior   = FALSE,
                        prior.plots = FALSE
                    )
                }
                count_data.bat <- t(bat_corrected)   # stay in CLR space for relab
                rownames(count_data.bat) <- rownames(count_data)
                end_time <- Sys.time()
                cat(c("bat runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.bat, paste0(output_root, "_bat.csv"), row.names = TRUE)
            }

            ### 2.9 rbe (removeBatchEffect on CLR; returns CLR-space output)
            if ('rbe' %in% used_methods && !file.exists(paste0(output_root, "_rbe.csv"))) {
                batch_info_rbe <- as.factor(setNames(as.character(metadata[, dataset]), metadata[[Sam_id]]))
                start_time <- Sys.time()
                if (is.null(covar)) {
                    rbe_design    <- matrix(1, nrow(count_data), 1)
                    count_data.rbe <- t(removeBatchEffect(
                        t(count_data.clr), batch = batch_info_rbe, design = rbe_design))
                } else {
                    rbe_design    <- stats::model.matrix(~ as.factor(metadata[, covar[1]]))
                    count_data.rbe <- t(removeBatchEffect(
                        t(count_data.clr), batch = batch_info_rbe, design = rbe_design))
                }
                rownames(count_data.rbe) <- rownames(count_data)
                end_time <- Sys.time()
                cat(c("rbe runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.rbe, paste0(output_root, "_rbe.csv"), row.names = TRUE)
            }

            ### 2.10 svd (deflate first PC of CLR cross-product; returns CLR-space output)
            if ('svd' %in% used_methods && !file.exists(paste0(output_root, "_svd.csv"))) {
                start_time <- Sys.time()
                Y_svd      <- as.matrix(count_data.clr)
                tmp.sd     <- apply(Y_svd, 2, sd)
                tmp.sd[tmp.sd == 0] <- 1
                tmp.mean   <- colMeans(Y_svd)
                Y_scaled   <- scale(Y_svd, center = TRUE, scale = TRUE)
                Y_cross    <- crossprod(Y_scaled)           # p x p
                svd.res    <- svd(Y_cross)
                svd.u1     <- svd.res$u[, 1]
                tmp.t1     <- Y_scaled %*% svd.u1 / drop(sqrt(crossprod(svd.u1)))
                tmp.c1     <- crossprod(Y_scaled, tmp.t1) / drop(crossprod(tmp.t1))
                Y_defl     <- Y_scaled - tmp.t1 %*% t(tmp.c1)
                count_data.svd <- sweep(sweep(Y_defl, 2, tmp.sd, `*`), 2, tmp.mean, `+`)
                rownames(count_data.svd) <- rownames(count_data)
                end_time <- Sys.time()
                cat(c("svd runtime", toString(difftime(end_time, start_time, unit="secs")), "seconds"), file=sink_file_name, append=TRUE)
                cat('\n', file=sink_file_name, append=TRUE)
                write.csv(count_data.svd, paste0(output_root, "_svd.csv"), row.names = TRUE)
            }

        }
        else {
            print("already done")
        }
    }

}


library(yaml)
current_path = getwd()
if (grepl("/benchmark", current_path)){
    Snakefile_path <- "../Snakefile"
} else{
    Snakefile_path <- "./Snakefile"
}

conn <- file(Snakefile_path,open="r")
linn <-readLines(conn)
for (i in 1:length(linn)){
    if (grepl("configfile", linn[i])){
        configfile_path = strsplit(linn[i], " ")[[1]][2]
        configfile_path = gsub("'", "", configfile_path)
        configfile_path = gsub("\"", "", configfile_path)
        break
    }
}
close(conn)

if (grepl("/benchmark", current_path)){
    config_object <- yaml.load_file(paste0("../", configfile_path))
} else{
    config_object <- yaml.load_file(configfile_path)
}

dataset_name <- config_object$dataset_name
batch_ref <- config_object$batch_ref
covar <- config_object$covar
used_methods <- config_object$used_R_methods
current_root <- paste0(config_object$src, '/cleaned_data/', dataset_name, '/', dataset_name)
output_root <- paste0(config_object$post_integration_outputs, '/', dataset_name, '/', dataset_name)

if(!dir.exists(paste0(config_object$post_integration_outputs, '/', dataset_name, '/'))){
    print(paste0(config_object$post_integration_outputs, '/', dataset_name, '/'))
    dir.create(paste0(config_object$post_integration_outputs, '/', dataset_name, '/'))
}

if(config_object$datatype=='count'){
    count = TRUE
} else {
    count = FALSE
}
print(used_methods)

if (grepl('microbiomeHD', dataset_name)){
    run_methods(paste0(current_root, "_count_data.csv"),
        paste0(current_root, "_meta_data.csv"),
        output_root,
        dataset = "Dataset",
        batch_ref = batch_ref,
        covar = covar,
        count = count,
        used_methods =  used_methods
    )
} else if (grepl('CMD', dataset_name)){
    run_methods(paste0(current_root, "_count_data.csv"),
        paste0(current_root, "_meta_data.csv"),
        output_root,
        dataset = "study_name",
        batch_ref = batch_ref,
        covar = covar,
        count = count,
        used_methods = used_methods
    )
}

# ── Nocovar second pass for leakage-free RF evaluation ──────────────────────
# Runs only when run_nocovar_for_rf: True in config (wml configs).
# Identical to the main run but with covar = NULL, saving to a separate folder.
if (!is.null(config_object$run_nocovar_for_rf) && isTRUE(config_object$run_nocovar_for_rf)) {
    nocovar_output_root <- paste0(config_object$post_integration_outputs_nocovar, '/', dataset_name, '/', dataset_name)
    if (!dir.exists(paste0(config_object$post_integration_outputs_nocovar, '/', dataset_name, '/'))) {
        dir.create(paste0(config_object$post_integration_outputs_nocovar, '/', dataset_name, '/'), recursive = TRUE)
    }
    print("Running nocovar pass for leakage-free RF evaluation...")
    if (grepl('microbiomeHD', dataset_name)){
        run_methods(paste0(current_root, "_count_data.csv"),
            paste0(current_root, "_meta_data.csv"),
            nocovar_output_root,
            dataset = "Dataset",
            batch_ref = batch_ref,
            covar = NULL,
            count = count,
            used_methods = used_methods
        )
    } else if (grepl('CMD', dataset_name)){
        run_methods(paste0(current_root, "_count_data.csv"),
            paste0(current_root, "_meta_data.csv"),
            nocovar_output_root,
            dataset = "study_name",
            batch_ref = batch_ref,
            covar = NULL,
            count = count,
            used_methods = used_methods
        )
    }
    print("Nocovar pass complete.")
}
