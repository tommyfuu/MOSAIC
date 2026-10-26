#!/usr/bin/env Rscript
# methods_benchmarking_sim_condsweep.R
#
# Simulation benchmarking for the condsweep scenario (OR=1.25 by default, cond_effect in
# {0.005,0.01,0.02,0.05,0.1}, batch_effect=0.5), count + relab, no/yes relation,
# covar and nocovar passes. Same run_methods() as methods_benchmarking_sim_052626.R.
#
# Usage: Rscript methods_benchmarking_sim_condsweep.R <iter> [combo] [OR]
#   combo: count_no, count_yes, relab_no, relab_yes or count_no_nocovar_only (default: all 4)

source(Find(file.exists, c("paths.R", "../paths.R", "../../paths.R")))
library(xtable)
library(sva)
library(limma)
library(vegan)
library(MMUPHin)
library(FDboost)
library(doParallel)
library(dplyr)
library(readr)
library(tibble)
library(stringr)
library(mixOmics)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
print(args)
GLOBAL_ITER <- if (length(args) >= 1) as.integer(args[1]) else {
    print("No iteration argument supplied; defaulting to iter = 1")
    1L
}
set.seed(GLOBAL_ITER)

# 2nd arg: datatype/relation combo (default "all")
COMBO_SELECT <- if (length(args) >= 2) args[2] else "all"
print(paste("COMBO_SELECT:", COMBO_SELECT))
# 3rd arg: odds ratio between batch and condition (default 1.25; output files are tagged with it)
CONDSWEEP_OR <- if (length(args) >= 3) as.numeric(args[3]) else 1.25
print(paste("CONDSWEEP_OR:", CONDSWEEP_OR))

# ConQuR
current_path <- getwd()
print("current path"); print(current_path)
conqur_path <- file.path(REPO_ROOT, "ConQuR")
source(paste0(conqur_path, "/ConQuR_help_functions.R"))
source(paste0(conqur_path, "/ConQuR_main_tune.R"))
source(paste0(conqur_path, "/ConQuR_help_functions_libsize_old.R"))
source(paste0(conqur_path, "/ConQuR_main_tune_libsize_old.R"))
source(paste0(conqur_path, "/ConQuR_help_functions_rel.R"))
source(paste0(conqur_path, "/ConQuR_main_tune_rel.R"))
source(paste0(conqur_path, "/supporting_functions.R"))

# Main benchmarking function
run_methods <- function(data_mat_path, meta_data_path, output_root,
                        batch_ref, dataset = "batchid", covar = c("cond"),
                        Sam_id = "Sam_id", transpose = TRUE, count = TRUE,
                        used_methods = c("combat_seq", "limma", "MMUPHin",
                                         "ConQuR", "ConQuR_libsize",
                                         "ruv3", "bmc", "bat", "rbe", "svd")) {

    print(output_root)
    sink_file_name <- paste0(output_root, "_runtime.txt")

    # create output directory if needed
    output_dir <- paste(head(strsplit(output_root, "/")[[1]], -1), collapse = "/")
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    # load data
    count_data <- if (transpose) {
        read.table(data_mat_path, sep = ",", header = TRUE, check.names = FALSE)
    } else {
        read.table(data_mat_path, sep = ",", header = TRUE, row.names = 1, check.names = FALSE)
    }
    metadata <- as.data.frame(read_csv(meta_data_path))
    num_core <- if (nrow(count_data) > 300) 10 else 2

    # batch factor + covariate matrix (needed by ConQuR and new methods)
    batchid  <- factor(metadata[, dataset])
    covar_df <- if (!is.null(covar) && length(covar) > 0) {
        if (length(covar) == 1) matrix(as.numeric(factor(metadata[, covar])))
        else apply(metadata[, covar], 2, function(x) as.numeric(factor(x)))
    } else NULL

    # needs_run: TRUE if any requested method is absent from runtime.txt
    needs_run <- function(methods_to_check) {
        if (length(methods_to_check) == 0) return(FALSE)
        rtf <- paste0(output_root, "_runtime.txt")
        if (!file.exists(rtf)) return(TRUE)
        lines <- readLines(rtf)
        !all(sapply(methods_to_check, function(m) any(grepl(m, lines))))
    }

    # write runtime header only once
    ensure_header <- function() {
        if (!file.exists(sink_file_name) ||
            !any(grepl("runtime documenting", readLines(sink_file_name))))
            cat("runtime documenting...\n", file = sink_file_name, append = TRUE)
    }

    ## CASE 1 - count data
    if (count) {
        count_data.clr <- logratio.transfo(count_data + 1, logratio = "CLR")

        count_methods_all      <- c("combat_seq", "limma", "MMUPHin",
                                    "ConQuR", "ConQuR_libsize",
                                    "Tune_ConQuR", "Tune_ConQuR_libsize",
                                    "bmc", "bat", "rbe", "svd")
        # ruv3 is not in this list: it runs below based only on whether its CSV exists
        count_methods_to_check <- intersect(used_methods, count_methods_all)

        if (needs_run(count_methods_to_check)) {
            ensure_header()

            ### 1.1 combat_seq
            if ("combat_seq" %in% used_methods &&
                !file.exists(paste0(output_root, "_combat_seq.csv"))) {
                batch_info <- as.factor(setNames(as.character(metadata[, dataset]), metadata[[Sam_id]]))
                start_time <- Sys.time()
                count_data.combat_seq <- if (is.null(covar)) {
                    t(ComBat_seq(t(count_data), batch = batch_info))
                } else if (length(covar) == 1) {
                    t(ComBat_seq(t(count_data), batch = batch_info,
                                  group = as.numeric(factor(metadata[, covar]))))
                } else {
                    t(ComBat_seq(t(count_data), batch = batch_info,
                                  covar_mod = apply(metadata[, covar], 2, function(x) as.numeric(factor(x)))))
                }
                end_time <- Sys.time()
                cat(c("combat_seq runtime", toString(difftime(end_time, start_time, unit = "secs")), "seconds"),
                    file = sink_file_name, append = TRUE); cat("\n", file = sink_file_name, append = TRUE)
                write.csv(count_data.combat_seq, paste0(output_root, "_combat_seq.csv"), row.names = TRUE)
            }

            ### 1.2 limma
            if ("limma" %in% used_methods &&
                !file.exists(paste0(output_root, "_limma.csv"))) {
                batch_info <- as.factor(setNames(as.character(metadata[, dataset]), metadata[[Sam_id]]))
                start_time <- Sys.time()
                count_data.limma <- if (is.null(covar)) {
                    t(removeBatchEffect(t(count_data.clr), batch = batch_info))
                } else if (length(covar) == 1) {
                    t(removeBatchEffect(t(count_data.clr), batch = batch_info,
                                         design = matrix(as.numeric(factor(metadata[, covar])))))
                } else {
                    t(removeBatchEffect(t(count_data.clr), batch = batch_info, design = covar_df))
                }
                count_data.limma <- clr(count_data.limma, inverse = TRUE)
                end_time <- Sys.time()
                cat(c("limma runtime", toString(difftime(end_time, start_time, unit = "secs")), "seconds"),
                    file = sink_file_name, append = TRUE); cat("\n", file = sink_file_name, append = TRUE)
                write.csv(count_data.limma, paste0(output_root, "_limma.csv"), row.names = TRUE)
            }

            ### 1.3 MMUPHin
            if ("MMUPHin" %in% used_methods &&
                !file.exists(paste0(output_root, "_MMUPHin.csv"))) {
                meta_mup <- metadata; row.names(meta_mup) <- metadata[[Sam_id]]
                feat_abd <- as.data.frame(t(count_data)); colnames(feat_abd) <- rownames(meta_mup)
                start_time <- Sys.time()
                fit <- adjust_batch(feature_abd = feat_abd, batch = dataset,
                                     covariates = covar, data = meta_mup,
                                     control = list(verbose = FALSE))
                count_data.MMUPHin <- t(fit$feature_abd_adj)
                end_time <- Sys.time()
                cat(c("MMUPHin runtime", toString(difftime(end_time, start_time, unit = "secs")), "seconds"),
                    file = sink_file_name, append = TRUE); cat("\n", file = sink_file_name, append = TRUE)
                write.csv(count_data.MMUPHin, paste0(output_root, "_MMUPHin.csv"), row.names = TRUE)
            }

            ### 1.4 ConQuR
            if ("ConQuR" %in% used_methods &&
                !file.exists(paste0(output_root, "_ConQuR.csv"))) {
                start_time <- Sys.time()
                count_data.ConQuR <- ConQuR(tax_tab = count_data, batchid = batchid,
                                             covariates = covar_df, batch_ref = batch_ref, simple_match = is.null(covar_df), num_core = num_core)
                end_time <- Sys.time()
                cat(c("ConQuR runtime", toString(difftime(end_time, start_time, unit = "secs")), "seconds"),
                    file = sink_file_name, append = TRUE); cat("\n", file = sink_file_name, append = TRUE)
                write.csv(count_data.ConQuR, paste0(output_root, "_ConQuR.csv"), row.names = TRUE)
            }

            ### 1.5 ConQuR_libsize
            if ("ConQuR_libsize" %in% used_methods &&
                !file.exists(paste0(output_root, "_ConQuR_libsize.csv"))) {
                start_time <- Sys.time()
                count_data.ConQuR_libsize <- ConQuR_libsize(tax_tab = count_data, batchid = batchid,
                                                             covariates = covar_df, batch_ref = batch_ref, simple_match = is.null(covar_df), num_core = num_core)
                end_time <- Sys.time()
                cat(c("ConQuR_libsize runtime", toString(difftime(end_time, start_time, unit = "secs")), "seconds"),
                    file = sink_file_name, append = TRUE); cat("\n", file = sink_file_name, append = TRUE)
                write.csv(count_data.ConQuR_libsize, paste0(output_root, "_ConQuR_libsize.csv"), row.names = TRUE)
            }

            ### 1.6 Tune_ConQuR
            if ("Tune_ConQuR" %in% used_methods &&
                !file.exists(paste0(output_root, "_Tune_ConQuR.csv"))) {
                start_time <- Sys.time()
                count_data.Tune_ConQuR <- Tune_ConQuR(
                    tax_tab = count_data, batchid = batchid, covariates = covar_df,
                    batch_ref_pool = batch_ref, logistic_lasso_pool = c(T, F),
                    quantile_type_pool = c("standard","lasso","composite"),
                    simple_match_pool = c(T, F), lambda_quantile_pool = c(NA,"2p/n","2p/logn"),
                    interplt_pool = c(T, F), frequencyL = 0, frequencyU = 1,
                    cutoff = 0.25, num_core = num_core)
                end_time <- Sys.time()
                cat(c("Tune_ConQuR runtime", toString(difftime(end_time, start_time, unit = "secs")), "seconds"),
                    file = sink_file_name, append = TRUE); cat("\n", file = sink_file_name, append = TRUE)
                write.csv(count_data.Tune_ConQuR$tax_final, paste0(output_root, "_Tune_ConQuR.csv"), row.names = TRUE)
            }

            ### 1.7 Tune_ConQuR_libsize
            if ("Tune_ConQuR_libsize" %in% used_methods &&
                !file.exists(paste0(output_root, "_Tune_ConQuR_libsize.csv"))) {
                start_time <- Sys.time()
                count_data.Tune_ConQuR_libsize <- Tune_ConQuR_libsize(
                    tax_tab = count_data, batchid = batchid, covariates = covar_df,
                    batch_ref_pool = batch_ref, logistic_lasso_pool = c(T, F),
                    quantile_type_pool = c("standard","lasso","composite"),
                    simple_match_pool = c(T, F), lambda_quantile_pool = c(NA,"2p/n","2p/logn"),
                    interplt_pool = c(T, F), frequencyL = 0, frequencyU = 1,
                    cutoff = 0.25, num_core = num_core)
                end_time <- Sys.time()
                cat(c("Tune_ConQuR_libsize runtime", toString(difftime(end_time, start_time, unit = "secs")), "seconds"),
                    file = sink_file_name, append = TRUE); cat("\n", file = sink_file_name, append = TRUE)
                write.csv(count_data.Tune_ConQuR_libsize$tax_final,
                           paste0(output_root, "_Tune_ConQuR_libsize.csv"), row.names = TRUE)
            }

            # Added methods (count): CLR in, inverse-CLR out. ruv3 runs after this block.

            ### 1.9 bmc
            if ("bmc" %in% used_methods &&
                !file.exists(paste0(output_root, "_bmc.csv"))) {
                start_time <- Sys.time()
                count_data.bmc <- as.matrix(count_data.clr)
                for (b in levels(batchid)) {
                    idx <- which(as.character(batchid) == b)
                    if (length(idx) == 0) next
                    count_data.bmc[idx, ] <- scale(count_data.bmc[idx, , drop = FALSE],
                                                    center = TRUE, scale = FALSE)
                }
                count_data.bmc <- clr(count_data.bmc, inverse = TRUE)
                rownames(count_data.bmc) <- rownames(count_data)
                end_time <- Sys.time()
                cat(c("bmc runtime", toString(difftime(end_time, start_time, unit = "secs")), "seconds"),
                    file = sink_file_name, append = TRUE); cat("\n", file = sink_file_name, append = TRUE)
                write.csv(count_data.bmc, paste0(output_root, "_bmc.csv"), row.names = TRUE)
            }

            ### 1.10 bat
            if ("bat" %in% used_methods &&
                !file.exists(paste0(output_root, "_bat.csv"))) {
                start_time <- Sys.time()
                bat_mod <- if (is.null(covar)) matrix(1, nrow(count_data), 1) else
                    stats::model.matrix(~ as.factor(metadata[, covar[1]]))
                bat_corrected <- sva::ComBat(dat = t(count_data.clr), batch = batchid,
                                              mod = bat_mod, par.prior = FALSE, prior.plots = FALSE)
                count_data.bat <- clr(t(bat_corrected), inverse = TRUE)
                rownames(count_data.bat) <- rownames(count_data)
                end_time <- Sys.time()
                cat(c("bat runtime", toString(difftime(end_time, start_time, unit = "secs")), "seconds"),
                    file = sink_file_name, append = TRUE); cat("\n", file = sink_file_name, append = TRUE)
                write.csv(count_data.bat, paste0(output_root, "_bat.csv"), row.names = TRUE)
            }

            ### 1.11 rbe
            if ("rbe" %in% used_methods &&
                !file.exists(paste0(output_root, "_rbe.csv"))) {
                batch_info_rbe <- as.factor(setNames(as.character(metadata[, dataset]), metadata[[Sam_id]]))
                rbe_design <- if (is.null(covar)) matrix(1, nrow(count_data), 1) else
                    stats::model.matrix(~ as.factor(metadata[, covar[1]]))
                start_time <- Sys.time()
                count_data.rbe <- clr(t(removeBatchEffect(t(count_data.clr),
                                                           batch = batch_info_rbe,
                                                           design = rbe_design)), inverse = TRUE)
                rownames(count_data.rbe) <- rownames(count_data)
                end_time <- Sys.time()
                cat(c("rbe runtime", toString(difftime(end_time, start_time, unit = "secs")), "seconds"),
                    file = sink_file_name, append = TRUE); cat("\n", file = sink_file_name, append = TRUE)
                write.csv(count_data.rbe, paste0(output_root, "_rbe.csv"), row.names = TRUE)
            }

            ### 1.12 svd
            if ("svd" %in% used_methods &&
                !file.exists(paste0(output_root, "_svd.csv"))) {
                start_time <- Sys.time()
                Y_svd    <- as.matrix(count_data.clr)
                tmp.sd   <- apply(Y_svd, 2, sd); tmp.sd[tmp.sd == 0] <- 1
                tmp.mean <- colMeans(Y_svd)
                Y_scaled <- scale(Y_svd, center = TRUE, scale = TRUE)
                svd.res  <- svd(crossprod(Y_scaled))
                svd.u1   <- svd.res$u[, 1]
                tmp.t1   <- Y_scaled %*% svd.u1 / drop(sqrt(crossprod(svd.u1)))
                tmp.c1   <- crossprod(Y_scaled, tmp.t1) / drop(crossprod(tmp.t1))
                Y_defl   <- Y_scaled - tmp.t1 %*% t(tmp.c1)
                count_data.svd <- clr(sweep(sweep(Y_defl, 2, tmp.sd, `*`), 2, tmp.mean, `+`),
                                       inverse = TRUE)
                rownames(count_data.svd) <- rownames(count_data)
                end_time <- Sys.time()
                cat(c("svd runtime", toString(difftime(end_time, start_time, unit = "secs")), "seconds"),
                    file = sink_file_name, append = TRUE); cat("\n", file = sink_file_name, append = TRUE)
                write.csv(count_data.svd, paste0(output_root, "_svd.csv"), row.names = TRUE)
            }

        } else {
            print("already done")
        }

        # ruv3 (count): skipped only if its CSV exists or ruv is not installed
        if ("ruv3" %in% used_methods && !file.exists(paste0(output_root, "_ruv3.csv"))) {
            if (!requireNamespace("ruv", quietly = TRUE)) {
                warning("ruv not installed; skipping ruv3.")
                ensure_header()
                cat("ruv3 SKIPPED (ruv package not installed)\n", file = sink_file_name, append = TRUE)
            } else {
                ensure_header()
                start_time <- Sys.time()
                Y_ruv   <- as.matrix(count_data.clr)
                ctl_ruv <- tryCatch({
                    if (!is.null(covar) && length(covar) > 0) {
                        pv <- sapply(seq_len(ncol(Y_ruv)), function(j) {
                            df_lm <- data.frame(y = Y_ruv[, j], metadata[, covar, drop = FALSE])
                            fit   <- lm(y ~ ., data = df_lm)
                            f     <- summary(fit)$fstatistic
                            if (is.null(f)) return(1.0)
                            pf(f[1], f[2], f[3], lower.tail = FALSE)
                        })
                        pv > 0.05
                    } else rep(TRUE, ncol(Y_ruv))
                }, error = function(e) rep(TRUE, ncol(Y_ruv)))
                if (sum(ctl_ruv) < 2) ctl_ruv <- rep(TRUE, ncol(Y_ruv))
                M_ruv <- if ("replicate" %in% colnames(metadata) &&
                              length(unique(metadata[["replicate"]])) < nrow(metadata)) {
                    ruv::replicate.matrix(metadata[["replicate"]])
                } else {
                    message("No technical replicates; using condition pseudo-replicates for RUV-III.")
                    if (!is.null(covar) && length(covar) > 0)
                        model.matrix(~ factor(metadata[, covar[1]]) - 1)
                    else matrix(1, nrow(Y_ruv), 1)
                }
                k_ruv <- min(max(1L, nlevels(batchid) - 1L),
                              max(1L, nrow(Y_ruv) - ncol(M_ruv) - 1L))
                Y_ruv[!is.finite(Y_ruv)] <- 0
                tryCatch({
                    ruv3_out        <- ruv::RUVIII(Y_ruv, M_ruv, ctl_ruv, k = k_ruv)
                    ruv3_corrected  <- if (is.list(ruv3_out)) ruv3_out$newY else ruv3_out
                    count_data.ruv3 <- clr(ruv3_corrected, inverse = TRUE)
                    rownames(count_data.ruv3) <- rownames(count_data)
                    end_time <- Sys.time()
                    cat(c("ruv3 runtime", toString(difftime(end_time, start_time, unit = "secs")), "seconds"),
                        file = sink_file_name, append = TRUE); cat("\n", file = sink_file_name, append = TRUE)
                    write.csv(count_data.ruv3, paste0(output_root, "_ruv3.csv"), row.names = TRUE)
                }, error = function(e) {
                    end_time <- Sys.time()
                    warning(paste("ruv3 SVD failed:", conditionMessage(e)))
                    cat(c("ruv3 runtime", toString(difftime(end_time, start_time, unit = "secs")), "seconds"),
                        file = sink_file_name, append = TRUE); cat("\n", file = sink_file_name, append = TRUE)
                    cat("ruv3 note: SVD failed; output is uncorrected (fallback to original data)\n",
                        file = sink_file_name, append = TRUE)
                    write.csv(count_data, paste0(output_root, "_ruv3.csv"), row.names = TRUE)
                })
            }
        }

    ## CASE 2 - relative abundance data
    } else {
        count_data.clr <- logratio.transfo(count_data + min(count_data[count_data > 0]),
                                            logratio = "CLR")
        print("NAs in CLR:");  print(sum(is.na(count_data.clr)))
        print("Infs in CLR:"); print(sum(is.infinite(count_data.clr)))

        relab_methods_all      <- c("combat", "limma", "MMUPHin", "ConQuR_rel",
                                    "Tune_ConQuR_rel", "bmc", "bat", "rbe", "svd")
        # ruv3 is not in this list: it runs below based only on whether its CSV exists
        relab_methods_to_check <- intersect(used_methods, relab_methods_all)

        if (needs_run(relab_methods_to_check)) {
            ensure_header()

            ### 2.1 combat
            if ("combat" %in% used_methods &&
                !file.exists(paste0(output_root, "_combat.csv"))) {
                batch_info <- as.factor(setNames(as.character(metadata[, dataset]), metadata[[Sam_id]]))
                start_time <- Sys.time()
                count_data.combat <- if (is.null(covar)) {
                    t(ComBat(t(count_data.clr), batch = batch_info, par.prior = FALSE, mod = NULL))
                } else if (length(covar) == 1) {
                    t(ComBat(t(count_data.clr), batch = batch_info, par.prior = FALSE,
                              mod = as.numeric(factor(metadata[, covar]))))
                } else {
                    t(ComBat(t(count_data.clr), batch = batch_info, par.prior = FALSE, mod = covar_df))
                }
                end_time <- Sys.time()
                cat(c("combat runtime", toString(difftime(end_time, start_time, unit = "secs")), "seconds"),
                    file = sink_file_name, append = TRUE); cat("\n", file = sink_file_name, append = TRUE)
                write.csv(count_data.combat, paste0(output_root, "_combat.csv"), row.names = TRUE)
            }

            ### 2.2 limma
            if ("limma" %in% used_methods &&
                !file.exists(paste0(output_root, "_limma.csv"))) {
                batch_info <- as.factor(setNames(as.character(metadata[, dataset]), metadata[[Sam_id]]))
                start_time <- Sys.time()
                count_data.limma <- if (is.null(covar) || length(covar) == 1) {
                    t(removeBatchEffect(t(count_data.clr), batch = batch_info))
                } else {
                    t(removeBatchEffect(t(count_data.clr), batch = batch_info, design = covar_df))
                }
                end_time <- Sys.time()
                cat(c("limma runtime", toString(difftime(end_time, start_time, unit = "secs")), "seconds"),
                    file = sink_file_name, append = TRUE); cat("\n", file = sink_file_name, append = TRUE)
                write.csv(count_data.limma, paste0(output_root, "_limma.csv"), row.names = TRUE)
            }

            ### 2.3 MMUPHin
            if ("MMUPHin" %in% used_methods &&
                !file.exists(paste0(output_root, "_MMUPHin.csv"))) {
                meta_mup <- metadata; row.names(meta_mup) <- metadata[[Sam_id]]
                feat_abd <- as.data.frame(t(count_data)); colnames(feat_abd) <- rownames(meta_mup)
                start_time <- Sys.time()
                fit <- adjust_batch(feature_abd = feat_abd, batch = dataset,
                                     covariates = covar, data = meta_mup,
                                     control = list(verbose = FALSE))
                count_data.MMUPHin <- t(fit$feature_abd_adj)
                end_time <- Sys.time()
                cat(c("MMUPHin runtime", toString(difftime(end_time, start_time, unit = "secs")), "seconds"),
                    file = sink_file_name, append = TRUE); cat("\n", file = sink_file_name, append = TRUE)
                write.csv(count_data.MMUPHin, paste0(output_root, "_MMUPHin.csv"), row.names = TRUE)
            }

            ### 2.4 ConQuR_rel
            if ("ConQuR_rel" %in% used_methods &&
                !file.exists(paste0(output_root, "_ConQuR_rel.csv"))) {
                start_time <- Sys.time()
                count_data.ConQuR_rel <- ConQuR_rel(tax_tab = count_data, batchid = batchid,
                                                     covariates = covar_df, batch_ref = batch_ref, simple_match = is.null(covar_df), num_core = num_core)
                end_time <- Sys.time()
                cat(c("ConQuR_rel runtime", toString(difftime(end_time, start_time, unit = "secs")), "seconds"),
                    file = sink_file_name, append = TRUE); cat("\n", file = sink_file_name, append = TRUE)
                write.csv(count_data.ConQuR_rel, paste0(output_root, "_ConQuR_rel.csv"), row.names = TRUE)
            }

            ### 2.5 Tune_ConQuR_rel
            if ("Tune_ConQuR_rel" %in% used_methods &&
                !file.exists(paste0(output_root, "_Tune_ConQuR_rel.csv"))) {
                start_time <- Sys.time()
                count_data.Tune_ConQuR_rel <- Tune_ConQuR_rel(
                    tax_tab = count_data, batchid = batchid, covariates = covar_df,
                    batch_ref_pool = batch_ref, logistic_lasso_pool = c(T, F),
                    quantile_type_pool = c("standard","lasso","composite"),
                    simple_match_pool = c(T, F), lambda_quantile_pool = c(NA,"2p/n","2p/logn"),
                    interplt_pool = c(T, F), frequencyL = 0, frequencyU = 1,
                    cutoff = 0.25, num_core = num_core)
                end_time <- Sys.time()
                cat(c("Tune_ConQuR_rel runtime", toString(difftime(end_time, start_time, unit = "secs")), "seconds"),
                    file = sink_file_name, append = TRUE); cat("\n", file = sink_file_name, append = TRUE)
                write.csv(count_data.Tune_ConQuR_rel$tax_final,
                           paste0(output_root, "_Tune_ConQuR_rel.csv"), row.names = TRUE)
            }

            # Added methods (relab): CLR in, CLR out. ruv3 runs after this block.

            ### 2.7 bmc
            if ("bmc" %in% used_methods &&
                !file.exists(paste0(output_root, "_bmc.csv"))) {
                start_time <- Sys.time()
                count_data.bmc <- as.matrix(count_data.clr)
                for (b in levels(batchid)) {
                    idx <- which(as.character(batchid) == b)
                    if (length(idx) == 0) next
                    count_data.bmc[idx, ] <- scale(count_data.bmc[idx, , drop = FALSE],
                                                    center = TRUE, scale = FALSE)
                }
                rownames(count_data.bmc) <- rownames(count_data)
                end_time <- Sys.time()
                cat(c("bmc runtime", toString(difftime(end_time, start_time, unit = "secs")), "seconds"),
                    file = sink_file_name, append = TRUE); cat("\n", file = sink_file_name, append = TRUE)
                write.csv(count_data.bmc, paste0(output_root, "_bmc.csv"), row.names = TRUE)
            }

            ### 2.8 bat
            if ("bat" %in% used_methods &&
                !file.exists(paste0(output_root, "_bat.csv"))) {
                start_time <- Sys.time()
                bat_mod <- if (is.null(covar)) matrix(1, nrow(count_data), 1) else
                    stats::model.matrix(~ as.factor(metadata[, covar[1]]))
                bat_corrected <- sva::ComBat(dat = t(count_data.clr), batch = batchid,
                                              mod = bat_mod, par.prior = FALSE, prior.plots = FALSE)
                count_data.bat <- t(bat_corrected)   # stay in CLR space
                rownames(count_data.bat) <- rownames(count_data)
                end_time <- Sys.time()
                cat(c("bat runtime", toString(difftime(end_time, start_time, unit = "secs")), "seconds"),
                    file = sink_file_name, append = TRUE); cat("\n", file = sink_file_name, append = TRUE)
                write.csv(count_data.bat, paste0(output_root, "_bat.csv"), row.names = TRUE)
            }

            ### 2.9 rbe
            if ("rbe" %in% used_methods &&
                !file.exists(paste0(output_root, "_rbe.csv"))) {
                batch_info_rbe <- as.factor(setNames(as.character(metadata[, dataset]), metadata[[Sam_id]]))
                rbe_design <- if (is.null(covar)) matrix(1, nrow(count_data), 1) else
                    stats::model.matrix(~ as.factor(metadata[, covar[1]]))
                start_time <- Sys.time()
                count_data.rbe <- t(removeBatchEffect(t(count_data.clr),
                                                       batch = batch_info_rbe, design = rbe_design))
                rownames(count_data.rbe) <- rownames(count_data)
                end_time <- Sys.time()
                cat(c("rbe runtime", toString(difftime(end_time, start_time, unit = "secs")), "seconds"),
                    file = sink_file_name, append = TRUE); cat("\n", file = sink_file_name, append = TRUE)
                write.csv(count_data.rbe, paste0(output_root, "_rbe.csv"), row.names = TRUE)
            }

            ### 2.10 svd
            if ("svd" %in% used_methods &&
                !file.exists(paste0(output_root, "_svd.csv"))) {
                start_time <- Sys.time()
                Y_svd    <- as.matrix(count_data.clr)
                tmp.sd   <- apply(Y_svd, 2, sd); tmp.sd[tmp.sd == 0] <- 1
                tmp.mean <- colMeans(Y_svd)
                Y_scaled <- scale(Y_svd, center = TRUE, scale = TRUE)
                svd.res  <- svd(crossprod(Y_scaled))
                svd.u1   <- svd.res$u[, 1]
                tmp.t1   <- Y_scaled %*% svd.u1 / drop(sqrt(crossprod(svd.u1)))
                tmp.c1   <- crossprod(Y_scaled, tmp.t1) / drop(crossprod(tmp.t1))
                Y_defl   <- Y_scaled - tmp.t1 %*% t(tmp.c1)
                count_data.svd <- sweep(sweep(Y_defl, 2, tmp.sd, `*`), 2, tmp.mean, `+`)
                rownames(count_data.svd) <- rownames(count_data)
                end_time <- Sys.time()
                cat(c("svd runtime", toString(difftime(end_time, start_time, unit = "secs")), "seconds"),
                    file = sink_file_name, append = TRUE); cat("\n", file = sink_file_name, append = TRUE)
                write.csv(count_data.svd, paste0(output_root, "_svd.csv"), row.names = TRUE)
            }

        } else {
            print("already done")
        }

        # ruv3 (relab): skipped only if its CSV exists or ruv is not installed
        if ("ruv3" %in% used_methods && !file.exists(paste0(output_root, "_ruv3.csv"))) {
            if (!requireNamespace("ruv", quietly = TRUE)) {
                warning("ruv not installed; skipping ruv3.")
                ensure_header()
                cat("ruv3 SKIPPED (ruv package not installed)\n", file = sink_file_name, append = TRUE)
            } else {
                ensure_header()
                start_time <- Sys.time()
                Y_ruv   <- as.matrix(count_data.clr)
                ctl_ruv <- tryCatch({
                    if (!is.null(covar) && length(covar) > 0) {
                        pv <- sapply(seq_len(ncol(Y_ruv)), function(j) {
                            df_lm <- data.frame(y = Y_ruv[, j], metadata[, covar, drop = FALSE])
                            fit   <- lm(y ~ ., data = df_lm)
                            f     <- summary(fit)$fstatistic
                            if (is.null(f)) return(1.0)
                            pf(f[1], f[2], f[3], lower.tail = FALSE)
                        })
                        pv > 0.05
                    } else rep(TRUE, ncol(Y_ruv))
                }, error = function(e) rep(TRUE, ncol(Y_ruv)))
                if (sum(ctl_ruv) < 2) ctl_ruv <- rep(TRUE, ncol(Y_ruv))
                M_ruv <- if ("replicate" %in% colnames(metadata) &&
                              length(unique(metadata[["replicate"]])) < nrow(metadata)) {
                    ruv::replicate.matrix(metadata[["replicate"]])
                } else {
                    message("No technical replicates; using condition pseudo-replicates for RUV-III.")
                    if (!is.null(covar) && length(covar) > 0)
                        model.matrix(~ factor(metadata[, covar[1]]) - 1)
                    else matrix(1, nrow(Y_ruv), 1)
                }
                k_ruv <- min(max(1L, nlevels(batchid) - 1L),
                              max(1L, nrow(Y_ruv) - ncol(M_ruv) - 1L))
                Y_ruv[!is.finite(Y_ruv)] <- 0
                tryCatch({
                    ruv3_out        <- ruv::RUVIII(Y_ruv, M_ruv, ctl_ruv, k = k_ruv)
                    ruv3_corrected  <- if (is.list(ruv3_out)) ruv3_out$newY else ruv3_out
                    count_data.ruv3 <- ruv3_corrected   # stay in CLR space for relab
                    rownames(count_data.ruv3) <- rownames(count_data)
                    end_time <- Sys.time()
                    cat(c("ruv3 runtime", toString(difftime(end_time, start_time, unit = "secs")), "seconds"),
                        file = sink_file_name, append = TRUE); cat("\n", file = sink_file_name, append = TRUE)
                    write.csv(count_data.ruv3, paste0(output_root, "_ruv3.csv"), row.names = TRUE)
                }, error = function(e) {
                    end_time <- Sys.time()
                    warning(paste("ruv3 SVD failed:", conditionMessage(e)))
                    cat(c("ruv3 runtime", toString(difftime(end_time, start_time, unit = "secs")), "seconds"),
                        file = sink_file_name, append = TRUE); cat("\n", file = sink_file_name, append = TRUE)
                    cat("ruv3 note: SVD failed; output is uncorrected (fallback to original data)\n",
                        file = sink_file_name, append = TRUE)
                    write.csv(count_data.clr, paste0(output_root, "_ruv3.csv"), row.names = TRUE)
                })
            }
        }
    }
}

# Helpers
mcsapply <- function(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
    FUN    <- match.fun(FUN)
    answer <- parallel::mclapply(X = X, FUN = FUN, ...)
    if (USE.NAMES && is.character(X) && is.null(names(answer))) names(answer) <- X
    if (!isFALSE(simplify) && length(answer))
        simplify2array(answer, higher = (simplify == "array"))
    else answer
}

# Slurm runner (one iteration per job)
# covar: pass c("cond") for standard run, NULL for nocovar/ML-safe run
scaled_slurm_methods_benchmarking <- function(output_dir, overall_path,
                                               method_l, or_l,
                                               cond_effect_val_l, batch_effect_val_l,
                                               iter, count = TRUE, covar = c("cond")) {
    for (or in or_l) {
        for (cond_effect_val in cond_effect_val_l) {
            for (batch_effect_val in batch_effect_val_l) {
                if (cond_effect_val + batch_effect_val <= 1) {
                    print(paste("or:", or, "| cond:", cond_effect_val,
                                "| batch:", batch_effect_val, "| iter:", iter))
                    data_file <- paste0(overall_path,
                        if (count) "/ibd_150_count_" else "/ibd_150_relab_",
                        or, "_", cond_effect_val, "_", batch_effect_val, "_iter_", iter, ".csv")
                    meta_file <- paste0(overall_path, "/ibd_150_meta_",
                                        or, "_", cond_effect_val, "_", batch_effect_val,
                                        "_iter_", iter, ".csv")
                    out_root  <- paste0(output_dir, "/out_",
                                        or, "_", cond_effect_val, "_", batch_effect_val,
                                        "_iter_", iter, "/ibd_",
                                        or, "_", cond_effect_val, "_", batch_effect_val,
                                        "_iter_", iter)
                    run_methods(data_file, meta_file, out_root,
                                 batch_ref    = "batch_0",
                                 dataset      = "batchid",
                                 covar        = covar,
                                 Sam_id       = "Sam_id",
                                 transpose    = TRUE,
                                 count        = count,
                                 used_methods = method_l)
                }
            }
        }
    }
}

# Parallel runner (all iterations in one call)
scaled_midas_methods_benchmarking <- function(output_dir, overall_path,
                                               method_l, or_l,
                                               cond_effect_val_l, batch_effect_val_l,
                                               num_iter, count = TRUE) {
    for (or in or_l) {
        for (cond_effect_val in cond_effect_val_l) {
            for (batch_effect_val in batch_effect_val_l) {
                if (cond_effect_val + batch_effect_val <= 1) {
                    mcsapply(seq_len(num_iter), function(iter)
                        scaled_slurm_methods_benchmarking(
                            output_dir, overall_path, method_l, or_l,
                            cond_effect_val_l, batch_effect_val_l, iter, count),
                        mc.cores = 5)
                }
            }
        }
    }
}

# Parameter grid
or_condsweep  <- c(CONDSWEEP_OR)
cond_sweep_l  <- c(0.005, 0.01, 0.02, 0.05, 0.1)
batch_single  <- c(0.5)

count_methods <- c("combat_seq", "limma", "MMUPHin", "ConQuR", "ConQuR_libsize",
                    "ruv3", "bmc", "bat", "rbe", "svd")
relab_methods <- c("combat", "limma", "MMUPHin", "ConQuR_rel",
                    "ruv3", "bmc", "bat", "rbe", "svd")

sim_base <- file.path(SIM_OUTPUT_ROOT, 'scenarios_simulation_data_052626')

# no-relation / count, covariate-free only (used for the OR=1 sweep, which only needs the leakage-free RF
# metrics). Outputs go to both the regular and _nocovar folders so evaluation finds them, which means the
# regular folder holds no covariate-guided corrections in this mode.
if (COMBO_SELECT == "count_no_nocovar_only") {
    for (od in c('/output_count_condsweep_no_relation_nocovar', '/output_count_condsweep_no_relation')) {
        print(paste("condsweep count no-relation - nocovar only ->", od))
        scaled_slurm_methods_benchmarking(
            paste0(sim_base, od),
            paste0(sim_base, '/count_condsweep_no_relation'),
            count_methods, or_condsweep, cond_sweep_l, batch_single, GLOBAL_ITER, count = TRUE, covar = NULL)
    }
}

# no-relation / count
if (COMBO_SELECT %in% c("all", "count_no")) {
    print("condsweep count no-relation")
    scaled_slurm_methods_benchmarking(
        paste0(sim_base, '/output_count_condsweep_no_relation'),
        paste0(sim_base, '/count_condsweep_no_relation'),
        count_methods, or_condsweep, cond_sweep_l, batch_single, GLOBAL_ITER, count = TRUE)

    print("condsweep count no-relation - nocovar")
    scaled_slurm_methods_benchmarking(
        paste0(sim_base, '/output_count_condsweep_no_relation_nocovar'),
        paste0(sim_base, '/count_condsweep_no_relation'),
        count_methods, or_condsweep, cond_sweep_l, batch_single, GLOBAL_ITER, count = TRUE, covar = NULL)
}

# no-relation / relab
if (COMBO_SELECT %in% c("all", "relab_no")) {
    print("condsweep relab no-relation")
    scaled_slurm_methods_benchmarking(
        paste0(sim_base, '/output_relab_condsweep_no_relation'),
        paste0(sim_base, '/count_condsweep_no_relation'),
        relab_methods, or_condsweep, cond_sweep_l, batch_single, GLOBAL_ITER, count = FALSE)

    print("condsweep relab no-relation - nocovar")
    scaled_slurm_methods_benchmarking(
        paste0(sim_base, '/output_relab_condsweep_no_relation_nocovar'),
        paste0(sim_base, '/count_condsweep_no_relation'),
        relab_methods, or_condsweep, cond_sweep_l, batch_single, GLOBAL_ITER, count = FALSE, covar = NULL)
}

# yes-relation / count
if (COMBO_SELECT %in% c("all", "count_yes")) {
    print("condsweep count yes-relation")
    scaled_slurm_methods_benchmarking(
        paste0(sim_base, '/output_count_condsweep_yes_relation'),
        paste0(sim_base, '/count_condsweep_yes_relation'),
        count_methods, or_condsweep, cond_sweep_l, batch_single, GLOBAL_ITER, count = TRUE)

    print("condsweep count yes-relation - nocovar")
    scaled_slurm_methods_benchmarking(
        paste0(sim_base, '/output_count_condsweep_yes_relation_nocovar'),
        paste0(sim_base, '/count_condsweep_yes_relation'),
        count_methods, or_condsweep, cond_sweep_l, batch_single, GLOBAL_ITER, count = TRUE, covar = NULL)
}

# yes-relation / relab
if (COMBO_SELECT %in% c("all", "relab_yes")) {
    print("condsweep relab yes-relation")
    scaled_slurm_methods_benchmarking(
        paste0(sim_base, '/output_relab_condsweep_yes_relation'),
        paste0(sim_base, '/count_condsweep_yes_relation'),
        relab_methods, or_condsweep, cond_sweep_l, batch_single, GLOBAL_ITER, count = FALSE)

    print("condsweep relab yes-relation - nocovar")
    scaled_slurm_methods_benchmarking(
        paste0(sim_base, '/output_relab_condsweep_yes_relation_nocovar'),
        paste0(sim_base, '/count_condsweep_yes_relation'),
        relab_methods, or_condsweep, cond_sweep_l, batch_single, GLOBAL_ITER, count = FALSE, covar = NULL)
}

print(paste("DONE - condsweep methods benchmarking complete for this iteration (combo:", COMBO_SELECT, ")"))
