#!/usr/bin/env Rscript
# Apply the R-based integration strategies to one replicate of a simulation scenario (all its settings).
# Called by Snakefile_simulation. Writes <out>/out_<or>_<cond>_<batch>_iter_<i>/ (condition used as a covariate)
# and <out>_nocovar/... (no condition covariate; used for the leakage-free random forest evaluation).
#
# Usage: Rscript run_sim_methods.R --data=<dir> --out=<dir> --iter=<i> --datatype=count|relab
#          --methods=combat_seq,limma,... [--or_l=1] [--cond_effect_val_l=0.5] [--batch_effect_val_l=0.5]
#          [--covar=both|yes|no]
# Same order of random draws as methods_benchmarking_sim_scenario.R: set.seed(iter), covariate pass, then
# no-covariate pass.

source(Find(file.exists, c("paths.R", "../paths.R", "../../paths.R")))
source(file.path(BENCH_DIR, "sim_lib", "sim_methods_lib.R"))

a <- list()
for (arg in commandArgs(trailingOnly = TRUE)) {
    kv <- regmatches(arg, regexec("^--([^=]+)=(.*)$", arg))[[1]]
    if (length(kv) != 3) stop(paste("bad argument:", arg))
    a[[kv[2]]] <- kv[3]
}
get <- function(key, default) if (is.null(a[[key]])) default else a[[key]]
num_list <- function(x) as.numeric(strsplit(x, ",")[[1]])

data_dir <- a[["data"]]; out_dir <- a[["out"]]; iter <- as.integer(a[["iter"]])
datatype <- get("datatype", "count")
methods  <- strsplit(a[["methods"]], ",")[[1]]
or_l     <- num_list(get("or_l", "1"))
cond_l   <- num_list(get("cond_effect_val_l", "0.5"))
batch_l  <- num_list(get("batch_effect_val_l", "0.5"))
covar    <- get("covar", "both")
if (is.null(data_dir) || is.null(out_dir) || is.na(iter)) stop("--data, --out and --iter are required")
is_count <- datatype == "count"

set.seed(iter)
if (covar %in% c("both", "yes"))
    scaled_slurm_methods_benchmarking(out_dir, data_dir, methods, or_l, cond_l, batch_l, iter, count = is_count)
if (covar %in% c("both", "no"))
    scaled_slurm_methods_benchmarking(paste0(out_dir, "_nocovar"), data_dir, methods, or_l, cond_l, batch_l, iter,
                                      count = is_count, covar = NULL)
print(paste("DONE - R integration strategies,", datatype, "iteration", iter))
