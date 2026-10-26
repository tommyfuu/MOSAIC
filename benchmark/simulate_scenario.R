#!/usr/bin/env Rscript
# Simulate one replicate of a simulation scenario (all its OR x condition x batch settings) with MIDASim.
# Called by Snakefile_simulation; every design parameter is an explicit argument (see config_simulation_scenarios.yml).
#
# Usage: Rscript simulate_scenario.R --out=<dir> --iter=<i> [--n_samples=450] [--resample_libsize=true]
#          [--p_cond=0.5] [--p_batch=0.5] [--or_l=1] [--cond_effect_val_l=0.5] [--batch_effect_val_l=0.5]
#          [--batch_libsize_related=false] [--batch_taxa_fraction=1]
# Lists are comma-separated. Replicate i uses set.seed(i) before drawing library sizes, as in
# generate_data_MIDASim_imbalance.R.

source(Find(file.exists, c("paths.R", "../paths.R", "../../paths.R")))
source(file.path(BENCH_DIR, "sim_lib", "midas_sim_lib.R"))

parse_args <- function(args) {
    out <- list()
    for (a in args) {
        kv <- regmatches(a, regexec("^--([^=]+)=(.*)$", a))[[1]]
        if (length(kv) != 3) stop(paste("bad argument:", a))
        out[[kv[2]]] <- kv[3]
    }
    out
}
a <- parse_args(commandArgs(trailingOnly = TRUE))
get <- function(key, default) if (is.null(a[[key]])) default else a[[key]]
num_list <- function(x) as.numeric(strsplit(x, ",")[[1]])
as_bool <- function(x) tolower(x) %in% c("true", "t", "1", "yes")

out_dir   <- a[["out"]]
iter      <- as.integer(a[["iter"]])
n_samples <- as.integer(get("n_samples", "450"))
resample  <- as_bool(get("resample_libsize", "true"))
p_cond    <- as.numeric(get("p_cond", "0.5"))
p_batch   <- as.numeric(get("p_batch", "0.5"))
or_l      <- num_list(get("or_l", "1"))
cond_l    <- num_list(get("cond_effect_val_l", "0.5"))
batch_l   <- num_list(get("batch_effect_val_l", "0.5"))
libsize_related <- as_bool(get("batch_libsize_related", "false"))
batch_frac <- as.numeric(get("batch_taxa_fraction", "1"))
if (is.null(out_dir) || is.na(iter)) stop("--out and --iter are required")

load(file.path(BENCH_DIR, "ibd_150.Rdata"))   # template: `otu` (taxa x samples), IBDMDB genus level
otu_original <- t(otu)
stopifnot(ncol(otu_original) == 301)
libsize_orig <- rowSums(otu_original)

set.seed(iter)
if (resample) {
    libsize_l <- sample(libsize_orig, n_samples, replace = TRUE)
} else {
    if (n_samples != length(libsize_orig)) stop("resample_libsize = false requires n_samples = ", length(libsize_orig))
    libsize_l <- libsize_orig
}
id_cond  <- condition_taxa(otu_original)
id_batch <- batch_taxa(otu_original, batch_frac)
print(paste("iter", iter, "| n", n_samples, "| resample", resample, "| p_cond", p_cond, "| p_batch", p_batch,
            "| OR", paste(or_l, collapse = ","), "| cond", paste(cond_l, collapse = ","),
            "| batch", paste(batch_l, collapse = ","), "| libsize_related", libsize_related,
            "| batch taxa", length(id_batch)))

simulate_grid(out_dir, otu_original, or_l, cond_l, batch_l, iter, libsize_l, id_cond, id_batch,
              batch_libsize_related = libsize_related, p_cond = p_cond, p_batch = p_batch)
print("simulation done")
