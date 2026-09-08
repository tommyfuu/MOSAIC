# generate_data_MIDASim_052626.R
#
# Generates ONLY the new reviewer-response scenarios.
# The original scenario (n=450, all taxa, p_batch=p_cond=0.5) has already
# been generated and does not need to be re-run.
#
# New scenarios:
#
#   1  Small cohort  n=150 (no replacement, use original lib sizes directly),
#                    all 301 taxa, p_batch=0.5, p_cond=0.5
#                    Full OR x effect grid; both confounding types
#
#   2  Half taxa     n=450, ~151 taxa (every-other by zero-inflation rank),
#                    OR=1, reduced effect grid {0,0.5}x{0,0.5},
#                    no lib-size confounding
#
#   3  40% batch     n=450, all 301 taxa, p_batch=0.4, p_cond=0.5,
#                    OR=1, cond=0.5, batch=0.5, no lib-size confounding
#
#   4  40% cond      n=450, all 301 taxa, p_cond=0.4, p_batch=0.5,
#                    OR=1, cond=0.5, batch=0.5, no lib-size confounding
#
# All scenarios: 1000 iterations (one per Slurm job via GLOBAL_ITER arg).
# Set base_output_root below before running.

# ============================================================
# LOAD DATA
# ============================================================
current_path <- getwd()
if (grepl("/benchmark", current_path)) {
    load("./ibd_150.Rdata")
} else {
    load("./benchmark/ibd_150.Rdata")
}

library(bindata)
library(MIDASim)
library(tibble)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
print(args)
GLOBAL_ITER <- if (length(args) == 1) {
    as.integer(args[1])
} else {
    print("No iteration argument supplied; defaulting to iter = 1")
    1L
}

# ============================================================
# BASE SETUP
# ============================================================
otu_original <- t(otu)
p            <- ncol(otu_original)   # 301 taxa
n_orig       <- nrow(otu_original)   # 150 samples in source data

stopifnot(p == 301)
print(paste("p =", p, "| n_orig =", n_orig))

libsize_orig <- rowSums(otu_original)

# Scenarios 2-4: n=450, sample with replacement (same as original)
libsize_450  <- sample(libsize_orig, n_orig * 3, replace = TRUE)

# Scenario 1: n=150, no replacement — use original 150 lib sizes directly
libsize_150  <- libsize_orig

# ============================================================
# HELPER: odds ratio -> binary correlation  (from ConQuR paper)
# ============================================================
bincorr <- function(OR, p1, p2) {
    if (OR == 1) {
        p11 <- p2 - p2 + p1 * p2
    } else {
        p11_1 <- p2 - (1/2/(1-OR) * (1-p1+OR*p1+p2-OR*p2 -
                   sqrt((-1+p1-OR*p1-p2+OR*p2)^2 - 4*(1-OR)*(p2-p1*p2))))
        p11_2 <- p2 - (1/2/(1-OR) * (1-p1+OR*p1+p2-OR*p2 -
                   sqrt((-1+p1-OR*p1-p2+OR*p2)^2) - 4*(1-OR)*(p2-p1*p2)))
        p11 <- if (p11_1 > 0 && p11_1 < p1 && p11_1 < p2) p11_1 else p11_2
    }
    (p11 - p1*p2) / sqrt(p1*(1-p1)*p2*(1-p2))
}

# ============================================================
# TAXA SELECTION
# ============================================================
zp <- colMeans(otu_original == 0)   # zero-inflation rate per taxon

# Condition taxa: 50 taxa spread evenly across zero-inflation spectrum
# (unchanged from original)
id_cond <- NULL
for (tau in 0:24/25) {
    id_cond <- c(id_cond, order(abs(zp - tau))[1:2])
}

# Batch taxa — scenarios 1, 3, 4: all 301 (same as original)
id_batch <- 1:p

# Batch taxa — scenario 2: every other taxon sorted by zero-inflation rate
# Selects 151 of 301 taxa covering the full prevalence spectrum
id_batch_half <- order(zp)[seq(1, p, by = 2)]

print(paste("id_cond length:",       length(id_cond)))
print(paste("id_batch length:",      length(id_batch)))
print(paste("id_batch_half length:", length(id_batch_half)))

# ============================================================
# CORE SIMULATION FUNCTION
# ============================================================
midas_simulate <- function(otu_original, n, or, cond_effect, batch_effect,
                            out_count, out_relab, out_meta,
                            p_cond = 0.5, p_batch = 0.5,
                            libsize_l = NULL, batch_libsize_related = FALSE,
                            id_batch_param = NULL) {

    if (is.null(libsize_l)) libsize_l <- rep(10000, n)

    id_batch_use <- if (!is.null(id_batch_param)) id_batch_param else id_batch

    print(paste("n:", n, "| or:", or,
                "| cond_effect:", cond_effect, "| batch_effect:", batch_effect,
                "| p_cond:", p_cond, "| p_batch:", p_batch,
                "| n_batch_taxa:", length(id_batch_use),
                "| libsize_related:", batch_libsize_related))
    print(paste("n_cond_taxa:", length(zp)))

    # Assign patients to (condition, batch) groups
    if (batch_libsize_related == TRUE) {
        p_batch_vec <- 1 / (1 + exp(-scale(libsize_l)))
        bin_corr    <- sapply(p_batch_vec, function(x) bincorr(or, p_cond, x))
        cond_batchid_vec <- matrix(ncol = 2, nrow = n)
        for (i in 1:n) {
            cond_batchid_vec[i, ] <- rmvbin(1, c(p_cond, p_batch_vec[i]),
                                             bincorr = (1 - bin_corr[i]) * diag(2) + bin_corr[i])
        }
    } else {
        bin_corr         <- bincorr(or, p_cond, p_batch)
        print(c("OR, p_cond, p_batch, bin_corr:", or, p_cond, p_batch, bin_corr))
        cond_batchid_vec <- rmvbin(n, c(p_cond, p_batch),
                                    bincorr = (1 - bin_corr) * diag(2) + bin_corr)
    }

    cond    <- cond_batchid_vec[, 1]
    batchid <- cond_batchid_vec[, 2]

    current_metadata <- data.frame(
        subjectid_text = sprintf("Subject_%d", 1:n),
        batchid        = ifelse(batchid == 0, "batch_0", "batch_1"),
        cond           = ifelse(cond    == 0, "cond_0",  "cond_1")
    )
    write.csv(current_metadata, file = out_meta, row.names = FALSE)

    print("MIDASim running...")
    fitted <- MIDASim.setup(otu_original, mode = 'nonparametric')

    perm_batch <- sample(id_batch_use, length(id_batch_use))
    perm_cond  <- sample(id_cond,      length(id_cond))

    rel0      <- fitted$mean.rel.abund
    rel_batch <- fitted$mean.rel.abund; rel_batch[id_batch_use] <- rel_batch[perm_batch]
    rel_cond  <- fitted$mean.rel.abund; rel_cond[id_cond]       <- rel_cond[perm_cond]

    rel10      <- unlist(lapply(fitted$rel.abund.1, mean))
    rel1_batch <- unlist(lapply(fitted$rel.abund.1, mean)); rel1_batch[id_batch_use] <- rel1_batch[perm_batch]
    rel1_cond  <- unlist(lapply(fitted$rel.abund.1, mean)); rel1_cond[id_cond]       <- rel1_cond[perm_cond]

    fitted_c0b0 <- MIDASim.modify(fitted, lib.size = libsize_l, SCAM = TRUE)
    fitted_c1b0 <- MIDASim.modify(fitted, lib.size = libsize_l, SCAM = TRUE,
        mean.rel.abund.1 = (1 - cond_effect) * rel10 + cond_effect * rel1_cond,
        mean.rel.abund   = (1 - cond_effect) * rel0  + cond_effect * rel_cond)
    fitted_c0b1 <- MIDASim.modify(fitted, lib.size = libsize_l, SCAM = TRUE,
        mean.rel.abund.1 = (1 - batch_effect) * rel10 + batch_effect * rel1_batch,
        mean.rel.abund   = (1 - batch_effect) * rel0  + batch_effect * rel_batch)
    fitted_c1b1 <- MIDASim.modify(fitted, lib.size = libsize_l, SCAM = TRUE,
        mean.rel.abund.1 = (1 - cond_effect - batch_effect) * rel10 + cond_effect * rel1_cond  + batch_effect * rel1_batch,
        mean.rel.abund   = (1 - cond_effect - batch_effect) * rel0  + cond_effect * rel_cond   + batch_effect * rel_batch)

    otu_c0b0 <- MIDASim(fitted_c0b0)$sim_count
    otu_c1b0 <- MIDASim(fitted_c1b0)$sim_count
    otu_c0b1 <- MIDASim(fitted_c0b1)$sim_count
    otu_c1b1 <- MIDASim(fitted_c1b1)$sim_count

    otu_sim <- matrix(NA, nrow = n, ncol = p)
    otu_sim[cond == 0 & batchid == 0, ] <- otu_c0b0[cond == 0 & batchid == 0, ]
    otu_sim[cond == 1 & batchid == 0, ] <- otu_c1b0[cond == 1 & batchid == 0, ]
    otu_sim[cond == 0 & batchid == 1, ] <- otu_c0b1[cond == 0 & batchid == 1, ]
    otu_sim[cond == 1 & batchid == 1, ] <- otu_c1b1[cond == 1 & batchid == 1, ]

    write.csv(otu_sim, file = out_count, row.names = FALSE)
    rela <- otu_sim / rowSums(otu_sim)
    write.csv(rela, file = out_relab, row.names = FALSE)
}

# ============================================================
# PER-ITERATION WRAPPER
# ============================================================
midas_generate_per_iter <- function(otu_original, output_root,
                                     or, cond_effect_val, batch_effect_val,
                                     iter, libsize_l, batch_libsize_related,
                                     p_cond = 0.5, p_batch = 0.5,
                                     id_batch_param = NULL) {
    out_count <- paste0(output_root, "/ibd_150_count_",
                        or, "_", cond_effect_val, "_", batch_effect_val, "_iter_", iter, ".csv")
    out_relab <- paste0(output_root, "/ibd_150_relab_",
                        or, "_", cond_effect_val, "_", batch_effect_val, "_iter_", iter, ".csv")
    out_meta  <- paste0(output_root, "/ibd_150_meta_",
                        or, "_", cond_effect_val, "_", batch_effect_val, "_iter_", iter, ".csv")

    if (file.exists(out_relab)) {
        print(paste(out_relab, "already exists, skipping"))
    } else {
        midas_simulate(otu_original, length(libsize_l),
                        or, cond_effect_val, batch_effect_val,
                        out_count, out_relab, out_meta,
                        p_cond = p_cond, p_batch = p_batch,
                        libsize_l = libsize_l,
                        batch_libsize_related = batch_libsize_related,
                        id_batch_param = id_batch_param)
    }
}

# ============================================================
# GRID RUNNER
# ============================================================
scaled_slurm_midas_data_generation <- function(output_root, otu_original,
                                                 or_l, cond_effect_val_l, batch_effect_val_l,
                                                 iter, libsize_l,
                                                 batch_libsize_related = FALSE,
                                                 p_cond = 0.5, p_batch = 0.5,
                                                 id_batch_param = NULL) {
    id_batch_log <- if (!is.null(id_batch_param)) id_batch_param else id_batch

    for (or in or_l) {
        for (cond_effect_val in cond_effect_val_l) {
            for (batch_effect_val in batch_effect_val_l) {
                if (cond_effect_val + batch_effect_val <= 1) {
                    sink(paste0(output_root, "/ibd_150_id_cond_",
                                or, "_", cond_effect_val, "_", batch_effect_val, "_iter_", iter, ".txt"))
                    cat(id_cond, "\n"); sink()

                    sink(paste0(output_root, "/ibd_150_id_batch_",
                                or, "_", cond_effect_val, "_", batch_effect_val, "_iter_", iter, ".txt"))
                    cat(id_batch_log, "\n"); sink()

                    midas_generate_per_iter(otu_original, output_root,
                                             or, cond_effect_val, batch_effect_val,
                                             iter, libsize_l, batch_libsize_related,
                                             p_cond = p_cond, p_batch = p_batch,
                                             id_batch_param = id_batch_param)
                }
            }
        }
    }
}

# ============================================================
# PARAMETER GRIDS
# Two selected combos based on original results: OR in {1, 1.25}, cond=0.5, batch=0.5.
# These represent the moderate-effect no-confounding and moderate-confounding cases
# that best differentiate methods in the original benchmark.
# ============================================================
or_one        <- c(1)
cond_single   <- c(0.5)
batch_single  <- c(0.5)

# ============================================================
# OUTPUT DIRECTORIES  — set base_output_root before running
# ============================================================
base_output_root <- '/athena/linglab/scratch/chf4012/simulation_outputs/scenarios_simulation_data_052626'

dir_s1_norel  <- file.path(base_output_root, 'count_n150_no_relation')
dir_s1_yesrel <- file.path(base_output_root, 'count_n150_yes_relation')
dir_s2        <- file.path(base_output_root, 'count_halftaxa')
dir_s3        <- file.path(base_output_root, 'count_batchfrac40')
dir_s4        <- file.path(base_output_root, 'count_condfrac40')

for (d in c(dir_s1_norel, dir_s1_yesrel, dir_s2, dir_s3, dir_s4)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ============================================================
# SCENARIO 1  Small cohort (n=150, no replacement, no lib-size confounding)
# OR in {1, 1.25}, cond=0.5, batch=0.5  — no-relation only
# ============================================================
scaled_slurm_midas_data_generation(
    dir_s1_norel, otu_original,
    or_one, cond_single, batch_single,
    iter = GLOBAL_ITER, libsize_l = libsize_150,
    batch_libsize_related = FALSE)

# ============================================================
# SCENARIO 2  Half taxa (every-other by zero-inflation rank)
# OR in {1, 1.25}, cond=0.5, batch=0.5
# ============================================================
scaled_slurm_midas_data_generation(
    dir_s2, otu_original,
    or_one, cond_single, batch_single,
    iter = GLOBAL_ITER, libsize_l = libsize_450,
    batch_libsize_related = FALSE,
    id_batch_param = id_batch_half)

# ============================================================
# SCENARIO 3  40% patients affected by batch  (p_batch=0.4)
# OR in {1, 1.25}, cond=0.5, batch=0.5
# ============================================================
scaled_slurm_midas_data_generation(
    dir_s3, otu_original,
    or_one, cond_single, batch_single,
    iter = GLOBAL_ITER, libsize_l = libsize_450,
    batch_libsize_related = FALSE,
    p_batch = 0.4)

# ============================================================
# SCENARIO 4  40% patients affected by condition  (p_cond=0.4)
# OR in {1, 1.25}, cond=0.5, batch=0.5
# ============================================================
scaled_slurm_midas_data_generation(
    dir_s4, otu_original,
    or_one, cond_single, batch_single,
    iter = GLOBAL_ITER, libsize_l = libsize_450,
    batch_libsize_related = FALSE,
    p_cond = 0.4)

print("All new scenarios complete for this iteration.")
