# generate_data_MIDASim_condsweep.R
#
# Condsweep scenario: OR=1.25, cond_effect in {0.005,0.01,0.02,0.05,0.1}, batch_effect=0.5,
# n=450, all 301 taxa, p_cond=p_batch=0.5; both confounding types. One iteration per job.

# Load data
source(Find(file.exists, c("paths.R", "../paths.R", "../../paths.R")))
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
GLOBAL_ITER <- if (length(args) >= 1) {
    as.integer(args[1])
} else {
    print("No iteration argument supplied; defaulting to iter = 1")
    1L
}

# Base setup
otu_original <- t(otu)
p            <- ncol(otu_original)   # 301 taxa
n_orig       <- nrow(otu_original)   # 150 samples in source data

stopifnot(p == 301)
print(paste("p =", p, "| n_orig =", n_orig))

libsize_orig <- rowSums(otu_original)

# Same as original main simulation: n=450, sample with replacement
libsize_450  <- sample(libsize_orig, n_orig * 3, replace = TRUE)

# Odds ratio -> binary correlation (from ConQuR paper)
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

# Taxa selection
zp <- colMeans(otu_original == 0)   # zero-inflation rate per taxon

# Condition taxa: 50 taxa spread evenly across zero-inflation spectrum
# (unchanged from original)
id_cond <- NULL
for (tau in 0:24/25) {
    id_cond <- c(id_cond, order(abs(zp - tau))[1:2])
}

# Batch taxa: all 301 (same as original)
id_batch <- 1:p

print(paste("id_cond length:",       length(id_cond)))
print(paste("id_batch length:",      length(id_batch)))

# Core simulation function
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

# Per-iteration wrapper
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

# Grid runner
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

# Scenario 5: cond_effect sweep at OR=1.25, batch=0.5, both confounding types
# optional 2nd arg: odds ratio (default 1.25); optional 3rd arg: relation to generate (no | yes | both, default both)
or_condsweep   <- if (length(args) >= 2) c(as.numeric(args[2])) else c(1.25)
REL_SELECT     <- if (length(args) >= 3) args[3] else "both"
cond_sweep_l   <- c(0.005, 0.01, 0.02, 0.05, 0.1)
batch_single   <- c(0.5)

base_output_root <- file.path(SIM_OUTPUT_ROOT, 'scenarios_simulation_data_052626')
dir_condsweep_norel  <- file.path(base_output_root, 'count_condsweep_no_relation')
dir_condsweep_yesrel <- file.path(base_output_root, 'count_condsweep_yes_relation')

for (d in c(dir_condsweep_norel, dir_condsweep_yesrel)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

if (REL_SELECT %in% c("no", "both")) scaled_slurm_midas_data_generation(
    dir_condsweep_norel, otu_original,
    or_condsweep, cond_sweep_l, batch_single,
    iter = GLOBAL_ITER, libsize_l = libsize_450,
    batch_libsize_related = FALSE)

if (REL_SELECT %in% c("yes", "both")) scaled_slurm_midas_data_generation(
    dir_condsweep_yesrel, otu_original,
    or_condsweep, cond_sweep_l, batch_single,
    iter = GLOBAL_ITER, libsize_l = libsize_450,
    batch_libsize_related = TRUE)

print("condsweep scenario complete for this iteration.")
