# Shared MIDASim simulation functions for the scenario workflow (simulate_scenario.R).
# Same simulation as generate_data_MIDASim*.R (identical sequence of random draws), with the condition and
# batch taxa passed in instead of read from globals.
#
# Design (see Methods): condition C_i and batch B_i are drawn from a joint Bernoulli distribution with
# P(C=1) = p_cond, P(B=1) = p_batch and odds ratio `or`; with batch_libsize_related = TRUE, P(B_i = 1) =
# 1 / (1 + exp(-z_i)) with z_i the standardized library size. For a sample in condition c and batch b,
# MIDASim is run with v(c, b) = (1 - s_C c - s_B b) v0 + s_C c v_C + s_B b v_B for both mean vectors
# (mean relative abundance, and mean relative abundance among non-zero entries), where v_C permutes v0 among
# the condition taxa and v_B among the batch taxa.

library(bindata)
library(MIDASim)

# odds ratio -> binary correlation (from the ConQuR paper)
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

# 50 condition taxa: the two taxa closest to each of 25 evenly spaced zero proportions
condition_taxa <- function(otu_original) {
    zp <- colMeans(otu_original == 0)
    id_cond <- NULL
    for (tau in 0:24/25) id_cond <- c(id_cond, order(abs(zp - tau))[1:2])
    id_cond
}

# batch taxa: a fraction of all taxa, spread evenly across zero-proportion ranks (fraction 1 = all taxa;
# fraction 0.5 = every other taxon, i.e. 151 of 301)
batch_taxa <- function(otu_original, fraction = 1) {
    p <- ncol(otu_original)
    if (fraction >= 1) return(1:p)
    zp <- colMeans(otu_original == 0)
    k <- ceiling(fraction * p)
    order(zp)[unique(round(seq(1, p, length.out = k)))]
}

midas_simulate <- function(otu_original, n, or, cond_effect, batch_effect,
                           out_count, out_relab, out_meta, id_cond, id_batch,
                           p_cond = 0.5, p_batch = 0.5,
                           libsize_l = NULL, batch_libsize_related = FALSE) {
    p <- ncol(otu_original)
    if (is.null(libsize_l)) libsize_l <- rep(10000, n)
    print(paste("n:", n, "| or:", or, "| cond_effect:", cond_effect, "| batch_effect:", batch_effect,
                "| p_cond:", p_cond, "| p_batch:", p_batch, "| n_batch_taxa:", length(id_batch),
                "| libsize_related:", batch_libsize_related))

    # assign samples to (condition, batch)
    if (batch_libsize_related) {
        p_batch_vec <- 1 / (1 + exp(-scale(libsize_l)))
        bin_corr    <- sapply(p_batch_vec, function(x) bincorr(or, p_cond, x))
        cond_batchid_vec <- matrix(ncol = 2, nrow = n)
        for (i in 1:n) {
            cond_batchid_vec[i, ] <- rmvbin(1, c(p_cond, p_batch_vec[i]),
                                            bincorr = (1 - bin_corr[i]) * diag(2) + bin_corr[i])
        }
    } else {
        bin_corr         <- bincorr(or, p_cond, p_batch)
        cond_batchid_vec <- rmvbin(n, c(p_cond, p_batch), bincorr = (1 - bin_corr) * diag(2) + bin_corr)
    }
    cond    <- cond_batchid_vec[, 1]
    batchid <- cond_batchid_vec[, 2]
    write.csv(data.frame(subjectid_text = sprintf("Subject_%d", 1:n),
                         batchid = ifelse(batchid == 0, "batch_0", "batch_1"),
                         cond    = ifelse(cond    == 0, "cond_0",  "cond_1")),
              file = out_meta, row.names = FALSE)

    fitted <- MIDASim.setup(otu_original, mode = 'nonparametric')
    perm_batch <- sample(id_batch, length(id_batch))
    perm_cond  <- sample(id_cond,  length(id_cond))

    rel0      <- fitted$mean.rel.abund
    rel_batch <- fitted$mean.rel.abund; rel_batch[id_batch] <- rel_batch[perm_batch]
    rel_cond  <- fitted$mean.rel.abund; rel_cond[id_cond]   <- rel_cond[perm_cond]
    rel10      <- unlist(lapply(fitted$rel.abund.1, mean))
    rel1_batch <- rel10; rel1_batch[id_batch] <- rel1_batch[perm_batch]
    rel1_cond  <- rel10; rel1_cond[id_cond]   <- rel1_cond[perm_cond]

    fitted_c0b0 <- MIDASim.modify(fitted, lib.size = libsize_l, SCAM = TRUE)
    fitted_c1b0 <- MIDASim.modify(fitted, lib.size = libsize_l, SCAM = TRUE,
        mean.rel.abund.1 = (1 - cond_effect) * rel10 + cond_effect * rel1_cond,
        mean.rel.abund   = (1 - cond_effect) * rel0  + cond_effect * rel_cond)
    fitted_c0b1 <- MIDASim.modify(fitted, lib.size = libsize_l, SCAM = TRUE,
        mean.rel.abund.1 = (1 - batch_effect) * rel10 + batch_effect * rel1_batch,
        mean.rel.abund   = (1 - batch_effect) * rel0  + batch_effect * rel_batch)
    fitted_c1b1 <- MIDASim.modify(fitted, lib.size = libsize_l, SCAM = TRUE,
        mean.rel.abund.1 = (1 - cond_effect - batch_effect) * rel10 + cond_effect * rel1_cond + batch_effect * rel1_batch,
        mean.rel.abund   = (1 - cond_effect - batch_effect) * rel0  + cond_effect * rel_cond  + batch_effect * rel_batch)

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
    write.csv(otu_sim / rowSums(otu_sim), file = out_relab, row.names = FALSE)
}

# one replicate of every (or, cond, batch) setting with cond + batch <= 1; file names as in the original scripts
simulate_grid <- function(output_root, otu_original, or_l, cond_effect_val_l, batch_effect_val_l, iter,
                          libsize_l, id_cond, id_batch, batch_libsize_related = FALSE,
                          p_cond = 0.5, p_batch = 0.5) {
    dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
    for (or in or_l) for (cond_effect_val in cond_effect_val_l) for (batch_effect_val in batch_effect_val_l) {
        if (cond_effect_val + batch_effect_val > 1) next
        tag <- paste0(or, "_", cond_effect_val, "_", batch_effect_val, "_iter_", iter)
        sink(paste0(output_root, "/ibd_150_id_cond_", tag, ".txt"));  cat(id_cond, "\n");  sink()
        sink(paste0(output_root, "/ibd_150_id_batch_", tag, ".txt")); cat(id_batch, "\n"); sink()
        out_relab <- paste0(output_root, "/ibd_150_relab_", tag, ".csv")
        if (file.exists(out_relab)) { print(paste(out_relab, "already exists, skipping")); next }
        midas_simulate(otu_original, length(libsize_l), or, cond_effect_val, batch_effect_val,
                       paste0(output_root, "/ibd_150_count_", tag, ".csv"), out_relab,
                       paste0(output_root, "/ibd_150_meta_", tag, ".csv"), id_cond, id_batch,
                       p_cond = p_cond, p_batch = p_batch, libsize_l = libsize_l,
                       batch_libsize_related = batch_libsize_related)
    }
}
