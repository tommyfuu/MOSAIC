## set r cran mirror
r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)

## install the required libraries
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ExperimentHub")
BiocManager::install("mia", force=TRUE)
BiocManager::install("curatedMetagenomicData", force=TRUE)
BiocManager::install("phyloseq", force=TRUE)

library(dplyr)
library(DT)
library(mia)
library(curatedMetagenomicData)

## to access a single study and convert to phyloseq
brooksB2017_metaObj = sampleMetadata |>
    filter(study_name == "BrooksB_2017") |>
    select(where(~ !any(is.na(.x))))  |>
    returnSamples("relative_abundance", rownames = "short")

brooksB2017_Phyloseq <- brooksB2017_metaObj|> 
  mia::makePhyloseqFromTreeSummarizedExperiment(assay.type = "relative_abundance")
brooksB2017_Phyloseq  <- prune_taxa(taxa_sums(brooksB2017_Phyloseq) >0, brooksB2017_Phyloseq)


### to access multiple studies and make a phyloseq object
IBD_3sets_metaObj = sampleMetadata |>
    filter(study_name %in% c("HMP_2019_ibdmdb", "IjazUZ_2017", "LiJ_2014")) |>
    select(where(~ !any(is.na(.x))))  |>
    returnSamples("relative_abundance", rownames = "short")

IBD_3sets_Phyloseq <- IBD_3sets_metaObj|> 
  mia::makePhyloseqFromTreeSummarizedExperiment(abund_values="relative_abundance")
IBD_3sets_Phyloseq   <- prune_taxa(taxa_sums(IBD_3sets_Phyloseq ) >0, IBD_3sets_Phyloseq)