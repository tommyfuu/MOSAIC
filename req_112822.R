install.packages("xtable")
install.packages("vegan")
install.packages("doParallel")
install.packages("dplyr")
install.packages("DT")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("sva")
BiocManager::install("limma")
BiocManager::install("MMUPHin")
BiocManager::install("curatedMetagenomicData")
BiocManager::install("mia")
BiocManager::install("phyloseq")
BiocManager::install("mixOmics")

install.packages("remotes")
remote::install_github("wdl2459/ConQuR")

