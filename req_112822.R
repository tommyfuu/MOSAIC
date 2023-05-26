# open R in terminal, then copy and paste the following command one by one to install the necessary packages
install.packages("xtable")
install.packages("vegan")
install.packages("doParallel")
install.packages("dplyr")
install.packages("DT")
install.packages("gmp")
install.packages("Rmpfr")


install.packages("BiocManager")

BiocManager::install("sva")
BiocManager::install("limma")
BiocManager::install("MMUPHin")
BiocManager::install("curatedMetagenomicData")
BiocManager::install("mia")
BiocManager::install("phyloseq")
BiocManager::install("mixOmics")
BiocManager::install("sparseDOSSA")

install.packages("remotes")
remotes::install_github("wdl2459/ConQuR", force=TRUE)
