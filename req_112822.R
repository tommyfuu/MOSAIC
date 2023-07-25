# open R in terminal, then copy and paste the following command one by one to install the necessary packages
install.packages("xtable")
install.packages("vegan")
install.packages("doParallel")
install.packages("dplyr")
install.packages("DT")
install.packages("gmp")
install.packages("Rmpfr")
install.packages("coda.base")
install.packages("ade4")
install.packages("compositions")
install.packages("bindata") # need to install in R
install.packages("FDboost") # need to install in R

install.packages("BiocManager") # need to install in R

BiocManager::install("sva") # need to install in R
BiocManager::install("annotate") # need to install in R
BiocManager::install("limma")
BiocManager::install("MMUPHin")
BiocManager::install("curatedMetagenomicData")
BiocManager::install("mia")
BiocManager::install("phyloseq")
BiocManager::install("mixOmics")


install.packages("remotes")
remotes::install_github("r-lib/xml2") 
remotes::install_github("mengyu-he/MIDAS")
# remotes::install_github("wdl2459/ConQuR", force=TRUE)