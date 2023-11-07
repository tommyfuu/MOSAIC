## R functions to clean and prune the dataset
library("dplyr")
library("phyloseq")
library("readr")
# load dataset as a phyloseq object
# OTU_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/cdi_3_microbiomeHD/cdi_schubert_results/RDP/cdi_schubert.otu_table.100.denovo.rdp_assigned'
# metadata_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/cdi_3_microbiomeHD/cdi_schubert_results/cdi_schubert.metadata.txt'
load_microbiomeHD_phyloseq <- function(OTU_oath, metadata_path){
OTU_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/benchmark/benchmarked_data/cdi_3_microbiomeHD_count_data.csv'
metadata_path = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/cdi_3_microbiomeHD/cdi_schubert_results/cdi_schubert.metadata.txt'
otu_mat<- t(read.csv(OTU_path,header=T, row.names = 1))

# turn otu_mat's rownames into a dataframe by splitting with ';__'
a = lapply(rownames(otu_mat)), function(x) strsplit(x, "?.")[[1]]
tax_table <- t(as.data.frame(do.call(cbind, lapply(rownames(otu_mat), function(x) strsplit(x, ".")[[1]]))))
raw_taxa_classes <- lapply(rownames(otu_mat)[[1]], function(x) strsplit(x, '.')[[1]])
taxa_classes <- lapply(raw_taxa_classes[[1]], function(x) strsplit(x, '__')[[1]][[1]])
colnames(tax_table) = taxa_classes

# change otu_mat's rownames to the first column of tax_table
rownames(tax_table) = rownames(otu_mat)

write.csv(tax_table, '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/test/tax_table.csv')
write.csv(otu_mat, '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/test/otu_mat.csv')


otu_mat<- read.table(OTU_path, sep="\t",header=T,check.names = F)
otu_mat<- read_csv(OTU_path, delim = "\t")
samples_df <- read_csv(metadata_path)

    # generate taxonomy table from OTU table's colnames, assign colnames as otu, Domain, Supergroup, Phylum, Class, Order, Family, Genus, Species
    tax_table <- data.frame(t(apply(otu_mat, 2, function(x) strsplit(x, ";")[[1]])))
}

dataset_clean_prune <- function(phyloseq_dataset, libsize_threshold, relab_threshold){
    # sample-wise cleaning: remove samples with less than libsize_threshold reads
    phyloseq_dataset <- prune_samples(sample_sums(phyloseq_dataset) > libsize_threshold, phyloseq_dataset)
    # taxon-wise cleaning: remove taxa with less than relab_threshold relative abundance
    phyloseq_dataset <- prune_taxa(taxa_sums(phyloseq_dataset) > relab_threshold, phyloseq_dataset)
}