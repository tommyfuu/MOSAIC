library(dplyr)
library(DT)
library(mia)
library(curatedMetagenomicData)
library(phyloseq)

fetch_studies_from_CMG <- function (list_of_studies, list_of_conditions, abund_values, counts = FALSE, save = FALSE, save_to = NULL) {

    study_l_str <- paste(list_of_studies,collapse="|")
    study_l_abund_str <- paste("(", study_l_str, ").", abund_values, sep="")

    # metaObj <- sampleMetadata |>
    #     filter(study_name %in% study_l_abund_str) |>
    #     filter(disease %in% c("IBD", "healthy")) |>
    #     returnSamples("relative_abundance", rownames = "short")

    # current_phylo_obj <- metaObj |> 
    #     mia::makePhyloseqFromTreeSummarizedExperiment(abund_values="relative_abundance")
    # current_phylo_obj <-
    #     returnSamples("relative_abundance") |>
    #     curatedMetagenomicData(study_l_abund_str, dryrun = FALSE, counts) |>
    #     mergeData() |>
    #     makePhyloseqFromTreeSummarizedExperiment(abund_values = abund_values)

    current_phylo_obj <-
        filter(sampleMetadata, disease %in% list_of_conditions & study_name %in% list_of_studies) |>
        select(where(~ !all(is.na(.x)))) |>
        returnSamples("relative_abundance") |>
        makePhyloseqFromTreeSummarizedExperiment(abund_values = abund_values)

    print(sample_data(current_phylo_obj))
    if(save){
        # currently can only save otu table, metadata, and taxonomy info
        write.csv(otu_table(current_phylo_obj), paste(save_to, "/otu_table_", study_l_str, '.csv', sep=""), row.names = TRUE)
        write.csv(sample_data(current_phylo_obj), paste(save_to, "/sample_table_", study_l_str, '.csv', sep=""), row.names = TRUE)
        write.csv(tax_table(current_phylo_obj), paste(save_to, "/tax_table_", study_l_str, '.csv', sep=""), row.names = TRUE)
    }
    return(current_phylo_obj)
}

# want to save real counts data
# a = fetch_studies_from_CMG(c("HMP_2019_ibdmdb", "IjazUZ_2017", "LiJ_2014"), c("IBD", 'healthy'), "relative_abundance", counts = TRUE, save = TRUE, save_to = '/home/fuc/harmonicMic/data/ibd_3_CMD')
# # PROBLEMATIC # b = fetch_studies_from_CMG(c("FrankelAE_2017", "GopalakrishnanV_2018", "LeeKA_2022", "MatsonV_2018", "WindTT_2020"), c("melanoma", "metastases"), "relative_abundance", counts = TRUE, save = TRUE, save_to = '/home/fuc/harmonicMic/data/melanoma_5_CMD')
# b = fetch_studies_from_CMG(c("FengQ_2015", "HanniganGD_2017", "ThomasAM_2018a", "YachidaS_2019", "ZellerG_2014"), c("adenoma", "healthy"), "relative_abundance", counts = TRUE, save = TRUE, save_to = '/home/fuc/harmonicMic/data/adenoma_5_CMD')
# c = fetch_studies_from_CMG(c("Castro-NallarE_2015", "FengQ_2015", "GhensiP_2019", "HMP_2019_t2d", "KarlssonFH_2013", "LiJ_2014", "MetaCardis_2020_a", "QinJ_2012", "SankaranarayananK_2015", "ThomasAM_2018a", "YuJ_2015"), c("T2D", 'healthy'),"relative_abundance", counts = TRUE, save = TRUE, save_to = '/home/fuc/harmonicMic/data/T2D_10_CMD')
d = fetch_studies_from_CMG(c("FengQ_2015", "GuptaA_2019"), c("CRC", 'healthy'),"relative_abundance", counts = TRUE, save = TRUE, save_to = '/athena/linglab/scratch/chf4012/mic_bc_benchmark/data/test')