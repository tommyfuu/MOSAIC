library(dplyr)
library(DT)
library(mia)
library(curatedMetagenomicData)


fetch_studies_from_CMG <- function (list_of_studies, abund_values, save = FALSE, save_to = NULL) {

    study_l_str <- paste(list_of_studies,collapse="|")
    study_l_abund_str <- paste("(", study_l_str, ").", abund_values, sep="")

    current_phylo_obj <-
    curatedMetagenomicData(study_l_abund_str, dryrun = FALSE) |>
    mergeData() |>
    makePhyloseqFromTreeSummarizedExperiment(abund_values = abund_values)

    if(save){
        # currently can only save otu table, metadata, and taxonomy info
        write.csv(otu_table(save_to), paste(save_to, "/otu_table_", study_l_str, '.csv'), row.names = TRUE)
        write.csv(sample_table(save_to), paste(save_to, "/sample_table_", study_l_str, '.csv'), row.names = TRUE)
        write.csv(tax_table(save_to), paste(save_to, "/tax_table_", study_l_str, '.csv'), row.names = TRUE)
    }
    return(current_phylo_obj)
}