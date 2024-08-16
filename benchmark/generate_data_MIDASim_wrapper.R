# load from yml file
library(yaml)
current_path = getwd()
if (grepl("/benchmark", current_path)){
    Snakefile_path <- "../Snakefile"
    source('./generate_data_MIDASim.R')
} else{
    Snakefile_path <- "./Snakefile"
    source('./benchmark/generate_data_MIDASim.R')
}

conn <- file(Snakefile_path,open="r")
linn <-readLines(conn)
for (i in 1:length(linn)){
    if (grepl("configfile", linn[i])){
        configfile_path = strsplit(linn[i], " ")[[1]][2]
        configfile_path = gsub("'", "", configfile_path)
        configfile_path = gsub("\"", "", configfile_path)
        break
    }
}
# close connection
close(conn)

if (grepl("/benchmark", current_path)){
    config_object <- yaml.load_file(paste0("../", configfile_path))
} else{
    config_object <- yaml.load_file(configfile_path)
}



# load or_l, cond_effect_val_l, and batch_effect_val_l
or_l <- config_object$or_l
cond_effect_val_l <- config_object$cond_effect_val_l
batch_effect_val_l <- config_object$batch_effect_val_l
output_root <- config_object$sim_output_root
GLOBAL_ITER <- config_object$iter

output_root_simulate <- paste0(output_root, "/simulate")
dir.create(output_root_simulate, showWarnings = FALSE, recursive = TRUE)

# generate data
scaled_slurm_midas_data_generation(output_root_simulate, otu_original, n, or_l, cond_effect_val_l, batch_effect_val_l, iter=GLOBAL_ITER, batch_libsize_related = FALSE, libsize_l=sampled_libsize_l)