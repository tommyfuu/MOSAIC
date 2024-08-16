# load from yml file
library(yaml)
current_path = getwd()
if (grepl("/benchmark", current_path)){
    Snakefile_path <- "../Snakefile"
    source('./methods_benchmarking_sim.R')
} else{
    Snakefile_path <- "./Snakefile"
    source('./benchmark/methods_benchmarking_sim.R')
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
method_l_count <- config_object$used_R_methods_count
method_l_relab <- config_object$used_R_methods_relab


datatypes <- c('count', 'relab')
relations <- c('no', 'yes')

overall_path <- paste0(output_root, "/simulate")

# apply data integration methods
for (datatype in datatypes) {
  for (relation in relations) {
    output_dir <- paste0(output_root, "/benchmark/", datatype, "_", relation, "_relation")
    # make dir and parent dir if needed
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    if (datatype == 'count') {
      count = TRUE
        method_l = method_l_count
    } else {
       count = FALSE
        method_l = method_l_relab
    }
    scaled_slurm_methods_bencharking(output_dir, overall_path, method_l, or_l, cond_effect_val_l, batch_effect_val_l, GLOBAL_ITER, count = count)

  }
}

