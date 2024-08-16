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
used_R_methods <- config_object$used_R_methods

overall_path <- paste0(output_root, "/simulate")

datatypes <- c('count', 'relab')
relations <- c('no', 'yes')

# apply data integration methods
for (datatype in datatypes) {
  for (relation in relations) {
    output_dir <- paste0(overall_path, "/benchmark/", datatype, "_", relation, "_relation/out_", datatype, "_", relation, "_relation_iter", GLOBAL_ITER)
    dir.create(output_dir, showWarnings = FALSE)
    if (datatype == 'count') {
      count = TRUE
    } else {
       count = FALSE
    }
    scaled_slurm_methods_bencharking(output_dir, overall_path, method_l, or_l, cond_effect_val_l, batch_effect_val_l, GLOBAL_ITER, count = count)

  }
}

