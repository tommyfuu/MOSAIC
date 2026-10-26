# Paths shared by the pipeline scripts; machine-specific values live in paths.yml.
# Scripts load it with: source(Find(file.exists, c("paths.R", "../paths.R", "../../paths.R")))
.paths_yml <- Find(file.exists, c("paths.yml", "../paths.yml", "../../paths.yml"))
if (is.null(.paths_yml)) stop("paths.yml not found; run from the repository or one of its subfolders")
REPO_ROOT <- normalizePath(dirname(.paths_yml))

.read_paths_yml <- function(path) {
    values <- list()
    for (line in readLines(path, warn = FALSE)) {
        line <- trimws(sub("#.*", "", line))
        if (grepl(":", line)) {
            key <- trimws(sub(":.*", "", line))
            value <- gsub("^['\"]|['\"]$", "", trimws(sub("^[^:]*:", "", line)))
            values[[key]] <- value
        }
    }
    values
}

.paths_values <- .read_paths_yml(.paths_yml)
SIM_OUTPUT_ROOT <- sub("/+$", "", Sys.getenv("SIM_OUTPUT_ROOT", .paths_values[["sim_output_root"]]))
BENCH_DIR <- file.path(REPO_ROOT, "benchmark")
DATA_DIR <- file.path(REPO_ROOT, "data")
