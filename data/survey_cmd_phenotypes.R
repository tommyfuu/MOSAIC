source(Find(file.exists, c("paths.R", "../paths.R", "../../paths.R")))
library(curatedMetagenomicData)
library(dplyr)

# Survey all disease phenotypes in curatedMetagenomicData gut studies,
# reporting how many studies and samples exist per disease vs healthy comparison.

meta <- sampleMetadata

# Focus on gut microbiome (body site = stool)
gut <- meta %>% filter(body_site == "stool")

cat("=== Total gut samples:", nrow(gut), "===\n")
cat("=== Total studies:", length(unique(gut$study_name)), "===\n\n")

# For each disease (excluding healthy/control), count:
#   - number of studies with that disease AND healthy controls
#   - number of disease samples
#   - number of healthy samples in those studies
diseases <- gut %>%
    filter(!disease %in% c("healthy", "control", NA)) %>%
    pull(disease) %>%
    unique() %>%
    sort()

results <- lapply(diseases, function(d) {
    # studies that have both this disease and healthy samples
    studies_with_disease <- gut %>%
        filter(disease == d) %>%
        pull(study_name) %>%
        unique()

    studies_with_healthy <- gut %>%
        filter(study_name %in% studies_with_disease,
               disease %in% c("healthy", "control")) %>%
        pull(study_name) %>%
        unique()

    n_shared_studies <- length(studies_with_healthy)
    n_disease_samples <- gut %>%
        filter(disease == d, study_name %in% studies_with_healthy) %>%
        nrow()
    n_healthy_samples <- gut %>%
        filter(study_name %in% studies_with_healthy,
               disease %in% c("healthy", "control")) %>%
        nrow()

    data.frame(
        disease           = d,
        n_studies         = n_shared_studies,
        n_disease_samples = n_disease_samples,
        n_healthy_samples = n_healthy_samples,
        total_samples     = n_disease_samples + n_healthy_samples
    )
})

results_df <- do.call(rbind, results) %>%
    filter(n_studies >= 2) %>%          # at least 2 studies (meaningful for batch correction)
    arrange(desc(n_studies), desc(total_samples))

cat("=== Diseases with >= 2 studies having both disease + healthy samples ===\n")
print(results_df, row.names = FALSE)

cat("\n=== Full table (all diseases, sorted by n_studies) ===\n")
full_df <- do.call(rbind, results) %>% arrange(desc(n_studies), desc(total_samples))
print(full_df, row.names = FALSE)

write.csv(results_df, file.path(REPO_ROOT, "data/cmd_phenotype_survey.csv"), row.names = FALSE)
cat("\nSaved to data/cmd_phenotype_survey.csv\n")
