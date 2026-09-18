# Generic ADNI plots and tables
# =============================
#
# This script contains one main function, `run_adni_plots_and_tables()`.
# Unlike plots_and_tables.R, it does not load a particular .Rdata file or
# assume that a specific clustering run is being used.  The caller supplies
# the cluster assignments, while all ADNI data are loaded through the
# functions in utils/adni_data_loaders.R.
#
# Expected cluster input:
#   * a data frame with columns `RID` and `Cluster`, or
#   * a named vector whose names are RIDs and whose values are clusters.
#
# Example:
#   assignments <- data.frame(RID = adniLTC@RID, Cluster = adniLTC@Cluster)
#   results <- run_adni_plots_and_tables(assignments,
#                                        output_dir = "~/R/EDAP-data/plots/generic")
#
# The returned list contains the combined descriptive table, statistical
# tables, and ggplot objects.  If output_dir is supplied, CSV and PNG files
# are also written there.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
})

# Locate the loader whether this script is run from the repository root or
# from the cluster-analysis directory.
loader_path <- if (file.exists("utils/adni_data_loaders.R")) {
  "utils/adni_data_loaders.R"
} else {
  file.path("..", "utils", "adni_data_loaders.R")
}
source(loader_path)

# Convert either supported cluster-input format into a standard data frame.
as_cluster_table <- function(cluster_assignments) {
  if (is.data.frame(cluster_assignments)) {
    required <- c("RID", "Cluster")
    if (!all(required %in% names(cluster_assignments))) {
      stop("A cluster data frame must contain columns named RID and Cluster.")
    }
    # Preserve an optional time_shift column.  It is not required, but if it
    # is present it allows the onset-age plot to use model-derived time shifts.
    keep <- intersect(c("RID", "Cluster", "time_shift"), names(cluster_assignments))
    return(cluster_assignments %>% select(all_of(keep)) %>% distinct(RID, .keep_all = TRUE))
  }

  if (!is.null(names(cluster_assignments))) {
    return(tibble(RID = as.numeric(names(cluster_assignments)),
                  Cluster = unname(cluster_assignments)))
  }

  stop("cluster_assignments must be a data frame or a named vector.")
}

# Run a one-way ANOVA and Tukey post-hoc comparisons when there are enough
# observations.  Returning empty tables makes the rest of the report robust
# to a cohort with only one cluster or too little data.
anova_tukey <- function(data, response) {
  data <- data %>%
    select(Cluster, value = all_of(response)) %>%
    filter(!is.na(Cluster), !is.na(value))

  summary <- data %>%
    group_by(Cluster) %>%
    summarise(mean = mean(value), sd = sd(value), n = n(), .groups = "drop")

  if (n_distinct(data$Cluster) < 2) {
    return(list(summary = summary, tukey = tibble()))
  }

  fit <- aov(value ~ Cluster, data = data)
  tukey <- as.data.frame(TukeyHSD(fit)$Cluster) %>%
    tibble::rownames_to_column("comparison") %>%
    separate(comparison, c("Cluster1", "Cluster2"), sep = "-") %>%
    arrange(Cluster1, Cluster2)

  list(summary = summary, tukey = tukey)
}

# Summarise a categorical variable as a percentage and count in each cluster.
categorical_summary <- function(data, variable) {
  data %>%
    filter(!is.na(Cluster), !is.na(.data[[variable]])) %>%
    count(Cluster, value = .data[[variable]], name = "n") %>%
    group_by(Cluster) %>%
    mutate(percent = 100 * n / sum(n)) %>%
    ungroup()
}

# Create the generic report.  Each loader argument is exposed so that the
# caller can choose the longitudinal/normalisation settings for MRI data.
run_adni_plots_and_tables <- function(cluster_assignments,
                                      output_dir = NULL,
                                      only_vol = TRUE,
                                      filter_n = 0,
                                      normalize = TRUE) {
  clusters <- as_cluster_table(cluster_assignments)
  clusters$RID <- as.numeric(clusters$RID)

  # Load all source data using the central ADNI loader functions.
  mri <- ucsf_longitudinal_all(only_vol = only_vol,
                               filter_n = filter_n,
                               normalize = normalize) %>%
    inner_join(clusters, by = "RID")
  # The MRI loader provides elapsed time from the ADNI baseline, but not a
  # model-specific time shift.  Use zero when no shift was supplied.
  if (!"time_shift" %in% names(mri)) mri$time_shift <- 0
  demographics <- get_demographics() %>% inner_join(clusters, by = "RID")
  diagnoses <- get_diagnoses() %>% inner_join(clusters, by = "RID")
  amyloid <- get_ab_df() %>% inner_join(clusters, by = "RID")
  tau_pet <- get_tau_pet() %>% inner_join(clusters, by = "RID")

  # Keep one participant-level record for participant characteristics.
  participant_demo <- demographics %>%
    group_by(RID) %>%
    summarise(across(c(AGE, PTGENDER, PTEDUCAT, APOE4), first),
              Cluster = first(Cluster), .groups = "drop")

  # Estimate age of onset by subtracting the time shift from baseline age.
  onset <- mri %>%
    group_by(RID) %>%
    summarise(time_shift = first(time_shift[!is.na(time_shift)], default = 0),
              Cluster = first(Cluster), .groups = "drop") %>%
    left_join(participant_demo %>% select(RID, AGE), by = "RID") %>%
    mutate(onset_age = AGE - time_shift)

  # Basic participant-level summaries used in the main descriptive table.
  mri_count <- mri %>% count(RID, name = "n_mri") %>% inner_join(clusters, by = "RID")
  demo_summary <- participant_demo %>%
    group_by(Cluster) %>%
    summarise(n = n(),
              age_mean = mean(AGE, na.rm = TRUE),
              age_sd = sd(AGE, na.rm = TRUE),
              onset_mean = mean(onset$onset_age[match(RID, onset$RID)], na.rm = TRUE),
              education_mean = mean(PTEDUCAT, na.rm = TRUE),
              education_sd = sd(PTEDUCAT, na.rm = TRUE),
              male_percent = 100 * mean(PTGENDER == "Male", na.rm = TRUE),
              .groups = "drop")
  mri_summary <- mri_count %>% group_by(Cluster) %>%
    summarise(mri_mean = mean(n_mri), mri_sd = sd(n_mri), .groups = "drop")

  # Baseline and highest observed diagnosis for each participant.
  diagnosis_summary <- diagnoses %>%
    group_by(RID) %>%
    summarise(DX_bl = first(na.omit(DX.bl), default = NA),
              DX_highest = max(DIAGNOSIS, na.rm = TRUE),
              Cluster = first(Cluster), .groups = "drop") %>%
    mutate(CN2CI = if_else(DX_bl == "CN", DX_highest != "CN", NA),
           MCI2AD = if_else(DX_bl == "MCI", DX_highest == "Dementia", NA))
  diagnosis_counts <- categorical_summary(diagnosis_summary, "DX_bl")

  # Participant-level APOE and amyloid summaries.
  apoe_summary <- categorical_summary(participant_demo, "APOE4")
  amyloid_summary <- amyloid %>%
    group_by(RID) %>%
    summarise(A4240_ms = min(A4240.ms, na.rm = TRUE),
              A4240_re = min(A4240.re, na.rm = TRUE),
              AB_positive = any(AB_any, na.rm = TRUE),
              Cluster = first(Cluster), .groups = "drop") %>%
    mutate(across(c(A4240_ms, A4240_re), ~ ifelse(is.infinite(.x), NA, .x)))
  amyloid_cluster_summary <- amyloid_summary %>% group_by(Cluster) %>%
    summarise(across(c(A4240_ms, A4240_re),
                     list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE))),
              amyloid_positive_percent = 100 * mean(AB_positive, na.rm = TRUE),
              .groups = "drop")

  # Tau-PET summary: use the maximum observed value per participant.
  tau_summary <- tau_pet %>%
    group_by(RID) %>%
    summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)),
              Cluster = first(Cluster), .groups = "drop") %>%
    mutate(across(where(is.numeric), ~ ifelse(is.infinite(.x), NA, .x)))
  tau_cluster_summary <- tau_summary %>% group_by(Cluster) %>%
    summarise(across(ends_with("SUVR"),
                     list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE))),
              .groups = "drop")

  # Combine descriptive summaries into one easy-to-export table.
  descriptive_table <- demo_summary %>%
    left_join(mri_summary, by = "Cluster") %>%
    left_join(amyloid_cluster_summary, by = "Cluster") %>%
    left_join(tau_cluster_summary, by = "Cluster")

  # Continuous-variable post-hoc tests.
  age_test <- anova_tukey(participant_demo, "AGE")
  onset_test <- anova_tukey(onset, "onset_age")
  education_test <- anova_tukey(participant_demo, "PTEDUCAT")

  # A simple, generic chi-squared test for baseline diagnosis and APOE.
  categorical_test <- function(data, variable) {
    data <- data %>% filter(!is.na(Cluster), !is.na(.data[[variable]]))
    if (n_distinct(data$Cluster) < 2 || n_distinct(data[[variable]]) < 2) return(tibble())
    result <- chisq.test(table(data$Cluster, data[[variable]]))
    tibble(variable = variable, statistic = unname(result$statistic), p_value = result$p.value)
  }
  categorical_tests <- bind_rows(categorical_test(participant_demo, "APOE4"),
                                 categorical_test(diagnosis_summary, "DX_bl"),
                                 categorical_test(amyloid_summary, "AB_positive"))

  # Create the plots.  They are returned, so callers can customise or arrange
  # them in a manuscript-specific layout.
  p_age <- ggplot(participant_demo, aes(factor(Cluster), AGE, fill = factor(Cluster))) +
    geom_boxplot(outlier.shape = NA, alpha = 0.25) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.7) +
    theme_classic() + labs(x = "Cluster", y = "Age") + guides(fill = "none")
  p_onset <- ggplot(onset, aes(factor(Cluster), onset_age, fill = factor(Cluster))) +
    geom_boxplot(outlier.shape = NA, alpha = 0.25) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.7) +
    theme_classic() + labs(x = "Cluster", y = "Estimated age of onset") + guides(fill = "none")
  p_diagnosis <- ggplot(diagnosis_counts,
                        aes(factor(Cluster), percent, fill = value)) +
    geom_col(position = "stack") + theme_classic() +
    labs(x = "Cluster", y = "Participants (%)", fill = "Baseline diagnosis")
  p_mri <- ggplot(mri, aes(Years, Cluster, group = RID, colour = factor(Cluster))) +
    geom_line(alpha = 0.35) + geom_point(alpha = 0.6) + theme_classic() +
    labs(x = "Years since baseline", y = "Cluster", colour = "Cluster")

  plots <- list(age = p_age, onset = p_onset, diagnosis = p_diagnosis, mri = p_mri)

  # Optionally write machine-readable tables and plots to disk.
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(descriptive_table, file.path(output_dir, "descriptive_table.csv"), row.names = FALSE)
    write.csv(categorical_tests, file.path(output_dir, "categorical_tests.csv"), row.names = FALSE)
    write.csv(age_test$tukey, file.path(output_dir, "age_tukey.csv"), row.names = FALSE)
    write.csv(onset_test$tukey, file.path(output_dir, "onset_tukey.csv"), row.names = FALSE)
    write.csv(education_test$tukey, file.path(output_dir, "education_tukey.csv"), row.names = FALSE)
    ggsave(file.path(output_dir, "age.png"), p_age, width = 5, height = 4, dpi = 300)
    ggsave(file.path(output_dir, "onset.png"), p_onset, width = 5, height = 4, dpi = 300)
    ggsave(file.path(output_dir, "diagnosis.png"), p_diagnosis, width = 6, height = 4, dpi = 300)
    ggsave(file.path(output_dir, "mri_visits.png"), p_mri, width = 7, height = 5, dpi = 300)
  }

  list(data = list(mri = mri, demographics = demographics, diagnoses = diagnoses,
                   amyloid = amyloid, tau_pet = tau_pet),
       tables = list(descriptive = descriptive_table,
                     diagnosis = diagnosis_counts,
                     apoe = apoe_summary,
                     amyloid = amyloid_cluster_summary,
                     tau_pet = tau_cluster_summary,
                     categorical_tests = categorical_tests,
                     age_tukey = age_test$tukey,
                     onset_tukey = onset_test$tukey,
                     education_tukey = education_test$tukey),
       plots = plots)
}
