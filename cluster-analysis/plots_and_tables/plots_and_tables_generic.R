# Generic ADNI and multi-cohort plots and tables
# =============================================
#
# This script provides two functions:
#
#   run_adni_plots_and_tables()
#     Creates a report for ADNI-only cluster assignments.
#     Data are loaded from utils/adni_data_loaders.R.
#
#   run_multicohort_plots_and_tables()
#     Creates a report for a multiLTC experiment containing ADNI and OASIS
#     participants.  ADNI and OASIS are loaded with their respective loader
#     files and then combined using harmonised variable names.
#
# Neither function loads a particular .Rdata file.  The caller supplies the
# cluster assignments, which makes the functions reusable across experiments.
#
# ADNI-only input:
#   * data frame with columns `RID` and `Cluster`, or
#   * named vector whose names are numeric RIDs and whose values are clusters.
#
# Multi-cohort input:
#   * a `multiLTC` object with @RID and @Cluster slots,
#   * data frame with columns `RID` and `Cluster`, or
#   * named vector whose names are prefixed IDs.
#
# Multi-cohort IDs must identify their source cohort using one of these forms:
#   * ADNI_<RID>, for example ADNI_1234
#   * OASIS_<RID>, for example OASIS_0456
#
# Optional `time_shift` can be included in an input data frame.  It is used
# when calculating estimated age of onset; if it is omitted, the shift is 0.
#
# Examples:
#   source("cluster-analysis/plots_and_tables_generic.R")
#
#   # ADNI-only analysis
#   adni_assignments <- data.frame(RID = adniLTC@RID,
#                                  Cluster = adniLTC@Cluster)
#   adni_results <- run_adni_plots_and_tables(
#     adni_assignments,
#     output_dir = "~/R/EDAP-data/plots/generic_adni"
#   )
#
#   # Multi-cohort analysis using a multiLTC object
#   multi_results <- run_multicohort_plots_and_tables(
#     multiLTC,
#     output_dir = "~/R/EDAP-data/plots/generic_multicohort"
#   )
#
# Both functions return a list containing the loaded data, summary/statistical
# tables, and ggplot objects.  If output_dir is supplied, CSV and PNG files
# are also written there.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(ggalluvial)
})

# Locate the loader whether this script is run from the repository root or
# from the cluster-analysis directory.
loader_path <- if (file.exists("utils/adni_data_loaders.R")) {
  "utils/adni_data_loaders.R"
} else {
  file.path("..", "utils", "adni_data_loaders.R")
}
source(loader_path)

# Load the OASIS functions into their own environment.  Both loader files
# define functions called get_demographics(), get_diagnoses(), and
# get_tau_pet(), so separate environments prevent one cohort from overwriting
# the other cohort's functions.
oasis_loader_path <- if (file.exists("utils/oasis_data_loaders.R")) {
  "utils/oasis_data_loaders.R"
} else {
  file.path("..", "utils", "oasis_data_loaders.R")
}
oasis_loader <- new.env()
source(oasis_loader_path, local = oasis_loader)

adni_loader <- new.env()
source(loader_path, local = adni_loader)

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
  group_vars <- intersect(c("Cohort", "Cluster"), names(data))
  data %>%
    filter(!is.na(Cluster), !is.na(.data[[variable]])) %>%
    mutate(value = .data[[variable]]) %>%
    group_by(across(all_of(c(group_vars, "value")))) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(across(all_of(group_vars))) %>%
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
    summarise(DX.bl = first(na.omit(DX.bl), default = NA),
              DX.highest = max(DIAGNOSIS, na.rm = TRUE),
              Cluster = first(Cluster), .groups = "drop") %>%
    mutate(CN2CI = if_else(DX.bl == "CN", DX.highest != "CN", NA),
           MCI2AD = if_else(DX.bl == "MCI", DX.highest == "Dementia", NA))
  diagnosis_counts <- categorical_summary(diagnosis_summary, "DX.bl")

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
                                 categorical_test(diagnosis_summary, "DX.bl"),
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

# Convert a multi-cohort cluster result into a standard table while retaining
# the cohort encoded in IDs such as ADNI_1234 and OASIS_0456.
as_multicohort_cluster_table <- function(cluster_assignments) {
  # This also allows the function to be called directly with a multiLTC S4
  # object, as used by the existing multi-cohort analysis scripts.
  if (is(cluster_assignments, "clusterObject")) {
    cluster_assignments <- tibble(
      RID = cluster_assignments@RID,
      Cluster = cluster_assignments@Cluster
    )
  }

  if (is.data.frame(cluster_assignments)) {
    if (!all(c("RID", "Cluster") %in% names(cluster_assignments))) {
      stop("A multi-cohort cluster data frame must contain RID and Cluster.")
    }
    keep <- intersect(c("RID", "Cluster", "time_shift"), names(cluster_assignments))
    clusters <- cluster_assignments %>% select(all_of(keep))
  } else if (!is.null(names(cluster_assignments))) {
    clusters <- tibble(RID = names(cluster_assignments),
                       Cluster = unname(cluster_assignments))
  } else {
    stop("cluster_assignments must be a data frame, a multiLTC object, or a named vector.")
  }

  clusters <- clusters %>%
    mutate(RID = as.character(RID),
           Cohort = case_when(
             grepl("^ADNI_", RID) ~ "ADNI",
             grepl("^OASIS_", RID) ~ "OASIS",
             TRUE ~ NA_character_
           ))

  if (any(is.na(clusters$Cohort))) {
    bad_ids <- unique(clusters$RID[is.na(clusters$Cohort)])
    stop("Every RID must start with ADNI_ or OASIS_. Invalid IDs: ",
         paste(bad_ids, collapse = ", "))
  }

  clusters %>% distinct(RID, .keep_all = TRUE)
}

# Add the cohort prefix used by multiLTC to the IDs returned by a cohort
# loader.  Keeping this operation in one helper makes the joins below easier
# to audit.
prefix_cohort_ids <- function(data, cohort, id_column = "RID") {
  data %>% mutate(RID = paste0(cohort, "_", .data[[id_column]]),
                  Cohort = cohort)
}

# Create plots and tables for a multi-cohort experiment.  The function loads
# ADNI and OASIS separately, harmonises their common variables, and then joins
# the result to the supplied cluster assignments.
run_multicohort_plots_and_tables <- function(cluster_assignments,
                                             output_dir = NULL,
                                             only_vol = TRUE,
                                             filter_n = 0,
                                             normalize = TRUE,
                                             unified_norm = TRUE,
                                             use_centiloid = FALSE) {
  clusters <- as_multicohort_cluster_table(cluster_assignments)

  # Load longitudinal MRI from both cohorts.
  multi_cohort_df <- read.csv("~/R/EDAP-data/MULTI_COHORT_4.csv", header = TRUE)
  # Filter out NACC
  multi_cohort_df <- filter_out(multi_cohort_df, Cohort == "NACC") %>%
    left_join(clusters, by = c("RID", "Cohort"))

  # Harmonise participant demographics.  OASIS uses different source column
  # names, so the conversion is explicit instead of relying on column order.
  adni_demographics <- adni_loader$get_demographics() %>%
    transmute(RID = paste0("ADNI_", RID),
              AGE,
              PTGENDER,
              PTEDUCAT,
              APOE4,
              Cohort = "ADNI")
  oasis_demographics <- oasis_loader$get_demographics() %>%
    transmute(RID = paste0("OASIS_", RID),
              AGE = AgeatEntry,
              PTGENDER = recode(as.character(GENDER), `1` = "Male", `2` = "Female"),
              PTEDUCAT = EDUC,
              APOE4,
              Cohort = "OASIS")
  demographics <- bind_rows(adni_demographics, oasis_demographics) %>%
    inner_join(clusters, by = c("RID", "Cohort"))

  # Harmonise diagnoses to the common labels used in the summary table.
  adni_diagnoses <- adni_loader$get_diagnoses() %>%
    transmute(RID = paste0("ADNI_", RID),
              DX.bl = as.character(DX.bl),
              DX = as.character(DIAGNOSIS),
              Cohort = "ADNI")
  oasis_diagnoses <- oasis_loader$get_diagnoses() %>%
    transmute(RID = paste0("OASIS_", RID),
              DX.bl = as.character(DX.bl),
              DX = as.character(DX),
              Cohort = "OASIS")
  diagnoses <- bind_rows(adni_diagnoses, oasis_diagnoses) %>%
    inner_join(clusters, by = c("RID", "Cohort")) %>%
    mutate(DX = factor(DX, levels = c("CN", "Impaired", "MCI", "Dementia"), ordered = TRUE),
           DX.bl = factor(DX.bl, levels = c("CN", "Impaired", "MCI", "Dementia"), ordered = TRUE)) 
  diagnosis_summary <- diagnoses %>%
    group_by(RID, Cohort, Cluster) %>%
    summarise(DX.bl = first(na.omit(DX.bl), default = NA_character_),
              DX.highest = case_when(
                any(DX == "Dementia", na.rm = TRUE) ~ "Dementia",
                any(DX == "MCI", na.rm = TRUE) ~ "MCI",
                any(DX == "Impaired", na.rm = TRUE) ~ "Impaired",
                any(DX == "CN", na.rm = TRUE) ~ "CN",
                TRUE ~ NA_character_),
              .groups = "drop") %>%
    mutate(
      DX.highest = factor(DX.highest, levels = c("CN", "Impaired", "MCI", "Dementia"), ordered = TRUE),
      CN2CI = if_else(DX.bl == "CN", DX.highest != "CN", NA),
      MCI2AD = if_else(DX.bl == "MCI", DX.highest == "Dementia", NA))

  # Load amyloid measures from both cohorts.  The common PET variable is
  # called PET_SUVR; ADNI CSF variables remain available in adni_amyloid.
  adni_amyloid <- adni_loader$get_ab_df() %>%
    transmute(RID = paste0("ADNI_", RID),
              tracer = as.character(TRACER),
              PET_SUVR = SUMMARY_SUVR,
              AB_positive = AB_pos.pet,
              Cohort = "ADNI")
  oasis_amyloid <- oasis_loader$get_oasis_ab(use_centiloid = use_centiloid) %>%
    transmute(RID = paste0("OASIS_", gsub("^OAS", "", OASISID)),
              tracer = as.character(tracer),
              PET_SUVR = if (use_centiloid) Centiloid_SUVR else PET_SUVR,
              AB_positive = AB_pos,
              Cohort = "OASIS") %>%
    mutate(tracer = ifelse(tracer == "AV45", "FBP", tracer))
  amyloid <- bind_rows(adni_amyloid, oasis_amyloid) %>%
    filter(!is.na(PET_SUVR)) %>%
    inner_join(clusters, by = c("RID", "Cohort"))
  amyloid_participant <- amyloid %>%
    group_by(RID, Cohort, Cluster, tracer) %>%
    summarise(PET_SUVR = max(PET_SUVR, na.rm = TRUE),
              AB_positive = any(AB_positive, na.rm = TRUE),
              .groups = "drop") %>%
    filter(tracer %in% c("FBP", "PIB"))
  
  # Specific Cohort cut-offs
  cutoffs <- data.frame(
    Cohort = c("ADNI", "OASIS", "ADNI", "OASIS"),
    tracer = c("FBP", "FBP", "PIB", "PIB"),
    cutoff = if (use_centiloid) c(1.11, 20.6, 1.21, 16.4) else c(1.11, 1.19, 1.21, 1.42)
  )

  # Tau-PET is returned in a harmonised long-ish table.  Region names differ
  # between cohorts, so the summary only aggregates shared numeric SUVR fields.
  adni_tau <- adni_loader$get_tau_pet() %>%
    prefix_cohort_ids("ADNI")
  oasis_tau <- oasis_loader$get_tau_pet() %>%
    prefix_cohort_ids("OASIS", id_column = "OASISID")
  tau_pet <- bind_rows(adni_tau, oasis_tau) %>%
    inner_join(clusters, by = c("RID", "Cohort"))
  tau_fields <- names(tau_pet)[grepl("SUVR", names(tau_pet), ignore.case = TRUE)]
  tau_summary <- tau_pet %>%
    group_by(RID, Cohort, Cluster) %>%
    summarise(across(all_of(tau_fields), ~ max(.x, na.rm = TRUE)), .groups = "drop") %>%
    mutate(across(all_of(tau_fields), ~ ifelse(is.infinite(.x), NA, .x)))

  # Calculate participant-level descriptive summaries, retaining Cohort so
  # pooled and cohort-specific results can both be inspected.
  participant_demo <- demographics %>%
    group_by(RID, Cohort, Cluster) %>%
    summarise(across(c(AGE, PTGENDER, PTEDUCAT, APOE4), first), .groups = "drop")
  onset <- multi_cohort_df %>%
    group_by(RID, Cohort, Cluster) %>%
    summarise(time_shift = first(time_shift[!is.na(time_shift)], default = 0),
              .groups = "drop") %>%
    left_join(participant_demo %>% select(RID, AGE), by = "RID") %>%
    mutate(onset_age = AGE - time_shift)
  mri_count <- multi_cohort_df %>% count(RID, Cohort, Cluster, name = "n_mri")

  descriptive_table <- participant_demo %>%
    group_by(Cohort, Cluster) %>%
    summarise(n = n(),
              age_mean = mean(AGE, na.rm = TRUE),
              age_sd = sd(AGE, na.rm = TRUE),
              education_mean = mean(PTEDUCAT, na.rm = TRUE),
              education_sd = sd(PTEDUCAT, na.rm = TRUE),
              male_percent = 100 * mean(PTGENDER == "Male", na.rm = TRUE),
              .groups = "drop") %>%
    left_join(onset %>% group_by(Cohort, Cluster) %>%
                summarise(onset_mean = mean(onset_age, na.rm = TRUE),
                          onset_sd = sd(onset_age, na.rm = TRUE), .groups = "drop"),
              by = c("Cohort", "Cluster")) %>%
    left_join(mri_count %>% group_by(Cohort, Cluster) %>%
                summarise(mri_mean = mean(n_mri), mri_sd = sd(n_mri), .groups = "drop"),
              by = c("Cohort", "Cluster")) %>%
    left_join(amyloid_participant %>% group_by(Cohort, Cluster) %>%
                summarise(pet_mean = mean(PET_SUVR, na.rm = TRUE),
                          pet_sd = sd(PET_SUVR, na.rm = TRUE),
                          amyloid_positive_percent = 100 * mean(AB_positive, na.rm = TRUE),
                          .groups = "drop"),
              by = c("Cohort", "Cluster"))

  #diagnosis_counts <- categorical_summary(diagnosis_summary, "DX.bl")
  diagnosis_counts <- diagnosis_summary %>% 
    mutate(DX.bl = recode(DX.bl, "Dementia" = "AD"), 
           DX.highest = recode(DX.highest, "Dementia" = "AD")) %>%
    group_by(Cluster) %>%
    count(DX.bl, DX.highest, name = "Freq")
  apoe_counts <- categorical_summary(participant_demo, "APOE4")

  # Pooled post-hoc tests are useful for the multi-cohort experiment, while
  # the descriptive table and plots above retain the cohort breakdown.
  age_test <- anova_tukey(participant_demo, "AGE")
  onset_test <- anova_tukey(onset, "onset_age")
  education_test <- anova_tukey(participant_demo, "PTEDUCAT")

  categorical_test <- function(data, variable) {
    data <- data %>% filter(!is.na(Cluster), !is.na(.data[[variable]]))
    if (n_distinct(data$Cluster) < 2 || n_distinct(data[[variable]]) < 2) return(tibble())
    result <- chisq.test(table(data$Cluster, data[[variable]]))
    tibble(variable = variable, statistic = unname(result$statistic), p_value = result$p.value)
  }
  categorical_tests <- bind_rows(
    categorical_test(participant_demo, "APOE4"),
    categorical_test(diagnosis_summary, "DX.bl"),
    categorical_test(amyloid_participant, "AB_positive")
  )

  # Cohort-faceted plots make it possible to see whether a cluster pattern is
  # driven by one source cohort.
  p_age <- ggplot(participant_demo, aes(factor(Cluster), AGE, fill = Cluster)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.3) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.6, aes(color = Cohort)) +
    scale_color_paletteer_d("ggthemes::wsj_red_green") +
    scale_fill_paletteer_d("ggthemes::Tableau_10") +
    theme_classic() +
    labs(x = "Cluster", y = "Age")
  
  p_onset <- ggplot(onset, aes(factor(Cluster), onset_age, fill = Cluster)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.3) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.6, aes(color = Cohort)) +
    scale_color_paletteer_d("ggthemes::wsj_red_green") +
    scale_fill_paletteer_d("ggthemes::Tableau_10") +
    theme_classic() +
    labs(x = "Cluster", y = "Estimated age of onset")
  
  plot_of_ages <- ggarrange(p_age, p_onset, nrow = 1, common.legend = TRUE, legend = "right")
  
  p_diagnosis <- ggplot(diagnosis_counts,
                        aes(factor(Cluster), percent, fill = value)) +
    scale_fill_paletteer_d("ggthemes::Tableau_10") +
    geom_col() + facet_wrap(~ Cohort) + theme_classic() +
    labs(x = "Cluster", y = "Participants (%)", fill = "Baseline diagnosis")
  
  gp_alluv_list <- lapply(unique(clusters$Cluster), function(c) {
    gp_alluv <- diagnosis_counts %>% filter(Cluster == c) %>%
      ggplot(aes(axis1 = DX.bl, axis2 = DX.highest, y = Freq)) +
      geom_alluvium(aes(fill = DX.highest), width = 1/12) +
      geom_stratum(width = 1/8, fill = "grey90", color = "grey40") +
      geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
      scale_x_discrete(limits = c("Baseline", "Final diagnosis"), expand = c(.1, .05)) +
      labs(x = NULL, y = "Count", title = paste0("Cluster ", c)) +
      scale_fill_brewer(palette = "YlOrRd", name = "Diagnosis") +
      theme_classic()
  })
  
  gp_alluv <- ggarrange(plotlist=gp_alluv_list, ncol=2, nrow=2, common.legend = TRUE, legend="right")
  #plot(gp_alluv)
  
  p_amyloid <- ggplot(amyloid_participant,
                      aes(factor(Cluster), PET_SUVR, fill = Cluster)) +
    geom_violin(alpha = 0.35, position = position_dodge(width = 0.8)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
               alpha = 0.6) + facet_wrap(tracer ~ Cohort, scales = "free_y") +
    geom_hline(
      data = cutoffs,
      aes(yintercept = cutoff),
      linetype = "dashed",
      inherit.aes = FALSE
    ) +
    scale_fill_paletteer_d("ggthemes::Tableau_10") +
    theme_classic() + labs(x = "Cluster", y = "Amyloid PET measure")
  
  # Do not include this, not sure what it adds
  p_mri <- ggplot(multi_cohort_df, aes(Years, Cluster, group = RID, colour = Cohort)) +
    geom_line(alpha = 0.3) + geom_point(alpha = 0.5) + facet_wrap(~ Cohort) +
    theme_classic() + labs(x = "Years since baseline", y = "Cluster", colour = "Cohort")
  
  
  
  plots <- list(age = p_age, onset = p_onset, p_of_ages = plot_of_ages, diagnosis = gp_alluv,
                amyloid = p_amyloid)

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    #write.csv(descriptive_table, file.path(output_dir, "multicohort_descriptive_table.csv"), row.names = FALSE)
    #write.csv(diagnosis_counts, file.path(output_dir, "multicohort_diagnosis.csv"), row.names = FALSE)
    #write.csv(apoe_counts, file.path(output_dir, "multicohort_apoe.csv"), row.names = FALSE)
    #write.csv(categorical_tests, file.path(output_dir, "multicohort_categorical_tests.csv"), row.names = FALSE)
    #write.csv(age_test$tukey, file.path(output_dir, "multicohort_age_tukey.csv"), row.names = FALSE)
    #write.csv(onset_test$tukey, file.path(output_dir, "multicohort_onset_tukey.csv"), row.names = FALSE)
    #write.csv(education_test$tukey, file.path(output_dir, "multicohort_education_tukey.csv"), row.names = FALSE)
    #write.csv(tau_summary, file.path(output_dir, "multicohort_tau_summary.csv"), row.names = FALSE)
    ggsave(file.path(output_dir, "multicohort_age.png"), p_age, width = 5, height = 4, dpi = 300)
    ggsave(file.path(output_dir, "multicohort_onset.png"), p_onset, width = 5, height = 4, dpi = 300)
    ggsave(file.path(output_dir, "multicohort_diagnosis.png"), gp_alluv, width = 8, height = 6, dpi = 300)
    ggsave(file.path(output_dir, "multicohort_amyloid.png"), p_amyloid, width = 8, height = 6, dpi = 300)
    ggsave(file.path(output_dir, "multicohort_ages.png"), p_amyloid, width = 9, height = 4, dpi = 500)
  }

  list(data = list(mri = mri, demographics = demographics,
                   diagnoses = diagnoses, amyloid = amyloid,
                   tau_pet = tau_pet),
       tables = list(descriptive = descriptive_table,
                     diagnosis = diagnosis_counts,
                     apoe = apoe_counts,
                     amyloid = amyloid_participant,
                     tau_pet = tau_summary,
                     categorical_tests = categorical_tests,
                     age_tukey = age_test$tukey,
                     onset_tukey = onset_test$tukey,
                     education_tukey = education_test$tukey),
       plots = plots)
}


# Multi-cohort analysis using a multiLTC object
run <- "exp_km_ab_ao_2"
load(paste("~/R/EDAP-data/LTC_MC/new/", run, ".Rdata", sep = ""))

multi_results <- run_multicohort_plots_and_tables(
  multiLTC,
  output_dir = "~/R/EDAP-data/plots/LTC_MC"
  )




