# Pairwise statistical tests for ADNI / multi-cohort cluster analyses
# ====================================================================
#
# This script deliberately contains statistical tests only.  It creates one
# row for every pair of clusters and returns the following columns:
#
#   Cluster pair, Age, Onset, Education, APOE4, DX, AB-CSF, AB-PET
#
# Numeric variables (Age, Onset, and Education) are tested with Tukey's HSD
# following a one-way ANOVA.  Categorical variables (APOE4, DX, AB-CSF, and
# AB-PET) are tested with a separate chi-squared test for every cluster pair.
# Pairwise categorical p-values are Bonferroni-adjusted within each variable.
# Tukey p-values are already adjusted by the Tukey HSD procedure.
# Longitudinal outcomes can be tested with `pairwise_longitudinal_tests()`;
# this fits a mixed model with a cluster-by-time interaction and a random
# intercept/slope for each participant.
#
# The easiest way to add another test is to add one row to
# `default_stat_test_specs()`:
#
#   tibble::add_row(specs,
#                   name = "Sex", column = "Sex", type = "categorical")
#
# The input to `run_pairwise_stat_tests()` must be a participant-level data
# frame containing `Cluster` and the columns named in `test_specs`.  The
# helper `load_multicohort_stat_data()` builds this data frame from a multiLTC
# object (or an RID/Cluster table) using both ADNI and OASIS data loaders.
#
# Example:
#   source("cluster-analysis/plots_and_tables/stat_tests.R")
#   stat_data <- load_multicohort_stat_data(multiLTC)
#   stat_table <- run_pairwise_stat_tests(stat_data)
#   print(stat_table)
#
#   longitudinal_table <- pairwise_longitudinal_tests(
#     longitudinal_data,
#     response = "TOTAL_WMH",
#     time = "Time",
#     age = "AGE"
#   )

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})

# Keep the loader environments separate because both loader files define
# functions such as get_demographics() and get_diagnoses().
loader_path <- if (file.exists("utils/adni_data_loaders.R")) {
  "utils/adni_data_loaders.R"
} else {
  file.path("..", "..", "utils", "adni_data_loaders.R")
}
oasis_loader_path <- if (file.exists("utils/oasis_data_loaders.R")) {
  "utils/oasis_data_loaders.R"
} else {
  file.path("..", "..", "utils", "oasis_data_loaders.R")
}
adni_loader <- new.env()
oasis_loader <- new.env()
source(loader_path, local = adni_loader)
source(oasis_loader_path, local = oasis_loader)

# Return the default test specification.  `name` becomes an output column,
# `column` identifies the input column, and `type` selects the test.
default_stat_test_specs <- function() {
  tibble(
    name = c("Age", "Onset", "Education", "APOE4", "DX.bl"),
    column = c("Age", "Onset", "Education", "APOE4", "DX.bl"),
    type = c("numeric", "numeric", "numeric", "categorical", "categorical")
  )
}

# Convert a multiLTC object, a data frame, or a named vector into IDs/clusters.
as_stat_cluster_table <- function(cluster_assignments) {
  if (is(cluster_assignments, "clusterObject")) {
    cluster_assignments <- tibble(RID = cluster_assignments@RID,
                                  Cluster = cluster_assignments@Cluster)
  }
  if (is.data.frame(cluster_assignments)) {
    if (!all(c("RID", "Cluster") %in% names(cluster_assignments))) {
      stop("cluster_assignments must contain columns named RID and Cluster.")
    }
    keep <- intersect(c("RID", "Cluster", "time_shift"), names(cluster_assignments))
    clusters <- cluster_assignments %>% select(all_of(keep))
  } else if (!is.null(names(cluster_assignments))) {
    clusters <- tibble(RID = names(cluster_assignments),
                       Cluster = unname(cluster_assignments))
  } else {
    stop("cluster_assignments must be a clusterObject, data frame, or named vector.")
  }

  clusters <- clusters %>% mutate(RID = as.character(RID),
    Cohort = case_when(grepl("^ADNI_", RID) ~ "ADNI",
                       grepl("^OASIS_", RID) ~ "OASIS",
                       grepl("^[0-9]+$", RID) ~ "ADNI",
                       TRUE ~ NA_character_),
    RID = if_else(Cohort == "ADNI" & grepl("^[0-9]+$", RID),
                  paste0("ADNI_", RID), RID))
  if (any(is.na(clusters$Cohort))) {
    stop("Each RID must start with ADNI_, OASIS_, or be numeric.")
  }
  clusters %>% distinct(RID, .keep_all = TRUE)
}

# Build participant-level data for the requested tests using both loaders.
load_multicohort_stat_data <- function(cluster_assignments,
                                       use_centiloid = FALSE,
                                       multi_cohort_file = "~/R/EDAP-data/MULTI_COHORT_4.csv") {
  clusters <- as_stat_cluster_table(cluster_assignments)
  time_shifts <- read.csv(path.expand(multi_cohort_file)) %>%
    select(RID, time_shift) %>% distinct(RID, .keep_all = TRUE)

  adni_demo <- adni_loader$get_demographics() %>%
    transmute(RID = paste0("ADNI_", RID), Age = AGE,
              Education = PTEDUCAT, APOE4, Cohort = "ADNI")
  oasis_demo <- oasis_loader$get_demographics() %>%
    transmute(RID = paste0("OASIS_", RID), Age = AgeatEntry,
              Education = EDUC, APOE4, Cohort = "OASIS")
  demographics <- bind_rows(adni_demo, oasis_demo) %>%
    inner_join(clusters, by = c("RID", "Cohort")) %>%
    distinct(RID, .keep_all = TRUE)

  # Collapse repeated diagnoses to one categorical outcome per participant.
  adni_dx <- adni_loader$get_diagnoses() %>%
    transmute(RID = paste0("ADNI_", RID), DX_bl = as.character(DX.bl),
              DX = as.character(DIAGNOSIS), Cohort = "ADNI")
  oasis_dx <- oasis_loader$get_diagnoses() %>%
    transmute(RID = paste0("OASIS_", RID), DX_bl = as.character(DX.bl),
              DX = as.character(DX), Cohort = "OASIS")
  diagnosis <- bind_rows(adni_dx, oasis_dx) %>%
    inner_join(clusters, by = c("RID", "Cohort")) %>%
    group_by(RID, Cohort, Cluster) %>%
    summarise(DX = case_when(
      any(DX == "Dementia", na.rm = TRUE) ~ "Dementia",
      any(DX == "MCI", na.rm = TRUE) ~ "MCI",
      any(DX == "Impaired", na.rm = TRUE) ~ "Impaired",
      any(DX == "CN", na.rm = TRUE) ~ "CN",
      TRUE ~ first(na.omit(DX_bl), default = NA_character_)),
      .groups = "drop")

  onset <- clusters %>% left_join(time_shifts, by = "RID") %>%
    mutate(time_shift = coalesce(time_shift, 0)) %>%
    select(RID, Cohort, Cluster, time_shift)
  
  # Co-Pathologies
  adni_copath <- adni_loader$get_copath() %>% 
    mutate(RID = paste0("ADNI_", RID),
           #log_WMH = log1p(TOTAL_WMH),
           Cohort = "ADNI")
  oasis_copath <- oasis_loader$get_copath() %>%
    mutate(
      #TOTAL_WMH = WM.hypointensities_volume / 1e3,
      #log_WMH = log1p(TOTAL_WMH),
           RID = paste0("OASIS_", RID),
           Cohort = "OASIS")
  
  #wmh_df <- adni_copath %>% select(RID, Cohort, AGE, Years, TOTAL_WMH, log_WMH)
  #wmh_df <- oasis_copath %>% select(RID, Cohort, AGE, Years, TOTAL_WMH, log_WMH) %>% rbind(wmh_df)
  
  lbd_df <- adni_copath %>% distinct(RID, Cohort, SAA) %>% drop_na(SAA) %>%
    rename(lbd_evidence = SAA) 
  lbd_df <- oasis_copath %>% distinct(RID, Cohort, lbdis) %>% drop_na(lbdis) %>%
    rename(lbd_evidence = lbdis) %>% rbind(lbd_df)

  demographics %>%
    left_join(onset, by = c("RID", "Cohort", "Cluster")) %>%
    mutate(Onset = Age - time_shift) %>%
    left_join(diagnosis %>% select(RID, DX, DX.bl), by = "RID") %>%
    left_join(lbd_df %>% select(RID, lbd_evidence), by = "RID") %>%
    select(RID, Cohort, Cluster, Age, Onset, Education, APOE4, DX, DX.bl, lbd_evidence)
}



# Return a stable, readable label for a pair of clusters.
cluster_pair_label <- function(cluster_1, cluster_2) {
  paste(sort(c(as.character(cluster_1), as.character(cluster_2))), collapse = " - ")
}

# Numeric variables: one-way ANOVA followed by Tukey's HSD.
pairwise_tukey <- function(data, column, cluster_col = "Cluster") {
  test_data <- data %>%
    transmute(value = suppressWarnings(as.numeric(.data[[column]])),
              Cluster = as.factor(.data[[cluster_col]])) %>%
    filter(!is.na(value), !is.na(Cluster))
  if (n_distinct(test_data$Cluster) < 2) return(tibble())
  tukey <- as.data.frame(TukeyHSD(aov(value ~ Cluster, data = test_data), "Cluster")$Cluster) %>%
    rownames_to_column("comparison")
  tibble(`Cluster pair` = vapply(strsplit(tukey$comparison, "-", fixed = TRUE),
                                 function(x) cluster_pair_label(x[1], x[2]), character(1)),
         p_value = tukey$`p adj`)
}

# Categorical variables: separate chi-squared test for each pair of clusters.
pairwise_chisq <- function(data, column, cluster_col = "Cluster",
                           p_adjust = "bonferroni") {
  test_data <- data %>%
    transmute(value = as.factor(.data[[column]]),
              Cluster = as.factor(.data[[cluster_col]])) %>%
    filter(!is.na(value), !is.na(Cluster))
  clusters <- sort(unique(test_data$Cluster))
  if (length(clusters) < 2) return(tibble())
  result <- lapply(combn(clusters, 2, simplify = FALSE), function(pair) {
    # Drop unused factor levels after subsetting; otherwise table() can keep
    # a third, empty cluster and chisq.test() returns NaN.
    subset <- test_data %>%
      filter(Cluster %in% pair) %>%
      mutate(Cluster = droplevels(Cluster), value = droplevels(value))
    contingency_table <- table(subset$Cluster, subset$value)
    if (nrow(contingency_table) < 2 || ncol(contingency_table) < 2) {
      return(tibble(`Cluster pair` = cluster_pair_label(pair[1], pair[2]),
                    p_value = NA_real_))
    }
    p_value <- tryCatch(
      suppressWarnings(chisq.test(contingency_table, correct = FALSE)$p.value),
      error = function(e) NA_real_)
    tibble(`Cluster pair` = cluster_pair_label(pair[1], pair[2]),
           p_value = p_value)
  }) %>% bind_rows()
  result$p_value <- p.adjust(result$p_value, method = p_adjust)
  result
}

# Pairwise longitudinal test based on the WMH model used in
# plots_and_tables.R.  For each reference cluster, the model estimates:
#   * p (cluster): difference between clusters at Time = 0;
#   * p (interaction): difference in longitudinal slopes.
#
# The function returns one canonical row per cluster pair.  If an age column
# is supplied, it is included as a covariate.  The data must contain one row
# per participant and visit, with a participant identifier in `id`.
pairwise_longitudinal_tests <- function(data,
                                        response,
                                        time,
                                        age = NULL,
                                        id = "RID",
                                        cluster = "Cluster",
                                        bonferroni = TRUE,
                                        return_models = FALSE) {
  if (!requireNamespace("lmerTest", quietly = TRUE)) {
    stop("pairwise_longitudinal_tests() requires the lmerTest package.")
  }
  required <- c(response, time, id, cluster, age)
  missing_columns <- setdiff(required[!is.na(required)], names(data))
  if (length(missing_columns)) {
    stop("Missing columns in longitudinal data: ",
         paste(missing_columns, collapse = ", "))
  }

  model_data <- data %>%
    transmute(Value = suppressWarnings(as.numeric(.data[[response]])),
              Time = suppressWarnings(as.numeric(.data[[time]])),
              RID = .data[[id]],
              Cluster = factor(.data[[cluster]]),
              AGE = if (is.null(age)) 0 else suppressWarnings(as.numeric(.data[[age]]))) %>%
    filter(!is.na(Value), !is.na(Time), !is.na(RID), !is.na(Cluster),
           !is.na(AGE))

  cluster_levels <- levels(droplevels(model_data$Cluster))
  if (length(cluster_levels) < 2) stop("At least two clusters are required.")

  model_results <- list()
  fitted_models <- list()

  for (reference_cluster in cluster_levels) {
    model_data$Cluster <- relevel(model_data$Cluster, ref = reference_cluster)
    model_formula <- if (is.null(age)) {
      Value ~ Time * Cluster + (Time | RID)
    } else {
      Value ~ Time * Cluster + AGE + (Time | RID)
    }
    fit <- tryCatch(
      lmerTest::lmer(model_formula, data = model_data, REML = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    if (return_models) fitted_models[[reference_cluster]] <- fit

    coefficient_table <- as.data.frame(summary(fit)$coefficients) %>%
      rownames_to_column("term")
    p_column <- "Pr(>|t|)"

    cluster_rows <- coefficient_table %>% filter(grepl("^Cluster", term))
    interaction_rows <- coefficient_table %>%
      filter(grepl("^(Time:Cluster|Cluster.*:Time)", term))
    if (!nrow(cluster_rows) && !nrow(interaction_rows)) next

    targets <- union(
      sub("^Cluster", "", cluster_rows$term),
      sub("^Cluster", "", sub("^Time:", "", sub(":Time$", "", interaction_rows$term)))
    )
    targets <- targets[nzchar(targets)]

    for (target_cluster in targets) {
      cluster_row <- cluster_rows %>%
        filter(sub("^Cluster", "", term) == target_cluster)
      interaction_row <- interaction_rows %>%
        filter(sub("^Cluster", "", sub("^Time:", "", sub(":Time$", "", term))) == target_cluster)
      model_results[[length(model_results) + 1L]] <- tibble(
        `Cluster pair` = cluster_pair_label(reference_cluster, target_cluster),
        `p (cluster)` = if (nrow(cluster_row)) cluster_row[[p_column]][1] else NA_real_,
        `p (interaction)` = if (nrow(interaction_row)) interaction_row[[p_column]][1] else NA_real_
      )
    }
  }

  result <- bind_rows(model_results) %>% distinct(`Cluster pair`, .keep_all = TRUE)
  if (bonferroni && nrow(result)) {
    n_comparisons <- length(cluster_levels) - 1
    result <- result %>% mutate(
      `p (cluster)` = pmin(`p (cluster)` * n_comparisons, 1),
      `p (interaction)` = pmin(`p (interaction)` * n_comparisons, 1))
  }
  if (return_models) return(list(table = result, models = fitted_models))
  result
}

# Run all specifications and reshape the results into the requested wide
# table.  Add future tests by adding a row to `test_specs`.
run_pairwise_stat_tests <- function(data,
                                    test_specs = default_stat_test_specs(),
                                    output_file = NULL,
                                    cluster_col = "Cluster") {
  missing_columns <- setdiff(c(cluster_col, test_specs$column), names(data))
  if (length(missing_columns)) {
    stop("Missing columns in data: ", paste(missing_columns, collapse = ", "))
  }
  all_clusters <- sort(unique(na.omit(data[[cluster_col]])))
  if (length(all_clusters) < 2) stop("At least two clusters are required.")
  pairs <- combn(all_clusters, 2, simplify = FALSE)
  output <- tibble(`Cluster pair` = vapply(pairs,
                                           function(x) cluster_pair_label(x[1], x[2]),
                                           character(1)))

  for (i in seq_len(nrow(test_specs))) {
    spec <- test_specs[i, ]
    test_result <- switch(spec$type,
      numeric = pairwise_tukey(data, spec$column, cluster_col),
      categorical = pairwise_chisq(data, spec$column, cluster_col),
      stop("Unknown test type: ", spec$type,
           ". Use 'numeric' or 'categorical'."))
    test_result <- test_result %>% rename(!!spec$name := p_value)
    output <- output %>% left_join(test_result, by = "Cluster pair")
  }
  if (!is.null(output_file)) write.csv(output, output_file, row.names = FALSE)
  output
}

main <- function() {
  run <- "exp_km_ab_ao_2"
  load(paste("~/R/EDAP-data/LTC_MC/new/", run, ".Rdata", sep = ""))
  
  cluster_assignments = multiLTC
  
  stat_data <- load_multicohort_stat_data(multiLTC)
  stat_table <- run_pairwise_stat_tests(stat_data)
  print(stat_table)
}

