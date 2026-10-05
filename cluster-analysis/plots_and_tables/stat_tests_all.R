library(ADNIMERGE2)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(tidyverse)

library(aricode)
library(infotheo)
library(paletteer)

library(lme4)
library(stats)
library(reshape2)
library(tibble)
library(RColorBrewer)

library(lmerTest)
library(emmeans)

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

source("~/R/LTC/utils/analysis_utils.R")

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
    mutate(RID = paste0("ADNI_", RID),
           Cohort = "ADNI") %>%
    distinct(RID, Cohort, DX.bl, DX.highest)
  oasis_dx <- oasis_loader$get_diagnoses() %>%
    mutate(RID = paste0("OASIS_", RID), 
           Cohort = "OASIS") %>% 
    distinct(RID, Cohort, DX.bl, DX.highest)
  
  diagnosis <- adni_dx %>% mutate(
    DX.bl = factor(DX.bl, levels = levels(oasis_dx$DX.bl), ordered = TRUE),
    DX.highest = factor(DX.highest, levels = levels(oasis_dx$DX.highest), ordered = TRUE),
  ) %>% bind_rows(oasis_dx) %>%
    inner_join(clusters, by = c("RID", "Cohort")) 
  
  onset <- clusters %>% left_join(time_shifts, by = "RID") %>%
    mutate(time_shift = coalesce(time_shift, 0)) %>%
    select(RID, Cohort, Cluster, time_shift)
  
  # Co-Pathologies
  adni_copath <- adni_loader$get_copath() %>% 
    mutate(RID = paste0("ADNI_", RID),
           Cohort = "ADNI")
  oasis_copath <- oasis_loader$get_copath() %>%
    mutate(
      RID = paste0("OASIS_", RID),
      Cohort = "OASIS")
  
  lbd_df <- adni_copath %>% distinct(RID, Cohort, SAA) %>% drop_na(SAA) %>%
    rename(lbd_evidence = SAA) 
  lbd_df <- oasis_copath %>% distinct(RID, Cohort, lbdis) %>% drop_na(lbdis) %>%
    rename(lbd_evidence = lbdis) %>% rbind(lbd_df)
  tdp_df <- read_csv(multi_cohort_file) %>% group_by(RID) %>% filter(Months == max(Months)) %>% ungroup() %>%
    mutate(TDP=(LH_MIDDLETEMPORAL+RH_MIDDLETEMPORAL+LH_INFERIORTEMPORAL+RH_INFERIORTEMPORAL)/(LH_HIPPOCAMPUS+RH_HIPPOCAMPUS)) %>% 
    drop_na(TDP) %>% distinct(RID, TDP)
  
  demographics %>%
    left_join(onset, by = c("RID", "Cohort", "Cluster")) %>%
    mutate(Onset = Age - time_shift) %>%
    left_join(diagnosis, by = c("RID", "Cohort", "Cluster")) %>%
    left_join(lbd_df %>% select(RID, lbd_evidence), by = "RID") %>%
    left_join(tdp_df, by = "RID") %>%
    select(RID, Cohort, Cluster, Age, Onset, Education, APOE4, DX.highest, DX.bl, lbd_evidence, TDP)
}

load_wmhs <- function(cluster_assignments,
                      multi_cohort_file = "~/R/EDAP-data/MULTI_COHORT_4.csv") {
  
  clusters <- as_stat_cluster_table(cluster_assignments)
  time_shifts <- read.csv(path.expand(multi_cohort_file)) %>%
    select(RID, time_shift) %>% distinct(RID, .keep_all = TRUE)
  
  adni_copath <- adni_loader$get_copath() %>% 
    mutate(RID = paste0("ADNI_", RID),
           log_WMH = log1p(TOTAL_WMH),
           Cohort = "ADNI")
  oasis_copath <- oasis_loader$get_copath() %>%
    mutate(
      TOTAL_WMH = WM.hypointensities_volume / 1e3,
      log_WMH = log1p(TOTAL_WMH),
      RID = paste0("OASIS_", RID),
      Cohort = "OASIS")
  
  wmh_df <- adni_copath %>% select(RID, Cohort, AGE, Years, TOTAL_WMH, log_WMH)
  wmh_df <- oasis_copath %>% select(RID, Cohort, AGE, Years, TOTAL_WMH, log_WMH) %>% rbind(wmh_df)
  
  wmh_df <- left_join(wmh_df, clusters, by = c("RID", "Cohort")) %>% drop_na(Cluster) %>%
    left_join(time_shifts, by = "RID") %>% mutate(Time = Years + time_shift)
  
  return(wmh_df)
}

load_amyloid_data <- function(cluster_assignments,
                              use_centiloid = FALSE,
                              multi_cohort_file = "~/R/EDAP-data/MULTI_COHORT_4.csv") {
  clusters <- as_stat_cluster_table(cluster_assignments)
  time_shifts <- read.csv(path.expand(multi_cohort_file)) %>%
    select(RID, time_shift) %>% distinct(RID, .keep_all = TRUE)
  
  adni_amyloid <- adni_loader$get_ab_df() %>%
    mutate(RID = paste0("ADNI_", RID),
              tracer = as.character(TRACER),
              Cohort = "ADNI") %>%
    rename(
      PET_SUVR = SUMMARY_SUVR,
      AB_positive = AB_pos.pet,
      CSF_RE = A4240.re,
      CSF_MS = A4240.ms,
    )
  oasis_amyloid <- oasis_loader$get_oasis_ab(use_centiloid = use_centiloid) %>%
    mutate(RID = paste0("OASIS_", gsub("^OAS", "", OASISID)),
              tracer = as.character(tracer),
              PET_SUVR = if (use_centiloid) Centiloid_SUVR else PET_SUVR,
              Cohort = "OASIS") %>%
    rename(
      AB_positive = AB_pos,
    ) %>%
    mutate(tracer = ifelse(tracer == "AV45", "FBP", tracer),
           Years = days_to_visit/365.25)
  amyloid <- bind_rows(adni_amyloid, oasis_amyloid) 
  
  # Specific Cohort cut-offs
  cutoffs <- data.frame(
    Cohort = c("ADNI", "OASIS", "ADNI", "OASIS"),
    tracer = c("FBP", "FBP", "PIB", "PIB"),
    cutoff = if (use_centiloid) c(1.11, 20.6, 1.21, 16.4) else c(1.11, 1.19, 1.21, 1.42)
  )
  # Could also add CSF cutoffs?
  
  amyloid <- amyloid %>% left_join(distinct(clusters, RID, Cluster), by = "RID") %>% drop_na(Cluster) %>%
    left_join(time_shifts, by = "RID") %>% left_join(cutoffs, by = c("Cohort", "tracer")) %>%
    mutate(Time = Years + time_shift)
  
  adni_demo <- adni_loader$get_demographics() %>%
    mutate(RID = paste0("ADNI_", RID))
  oasis_demo <- oasis_loader$get_demographics() %>%
    transmute(RID = paste0("OASIS_", RID), 
              AGE = AgeatEntry)
  ages <- bind_rows(adni_demo, oasis_demo) %>%
    inner_join(clusters, by = "RID") %>%
    distinct(RID, AGE)
  
  amyloid <- left_join(amyloid, ages, by = "RID")
  
  
  return(amyloid)
}


run_tukey_and_chi_tests <- function() {
  run <- "exp_km_ab_ao_2"
  load(paste("~/R/EDAP-data/LTC_MC/new/", run, ".Rdata", sep = ""))
  
  stat_data <- load_multicohort_stat_data(multiLTC)
  
  age_res <- anova_tukey_test(stat_data, "Age")
  onset_res <- anova_tukey_test(stat_data, "Onset")
  educ_res <- anova_tukey_test(stat_data, "Education")
  tdp_res <- anova_tukey_test(stat_data, "TDP")
  
  apoe4_res <- cat_chi_test(stat_data, "APOE4")
  dx_res <- cat_chi_test(stat_data, "DX.bl")
  dxt_res <- cat_chi_test(stat_data, "DX.highest")
  
  lbd_res <- cat_chi_test(stat_data, "lbd_evidence")
  
  stats_table <- age_res$tukey %>% select(Cluster1, Cluster2, lwr, upr, p.adj) %>% mutate(Variable = "Age")
  stats_table <- onset_res$tukey %>% select(Cluster1, Cluster2, lwr, upr, p.adj) %>% mutate(Variable = "Onset") %>% bind_rows(stats_table)
  stats_table <- educ_res$tukey %>% select(Cluster1, Cluster2, lwr, upr, p.adj) %>% mutate(Variable = "Educat.") %>% bind_rows(stats_table)
  
  stats_table <- apoe4_res %>% select(Cluster1, Cluster2, p.adj) %>% mutate(Variable = "APOE4") %>% bind_rows(stats_table)
  stats_table <- dx_res %>% select(Cluster1, Cluster2, p.adj) %>% mutate(Variable = "DX.bl") %>% bind_rows(stats_table)
  stats_table <- dxt_res %>% select(Cluster1, Cluster2, p.adj) %>% mutate(Variable = "DX.tr") %>% bind_rows(stats_table)
  stats_table <- lbd_res %>% select(Cluster1, Cluster2, p.adj) %>% mutate(Variable = "LBD") %>% bind_rows(stats_table)
  stats_table <- tdp_res$tukey %>% select(Cluster1, Cluster2, lwr, upr, p.adj) %>% mutate(Variable = "TDP") %>% bind_rows(stats_table)
  
  return(stats_table)
}

run_lmer_tests <- function() {
  run <- "exp_km_ab_ao_2"
  load(paste("~/R/EDAP-data/LTC_MC/new/", run, ".Rdata", sep = ""))
  
  wmh_data <- load_wmhs(multiLTC)
  
  ab_data <- load_amyloid_data(multiLTC, use_centiloid = FALSE)
  
  #wmh_res <- pairwise_emm_from_lmm(wmh_data, "log_WMH", "Time", correct_for_age = TRUE, by_time = FALSE)
  wmh_res <- pairwise_emm_from_lmm(wmh_data, "TOTAL_WMH", "Time", correct_for_age = TRUE, by_time = FALSE) %>%
    mutate(Variable = "WMH")
  ab_res <- ab_data %>% filter(tracer == "FBP") %>% 
    pairwise_emm_from_lmm("PET_SUVR", "Time", correct_for_age = TRUE, by_time = FALSE) %>%
    mutate(Variable = "AB_PET (FBP)")
  ab_res <- ab_data %>% filter(tracer == "PIB") %>%
    pairwise_emm_from_lmm("PET_SUVR", "Time", correct_for_age = TRUE, by_time = FALSE) %>%
    mutate(Variable = "AB_PET (PIB)") %>% bind_rows(ab_res)
  
  ab_res <- ab_data %>% filter_out(is.na(CSF_RE)) %>%
    pairwise_emm_from_lmm("CSF_RE", "Time", correct_for_age = TRUE, by_time = FALSE) %>%
    mutate(Variable = "AB42/AB40 (IA)") %>% bind_rows(ab_res)
  ab_res <- ab_data %>% filter_out(is.na(CSF_MS)) %>%
    pairwise_emm_from_lmm("CSF_MS", "Time", correct_for_age = TRUE, by_time = FALSE) %>%
    mutate(Variable = "AB42/AB40 (MS)") %>% bind_rows(ab_res)
  
  bind_rows(wmh_res, ab_res)
}

# Run the pairwise tests and return LaTeX tables in the same compact layout as
# the manuscript tables. Tukey tests include their confidence interval and
# adjusted p-value; chi-squared and mixed-model tests show adjusted p-values.
# `run_tukey_and_chi_tests()` and `run_lmer_tests()` use the configured
# multi-cohort analysis in this script.
stat_tests_latex_tables <- function() {
  if (!requireNamespace("knitr", quietly = TRUE)) {
    stop("stat_tests_latex_tables() requires the knitr package.")
  }

  format_p <- function(x) {
    ifelse(is.na(x), "", sprintf("$%.*f$", 3, x))
  }
  format_ci <- function(lower, upper) {
    ifelse(is.na(lower) | is.na(upper), "",
           sprintf("($%.*f$, $%.*f$)", 2, lower, 2, upper))
  }
  pair_key <- function(x) paste(x$Cluster1, x$Cluster2, sep = "-")

  tests <- run_tukey_and_chi_tests()
  tukey_vars <- c("Age", "Onset", "TDP")
  p_only_vars <- setdiff(unique(tests$Variable), tukey_vars)
  pairs <- tests %>% distinct(Cluster1, Cluster2) %>% arrange(Cluster1, Cluster2)
  result <- tibble(Clusters = pair_key(pairs))

  add_result <- function(data, variable, label = variable, include_ci = FALSE) {
    selected <- data %>% filter(Variable == variable) %>%
      mutate(Clusters = pair_key(.))
    if (include_ci) {
      selected <- selected %>% mutate(CI = format_ci(lwr, upr), P = format_p(p.adj)) %>%
        select(Clusters, CI, P)
      names(selected)[2:3] <- paste0(label, c("_CI", "_p"))
    } else {
      selected <- selected %>% mutate(P = format_p(p.adj)) %>% select(Clusters, P)
      names(selected)[2] <- label
    }
    result <<- result %>% left_join(selected, by = "Clusters")
  }

  for (variable in tukey_vars) {
    if (variable %in% tests$Variable) add_result(tests, variable, variable, TRUE)
  }
  add_result(tests, "Educat.", "Educat.", FALSE)
  
  labels <- c(APOE4 = "APOE4", DX.bl = "DX (baseline)", DX.tr = "DX (last)", LBD = "LBD+")
  for (variable in p_only_vars) {
    label <- if (variable %in% names(labels)) unname(labels[[variable]]) else variable
    add_result(tests, variable, label, FALSE)
  }

  # A second table contains pairwise longitudinal outcomes. Their contrasts
  # are estimates at the model's reference time, so report p-values only.
  longitudinal <- run_lmer_tests()
  longitudinal_vars <- unique(longitudinal$Variable)
  long_table <- tibble(Clusters = character())
  if (length(longitudinal_vars)) {
    long_pairs <- longitudinal %>% distinct(Cluster1, Cluster2) %>%
      arrange(Cluster1, Cluster2)
    long_table <- tibble(Clusters = pair_key(long_pairs))
    for (variable in longitudinal_vars) {
      selected <- longitudinal %>% filter(Variable == variable) %>%
        mutate(Clusters = pair_key(.)) %>%
        transmute(Clusters, !!variable := format_p(p.value))
      long_table <- left_join(long_table, selected, by = "Clusters")
    }
  }

  names(result)[1] <- "Clusters"
  
  long_table <- result %>% select(Clusters, TDP_CI, TDP_p, `LBD+`) %>% right_join(long_table, by = "Clusters")
  long_table <- select(long_table, Clusters, `AB42/AB40 (IA)`, `AB42/AB40 (MS)`, `AB_PET (PIB)`, `AB_PET (FBP)`, everything())
  
  demo_table <- result %>% select(-c("TDP_CI", "TDP_p", "LBD+"))
  
  list(
    demographics = knitr::kable(demo_table, format = "latex", escape = FALSE,
                                booktabs = TRUE, row.names = FALSE),
    longitudinal = knitr::kable(long_table, format = "latex", escape = FALSE,
                                 booktabs = TRUE, row.names = FALSE)
  )
}




#### Plots ####
make_plots <- function() {
  run <- "exp_km_ab_ao_2"
  load(paste("~/R/EDAP-data/LTC_MC/new/", run, ".Rdata", sep = ""))
  
  wmh_data <- load_wmhs(multiLTC)
  stat_data <- load_multicohort_stat_data(multiLTC)
  
  wmh_data$Cluster <- factor(wmh_data$Cluster, levels = sort(levels(wmh_data$Cluster)))
  gp_wmh <- ggplot(wmh_data, aes(x=Time, y=TOTAL_WMH, colour = Cluster, fill=Cluster, group=RID)) +
    geom_point(shape=21, color="black") +
    geom_line() +
    facet_wrap(~Cluster) +
    labs(x = "Years (shifted)", y = "White Matter Hyperintensity", color = "Cluster", fill = "Cluster") +
    scale_color_paletteer_d("ggthemes::Tableau_10") +
    scale_fill_paletteer_d("ggthemes::Tableau_10") +
    theme_classic() 
  plot(gp_wmh)
  
  
  
  tp <- ggplot(stat_data, aes(x=Cluster, y = TDP, fill = Cluster)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.2) +
    #geom_violin(alpha=0.5) +
    geom_point(size = 0.8, position = position_jitter(width = 0.2, height = 0), aes(color=Cluster)) +
    labs(x = "Cluster", y = "TDP-43 MRI ratio") +
    scale_color_paletteer_d("ggthemes::Tableau_10") +
    scale_fill_paletteer_d("ggthemes::Tableau_10") +
    theme_classic()
  tp
  
  lbd_summary <- stat_data %>% group_by(Cluster) %>% 
    summarise(Ratio.evidence = mean(lbd_evidence, na.rm=TRUE)) %>% ungroup()
  lbdp <- ggplot(lbd_summary, aes(x=Cluster, y=Ratio.evidence, fill=Cluster)) +
    geom_col() +
    labs(x = "Cluster", y = "LBD+ (%)") + 
    scale_color_paletteer_d("ggthemes::Tableau_10") +
    scale_fill_paletteer_d("ggthemes::Tableau_10") +
    theme_classic()
  lbdp
  
  # Save plots
  ggsave(gp_wmh, filename = paste("~/R/EDAP-data/plots/LTC_MC/", run, "_wmh.png", sep=""), width = 7, height = 5, dpi = 300)
  ggsave(tp, filename = paste("~/R/EDAP-data/plots/LTC_MC/", run, "_tdp.png", sep=""), width = 7, height = 5, dpi = 300)
  ggsave(lbdp, filename = paste("~/R/EDAP-data/plots/LTC_MC/", run, "_saa.png", sep=""), width = 6.5, height = 5, dpi = 300)
}

latex_res <- stat_tests_latex_tables()  

print(latex_res$demographics)
print(latex_res$longitudinal)

make_plots()

