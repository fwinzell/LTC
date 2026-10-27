# Generic BioFINDER cluster plots and tables
# =========================================
#
# Source this file, then pass a clusterObject or a data frame with `sid` and
# `Cluster` columns to run_biofinder_plots_and_tables(). The function returns
# its input data, descriptive/statistical tables, LaTeX table data, and plots.
# It only analyzes MRI-derived values for the TDP-43 signature; MRI scan count
# and follow-up are retained as cohort descriptors.
#
# Example:
#   assignments <- load_biofinder_cluster_assignments()
#   bf_results <- run_biofinder_plots_and_tables(assignments)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(ggalluvial)
  library(paletteer)
  library(knitr)
  library(purrr)
})

bf_loader_path <- if (file.exists("utils/biofinder_data_loaders.R")) {
  "utils/biofinder_data_loaders.R"
} else if (file.exists("../../utils/biofinder_data_loaders.R")) {
  "../../utils/biofinder_data_loaders.R"
} else {
  "~/R/LTC/utils/biofinder_data_loaders.R"
}
bf_loader <- new.env(parent = globalenv())
source(path.expand(bf_loader_path), local = bf_loader)

source("~/R/LTC/utils/analysis_utils.R")

as_biofinder_cluster_table <- function(cluster_assignments) {
  if (methods::is(cluster_assignments, "clusterObject")) {
    cluster_assignments <- data.frame(sid = cluster_assignments@RID,
                                      Cluster = cluster_assignments@Cluster)
  }
  if (is.data.frame(cluster_assignments)) {
    if (!all(c("sid", "Cluster") %in% names(cluster_assignments))) {
      stop("A BioFINDER cluster data frame must contain `sid` and `Cluster` columns.")
    }
    clusters <- cluster_assignments %>% select(sid, Cluster)
  } else if (!is.null(names(cluster_assignments))) {
    clusters <- tibble(sid = names(cluster_assignments),
                       Cluster = unname(cluster_assignments))
  } else {
    stop("cluster_assignments must be a clusterObject, data frame, or named vector.")
  }
  clusters %>% mutate(sid = as.character(sid), Cluster = as.character(Cluster)) %>%
    filter(!is.na(sid), !is.na(Cluster)) %>% distinct(sid, .keep_all = TRUE)
}

load_biofinder_cluster_assignments <- function(
    run = "exp_km_ab_bf", model_dir = "~/R/EDAP-data/LTC_MC/new") {
  model_env <- new.env(parent = globalenv())
  load(file.path(path.expand(model_dir), paste0(run, ".Rdata")), envir = model_env)
  if (!exists("bFLTC", envir = model_env, inherits = FALSE)) {
    stop("The model file does not contain an object named `bFLTC`.")
  }
  model_env$bFLTC
}

# Summarise a categorical variable as a percentage and count in each cluster.
categorical_summary <- function(data, variable, group_vars = c("Cohort", "Cluster")) {
  group_vars <- intersect(group_vars, names(data))
  data %>%
    filter(!is.na(Cluster), !is.na(.data[[variable]])) %>%
    mutate(value = .data[[variable]]) %>%
    group_by(across(all_of(c(group_vars, "value")))) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(across(all_of(group_vars))) %>%
    mutate(percent = 100 * n / sum(n)) %>%
    ungroup()
}


run_biofinder_plots_and_tables <- function(
    output_dir = "~/R/EDAP-data/plots/LTC_MC/BioFINDER",
    run = "exp_km_ab_bf",
    dpm_path = "~/R/EDAP-data/BioFINDER/DPM_BioFINDER.csv",
    save_plots = TRUE,
    print_tables = TRUE) {
  
  clusters <- as_biofinder_cluster_table(load_biofinder_cluster_assignments(run = run))
  cluster_levels <- sort(unique(clusters$Cluster))

  # Load longitudinal MRI.
  mri_df <- bf_loader$get_mri_data_updated(normalize = TRUE) %>%
    mutate(sid = as.character(sid)) %>% inner_join(clusters, by = "sid")
  mri_count <- mri_df %>% group_by(sid, Cluster) %>% summarise(
    n_mri = n(), follow_up = max(Years, na.rm = TRUE) - min(Years, na.rm = TRUE),
    .groups = "drop")
  
  # Diagnoses
  bf_dx <- bf_loader$get_diagnoses() %>%
    group_by(sid) %>%
    mutate(diag.highest = ifelse(
      diag.bl == "Dementia" | 
        converted_dementia_dich_baseline_variable %in% 1,
      "Dementia",
      ifelse(diag.bl == "MCI" | 
               converted_MCI_dich_baseline_variable %in% 1,
             "MCI", 
             as.character(diag.bl))
    )) %>% ungroup() %>%
    mutate(
      diag.bl = factor(diag.bl, levels = c("Normal", "SCD", "MCI", "Dementia"), ordered = TRUE),
      diag.highest =  factor(diag.highest, levels = c("Normal", "SCD", "MCI", "Dementia"), ordered = TRUE))
  
  dx_df <- bf_dx %>% mutate(sid = as.character(sid)) %>%
    left_join(clusters, by = "sid") %>% filter(!is.na(Cluster)) %>%
    group_by(sid) %>% arrange(Visit, .by_group = TRUE) %>%
    mutate(CN2CI = ifelse(diag.bl %in% c("Normal", "SCD"), diag.highest %in% c("MCI", "Dementia"), NA),
           MCI2AD = ifelse(diag.bl == "MCI", diag.highest == "Dementia", NA)) %>%
    ungroup()

  # Demographics
  demographics <- bf_loader$get_demographics() %>%
    transmute(sid = as.character(sid), AGE = age_bl,
              EDUC = education_level_years_baseline_variable,
              GENDER = ifelse(gender_baseline_variable == 1, "Female", "Male"), 
              APOE4) %>%
    distinct(sid, .keep_all = TRUE) %>% inner_join(clusters, by = "sid")
  demographics <- demographics %>% mutate(Cluster = factor(Cluster, levels = cluster_levels))
  
  dpm <- read.csv(path.expand(dpm_path)) %>%
    transmute(sid = as.character(sid), time_shift) %>% distinct(sid, .keep_all = TRUE)
  onset <- demographics %>% left_join(dpm, by = "sid") %>%
    mutate(onset_age = AGE - time_shift) %>% filter(is.finite(onset_age))

  # Amyloid status and measured PET/CSF values, summarized at participant level.
  amyloid <- bf_loader$get_ab_df() %>% mutate(sid = as.character(sid)) %>%
    inner_join(clusters, by = "sid") %>%
    rename(PET_SUVR = fnc_ber_com_composite,
           CSF_RE = CSF_Ab42_Ab40_ratio_imputed_Elecsys_2020_2022)
  
  amyloid_pet_participant <- amyloid %>%
    group_by(sid, Cluster) %>%
    summarise(PET_SUVR = max(PET_SUVR, na.rm = TRUE),
              AB_PET = any(AB_PET, na.rm = TRUE),
              .groups = "drop") 
  
  amyloid_csf_participant <- amyloid %>% 
    group_by(sid, Cluster) %>%
    summarise(CSF_RE = na_if(min(CSF_RE, na.rm = TRUE), Inf),
              Ab_CSF = any(AB_CSF, na.rm=TRUE),
              .groups = "drop")
  
  # Specific Cohort cut-offs
  cutoffs <- list(
    PET = 1.03,
    CSF = 0.080
  )


  descriptive_table <- demographics %>%
    group_by(Cluster) %>%
    summarise(n = n(),
              age_mean = mean(AGE, na.rm = TRUE),
              age_sd = sd(AGE, na.rm = TRUE),
              education_mean = mean(EDUC, na.rm = TRUE),
              education_sd = sd(EDUC, na.rm = TRUE),
              male_percent = 100 * mean(GENDER == "Male", na.rm = TRUE),
              fem_percent = 100 * mean(GENDER == "Female", na.rm = TRUE),
              .groups = "drop") %>%
    left_join(onset %>% group_by(Cluster) %>%
                summarise(onset_mean = mean(onset_age, na.rm = TRUE),
                          onset_sd = sd(onset_age, na.rm = TRUE), 
                          timeshift_mean = mean(time_shift, na.rm=TRUE),
                          timeshift_sd = sd(time_shift, na.rm=TRUE),
                          .groups = "drop"),
              by = c("Cluster")) %>%
    left_join(mri_count %>% group_by(Cluster) %>%
                summarise(mri_mean = mean(n_mri), mri_sd = sd(n_mri), .groups = "drop"),
              by = c("Cluster")) 

  diagnosis_counts <- dx_df %>% distinct(sid, diag.bl, diag.highest, Cluster) %>%
    mutate(diag.bl = recode(diag.bl, "Dementia" = "AD", "Normal" = "CN"), 
           diag.highest = recode(diag.highest, "Dementia" = "AD", "Normal" = "CN")) %>%
    group_by(Cluster) %>%
    count(diag.bl, diag.highest, name = "Freq")
  
  
  apoe_counts <- categorical_summary(demographics, "APOE4", group_vars = c("Cluster"))

  # Cohort-faceted plots make it possible to see whether a cluster pattern is
  # driven by one source cohort.
  p_age <- ggplot(demographics, aes(factor(Cluster), AGE, fill = Cluster)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.3) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.6) +
    scale_fill_paletteer_d("ggthemes::Tableau_10") +
    theme_classic() +
    labs(x = "Cluster", y = "Age")
  
  p_onset <- ggplot(onset, aes(factor(Cluster), onset_age, fill = Cluster)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.3) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.6) +
    scale_fill_paletteer_d("ggthemes::Tableau_10") +
    theme_classic() +
    labs(x = "Cluster", y = "Estimated age of onset")
  
  p_timeshift <- ggplot(onset, aes(factor(Cluster), time_shift, fill = Cluster)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.3) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.6) +
    scale_color_paletteer_d("ggthemes::wsj_red_green") +
    scale_fill_paletteer_d("ggthemes::Tableau_10") +
    theme_classic() +
    labs(x = "Cluster", y = "MCDP time-shift")

  
  plot_of_ages <- ggarrange(p_age, p_onset, nrow = 1, common.legend = TRUE, legend = "right")
  
  gp_alluv_list <- lapply(cluster_levels, function(c) {
    gp_alluv <- diagnosis_counts %>% filter(Cluster == c) %>%
      ggplot(aes(axis1 = diag.bl, axis2 = diag.highest, y = Freq)) +
      geom_alluvium(aes(fill = diag.highest), width = 1/12) +
      geom_stratum(width = 1/8, fill = "grey90", color = "grey40") +
      geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
      scale_x_discrete(limits = c("Baseline", "Final diagnosis"), expand = c(.1, .05)) +
      labs(x = NULL, y = "Count", title = paste0("Cluster ", c)) +
      scale_fill_brewer(palette = "YlOrRd", name = "Diagnosis") +
      theme_classic()
  })
  
  gp_alluv <- ggarrange(plotlist=gp_alluv_list, ncol=2, nrow=1, common.legend = TRUE, legend="right")
  plot(gp_alluv)
  
  p_ab_pet <- ggplot(amyloid_pet_participant,
                      aes(factor(Cluster), PET_SUVR, fill = Cluster)) +
    geom_violin(alpha = 0.35, position = position_dodge(width = 0.8)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
               alpha = 0.6) +
    geom_hline(
      aes(yintercept = cutoffs$PET),
      linetype = "dashed"
    ) +
    scale_fill_paletteer_d("ggthemes::Tableau_10") +
    theme_classic() + labs(x = "Cluster", y = "Amyloid-PET")
  
  p_ab_csf <- ggplot(amyloid_csf_participant,
                     aes(factor(Cluster), CSF_RE, fill = Cluster)) +
    geom_violin(alpha = 0.35, position = position_dodge(width = 0.8)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
               alpha = 0.6) +
    geom_hline(
      aes(yintercept = cutoffs$CSF),
      linetype = "dashed"
    ) +
    scale_fill_paletteer_d("ggthemes::Tableau_10") +
    theme_classic() + labs(x = "Cluster", y = "AB42/AB40 CSF")
  
  p_amyloid <- ggarrange(p_ab_pet, p_ab_csf, ncol=2, common.legend = TRUE, legend = "right")
  
  plots <- list(age = p_age, onset = p_onset, p_of_ages = plot_of_ages, diagnosis = gp_alluv,
                amyloid = p_amyloid, timeshift = p_timeshift)

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(file.path(output_dir, "bf_age.png"), p_age, width = 5, height = 4, dpi = 300)
    ggsave(file.path(output_dir, "bf_onset.png"), p_onset, width = 5, height = 4, dpi = 300)
    ggsave(file.path(output_dir, "bf_timeshift.png"), p_timeshift, width = 5, height = 4, dpi = 300)
    ggsave(file.path(output_dir, "bf_diagnosis.png"), gp_alluv, width = 8, height = 6, dpi = 300)
    ggsave(file.path(output_dir, "bf_amyloid.png"), p_amyloid, width = 9, height = 4, dpi = 300)
    ggsave(file.path(output_dir, "bf_ages.png"), plot_of_ages, width = 9, height = 4, dpi = 500)
  }

  list(data = list(demographics = demographics,
                   onset = onset,
                   diagnoses = dx_df, amyloid = amyloid),
       tables = list(descriptive = descriptive_table,
                     diagnosis = dx_df %>% distinct(sid, Cluster, diag.bl, diag.highest, CN2CI, MCI2AD),
                     diagnosis_counts = diagnosis_counts,
                     apoe = apoe_counts,
                     amyloid_pet = amyloid_pet_participant,
                     amyloid_csf = amyloid_csf_participant
                     ),
       plots = plots)
}


# Stats table 
run_biofinder_stats <- function(results, save_plots = FALSE) {
  demographics <- results$data$demographics
  onset <- results$data$onset
  amyloid <- results$data$amyloid
  diagnoses <- results$data$diagnoses
  clusters <- demographics %>% select(sid, Cluster) %>% distinct()

  # Keep the test inputs participant- or visit-level as appropriate, and use
  # the same ANOVA/Tukey and pairwise mixed-model helpers as stat_tests_all.R.
  age_res <- anova_tukey_test(demographics %>% filter(!is.na(AGE)), "AGE")
  onset_res <- anova_tukey_test(onset %>% filter(is.finite(onset_age)), "onset_age")
  education_res <- anova_tukey_test(demographics %>% filter(!is.na(EDUC)), "EDUC")

  # Chi-squared tests compare cluster pairs for categorical cohort outcomes.
  apoe_res <- cat_chi_test(demographics %>% filter(!is.na(APOE4)), "APOE4")
  dx_bl_res <- cat_chi_test(diagnoses %>% distinct(sid, Cluster, diag.bl) %>%
                              filter(!is.na(diag.bl)), "diag.bl")
  dx_highest_res <- cat_chi_test(diagnoses %>% distinct(sid, Cluster, diag.highest) %>%
                                   filter(!is.na(diag.highest)), "diag.highest")
  

  # Longitudinal MRI and amyloid measurements use participant random effects;
  # attach baseline age because the shared helper adjusts for age by default.
  copaths <- bf_loader$get_copath() %>%
    mutate(sid = as.character(sid)) %>%
    inner_join(clusters, by = "sid") %>%
    left_join(demographics %>% select(sid, AGE), by = "sid") %>%
    mutate(RID = sid)
  tdp_df <- copaths %>% group_by(sid, Cluster) %>%
    slice_max(Years, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(sid, Cluster, TDP = tdp)
  
  tdp_res <- anova_tukey_test(tdp_df, "TDP")
  mri_res <- pairwise_emm_from_lmm(copaths, "tdp", "Years", correct_for_age = TRUE,
                                   by_time = FALSE) %>% mutate(Variable = "TDP")
  wmh_res <- pairwise_emm_from_lmm(copaths, "wmhs_total", "Years", correct_for_age = TRUE,
                                   by_time = FALSE) %>% mutate(Variable = "WMH")

  amyloid_long <- amyloid %>%
    left_join(demographics %>% select(sid, AGE), by = "sid") %>%
    mutate(RID = as.character(sid))
  pet_res <- amyloid_long %>% filter(is.finite(PET_SUVR)) %>%
    pairwise_emm_from_lmm("PET_SUVR", "Years", correct_for_age = TRUE,
                           by_time = FALSE) %>% mutate(Variable = "AB_PET")
  csf_res <- amyloid_long %>% anova_tukey_test("CSF_RE")
  
  ## make plots
  copaths <- left_join(copaths, distinct(onset, sid, time_shift), by="sid") %>%
    mutate(Time = Years + time_shift)
  
  gp_wmh <- ggplot(copaths, aes(x=Time, y=wmhs_total, colour = Cluster, fill=Cluster, group=sid)) +
    geom_point(shape=21, color="black") +
    geom_line() +
    facet_wrap(~Cluster) +
    labs(x = "Years (shifted)", y = "White Matter Hyperintensity", color = "Cluster", fill = "Cluster") +
    scale_color_paletteer_d("ggthemes::Tableau_10") +
    scale_fill_paletteer_d("ggthemes::Tableau_10") +
    theme_classic() 
  plot(gp_wmh)
  
  tp <- ggplot(tdp_df, aes(x=Cluster, y = TDP, fill = Cluster)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.2) +
    #geom_violin(alpha=0.5) +
    geom_point(size = 0.8, position = position_jitter(width = 0.2, height = 0), aes(color=Cluster)) +
    labs(x = "Cluster", y = "TDP-43 MRI ratio") +
    scale_color_paletteer_d("ggthemes::Tableau_10") +
    scale_fill_paletteer_d("ggthemes::Tableau_10") +
    theme_classic()
  plot(tp)
  
  if (save_plots) {
    ggsave(gp_wmh, filename = paste("~/R/EDAP-data/plots/LTC_MC/BioFINDER/bf_wmh.png", sep=""), width = 5, height = 4, dpi = 300)
    ggsave(tp, filename = paste("~/R/EDAP-data/plots/LTC_MC/BioFINDER/bf_tdp.png", sep=""), width = 5, height = 4, dpi = 300)
  }
  
  list(
    anova_tukey = list(Age = age_res, Onset = onset_res, Education = education_res,
                       TDP = tdp_res, AB_CSF = csf_res),
    cat_chi_sq = list(APOE4 = apoe_res, diag.bl = dx_bl_res,
                      diag.highest = dx_highest_res),
    pairwise_emm = bind_rows(wmh_res, pet_res)
  )
}

# Useful functions for latex table construction
mean_sd_latex <- function(mean_value, sd_value, digits=2) {
  sprintf(paste0("$%.", digits, "f \\pm %.", digits, "f$"), mean_value, sd_value)
}
percent_latex <- function(p) {
  if (!p) return("--")
  sprintf("$%.1f\\%%$", p)
}
count_latex <- function(x) sprintf("$%d$", as.integer(x))

print_latex_table <- function(results) {
  # Convert to latex formatted nice table 
  desc <- results$tables$descriptive
  apoe <- results$tables$apoe
  diag <- results$tables$diagnosis %>% group_by(Cluster) %>%
    summarise(n=n(),
              n_cn=sum(diag.bl == "Normal"),
              n_scd = sum(diag.bl == "SCD"),
              n_mci=sum(diag.bl == "MCI"),
              n_dem=sum(diag.bl == "Dementia"),
              n_cn_ci=sum(CN2CI, na.rm = TRUE),
              n_mci_ad=sum(MCI2AD, na.rm = TRUE))
  
  ab_pet <- results$tables$amyloid_pet %>%
    filter(is.finite(PET_SUVR)) %>%
    group_by(Cluster) %>%
    summarise(Mean = mean(PET_SUVR, na.rm=TRUE),
              SD = sd(PET_SUVR, na.rm=TRUE))
  ab_csf <- results$tables$amyloid_csf %>%
    group_by(Cluster) %>% 
    summarise(CSF_RE_mean = mean(CSF_RE, na.rm=TRUE),
              CSF_RE_sd = sd(CSF_RE, na.rm=TRUE),
    )
  
  clusterList <- list()
  for (cl in desc$Cluster) { 
    row <- filter(desc, Cluster == cl)
    d <- filter(diag, Cluster == cl)
    
    metric_rows <- tibble(
      Variable=c("N", "Age", "Time-shift", "Onset", "Educat.", "Num. MRI", "Female",
                 "APOE e4 (0)", "APOE e4 (1)", "APOE e4 (2)", 
                 "CN", "SCD", "MCI", "AD", "CN - MCI/AD", "MCI - AD",
                 "AB42/AB40 (IA)", "AB-PET (FMM)"),
      Values=c(
        count_latex(row$n),
        mean_sd_latex(row$age_mean, row$age_sd),
        mean_sd_latex(row$timeshift_mean, row$timeshift_sd),
        mean_sd_latex(row$onset_mean, row$onset_sd),
        mean_sd_latex(row$education_mean, row$education_sd),
        mean_sd_latex(row$mri_mean, row$mri_sd),
        percent_latex(row$fem_percent),
        
        percent_latex(filter(apoe, Cluster == cl & value == 0)$percent),
        percent_latex(filter(apoe, Cluster == cl & value == 1)$percent),
        percent_latex(filter(apoe, Cluster == cl & value == 2)$percent),
        
        percent_latex(100*d$n_cn/d$n),
        percent_latex(100*d$n_scd/d$n),
        percent_latex(100*d$n_mci/d$n),
        percent_latex(100*d$n_dem/d$n),
        percent_latex(100*d$n_cn_ci/d$n_cn),
        percent_latex(100*d$n_mci_ad/d$n_mci),
        
        do.call(mean_sd_latex, c(unname(filter(ab_csf, Cluster == cl)[c("CSF_RE_mean", "CSF_RE_sd")]), 
                                 list(digits = 3))),
        
        do.call(mean_sd_latex, c(unname(filter(ab_pet, Cluster == cl)[c("Mean", "SD")]), 
                                 list(digits = 3)))
      )
    ) %>%
      rename(!!cl := Values)
    
    clusterList[[cl]] <- metric_rows
  }
  
  latex_out <- purrr::reduce(clusterList, left_join, by="Variable")
  knitr::kable(latex_out, format="latex", escape=FALSE, booktabs=TRUE,
               align=c("l", "l", "l", "l", "l"))
}

print_latex_stats <- function(stats) {
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
  
  tukey_table <- purrr::imap_dfr(
    stats$anova_tukey,
    ~ .x$tukey %>% mutate(Variable = .y)
  )
  
  chi_table <- purrr::imap_dfr(
    stats$cat_chi_sq,
    ~ .x %>% mutate(Variable = .y)
  )
  pair_key <- function(x) paste(x$Cluster1, x$Cluster2, sep = "-")
  
  pairs <- tukey_table %>% distinct(Cluster1, Cluster2) %>% arrange(Cluster1, Cluster2)
  stats_table <- tibble(Clusters = pair_key(pairs))
  
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
    stats_table <<- stats_table %>% left_join(selected, by = "Clusters")
  }
  
  for (variable in tukey_table$Variable) {
    add_result(tukey_table, variable, variable, TRUE)
  }
  
  for (variable in chi_table$Variable) {
    add_result(chi_table, variable, variable, FALSE)
  }
  
  long_table <- stats$pairwise_emm %>% mutate(Clusters = pair_key(.)) %>%
    select(Clusters, Variable, p.value) %>% 
    pivot_wider(id_cols = Clusters, names_from = Variable, values_from = p.value)
    
  
  demo_table <- stats_table %>% select(Clusters, Age_CI, Age_p, Onset_CI, Onset_p,
                                  Education_p, diag.highest, diag.bl, APOE4)
  
  long_table <- stats_table %>% select(Clusters, AB_CSF_CI, AB_CSF_p, TDP_CI, TDP_p) %>%
    right_join(long_table, by = "Clusters") %>% select(Clusters, AB_PET, everything())
  
  list(
    demographics = knitr::kable(demo_table, format = "latex", escape = FALSE,
                                booktabs = TRUE, row.names = FALSE),
    longitudinal = knitr::kable(long_table, format = "latex", escape = FALSE,
                                booktabs = TRUE, row.names = FALSE)
  )

}

# Get tables and plot plots
results <- run_biofinder_plots_and_tables(save_plots = FALSE)

table <- print_latex_table(results)
print(table)

# Run statistical tests
stats <- run_biofinder_stats(results, save_plots = TRUE)

stats_tables <- print_latex_stats(stats)

print(stats_tables$demographics)
print(stats_tables$longitudinal)
