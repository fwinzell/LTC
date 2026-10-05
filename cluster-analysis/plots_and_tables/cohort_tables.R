# Cohort-level baseline summaries split by whether a participant ever received
# an Alzheimer's disease diagnosis. Run from any working directory with the
# project dependencies and cohort data paths configured as in the loaders.
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)

multi_cohort_df <- read.csv("~/R/EDAP-data/MULTI_COHORT_4.csv", header = TRUE)
all.vars <- c(grepv("^(RH_|LH_|CC_)", colnames(multi_cohort_df)), "BRAINSTEM")

adni_dl <- new.env()
source("~/R/LTC/utils/adni_data_loaders.R", local=adni_dl)

oasis_dl <- new.env()
source("~/R/LTC/utils/oasis_data_loaders.R", local=oasis_dl)

bf_dl <- new.env()
source("~/R/LTC/utils/biofinder_data_loaders.R", local=bf_dl)


source("~/R/LTC/utils/analysis_utils.R")

make_table <- function(dx_df, demo_df, ab_df, rids, id_col="RID", visit_col="VISCODE2") {
  table_df <- left_join(dx_df, demo_df, by=id_col) 
  table_df <- left_join(table_df, ab_df, by=c(id_col, visit_col)) 
}

# ADNI: baseline demographics and diagnosis history; MRI and biomarker values
# are summarized per person when the relevant loader data are available.
adni_dx <- adni_dl$get_diagnoses()
adni_demo <- adni_dl$get_demographics() %>% distinct(RID, .keep_all = TRUE)
adni_ab <- adni_dl$get_ab_df() %>% select(-CI)

adni_mri <- multi_cohort_df %>% filter(Cohort == "ADNI") %>% mutate(RID = as.numeric(gsub("ADNI_", "", RID))) 
adni_rids <- adni_mri %>% distinct(RID) %>% unlist() %>% unname()

adni_table <- left_join(adni_dx, adni_demo, by="RID") %>%
  left_join(adni_ab, by=c("RID", "VISCODE2")) %>% 
  filter(RID %in% adni_rids)

adni_follow_up <- adni_mri %>% group_by(RID) %>% 
  summarise(n_mri=n(),
            follow_up=max(Years, na.rm=TRUE) - min(Years, na.rm=TRUE)) %>% ungroup()

adni_table <- left_join(adni_table, adni_follow_up, by="RID") %>%
  rename(MMSE = MMSCORE,
         GENDER = PTGENDER,
         EDUC = PTEDUCAT) %>%
  mutate(MMSE = as.integer(MMSE))

# OASIS: baseline diagnosis plus all available diagnosis visits to determine
# whether AD/dementia was ever recorded.
oasis_dx <- oasis_dl$get_diagnoses() %>% select(-RID)

oasis_demo <- oasis_dl$get_demographics() %>% select(-RID)

oasis_ab <- oasis_dl$get_oasis_ab(use_centiloid = FALSE)

oasis_mri <- multi_cohort_df %>% filter(Cohort == "OASIS") %>% mutate(RID = as.numeric(gsub("OASIS_", "", RID))) 
oasis_rids <- oasis_mri %>% distinct(RID) %>% unlist() %>% unname()

oasis_table <- left_join(oasis_dx, oasis_demo, by="OASISID") %>%
  full_join(oasis_ab, by=c("OASISID", "days_to_visit")) %>% 
  mutate(RID = as.numeric(gsub("^OAS", "", OASISID))) %>%
  filter(RID %in% oasis_rids)

oasis_table <- oasis_table %>% group_by(RID) %>% 
  fill(DX.bl, .direction = "downup") %>%
  fill(GENDER, .direction = "downup") %>%
  ungroup()

oasis_follow_up <- oasis_mri %>% group_by(RID) %>% 
  summarise(n_mri=n(),
            follow_up=max(Years, na.rm=TRUE) - min(Years, na.rm=TRUE)) %>% ungroup()

oasis_table <- left_join(oasis_table, oasis_follow_up, by="RID") %>%
  rename(AGE = AgeatEntry,
         TRACER = tracer) %>%
  mutate(GENDER = ifelse(GENDER == 1, "Male", "Female"))


# BioFINDER
bf_dx <- bf_dl$get_diagnoses()

bf_demo <- bf_dl$get_demographics()

bf_ab <- bf_dl$get_ab_df() %>% select(-c("cognitive_status_baseline_variable", 
                                         "mmse_score", "cdr_global_clinical"))

bf_mri <- bf_dl$get_mri_data_updated(normalize = TRUE) 
ab_pos_ids <- bf_ab %>% group_by(sid) %>% mutate(AB_any = any(AB)) %>% ungroup() %>%
  filter(AB_any) %>% select(sid) %>% unlist()

bf_mri <- filter(bf_mri, sid %in% ab_pos_ids)
bf_dpm <- read.csv("~/R/EDAP-data/BioFINDER/DPM_BioFINDER.csv") %>% distinct(sid, time_shift)
bf_mri <- left_join(bf_mri, bf_dpm, by="sid") %>% filter(!is.na(time_shift)) %>%
  mutate(Time = Years + time_shift) 
  #rename(RID = sid)

bf_sids <- unique(bf_mri$sid)

bf_table <- left_join(bf_dx, bf_demo, by="sid") %>% 
  left_join(bf_ab, by=c("sid", "Visit")) %>%
  filter(sid %in% bf_sids)

bf_follow_up <- bf_mri %>% group_by(sid) %>% 
  summarise(n_mri=n(),
            follow_up=max(Years, na.rm=TRUE) - min(Years, na.rm=TRUE)) %>% ungroup()

bf_table <- left_join(bf_table, bf_follow_up, by="sid") %>%
  rename(MMSE = mmse_score, 
         AGE = age_bl,
         EDUC = education_level_years_baseline_variable,
         GENDER = gender_baseline_variable,
         DX.bl = diag.bl,
         AB_PET_SUVR = fnc_ber_com_composite) %>%
  mutate(GENDER = ifelse(GENDER == 0, "Male", "Female"),
         TRACER = ifelse(!is.na(AB_PET_SUVR), "FMM", NA))


table(adni_table$DX.bl, useNA = "ifany")
table(oasis_table$DX.bl, useNA = "ifany")
table(bf_table$DX.bl, useNA = "ifany")

table(adni_table$GENDER, useNA = "ifany")
table(oasis_table$GENDER, useNA = "ifany")
table(bf_table$GENDER, useNA = "ifany")

summarize_cohort <- function(table_df, id_col = "RID") {
  # These inputs contain repeated visits; reduce to one row per participant
  # before summarizing demographics and MRI follow-up.
  subjects <- table_df %>%
    filter(!is.na(.data[[id_col]])) %>%
    distinct(.data[[id_col]], .keep_all=TRUE)

  as_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))
  mean_sd_latex <- function(x, digits=2) {
    x <- as_numeric(x)
    x <- x[is.finite(x)]
    if (!length(x)) return("--")
    sd_value <- if (length(x) > 1) stats::sd(x) else 0
    sprintf(paste0("$%.", digits, "f \\pm %.", digits, "f$"), mean(x), sd_value)
  }
  percent_latex <- function(n, denominator) {
    if (!denominator) return("--")
    sprintf("$%.1f\\%%$", 100*n/denominator)
  }
  count_latex <- function(x) sprintf("$%d$", as.integer(x))
  visit_col <- intersect(c("VISCODE2", "days_to_visit", "Visit"), names(table_df))[1]

  # Keep repeated MMSE observations. Remove only duplicate copies of the same
  # subject/visit/value that can result from joining PET data onto visit rows.
  if ("MMSE" %in% names(table_df)) {
    mmse_observations <- table_df %>%
      transmute(Subject=as.character(.data[[id_col]]),
                VisitKey=if (!is.na(visit_col)) as.character(.data[[visit_col]]) else NA_character_,
                MMSE=as_numeric(MMSE)) %>%
      filter(!is.na(Subject), is.finite(MMSE)) %>%
      distinct(Subject, VisitKey, MMSE) %>%
      pull(MMSE)
  } else {
    mmse_observations <- numeric()
  }

  metric_rows <- tibble(
    Variable=c("N", "Age (Y, at baseline)", "Education (Y)", "Female (\\%)",
               "Num. MRIs", "Follow up (Y)",
               "APOE-e4 (0)", "APOE-e4 (1)", "APOE-e4 (2)", 
               "MMSE"),
    Values=c(
      count_latex(n_distinct(subjects[[id_col]])),
      mean_sd_latex(subjects$AGE),
      mean_sd_latex(subjects$EDUC),
      percent_latex(sum(as.character(subjects$GENDER) == "Female", na.rm=TRUE),
                    sum(!is.na(subjects$GENDER))),
      mean_sd_latex(subjects$n_mri),
      mean_sd_latex(subjects$follow_up),
      percent_latex(sum(as_numeric(subjects$APOE4) == 0, na.rm=TRUE),
                    sum(!is.na(subjects$APOE4))),
      percent_latex(sum(as_numeric(subjects$APOE4) == 1, na.rm=TRUE),
                    sum(!is.na(subjects$APOE4))),
      percent_latex(sum(as_numeric(subjects$APOE4) == 2, na.rm=TRUE),
                    sum(!is.na(subjects$APOE4))),
      mean_sd_latex(mmse_observations)
    )
  )

  # Use all distinct visit-level PET measurements for population-level means
  # and SDs; remove only exact duplicates introduced by table joins.
  tracer_col <- intersect(c("TRACER", "tracer"), names(table_df))[1]
  value_col <- intersect(c("SUMMARY_SUVR", "PET_SUVR", "AB_PET_SUVR"), names(table_df))[1]
  pet_rows <- tibble(Variable=c("A$\\beta$-PET (FBP)", "A$\\beta$-PET (FMM)",
                                "A$\\beta$-PET (PiB)"), Values=rep("--", 3))
  if (!is.na(tracer_col) && !is.na(value_col)) {
    pet_summary <- table_df %>%
      transmute(Subject=as.character(.data[[id_col]]),
                VisitKey=if (!is.na(visit_col)) as.character(.data[[visit_col]]) else NA_character_,
                Tracer=toupper(as.character(.data[[tracer_col]])),
                SUVR=as_numeric(.data[[value_col]])) %>%
      mutate(Tracer=case_when(
        Tracer %in% c("AV45", "FLORBETAPIR", "FBP") ~ "FBP",
        Tracer %in% c("FLUTEMETAMOL", "FMM") ~ "FMM",
        Tracer %in% c("PIB", "PITTSBURGH COMPOUND B") ~ "PIB",
        TRUE ~ Tracer
      )) %>%
      filter(!is.na(Subject), Tracer %in% c("FBP", "FMM", "PIB"), is.finite(SUVR)) %>%
      distinct(Subject, VisitKey, Tracer, SUVR) %>%
      group_by(Tracer) %>%
      summarise(Value=mean_sd_latex(SUVR, digits=3), .groups="drop")
    tracer_labels <- c("FBP"="A$\\beta$-PET (FBP)", "FMM"="A$\\beta$-PET (FMM)",
                       "PIB"="A$\\beta$-PET (PiB)")
    if (nrow(pet_summary)) {
      pet_rows <- pet_rows %>% rows_update(
        pet_summary %>% transmute(Variable=unname(tracer_labels[Tracer]), Values=Value),
        by="Variable")
    }
  }

  bind_rows(metric_rows, pet_rows) %>% select(Variable, Values)
}

adni_summary <- summarize_cohort(adni_table, id_col="RID") %>% rename(ADNI=Values)
oasis_summary <- summarize_cohort(oasis_table, id_col="RID") %>% rename(OASIS=Values)
bf_summary <- summarize_cohort(bf_table, id_col="sid") %>% rename(BioFINDER=Values)

cohort_summary <- adni_summary %>%
  full_join(oasis_summary, by="Variable") %>%
  full_join(bf_summary, by="Variable") %>%
  select(Variable, ADNI, OASIS, BioFINDER)

knitr::kable(cohort_summary, format="latex", escape=FALSE, booktabs=TRUE,
             align=c("l", "c", "c", "c"))


# CN
adni_cn <- adni_table %>% filter(DX.bl == "CN") %>% summarize_cohort(id_col="RID") %>% rename(ADNI=Values)
oasis_cn <- oasis_table %>% filter(DX.bl == "CN") %>% summarize_cohort(id_col="RID") %>% rename(OASIS=Values)
bf_cn <- bf_table %>% filter(DX.bl %in%  c("Normal", "SCD")) %>% summarize_cohort(id_col="sid") %>% rename(BioFINDER=Values)

cohort_summary <- adni_cn %>%
  full_join(oasis_cn, by="Variable") %>%
  full_join(bf_cn, by="Variable") %>%
  select(Variable, ADNI, OASIS, BioFINDER)

knitr::kable(cohort_summary, format="latex", escape=FALSE, booktabs=TRUE,
             align=c("l", "l", "l", "l"))


# MCI
adni_mci <- adni_table %>% filter(DX.bl == "MCI") %>% summarize_cohort(id_col="RID") %>% rename(ADNI=Values)
oasis_mci <- oasis_table %>% filter(DX.bl %in% c("Impaired", "MCI")) %>% summarize_cohort(id_col="RID") %>% rename(OASIS=Values)
bf_mci <- bf_table %>% filter(DX.bl == "MCI") %>% summarize_cohort(id_col="sid") %>% rename(BioFINDER=Values)

cohort_summary <- adni_mci %>%
  full_join(oasis_mci, by="Variable") %>%
  full_join(bf_mci, by="Variable") %>%
  select(Variable, ADNI, OASIS, BioFINDER)

knitr::kable(cohort_summary, format="latex", escape=FALSE, booktabs=TRUE,
             align=c("l", "l", "l", "l"))


# Dementia
adni_dementia <- adni_table %>% filter(DX.bl == "Dementia") %>% summarize_cohort(id_col="RID") %>% rename(ADNI=Values)
oasis_dementia <- oasis_table %>% filter(DX.bl == "Dementia") %>% summarize_cohort(id_col="RID") %>% rename(OASIS=Values)
bf_dementia <- bf_table %>% filter(DX.bl == "Dementia") %>% summarize_cohort(id_col="sid") %>% rename(BioFINDER=Values)

cohort_summary <- adni_dementia %>%
  full_join(oasis_dementia, by="Variable") %>%
  full_join(bf_dementia, by="Variable") %>%
  select(Variable, ADNI, OASIS, BioFINDER)

knitr::kable(cohort_summary, format="latex", escape=FALSE, booktabs=TRUE,
             align=c("l", "l", "l", "l"))


