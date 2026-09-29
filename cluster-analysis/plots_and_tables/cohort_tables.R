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


# Summarize continuous measures as mean (SD), categorical measures as n (%),
# and retain a participant count for each diagnosis group.
summarize_cohort <- function(data, cohort, id, group, continuous, categorical = character()) {
  
}

# ADNI: baseline demographics and diagnosis history; MRI and biomarker values
# are summarized per person when the relevant loader data are available.
adni_dx <- adni_dl$get_diagnoses()

adni_demo <- adni_dl$get_demographics() %>% distinct(RID, .keep_all = TRUE)

adni_ab <- adni_dl$get_ab_df()
adni_tau <- adni_dl$get_tau_pet()

adni_rids <- multi_cohort_df %>% filter(Cohort == "ADNI") %>% mutate(RID = as.numeric(gsub("ADNI_", "", RID))) %>%
  distinct(RID) %>% unlist()

# OASIS: baseline diagnosis plus all available diagnosis visits to determine
# whether AD/dementia was ever recorded.
oasis_dx <- oasis_dl$get_diagnoses()

oasis_demo <- oasis_dl$get_demographics()

oasis_ab <- oasis_dl$get_oasis_ab()
oasis_tau <- oasis_dl$get_tau_pet(get_braak = TRUE)

oasis_rids <- multi_cohort_df %>% filter(Cohort == "OASIS") %>% mutate(RID = as.numeric(gsub("OASIS_", "", RID))) %>%
  distinct(RID) %>% unlist()

# BioFINDER
bf_dx <- bf_dl$get_diagnoses()

bf_demo <- bf_dl$get_demographics()

bf_ab <- bf_dl$get_ab_df()
bf_tau <- bf_dl$get_tau_data()

bf_mri <- bf_dl$get_mri_data_updated(normalize = TRUE) 
ab_pos_ids <- bf_ab %>% group_by(sid) %>% mutate(AB_any = any(AB)) %>% ungroup() %>%
  filter(AB_any) %>% select(sid) %>% unlist()

bf_mri <- filter(bf_mri, sid %in% ab_pos_ids)
bf_dpm <- read.csv("~/R/EDAP-data/BioFINDER/DPM_BioFINDER.csv") %>% distinct(sid, time_shift)
bf_mri <- left_join(bf_mri, bf_dpm, by="sid") %>% filter(!is.na(time_shift)) %>%
  mutate(Time = Years + time_shift) %>%
  rename(RID = sid)

bf_rids <- unique(bf_mri$RID)











