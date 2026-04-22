library(dplyr)
library(ggplot2)
library(stringr)
OASIS_DIR = "~/R/EDAP-data/OASIS/OASIS3_data_files/scans/"

load_oasis_mri_data <- function(unified_norm=TRUE) {
  OASIS_DIR = "~/R/EDAP-data/OASIS/OASIS3_data_files/scans/"
  
  oasis_fs_data <- read.csv(paste0(OASIS_DIR, "FS-Freesurfer_output/resources/csv/files/OASIS3_Freesurfer_output.csv")) %>%
    rename(OASISID = Subject) %>%
    mutate(days_to_visit = as.numeric(str_extract(MR_session, "\\d+$"))) %>%
    mutate(Months = round(days_to_visit/30.5), Years = Months/12)
  
  mri_cols <- c(grepv("^[lr]h_[a-z]+_volume$", colnames(oasis_fs_data)), 
                grepv("^(Left|Right)\\.[a-zA-Z._]+_volume$", colnames(oasis_fs_data)))
  mri_cols <- setdiff(mri_cols, grepv("vessel|Vent|WM|choroid", mri_cols))
  
  oasis_dx <- read.csv(paste0(OASIS_DIR, "UDSd1-Form_D1__Clinician_Diagnosis___Cognitive_Status_and_Dementia/resources/csv/files/OASIS3_UDSd1_diagnoses.csv")) %>%
    select(OASISID, days_to_visit, age.at.visit, NORMCOG, DEMENTED, MCIAMEM, MCIAPLUS, MCINON1, MCINON2, IMPNOMCI, alzdis) %>%
    filter_out(DEMENTED == 1 & alzdis == 0) 
  
  oasis_bl <- filter(oasis_dx, days_to_visit == 0) %>% filter(if_any(NORMCOG:IMPNOMCI, ~ !is.na(.))) %>%
    mutate(across(NORMCOG:IMPNOMCI, ~ ifelse(is.na(.), 0, .)),
           across(c(NORMCOG, MCIAMEM, MCIAPLUS, MCINON1, MCINON2, IMPNOMCI), ~ ifelse(DEMENTED == 1, 0, .)), ) %>% 
    rowwise() %>% mutate(MCI = sum(MCIAMEM, MCIAPLUS, MCINON1, MCINON2),
                         hasdiag = sum(NORMCOG, DEMENTED, MCIAMEM, MCIAPLUS, MCINON1, MCINON2, IMPNOMCI)) %>%
    filter(hasdiag != 0) %>%
    rename(CN = NORMCOG, AD = DEMENTED, PMCI = IMPNOMCI) %>%
    select(OASISID, CN, PMCI, MCI, AD) %>%
    mutate(DX.bl = case_when(
      CN == 1 ~ "CN",
      PMCI == 1 ~ "Impaired",
      MCI == 1 ~ "MCI",
      AD == 1 ~ "Dementia"
    ))
  
  oasis_fs_data <- filter(oasis_fs_data, OASISID %in% oasis_bl$OASISID)
  
  
  if (unified_norm) {
    b_weight = 1.513e-5 # Regression coefficient from data dictionary
    oasis_fs_data <- select(oasis_fs_data, OASISID, days_to_visit, "IntraCranialVol", all_of(mri_cols)) %>%
      group_by(OASISID) %>% 
      mutate(meanICV = mean(IntraCranialVol)) %>% ungroup() %>%
      mutate_at(mri_cols, ~ . - b_weight*(IntraCranialVol-meanICV))
  } else {
    oasis_fs_data <- select(oasis_fs_data, OASISID, days_to_visit, "IntraCranialVol", all_of(mri_cols)) %>%
      mutate_at(mri_cols, ~ . /IntraCranialVol)
  }
  
  oasis_fs_data <- inner_join(oasis_fs_data, oasis_bl)
  
  oasis_fs_data$RID <- as.numeric(gsub("^OAS", "", oasis_fs_data$OASISID))
  oasis_fs_data <- mutate(oasis_fs_data, Months = round(days_to_visit/30.5), Years = Months/12)
    
  return(oasis_fs_data)
}

get_oasis_ab <- function(use_centiloid = TRUE) {
  OASIS_DIR = "~/R/EDAP-data/OASIS/OASIS3_data_files/scans/"
  
  if (use_centiloid) {
    oasis_pet <- read.csv(paste0(OASIS_DIR, "Centiloid-Amyloid_Centiloid_Values/resources/csv/files/OASIS3_amyloid_centiloid.csv")) %>%
      mutate(days_to_visit = as.numeric(str_extract(oasis_session_id, "\\d+$"))) %>%
      rename(OASISID = subject_id,
             Centiloid_SUVR = Centiloid_fSUVR_rsf_TOT_CORTMEAN) %>% select(OASISID, days_to_visit, Centiloid_SUVR, tracer) %>%
      mutate(AB_pos = ifelse(tracer=="PIB", Centiloid_SUVR > 16.4, Centiloid_SUVR > 20.6))
    
  } else {
    oasis_pet <- read.csv(paste0(OASIS_DIR, "PUP-PUP_output/resources/csv/files/OASIS3_PUP.csv")) %>%
      select(PUP_PUPTIMECOURSEDATA.ID, tracer, FSId, MRId, Centil_fSUVR_rsf_TOT_CORTMEAN, PET_fSUVR_rsf_TOT_CORTMEAN) %>%
      mutate(days_to_visit = as.numeric(str_extract(PUP_PUPTIMECOURSEDATA.ID, "\\d+$")),
             OASISID = str_extract(PUP_PUPTIMECOURSEDATA.ID, "^OAS\\d+")) %>%
      select(OASISID, tracer, days_to_visit, PET_fSUVR_rsf_TOT_CORTMEAN) %>%
      rename(PET_SUVR = PET_fSUVR_rsf_TOT_CORTMEAN) %>%
      mutate(AB_pos = ifelse(tracer=="PIB", PET_SUVR > 1.42, PET_SUVR > 1.19))
      #pivot_wider(names_from = tracer, values_from = PET_fSUVR_rsf_TOT_CORTMEAN)
  }
  oasis_pet$tracer <- factor(oasis_pet$tracer)
  return(oasis_pet)
}

get_ab_pos_ids <- function() {
  ab_pos_ids <- get_oasis_ab(use_centiloid = TRUE) %>% filter(AB_pos) %>% select(OASISID) 
  ab_pos_ids <- get_oasis_ab(use_centiloid = FALSE) %>% filter(AB_pos) %>% select(OASISID) %>%
    rbind(ab_pos_ids) %>% distinct() %>% unlist()
  return(ab_pos_ids)
}


get_dpm_data <- function(use_centiloid=TRUE) {
  OASIS_DIR = "~/R/EDAP-data/OASIS/OASIS3_data_files/scans/"
  
  oasis_dx <- read.csv(paste0(OASIS_DIR, "UDSd1-Form_D1__Clinician_Diagnosis___Cognitive_Status_and_Dementia/resources/csv/files/OASIS3_UDSd1_diagnoses.csv")) %>%
    select(OASISID, days_to_visit, age.at.visit, NORMCOG, DEMENTED, MCIAMEM, MCIAPLUS, MCINON1, MCINON2, IMPNOMCI, alzdis) %>%
    filter_out(DEMENTED == 1 & alzdis == 0) 
  
  oasis_bl <- filter(oasis_dx, days_to_visit == 0) %>% filter(if_any(NORMCOG:IMPNOMCI, ~ !is.na(.))) %>%
    mutate(across(NORMCOG:IMPNOMCI, ~ ifelse(is.na(.), 0, .)),
           across(c(NORMCOG, MCIAMEM, MCIAPLUS, MCINON1, MCINON2, IMPNOMCI), ~ ifelse(DEMENTED == 1, 0, .)), ) %>% 
    rowwise() %>% mutate(MCI = sum(MCIAMEM, MCIAPLUS, MCINON1, MCINON2),
                         hasdiag = sum(NORMCOG, DEMENTED, MCIAMEM, MCIAPLUS, MCINON1, MCINON2, IMPNOMCI)) %>%
    filter(hasdiag != 0) %>%
    rename(CN = NORMCOG, AD = DEMENTED, PMCI = IMPNOMCI) %>%
    select(OASISID, CN, PMCI, MCI, AD) %>%
    mutate(DX.bl = case_when(
      CN == 1 ~ "CN",
      PMCI == 1 ~ "Impaired",
      MCI == 1 ~ "MCI",
      AD == 1 ~ "Dementia"
    ))
  
  
  oasis_cog <- read.csv(paste0(OASIS_DIR, "UDSb4-Form_B4__Global_Staging__CDR__Standard_and_Supplemental/resources/csv/files/OASIS3_UDSb4_cdr.csv")) %>%
    select(OASISID, days_to_visit, MMSE, CDRSUM) %>% filter(OASISID %in% oasis_bl$OASISID)
  
  if (use_centiloid) {
    oasis_pet <- read.csv(paste0(OASIS_DIR, "Centiloid-Amyloid_Centiloid_Values/resources/csv/files/OASIS3_amyloid_centiloid.csv")) %>%
      mutate(days_to_visit = as.numeric(str_extract(oasis_session_id, "\\d+$"))) %>%
      rename(OASISID = subject_id,
             Centiloid_SUVR = Centiloid_fSUVR_rsf_TOT_CORTMEAN) %>% select(OASISID, days_to_visit, Centiloid_SUVR) %>%
      filter(OASISID %in% oasis_bl$OASISID)
  } else {
    oasis_pet <- read.csv(paste0(OASIS_DIR, "PUP-PUP_output/resources/csv/files/OASIS3_PUP.csv")) %>%
      select(PUP_PUPTIMECOURSEDATA.ID, tracer, FSId, MRId, Centil_fSUVR_rsf_TOT_CORTMEAN, PET_fSUVR_rsf_TOT_CORTMEAN) %>%
      mutate(days_to_visit = as.numeric(str_extract(PUP_PUPTIMECOURSEDATA.ID, "\\d+$")),
             OASISID = str_extract(PUP_PUPTIMECOURSEDATA.ID, "^OAS\\d+")) %>%
      select(OASISID, tracer, days_to_visit, PET_fSUVR_rsf_TOT_CORTMEAN) %>%
      pivot_wider(names_from = tracer, values_from = PET_fSUVR_rsf_TOT_CORTMEAN) %>%
      filter(OASISID %in% oasis_bl$OASISID)
  }
  
  oasis_dpm <- full_join(oasis_cog, oasis_pet, by=c("OASISID", "days_to_visit")) %>%
    left_join(oasis_bl, by="OASISID") %>%
    mutate(Months = round(days_to_visit/30.5), 
           Years = Months/12) 
}

get_demographics <- function() {
  oasis_demog <- read.csv(paste0(OASIS_DIR, "/demo-demographics/resources/csv/files/OASIS3_demographics.csv")) %>%
    mutate(APOE = as.numeric(APOE)) %>%
    mutate(APOE4 = case_when(
      APOE == 44 ~ 2,
      APOE == 34 ~ 1,
      APOE == 24 ~ 1,
      APOE == 33 ~ 0,
      APOE < 24 ~ 0
    ))
  oasis_demog$RID <- as.numeric(gsub("^OAS", "", oasis_demog$OASISID))
  return(oasis_demog)
}

get_diagnoses <- function() {
  oasis_dx <- read.csv(paste0(OASIS_DIR, "UDSd1-Form_D1__Clinician_Diagnosis___Cognitive_Status_and_Dementia/resources/csv/files/OASIS3_UDSd1_diagnoses.csv")) %>%
    select(OASISID, days_to_visit, age.at.visit, NORMCOG, DEMENTED, MCIAMEM, MCIAPLUS, MCINON1, MCINON2, IMPNOMCI, alzdis) %>%
    filter_out(DEMENTED == 1 & alzdis == 0) %>% filter(if_any(NORMCOG:IMPNOMCI, ~ !is.na(.))) %>%
    mutate(across(NORMCOG:IMPNOMCI, ~ ifelse(is.na(.), 0, .)),
           across(c(NORMCOG, MCIAMEM, MCIAPLUS, MCINON1, MCINON2, IMPNOMCI), ~ ifelse(DEMENTED == 1, 0, .)), ) %>% 
    rowwise() %>% mutate(MCI = sum(MCIAMEM, MCIAPLUS, MCINON1, MCINON2),
                         hasdiag = sum(NORMCOG, DEMENTED, MCIAMEM, MCIAPLUS, MCINON1, MCINON2, IMPNOMCI)) %>%
    filter(hasdiag != 0) %>%
    mutate(DX = case_when(
      NORMCOG == 1 ~ "CN",
      IMPNOMCI == 1 ~ "Impaired",
      MCI == 1 ~ "MCI",
      DEMENTED == 1 ~ "Dementia"
    ))
  
  oasis_dx <- oasis_dx %>% group_by(OASISID) %>% mutate(DX.bl = ifelse(min(days_to_visit) == 0, DX, NA)) %>% ungroup()
  
  oasis_dx$RID <- as.numeric(gsub("^OAS", "", oasis_dx$OASISID))
  
  return(oasis_dx)
}

get_tau_pet <- function(get_braak=FALSE) {
  if (get_braak) {
    oasis_tau <- read.csv(paste0(OASIS_DIR, "BRAAK-BRAAK_BRAAK_TAUOPATHY/resources/csv/files/OASIS3_AV1451_braak_tauopathy.csv")) %>%
      rename(OASISID = OASIS_ID) 
  } else {
    oasis_tau <- read.csv(paste0(OASIS_DIR, "PUP-AV1451-PUP_output/resources/csv/files/OASIS3_AV1451_PUP.csv")) %>%
      mutate(days_to_visit = as.numeric(str_extract(PUP_PUPTIMECOURSEDATA.ID, "\\d+$")),
             OASISID = str_extract(PUP_PUPTIMECOURSEDATA.ID, "^OAS\\d+")) 
    
    cols <- grepv("PET_fSUVR_rsf", colnames(oasis_tau))
    oasis_tau <- select(oasis_tau, OASISID, tracer, days_to_visit, all_of(cols))
  }
  
  oasis_tau$RID <- as.numeric(gsub("^OAS", "", oasis_tau$OASISID))
  
  return(oasis_tau)
}

get_copath <- function() {
  oasis_bl <- get_diagnoses() %>% distinct(OASISID, RID, DX.bl) 
  
  oasis_cp <- read.csv(paste0(OASIS_DIR, "UDSd1-Form_D1__Clinician_Diagnosis___Cognitive_Status_and_Dementia/resources/csv/files/OASIS3_UDSd1_diagnoses.csv")) %>%
    distinct(OASISID, amndem, lbdis, cvd) %>% group_by(OASISID) %>% mutate(
      amndem = ifelse(any(!is.na(amndem)), max(amndem, na.rm=TRUE), NA),
      lbdis = ifelse(any(!is.na(lbdis)), max(lbdis, na.rm=TRUE), NA),
      cvd = ifelse(any(!is.na(cvd)), max(cvd, na.rm=TRUE), NA)
    ) %>% ungroup() %>% 
    left_join(oasis_bl, by="OASISID") %>% distinct()
  
  return(oasis_cp)
  
}
