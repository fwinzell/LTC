library(dplyr)
library(lubridate)
library(stringr)
library(readxl)


# Amyloid-Beta status
get_ab_df <- function(path = "~/R/EDAP-data/BioFINDER/data_for_Filip/data_filip.csv") {
  biofinder <- read.csv(path)
  
  ab_df <- select(biofinder, sid, Visit, cognitive_status_baseline_variable,
                  mmse_score, cdr_global_clinical,
                  ab_status, ab_binary_PET, ab_binary_CSF, Dx_ab, 
                  fnc_ber_com_composite, 
                  CSF_Ab42_Ab40_ratio_imputed_Elecsys_2020_2022
  ) %>%
    filter(!if_all(c(ab_status, ab_binary_PET, ab_binary_CSF, fnc_ber_com_composite, CSF_Ab42_Ab40_ratio_imputed_Elecsys_2020_2022), is.na)) %>%
    mutate(AB_CSF = CSF_Ab42_Ab40_ratio_imputed_Elecsys_2020_2022 < 0.080,
           AB_PET = fnc_ber_com_composite > 1.03) %>%
    rowwise() %>%
    mutate(CI = any(cognitive_status_baseline_variable %in% c("MCI", "Dementia"), 
                    cdr_global_clinical > 0.0, 
                    mmse_score < 26, 
                    na.rm=TRUE))
  
  ab_df <- ab_df %>% rowwise() %>% mutate(AB = ifelse(CI, # Cognitively impaired?
                                                      any(AB_PET, na.rm=TRUE), # yes - require AB PET +
                                                      any(AB_PET, AB_CSF, na.rm=TRUE))) # no - any biomarker
  
  #table(ab_df[c('ab_status', 'AB')])
  return(ab_df)
  
}

# Disease Progression Modeling data
get_biofinder_dpm <- function(path = "~/R/EDAP-data/BioFINDER/data_for_Filip/data_filip.csv") {
  biofinder <- read.csv(path)
  
  ab_df <- get_ab_df(path)
  
  dpm_df <- biofinder %>% select(uid, sid, subject_id, Visit, visit_date, 
                                 mmse_date, mmse_score, 
                                 adas_delayed_word_recall, adas_immediate_word_recall_average,
                                 fnc_ber_com_composite,
                                 diagnosis_baseline_variable, cognitive_status_baseline_variable,
                                 cdr_global_clinical, cdr_sum_of_boxes_clinical) %>%
    mutate(adas_score = adas_delayed_word_recall + adas_immediate_word_recall_average) %>%
    # Require at least one measurement of ADAS, MMSE or AB-PET
    filter(!if_all(c(adas_score, mmse_score, fnc_ber_com_composite), is.na)) %>%
    # Impute missing diagnosis
    mutate(diag = ifelse(
      cognitive_status_baseline_variable == "TBD",
        ifelse(cdr_global_clinical >= 1.0, "Dementia", 
             ifelse(cdr_global_clinical >= 0.5, "MCI", "Normal")),
      cognitive_status_baseline_variable
      ),
      diag = factor(diag, levels = c("Normal", "SCD", "MCI", "Dementia"))
    ) %>% drop_na(diag) 
  
  # Ensure everyone has a baseline visit and visit_date
  dpm_df <- dpm_df %>% group_by(sid) %>% 
    mutate(hasbl = any(Visit == 0)) %>% ungroup() %>% filter(hasbl) %>% select(-hasbl) %>%
    mutate(visit_date_filled = ifelse(trimws(visit_date) == "", mmse_date, visit_date)) 
  
  dpm_df <- dpm_df %>% arrange(sid, Visit) %>% group_by(sid) %>% 
    mutate(baseline_date = .data$visit_date[1],
           Years = interval(baseline_date, visit_date_filled) / years(1),
           Years = ifelse(is.na(Years), Visit, Years)) %>% ungroup() 
  
  missing_dx <- filter(dpm_df, cognitive_status_baseline_variable == "TBD")
  
  dpm_df <- ab_df %>% select(sid, Visit, AB) %>% mutate(ab_status = as.numeric(AB)) %>% 
    select(-AB) %>% right_join(dpm_df, by=c('sid', 'Visit'))
  
  # Impute any missing AB status by carrying negatives backwards, and positive forwards
  # assuming monoticity
  missing_ab <- dpm_df %>% filter(is.na(ab_status)) %>% distinct(sid) %>% unlist()
  
  for (id in missing_ab) {
    subj <- dpm_df %>% filter(sid == id) %>% arrange(Visit)
    idxs <- which(!is.na(subj$ab_status))
    if (length(idxs) > 0) {
      last_i = 1
      next_i = c(idxs, nrow(subj)+1)
      count = 1
      for (i in idxs) {
        count = count + 1
        if (subj$ab_status[i] == 0) {
          # Carry backwards
          subj[last_i:i, "ab_status"] <- subj$ab_status[i]
          last_i <- i
        } else {
          # Carry forwards
          subj[i:(next_i[count]-1), "ab_status"] <- subj$ab_status[i]
          last_i <- next_i[count]
        }
      }
      dpm_df <- rows_update(dpm_df, subj, by = c("sid", "Visit"))
    }
  }
  
  length(missing_ab)
  missing_ab <- dpm_df %>% filter(is.na(ab_status)) %>% distinct(sid) %>% unlist()
  length(missing_ab)
  
  dpm_df <- rowwise(dpm_df) %>% 
    mutate(DXCI = diag %in% c("MCI", "Dementia"), 
           CI = any(DXCI, cdr_global_clinical > 0.0, mmse_score < 26, na.rm=TRUE),
           invalid = !ab_status & CI) %>%
    filter(invalid == FALSE | is.na(invalid))
  
  # Indicator for negative AB baseline status
  ab.bl <- dpm_df %>% select(sid, Visit, ab_status) %>% drop_na(ab_status) %>% filter(Visit == 0) %>% 
    mutate(negAB.bl = 1-ab_status) %>% select(sid, negAB.bl) %>% distinct()
  
  dpm_df <- left_join(dpm_df, ab.bl, by = "sid") %>%
    drop_na(negAB.bl)
  
  dpm_df <- mutate(dpm_df,
    CU = ifelse(diag == "Normal", 1, 0),
    SCD = ifelse(diag == "SCD", 1, 0),
    MCI = ifelse(diag == "MCI", 1, 0),
    AD = ifelse(diag == "Dementia", 1, 0)) %>%
    rename(adas = adas_score,
           mmse = mmse_score,
           fnc_suvr = fnc_ber_com_composite)
  

  
  # Convert to long format
  cols <- c("sid", "Years", "diag", "CU", "SCD", "MCI", "AD", "negAB.bl")
  dpm_long <- rbind(dpm_df[c("adas", cols)] |> mutate(scale = "adas") |> rename(value = adas), 
                    dpm_df[c("mmse", cols)] |> mutate(scale = "mmse") |> rename(value = mmse),
                    dpm_df[c("fnc_suvr", cols)] |> mutate(scale = "fnc_suvr") |> rename(value = fnc_suvr))
  dpm_long <- drop_na(dpm_long, value)
  
  dpm_long <- dpm_long %>% mutate(
    scale = as.factor(scale),
    diag = as.factor(diag)
  )
  
  dpm_long$t <- dpm_long$Years
  
  return(dpm_df)
}


## MRI vars
get_mri_data <- function(normalize=TRUE) {
  biofinder <- read.csv("~/R/EDAP-data/BioFINDER/data_for_Filip/data_filip.csv")
  
  dx_df <- biofinder %>% select(sid, Visit, visit_date, 
                                 diagnosis_baseline_variable, cognitive_status_baseline_variable,
                                 cdr_global_clinical) %>%
    # Impute missing diagnosis
    mutate(diag.bl = ifelse(
      cognitive_status_baseline_variable == "TBD",
      ifelse(cdr_global_clinical >= 1.0, "Dementia", 
             ifelse(cdr_global_clinical >= 0.5, "MCI", "Normal")),
      cognitive_status_baseline_variable
    ),
    diag.bl = factor(diag.bl, levels = c("Normal", "SCD", "MCI", "Dementia"))
    ) %>% drop_na(diag.bl) %>% distinct(sid, diag.bl) 
  
  biofinder <- mutate(biofinder, visit_date_filled = ifelse(trimws(visit_date) == "", mmse_date, visit_date))
  bl_dates <- biofinder %>% 
    distinct(sid, Visit, visit_date_filled) %>% filter(Visit == 0) %>%
    select(-Visit) %>% rename(baseline_date = visit_date_filled)
  
  cortical_vols <- readxl::read_excel("~/R/EDAP-data/BioFINDER/data_for_Filip/cortical_vols_filip.xlsx") 
  cortical_vols <- filter(cortical_vols, !if_all(colnames(cortical_vols), is.na)) %>%
    mutate(mri_date = str_extract(csv_icv__index, "\\d{8}"),
           mri_date = ymd(mri_date)) 
  
  subcortical_vols <- grepv("samseg_vols", colnames(biofinder))
  subcortical_df <- biofinder %>% select(sid, Visit, cognitive_status_baseline_variable, diagnosis_baseline_variable,
                                           all_of(subcortical_vols), icv_mm3) %>%
    mutate(mri_date = str_extract(csv_samseg_vols__index, "\\d{8}"),
             mri_date = ymd(mri_date)) 
      
  subcortical_df <- filter(subcortical_df, !if_all(all_of(subcortical_vols[subcortical_vols != "csv_samseg_vols__index"]), is.na)) #%>%
    
  
  mri_df <- full_join(cortical_vols, subcortical_df, by=intersect(colnames(subcortical_df), colnames(cortical_vols))) %>%
    left_join(bl_dates, by="sid") %>%
    mutate(Years = interval(baseline_date, mri_date) / years(1))
  
  mri_df <- left_join(mri_df, dx_df, by='sid') 
  
  mri_df <- select(mri_df, 
                   -c("samseg_vols_3rd_Ventricle", "samseg_vols_4th_Ventricle", "samseg_vols_5th_Ventricle", 
                      "samseg_vols_CSF", "samseg_vols_Fluid_Inside_Eyes", "samseg_vols_Left_vessel", "samseg_vols_Right_vessel", 
                      "samseg_vols_Right_Inf_Lat_Vent", "samseg_vols_Right_Lateral_Ventricle", 
                      "samseg_vols_Left_Inf_Lat_Vent", "samseg_vols_Left_Lateral_Ventricle",
                      "samseg_vols_Soft_Nonbrain_Tissue", "samseg_vols_Unknown", "samseg_vols_Skull",                        
                      "samseg_vols_WM_hypointensities", "samseg_vols_non_WM_hypointensities", 
                      "samseg_vols_Right_Cerebral_Cortex", "samseg_vols_Right_Cerebral_White_Matter",
                      "samseg_vols_Left_Cerebral_Cortex", "samseg_vols_Left_Cerebral_White_Matter",
                      "samseg_vols_scan", "samseg_vols_total_scans"))
  
  names(mri_df) <- gsub("^aparc_grayvol_", "", names(mri_df))
  names(mri_df) <- gsub("^samseg_vols_", "", names(mri_df))
  
  convert_name <- function(x) {
    words <- regmatches(x, gregexpr("[A-Za-z]*", x))[[1]]
    
    if (length(words) == 0) {
      return(toupper(x))
    }
    
    if (words[1] == "Right") {
      region <- toupper(paste0(words[-1], collapse = ""))
      paste0("RH_", region)
    } else if (words[1] == "Left") {
      region <- toupper(paste0(words[-1], collapse = ""))
      paste0("LH_", region)
    } else if (words[length(words)] == "R") {
      region <- toupper(paste0(words[-length(words)], collapse = ""))
      paste0("RH_", region)
    } else if (words[length(words)] == "L") {
      region <- toupper(paste0(words[-length(words)], collapse = ""))
      paste0("LH_", region)
    } else {
      #toupper(paste0(words, collapse = ""))
      x
    }
  }
  
  test <- sapply(colnames(mri_df), convert_name)
  names(mri_df) <- sapply(names(mri_df), convert_name)
  
  mri_df <- rename(mri_df, 
                   OPTICCHIASM = Optic_Chiasm,
                   BRAINSTEM = Brain_Stem)
  
  mri_vars <- c(grepv("^(RH_|LH_)", colnames(mri_df)), "BRAINSTEM", "OPTICCHIASM")
  
  # Normalize by ICV
  if (normalize) {
    for (varname in mri_vars) {
      mri_df[[varname]] <- mri_df[[varname]] / mri_df$icv_mm3
    }
  }
  
  
  return(mri_df)
  
}

test_things <- function() {
  biofinder <- read.csv("~/R/EDAP-data/BioFINDER/data_for_Filip/data_filip.csv")
  
  cortical_vols <- grepv("aparc_ct", colnames(biofinder))
  subcortical_vols <- grepv("samseg_vols", colnames(biofinder))
  
  avg_scn <- biofinder$aparc_ct_avg_scan
  table(avg_scn, useNA= "ifany")
  skull <- biofinder$samseg_vols_Skull
  
  
  mmse_cols <- grepv("mmse", colnames(biofinder))
  adas_cols <- grepv("adas", colnames(biofinder))
  cdr_cols <- grepv("cdr", colnames(biofinder))
  
  date_cols <- grepv("date", colnames(biofinder))
  
  
  
  ab_pet_cols <- grepv("fnc", colnames(biofinder))
  
  ab_pet <- select(biofinder, subject_id, all_of(ab_pet_cols), ab_status, cognitive_status_baseline_variable) %>%
    filter(has_proc_csv_fncbw_sr_mr_fs == 1) %>%
    mutate(ab_status = as.factor(ab_status), 
           diag = as.factor(cognitive_status_baseline_variable))
  
  table(ab_pet$ab_status)
  
  library(ggplot2) 
  
  ggplot(ab_pet, aes(y=fnc_ber_com_composite, fill = diag, group = diag)) +
    geom_boxplot() +
    geom_hline(yintercept = 1.03, linetype = "dashed")
  
  
  
  # Amyloid-Beta stuff
  ab_df <- select(biofinder, sid, Visit, diagnosis_baseline_variable, ab_status, ab_binary_PET, ab_binary_CSF, Dx_ab, 
                  fnc_ber_com_composite, 
                  CSF_Ab42_Ab40_ratio_imputed_Elecsys_2020_2022
  ) %>%
    filter(!if_all(c(ab_status, ab_binary_PET, ab_binary_CSF, fnc_ber_com_composite, CSF_Ab42_Ab40_ratio_imputed_Elecsys_2020_2022), is.na)) %>%
    mutate(AB_CSF = CSF_Ab42_Ab40_ratio_imputed_Elecsys_2020_2022 < 0.080,
           AB_PET = fnc_ber_com_composite > 1.03) 
  
  
  table(ab_df[c("ab_binary_PET", "AB_PET")])
  
  wat <- filter(ab_df, ab_binary_PET == 0 & AB_PET == TRUE)
  
  table(ab_df[c("ab_binary_CSF", "AB_CSF")])
  
  wat <- filter(ab_df, (ab_binary_CSF == 0 & AB_CSF == TRUE) | (ab_binary_CSF == 1 & AB_CSF == FALSE))
  
  ggplot(ab_df, aes(y=fnc_ber_com_composite, fill = AB_PET, group = AB_PET)) +
    geom_boxplot() +
    geom_hline(yintercept = 1.03, linetype = "dashed")
  
  ggplot(ab_df, aes(y=CSF_Ab42_Ab40_ratio_imputed_Elecsys_2020_2022, fill = AB_CSF, group = AB_CSF)) +
    geom_boxplot() +
    geom_hline(yintercept = 0.080, linetype = "dashed")
  
  table(biofinder$Dx_ab)
  
  
  ### DX ###
  dx_df <- select(biofinder, sid, Visit,
                  diagnosis_baseline_variable, underlying_etiology_text_baseline_variable, 
                  etiology_genetic_variant_baseline_variable, etiology_extra_comment_baseline_variable,
                  cognitive_status_baseline_variable, converted_dementia_date_baseline_variable,
                  converted_dementia_dich_baseline_variable, last_visit_date_dementia_nonconverted_baseline_variable,
                  converted_MCI_date_baseline_variable, converted_MCI_dich_baseline_variable, 
                  last_visit_date_MCI_nonconverted_baseline_variable, PDrel_susp_dis_dich_baseline_baseline_variable, neuropathology_done_baseline_variable,                   
                  gds_clinical, cdr_global_clinical, cdr_sum_of_boxes_clinical, last_visit_date_Pdrel_dis_nonconverted_baseline_variable
                  ) 
  
  tbds <- filter(dx_df, cognitive_status_baseline_variable == "TBD") %>% distinct(sid) %>% unlist()
  
  dx_df <- filter(dx_df, sid %in% tbds)
  
  
  
  
}

