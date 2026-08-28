library(dplyr)
library(ADNIMERGE2)
library(lubridate)
library(stringr)

ucsf_longitudinal_all <- function(only_vol=FALSE, filter_n=1, normalize=TRUE) {
  # Update: subjects with PTID: 381_S_#### should be removed, ADNIMERGE2 is not updated ...
  remove_rids <- filter(ADNIMERGE2::REGISTRY, grepl("^381_S_.*", ADNIMERGE2::REGISTRY$PTID)) %>% select(RID, VISCODE2)
  
  # Load all Longitudinal UCSF datasets
  ucsf_data1 <- ADNIMERGE2::UCSFFSL51 # Final run w all 2022
  ucsf_data2 <- ADNIMERGE2::UCSFFSL # 2016
  
  cols <- intersect(colnames(ucsf_data1), colnames(ucsf_data2))
  ucsf_data1 <- ucsf_data1 %>% select(all_of(cols))
  ucsf_data2 <- ucsf_data2 %>% select(all_of(cols))
  
  ucsf_data <- rbind(ucsf_data1, ucsf_data2) #%>% distinct(RID, VISCODE2, .keep_all = TRUE)
  
  ucsf_data3 <- ADNIMERGE2::UCSFFSL51ALL # 2022 base images
  cols <- intersect(colnames(ucsf_data), colnames(ucsf_data3))
  ucsf_data3 <- select(ucsf_data3, all_of(cols)) 
  ucsf_data <- ucsf_data %>% rbind(ucsf_data3) #%>% distinct(RID, VISCODE2, .keep_all = TRUE)
  
  ucsf_data4 <- ADNIMERGE2::UCSFFSL51Y1 # 2016 base images
  cols <- intersect(colnames(ucsf_data), colnames(ucsf_data4))
  ucsf_data4 <- select(ucsf_data4, all_of(cols)) 
  ucsf_data <- ucsf_data %>% rbind(ucsf_data4) #%>% distinct(RID, VISCODE2, .keep_all = TRUE)
  
  ucsf_data$RID <- as.numeric(ucsf_data$RID)
  
  ucsf_data <- anti_join(ucsf_data, remove_rids, by=c("RID", "VISCODE2"))
  
  ucsf_data <- ADNIMERGE2::REGISTRY %>% filter(VISCODE2 == "bl") %>% select(RID, EXAMDATE) %>%
    rename(EXAMDATE.bl = EXAMDATE) %>% right_join(ucsf_data, by="RID") %>% drop_na(EXAMDATE.bl)
  
  ucsf_data <- mutate(ucsf_data, Months = interval(EXAMDATE.bl, EXAMDATE) %/% months(1))
  
  dx.df <- ADNIMERGE2::DXSUM %>% filter(VISCODE2 %in% c("bl", "sc")) %>% select(RID, VISCODE2, DIAGNOSIS) %>% 
    mutate(priority = case_when(
      VISCODE2 == "bl" ~ 1,
      VISCODE2 == "sc" ~ 2,
      TRUE ~ 3
    )) %>%
    arrange(RID, priority) %>%
    group_by(RID) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    select(RID, DIAGNOSIS, VISCODE2) %>% rename(DX.bl = DIAGNOSIS)
  
  ucsf_data <- dx.df %>% select(-VISCODE2) %>% right_join(ucsf_data, by="RID") %>% drop_na(DX.bl)
  
  
  # Drop the following variables
  #ST128SV - WMHypoIntensities
  #ST125SV - RightVessel
  #ST68SV - NonWMHypoIntensities
  #ST66SV - LeftVessel
  #ST7SV - CSF
  #ST80SV - RightChoroidPlexus
  #ST21SV - LeftChoroidPlexus
  #Ventricles
  #ST8SV
  #ST96SV
  #ST37SV
  #ST9SV
  #ST127SV
  #ST30SV
  
  ucsf_data <- ucsf_data |> select(-ST8SV, -ST128SV, -ST125SV, -ST68SV, 
                                   -ST66SV, -ST7SV, -ST80SV, -ST21SV, 
                                   -ST96SV, -ST37SV, -ST9SV, -ST127SV, 
                                   -ST30SV, -ST89SV)
  
  # Filter Bad Quality MRI
  # Pass: All regions OK
  # Patial: Some region have failed, remove failed regions
  ucsf_partial <- filter(ucsf_data, OVERALLQC == "Partial")
  
  # Remove failed regions from partial
  datadic <- datadict_as_tibble(ADNIMERGE2::DATADIC) %>% filter(TBLNAME == "UCSFFSL") 
  qctable <- datadic %>% filter(str_detect(FLDNAME, "QC")) %>% select(FLDNAME, TEXT) %>%
    mutate(ST_codes = str_extract_all(TEXT, "ST\\d+")) %>%
    tidyr::unnest(ST_codes)
  qcList <- c("BGQC", "TEMPQC", "CWMQC", "FRONTQC", "INSULAQC", "OCCQC", "PARQC", "VENTQC")
  for (qc_var in qcList) {
    codes <- qctable %>% filter(FLDNAME == qc_var) %>% select(ST_codes) %>% unlist()
    ucsf_partial <- ucsf_partial %>% mutate(across(starts_with(codes),
                                                   ~ if_else(.data[[qc_var]] == "Fail", NA, .)))
  }
  
  ucsf_data <- filter(ucsf_data, OVERALLQC == "Pass") %>% rbind(ucsf_partial) %>%
    arrange(RID, Months)
  
  # Filter duplicated visits by taking the most recent processing (2022 vs 2016)
  ucsf_data <- ucsf_data %>% group_by(RID, EXAMDATE) %>%
    slice_max(order_by = update_stamp, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  ucsf_data |> filter(!is.na(Months)) -> ucsf_data
  
  icv <- ucsf_data %>% select(RID, VISCODE2, ST10CV)
  ucsf_data <- select(ucsf_data, -ST10CV)
  
  vol.idxs <- grepl("^ST\\d+(CV|SV)$", colnames(ucsf_data))
  vol.vars <- colnames(ucsf_data)[vol.idxs == 1]
  all.idxs <- grepl("^ST\\d+(CV|SV|SA|TA)$", colnames(ucsf_data))
  all.vars <- colnames(ucsf_data)[all.idxs == 1]
  
  # Normalize by ICV
  if (normalize) {
    for (varname in vol.vars) {
      ucsf_data[[varname]] <- ucsf_data[[varname]] / icv$ST10CV
    }
  }
  
  if (only_vol) {
    ucsf_data <- ucsf_data %>% select(RID, DX.bl, Months, all_of(vol.vars)) %>%
      filter(!if_all(all_of(vol.vars), is.na))
  } else {
    ucsf_data <- ucsf_data %>% select(RID, DX.bl, Months, all_of(all.vars)) %>%
      filter(!if_all(all_of(all.vars), is.na))
  }
  
  filtered_ids <- ucsf_data |>
    group_by(RID) |>
    tally() |>
    filter(n >= filter_n) |>
    pull(RID)
  
  ucsf_data |> filter(RID %in% filtered_ids) -> ucsf_data
  
  ucsf_data$Years <- ucsf_data$Months/12
  return(ucsf_data)
}

ucsf_xs_all <- function(only_vol=FALSE, filter_n=0) {
  ucsf_data <- ADNIMERGE2::UCSFFSX7
  
  ucsf_data <- ucsf_data |> select(-ST8SV, -ST128SV, -ST125SV, -ST68SV, 
                                   -ST66SV, -ST7SV, -ST80SV, -ST21SV, 
                                   -ST96SV, -ST37SV, -ST9SV, -ST127SV, 
                                   -ST30SV, -ST89SV)
  
  ucsf_data <- ucsf_data %>%
    mutate(
      M = VISCODE2 %>%
        gsub("scmri|bl|sc|m", "0", .) %>%
        as.numeric()
    ) 
  
  # Calculate time since baseline
  # Remove participants without baseline visit (only screening not sufficient)
  
  #missing_m <- ucsf_data %>% filter(is.na(Months)) %>% select(RID, EXAMDATE)
  ucsf_data <- ADNIMERGE2::REGISTRY %>% filter(VISCODE2 == "bl") %>% select(RID, EXAMDATE) %>%
    rename(EXAMDATE.bl = EXAMDATE) %>% right_join(ucsf_data, by="RID") %>% drop_na(EXAMDATE.bl)
  
  ucsf_data <- mutate(ucsf_data, Months = interval(EXAMDATE.bl, EXAMDATE) %/% months(1))
  
  dx.df <- ADNIMERGE2::DXSUM %>% filter(VISCODE2 %in% c("bl", "sc")) %>% select(RID, VISCODE2, DIAGNOSIS) %>% 
    mutate(priority = case_when(
      VISCODE2 == "bl" ~ 1,
      VISCODE2 == "sc" ~ 2,
      TRUE ~ 3
    )) %>%
    arrange(RID, priority) %>%
    group_by(RID) %>%
    slice(1) %>%
    ungroup() %>%
    select(RID, DIAGNOSIS, VISCODE2) %>% rename(DX.bl = DIAGNOSIS)
  
  ucsf_data <- dx.df %>% select(-VISCODE2) %>% right_join(ucsf_data, by="RID") %>% drop_na(DX.bl)
  
  ucsf_data |> filter(OVERALLQC != "Fail" | is.na(OVERALLQC)) -> ucsf_data
  
  icv <- ucsf_data %>% select(RID, VISCODE2, ST10CV)
  ucsf_data <- select(ucsf_data, -ST10CV)
  
  vol.idxs <- grepl("^ST\\d+(CV|SV)$", colnames(ucsf_data))
  vol.vars <- colnames(ucsf_data)[vol.idxs == 1]
  all.idxs <- grepl("^ST\\d+(CV|SV|SA|TA)$", colnames(ucsf_data))
  all.vars <- colnames(ucsf_data)[all.idxs == 1]
  
  # Normalize by ICV
  for (varname in vol.vars) {
    ucsf_data[[varname]] <- ucsf_data[[varname]] / icv$ST10CV
  }
  
  if (only_vol) {
    ucsf_data <- ucsf_data %>% select(RID, M, DX.bl, Months, all_of(vol.vars)) %>%
      filter(!if_all(all_of(vol.vars), is.na))
  } else {
    ucsf_data <- ucsf_data %>% select(RID, M, DX.bl, Months, all_of(all.vars)) %>%
      filter(!if_all(all_of(all.vars), is.na))
  }
  
  filtered_ids <- ucsf_data |>
    group_by(RID) |>
    tally() |>
    filter(n >= filter_n) |>
    pull(RID)
  
  ucsf_data |> filter(RID %in% filtered_ids) -> ucsf_data
  
  ucsf_data$Years <- ucsf_data$Months/12
  return(ucsf_data)
}

get_ab_df <- function() {
  # Update: subjects with PTID: 381_S_#### should be removed, ADNIMERGE2 is not updated ...
  remove_rids <- filter(ADNIMERGE2::REGISTRY, grepl("^381_S_.*", ADNIMERGE2::REGISTRY$PTID)) %>% select(RID, VISCODE2)
  
  
  bl.df <- ADNIMERGE2::REGISTRY %>% filter(VISCODE2 == "bl") %>% select(RID, EXAMDATE) %>%
    rename(EXAMDATE.bl = EXAMDATE) %>% distinct(RID, .keep_all = TRUE)
  
  # MASS SPECTROMETRY
  # Only baselines 
  upennms1 <- read.csv("~/R/EDAP-data/ADNI/UPENNMS/UPENNMSMSABETA_10Mar2026.csv") %>% 
    select(RID, VISCODE, DRAWDATE, ABETA42, ABETA40) %>% 
    rename(VISCODE2 = VISCODE) %>%
    mutate(EXAMDATE = parse_date_time(DRAWDATE, orders=c("Ymd")), data = "upennms1") %>% 
    select(-DRAWDATE) 
  
  upennms <- UPENNMSMSABETA2CRM %>% select(RID, VISCODE2, EXAMDATE, ABETA42CRM, ABETA40) %>% 
    mutate(EXAMDATE = parse_date_time(EXAMDATE, orders=c("Ymd")), data = "upennms_crm") %>% 
    rename(ABETA42 = ABETA42CRM) %>%
    rbind(upennms1) %>% mutate(A4240 = ABETA42/ABETA40,) %>% group_by(RID, VISCODE2) %>%
    slice_max(data == "upennms_crm", n=1, with_ties=FALSE) %>% ungroup()
  
  
  # ROCHE ELECSYS
  upenn_new <- ADNIMERGE2::UPENNBIOMKADNIDIAN2017 %>% 
    select(RID, EXAMDATE, VISCODE2, ABETA, AB40, A4240) %>% rename(ABETA42 = ABETA,
                                                                   ABETA40 = AB40) %>% 
    mutate(data = "upenn_new") %>%
    mutate(EXAMDATE = parse_date_time(EXAMDATE, orders=c("Ymd")))
  
  upenn3 <- ADNIMERGE2::UPENNBIOMK_ROCHE_ELECSYS %>%
    select(RID, VISCODE2, EXAMDATE, ABETA40, ABETA42) %>% drop_na(ABETA40) %>%
    mutate(A4240 = ABETA42/ABETA40, data="upenn") %>%
    mutate(EXAMDATE = parse_date_time(EXAMDATE, orders=c("Ymd")))
  
  upennabeta <- rbind(upenn3, upenn_new) %>% group_by(RID, VISCODE2) %>%
    slice_max(data == "upenn_new", n=1, with_ties=FALSE) %>% ungroup()
  
  # AB-PET
  ab_pet_ucb <- ADNIMERGE2::UCBERKELEY_AMY_6MM
  ab_pet_pib <- ADNIMERGE2::PIBPETSUVR %>% rowwise() %>% mutate(SUMMARY_SUVR = mean(ACG, FRC, PAR, PRC)) %>%
    mutate(AMYLOID_STATUS = SUMMARY_SUVR > 1.21)
  
  ab_pet <- ab_pet_ucb %>% select(RID, VISCODE2, TRACER, SUMMARY_SUVR, AMYLOID_STATUS)
  
  ab_pet <- ab_pet_pib %>% mutate(TRACER = "PIB") %>% rename(VISCODE2 = VISCODE) %>%
    select(RID, VISCODE2, TRACER, SUMMARY_SUVR, AMYLOID_STATUS) %>% rbind(ab_pet)
  
  # Diagnoses
  # CI = cognitively impaired, by diagnosis 
  # Not sure if we should only sue the cognitive tests for this one
  
  dx.df <- ADNIMERGE2::DXSUM %>% distinct(RID, VISCODE2, DIAGNOSIS)
  dx.df <- ADNIMERGE2::CDR %>% select(RID, VISCODE2, CDGLOBAL) %>% right_join(dx.df, by =c("RID", "VISCODE2"))
  dx.df <- ADNIMERGE2::MMSE %>% distinct(RID, VISCODE2, .keep_all=TRUE) %>% select(RID, VISCODE2, MMSCORE) %>%
    right_join(dx.df, by =c("RID", "VISCODE2")) %>% rowwise() %>%
    mutate(CI = any(DIAGNOSIS %in% c("MCI", "Dementia"), CDGLOBAL > 0.0, MMSCORE < 26, na.rm=TRUE)) %>% arrange(RID, VISCODE2)
  
  # AB status definition
  upennabeta <- upennabeta %>% mutate(AB_pos.re = A4240 < 0.0666,
                                      method = "Roche Elecsys") %>%
    select(RID, VISCODE2, A4240, AB_pos.re) 
  upennms <- upennms %>% mutate(AB_pos.ms = A4240 < 0.138,
                                method = "2D-UPLC",
                                VISCODE2 = as.character(VISCODE2), 
                                ABETA42 = as.numeric(ABETA42),
                                ABETA40 = as.numeric(ABETA40),
                                A4240 = as.numeric(A4240)) %>%
    select(RID, VISCODE2, A4240, AB_pos.ms)
  
  # Join data frames
  ab_csf <- full_join(upennabeta, upennms, by = c("RID", "VISCODE2"), suffix=c(".re", ".ms")) %>%
    rowwise() %>% mutate(AB_pos.csf = any(AB_pos.re, AB_pos.ms, na.rm=TRUE))
  
  ab_df <- ab_csf %>% select(RID, VISCODE2, A4240.re, A4240.ms, AB_pos.csf) %>% 
    full_join(ab_pet, by=c("RID", "VISCODE2")) %>% mutate(AB_pos.pet = as.logical(AMYLOID_STATUS)) %>%
    select(-AMYLOID_STATUS)
  
  ab_df <- dx.df %>% select(RID, VISCODE2, CI) %>% unique() %>% right_join(ab_df, by=c("RID", "VISCODE2"))
  
  ab_df <- ab_df %>% rowwise %>% mutate(AB = ifelse(CI, # Cognitively impaired?
                                                    any(AB_pos.pet, na.rm=TRUE), # yes - require AB PET +
                                                    any(AB_pos.pet, AB_pos.csf, na.rm=TRUE))) # no - any biomarker
  
  ab_df <- ab_df %>% rowwise %>% mutate(AB_any = any(AB_pos.pet, AB_pos.csf, na.rm=TRUE)) # any biomarker
  
  
  ab_df <- distinct(ab_df, RID, VISCODE2, .keep_all = TRUE)
  
  # Remove invalid subjects
  ab_df <- anti_join(ab_df, remove_rids, by=c("RID", "VISCODE2"))
  
  return(ab_df)
  
}

get_ab_pos_ids <- function() {
  ab_df <- get_ab_df()
  # Any positive Abeta biomarker at any time point
  ab_df %>% select(RID, AB_any) %>% group_by(RID) %>% 
    mutate(AB_any = if_any(AB_any)) %>% distinct() %>% filter(AB_any) %>% 
    select(RID) %>% unlist() %>% unname() -> ab_pos_rids_any
  return(ab_pos_rids_any)
}

get_tau_pet <- function() {
  # TAU PET
  ## Partial Volume Correction (PVC)
  # Tau data corrected for partial volume effects using the Geometric Transfer Matrix (GTM) approach. 
  # Using the MRI closest in time to tau scan, the GTM approach models all FreeSurfer-defined ROIs as well as regions in which offtarget
  # binding is common (e.g., choroid plexus in FTP; meninges in MK6240) to reduce
  # contamination from these regions into neighboring regions of interest. 
  require(ADNIMERGE2)
  library(epitools)
  library(lubridate)
  
  tau.pet.pvc <- ADNIMERGE2::UCBERKELEY_TAUPVC_6MM %>% 
    filter(TRACER == "FTP")
  
  bl.df <- ADNIMERGE2::REGISTRY %>% filter(VISCODE2 == "bl") %>% select(RID, EXAMDATE) %>%
    rename(EXAMDATE.bl = EXAMDATE) %>% distinct(RID, .keep_all = TRUE)
  
  tau.pet.pvc <- left_join(tau.pet.pvc, bl.df, by = "RID") %>% 
    mutate(Months = interval(EXAMDATE.bl, SCANDATE) %/% months(1),
           Years = Months/12)
  
  return(tau.pet.pvc)
}


get_diagnoses <- function() {
  dx.df <- ADNIMERGE2::DXSUM %>% distinct(RID, VISCODE2, DIAGNOSIS)
  dx.df <- ADNIMERGE2::CDR %>% select(RID, VISCODE2, CDGLOBAL) %>% right_join(dx.df, by =c("RID", "VISCODE2"))
  dx.df <- ADNIMERGE2::MMSE %>% distinct(RID, VISCODE2, .keep_all=TRUE) %>% select(RID, VISCODE2, MMSCORE) %>%
    right_join(dx.df, by =c("RID", "VISCODE2")) %>% rowwise() %>%
    mutate(CI = any(DIAGNOSIS %in% c("MCI", "Dementia"), CDGLOBAL > 0.0, MMSCORE < 26, na.rm=TRUE)) %>% arrange(RID, VISCODE2)
  
  dx.df <- dx.df %>% drop_na(DIAGNOSIS) %>%
    mutate(DX.bl = ifelse(VISCODE2 == "bl", DIAGNOSIS, NA)) %>%
    group_by(RID) %>% mutate(
      DX.bl = first(DX.bl, na.rm=TRUE, default = NA)
    ) %>% ungroup() %>%
    mutate(DIAGNOSIS = factor(DIAGNOSIS, levels = c("CN", "MCI", "Dementia"), ordered = TRUE),
           DX.bl = factor(DX.bl, levels = c("CN", "MCI", "Dementia"), ordered = TRUE)) %>%
  
  return(dx.df)
}

get_demographics <- function() {
  adni_demog <- ADNIMERGE2::PTDEMOG %>% 
    mutate(VISDATE = parse_date_time(VISDATE, orders=c("Ymd")),
           PTDOB = parse_date_time(PTDOB, orders=c("m/y")), 
           AGE = interval(PTDOB, VISDATE) %/% years(1)) %>% group_by(RID) %>%
    mutate(AGE = min(AGE, na.rm=TRUE),
           PTEDUCAT = max(PTEDUCAT, na.rm=TRUE)) %>% ungroup() %>%
    select(RID, AGE, PTGENDER, PTEDUCAT) %>% na.omit() %>% distinct() %>% arrange(RID, AGE)
  
  apoe.df <- ADNIMERGE2::APOERES %>% 
    select(RID, GENOTYPE) %>% 
    mutate(APOE4 = str_count(GENOTYPE, "4"))
  
  adni_demog <- left_join(adni_demog, apoe.df, by="RID")
  return(adni_demog)
}

