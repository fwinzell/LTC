library(data.table)
library(lubridate)
library(ggplot2)
library(tidyr)

plot_histograms <- function(df, cols) {
  # Usage:
  #plot_histograms(data.frame(nacc_mri_scan), c("RIGHT_HIPPOCAMPUS", "RIGHT_HIPPO"))
  df_long <- pivot_longer(df[cols], cols = everything(), 
                          names_to = "variable", values_to = "value")
  
  ggplot(df_long, aes(x = value, fill = variable)) +
    geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
    theme_minimal()
}


load_nacc_mri <- function(ab_positives=TRUE) {
  nacc_ft <- fread("~/R/EDAP-data/NACC/investigator_ftldlbd_nacc72.csv")
  
  nacc_mri_scan <- fread("~/R/EDAP-data/NACC/investigator_scan_mri_nacc72/investigator_scan_mrisbm_nacc72.csv")
  
  nacc_bl <- nacc_ft %>% select(NACCID, VISITDAY, VISITMO, VISITYR, NACCFDYS) %>% filter(NACCFDYS == 0) %>% 
    distinct() %>% mutate(DATE.bl = parse_date_time(paste(VISITYR, VISITMO, VISITDAY, sep="-"), orders=c("Y-m-d")))
  
  nacc_mri_scan <- nacc_mri_scan %>% mutate(DATE = parse_date_time(SCANDT, orders=c("Y-m-d")))
  nacc_mri_scan <- nacc_bl %>% select(NACCID, DATE.bl) %>% right_join(nacc_mri_scan, by="NACCID") %>%
    mutate(Months = interval(DATE.bl, DATE) %/% months(1),
           Years = Months/12)
  
  mri_cols <- setdiff(c(colnames(nacc_mri_scan)[grep("^(LEFT|RIGHT|CC)_", colnames(nacc_mri_scan))], 
                colnames(nacc_mri_scan)[grep("^[LR]H_[^_]*_GVOL", colnames(nacc_mri_scan))],
                "BRAIN_STEM"), c("LEFT_HIPPO", "RIGHT_HIPPO"))
  nacc_mri_scan <- select(nacc_mri_scan, NACCID, DATE, DATE.bl, Months, Years, ESTIMATEDTOTALINTRACRANIALVOL, all_of(mri_cols)) %>%
    #mutate(LEFT_HIPPO = LEFT_HIPPO*1000,
    #       RIGHT_HIPPO = RIGHT_HIPPO*1000) %>%
    mutate_at(mri_cols, ~ . /ESTIMATEDTOTALINTRACRANIALVOL)
  
  nacc_bl <- nacc_ft %>% filter(NACCFDYS == 0) %>% select(NACCID, NACCUDSD) %>%
    mutate(NACCDXBL = factor(NACCUDSD, labels = c("CN", "Impaired", "MCI", "Dementia"))) %>% select(-NACCUDSD)
  nacc_mri_scan <- left_join(nacc_mri_scan, nacc_bl, by="NACCID")           
  
  if (ab_positives) {
    amy_pos_id <- nacc_ft %>% filter(AMYLPET == 1 | AMYLCSF == 1) %>% select(NACCID) %>% unique() %>% unlist()
    nacc_mri_scan <- filter(nacc_mri_scan, NACCID %in% amy_pos_id) 
  }
  
  return(nacc_mri_scan)
}

load_nacc_mri_mp <- function() {
  # Mixed protocols, do not use?
  nacc_ft <- fread("~/R/EDAP-data/NACC/investigator_ftldlbd_nacc72.csv")
  
  nacc_bl <- nacc_ft %>% select(NACCID, VISITDAY, VISITMO, VISITYR, NACCFDYS) %>% filter(NACCFDYS == 0) %>% 
    distinct() %>% mutate(DATE.bl = parse_date_time(paste(VISITYR, VISITMO, VISITDAY, sep="-"), orders=c("Y-m-d")))
  
  nacc_mri_mp <- fread("~/R/EDAP-data/NACC/investigator_mri_nacc72.csv") %>%
    filter(!(NACCICV %in% c(8888.888, 9999.999)))
  
  nacc_mri_mp <- mutate(nacc_mri_mp, DATE = parse_date_time(paste(MRIYR, MRIMO, MRIDY, sep="-"), orders=c("Y-m-d")))
  nacc_mri_mp <- nacc_bl %>% select(NACCID, DATE.bl) %>% right_join(nacc_mri_mp, by="NACCID") %>%
    mutate(Months = interval(DATE.bl, DATE) %/% months(1),
           Years = Months/12)
}

plot_nacc_pet <- function() {
  nacc_pet_scan <- fread("~/R/EDAP-data/NACC/investigator_scan_pet_nacc72/investigator_scan_amyloidpetgaain_nacc72.csv") %>%
    select(NACCID, NACCADC, LONIUID, SCANDATE, TRACER, AMYLOID_STATUS, GAAIN_SUMMARY_SUVR, GAAIN_WHOLECEREBELLUM_SUVR, 
           GAAIN_COMPOSITE_REF_SUVR, GAAIN_CEREBELLUM_CORTEX)
  nacc_pet_scan_mp <- fread("~/R/EDAP-data/NACC/investigator_scan_pet_nacc72/investigator_scan_mp_amyloidpetgaain_nacc72.csv") %>%
    select(NACCID, NACCADC, LONIUID, SCANDATE, TRACER, AMYLOID_STATUS, GAAIN_SUMMARY_SUVR, GAAIN_WHOLECEREBELLUM_SUVR, 
           GAAIN_COMPOSITE_REF_SUVR, GAAIN_CEREBELLUM_CORTEX)
  
  nacc_pet <- rbind(nacc_pet_scan, nacc_pet_scan_mp) %>% arrange(NACCID, SCANDATE) %>%
    mutate(DATE = parse_date_time(SCANDATE, "ymd")) %>% select(-SCANDATE) %>% 
    drop_na(GAAIN_SUMMARY_SUVR) %>% mutate(TRACER = factor(TRACER, 
                                                           labels = c("PiB", "FBP", "FBB", "NAV")))
  
  ggplot(nacc_pet, aes(x = TRACER, y = GAAIN_SUMMARY_SUVR, group = TRACER, color=TRACER)) +
    geom_violin(trim = FALSE, alpha = 0.4) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    geom_jitter(width = 0.12, alpha = 0.35, size = 1) +
    labs(
      x = "Amyloid PET tracer",
      y = "GAAIN / Centiloid value",
      title = "Distribution of amyloid PET values by tracer"
    ) +
    theme_minimal()
  
  ggplot(nacc_pet, aes(x = GAAIN_SUMMARY_SUVR, fill = TRACER)) +
    geom_density(alpha = 0.35) +
    labs(
      x = "GAAIN / Centiloid value",
      y = "Density",
      title = "Distribution of amyloid PET values by tracer"
    ) +
    theme_minimal()
}

get_dpm_data <- function () {
  nacc_ft <- fread("~/R/EDAP-data/NACC/investigator_ftldlbd_nacc72.csv")
  
  nacc_dx <- select(nacc_ft, NACCID, VISITDAY, VISITMO, VISITYR, NACCFDYS, NORMCOG, DEMENTED, NACCUDSD, AMNDEM, NACCLBDS, 
                    AMYLPET, AMYLCSF) %>% arrange(NACCID, NACCFDYS)
  # Diagnosis in NACCUDSD
  # 1 - Normal
  # 2 - Impaired not-MCI
  # 3 - MCI
  # 4 - Dementia
  
  # Abnormal AB status
  # 0 - No
  # 1 - Yes
  # 8 - unknown
  # -4 - not applicable
  
  # There are more variables that might be of interest here later
  
  # Filter out subjects with CI being AB-PET negative
  invalid_ids <- nacc_dx %>% select(NACCID, VISITDAY, VISITMO, VISITYR, NACCUDSD, AMYLPET) %>%
    filter(NACCUDSD > 2 & AMYLPET == 0)
  
  nacc_bl <- nacc_dx %>% filter(NACCFDYS == 0) %>% mutate(NACCDXBL = NACCUDSD, 
                                                          AMYLPET = ifelse(AMYLPET %in% c(8, -4), NA, AMYLPET),
                                                          AMYLCSF = ifelse(AMYLCSF %in% c(8, -4), NA, AMYLCSF)) %>%
    filter(!(is.na(AMYLPET) & is.na(AMYLCSF))) %>% rowwise() %>% mutate(AMYLBL = pmin(sum(AMYLPET, AMYLCSF, na.rm=TRUE), 1)) %>%
    mutate(DATE.bl = parse_date_time(paste0(VISITYR, "-", VISITMO, "-", VISITDAY), "ymd")) %>%
    select(NACCID, DATE.bl, NACCDXBL, AMYLBL)
  
  #nacc_dates <- nacc_dx %>% select(NACCID, VISITDAY, VISITMO, VISITYR, NACCFDYS) %>%
  #  mutate(DATE = parse_date_time(paste0(VISITYR, "-", VISITMO, "-", VISITDAY), "ymd")) %>% 
  #  select(-c(VISITDAY, VISITMO, VISITYR))
  
  nacc_cog <- nacc_ft %>% select(NACCID, VISITDAY, VISITMO, VISITYR, CDRSUM, CDRGLOB, NACCMMSE, NACCMOCA) %>%
    mutate(NACCMMSE = ifelse(NACCMMSE > 30 | NACCMMSE < 0, NA, NACCMMSE), 
           NACCMOCA = ifelse(NACCMOCA > 30 | NACCMOCA < 0, NA, NACCMOCA),
           CDRSUM = ifelse(CDRSUM > 18, NA, CDRSUM),
           DATE = parse_date_time(paste0(VISITYR, "-", VISITMO, "-", VISITDAY), "ymd"))
  
  # PET
  nacc_pet_scan <- fread("~/R/EDAP-data/NACC/investigator_scan_pet_nacc72/investigator_scan_amyloidpetgaain_nacc72.csv") %>%
    select(NACCID, NACCADC, LONIUID, SCANDATE, TRACER, AMYLOID_STATUS, GAAIN_SUMMARY_SUVR, GAAIN_WHOLECEREBELLUM_SUVR, 
           GAAIN_COMPOSITE_REF_SUVR, GAAIN_CEREBELLUM_CORTEX)
  nacc_pet_scan_mp <- fread("~/R/EDAP-data/NACC/investigator_scan_pet_nacc72/investigator_scan_mp_amyloidpetgaain_nacc72.csv") %>%
    select(NACCID, NACCADC, LONIUID, SCANDATE, TRACER, AMYLOID_STATUS, GAAIN_SUMMARY_SUVR, GAAIN_WHOLECEREBELLUM_SUVR, 
           GAAIN_COMPOSITE_REF_SUVR, GAAIN_CEREBELLUM_CORTEX)
  
  # TRACERS
  # 2. PiB - 2332
  # 3. FBP - 1188
  # 4. FBB - 1538
  # 5. NAV - 170
  # Use GAAIN values to include all?
  
  nacc_pet <- rbind(nacc_pet_scan, nacc_pet_scan_mp) %>% arrange(NACCID, SCANDATE) %>%
    select(NACCID, SCANDATE, GAAIN_SUMMARY_SUVR) %>% rename(GAAIN_SUVR = GAAIN_SUMMARY_SUVR) %>%
    mutate(DATE = parse_date_time(SCANDATE, "ymd")) %>% select(-SCANDATE)
  
  nacc_dpm <- select(nacc_cog, NACCID, DATE, CDRSUM, NACCMMSE, NACCMOCA) %>%
    full_join(nacc_pet, by=c("NACCID", "DATE")) %>%
    left_join(nacc_bl, by="NACCID") %>% drop_na(AMYLBL) %>%
    filter(!(NACCID %in% invalid_ids$NACCID)) %>%
    mutate(CN = ifelse(NACCDXBL == 1, 1, 0),
           PMCI = ifelse(NACCDXBL == 2, 1, 0),
           MCI = ifelse(NACCDXBL == 3, 1, 0),
           AD = ifelse(NACCDXBL == 4, 1, 0))
  
  nacc_dpm <- mutate(nacc_dpm, Months = interval(DATE.bl, DATE) %/% months(1)) %>%
    mutate(Years = Months/12, 
           negAB.bl = AMYLBL == 0)
  
  return(nacc_dpm)
  
}

