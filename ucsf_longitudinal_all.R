ucsf_longitudinal_all <- function(only_vol=FALSE, filter_n=1) {
  # Load all Longitudinal UCSF datasets
  ucsf_data1 <- ADNIMERGE2::UCSFFSL
  ucsf_data2 <- ADNIMERGE2::UCSFFSL51
  
  cols <- intersect(colnames(ucsf_data1), colnames(ucsf_data2))
  ucsf_data1 <- ucsf_data1 %>% select(all_of(cols))
  ucsf_data2 <- ucsf_data2 %>% select(all_of(cols))
  
  ucsf_data <- rbind(ucsf_data1, ucsf_data2)
  
  ucsf_data3 <- ADNIMERGE2::UCSFFSL51ALL
  cols <- intersect(colnames(ucsf_data), colnames(ucsf_data3))
  ucsf_data <- ucsf_data3 %>% select(all_of(cols)) %>% rbind(ucsf_data)
  
  ucsf_data4 <- ADNIMERGE2::UCSFFSL51Y1
  cols <- intersect(colnames(ucsf_data), colnames(ucsf_data4))
  ucsf_data <- ucsf_data4 %>% select(all_of(cols)) %>% rbind(ucsf_data)
  
  ucsf_data$RID <- as.numeric(ucsf_data$RID)
  
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
  
  ucsf_data |> filter(!is.na(Months)) -> ucsf_data
  
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