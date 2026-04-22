get_ab_df <- function() {
  # Make sure adnimerge is loaded
  #require(ADNIMERGE)
  require(ADNIMERGE2)
  library(lubridate)
  
  # Update: subjects with PTID: 381_S_#### should be removed, ADNIMERGE2 is not updated ...
  remove_rids <- filter(ADNIMERGE2::REGISTRY, grepl("^381_S_.*", ADNIMERGE2::REGISTRY$PTID)) %>% select(RID, VISCODE2)
  
  # Baseline examdates
  #bl.df <- adnimerge %>% select(RID, EXAMDATE.bl) %>% distinct() %>%
  #  mutate(RID = as.numeric(RID),
  #         EXAMDATE.bl = parse_date_time(EXAMDATE.bl, orders=c("Ymd")))
  
  bl.df <- ADNIMERGE2::REGISTRY %>% filter(VISCODE2 == "bl") %>% select(RID, EXAMDATE) %>%
    rename(EXAMDATE.bl = EXAMDATE) %>% distinct(RID, .keep_all = TRUE)
  
  # MASS SPECTROMY 
  #upennms2 <- upennmsmsabeta2 %>% select(RID, VISCODE, EXAMDATE, ABETA42, ABETA40) %>%
  #  mutate(EXAMDATE = parse_date_time(EXAMDATE, orders=c("Ymd")),
  #         data = "upennms2") %>% 
  #  left_join(bl.df, by = "RID") %>% 
  #  mutate(M = interval(EXAMDATE.bl, EXAMDATE) %/% months(1)) %>%
  #  mutate(M = round(M/6)*6) %>% mutate(VISCODE = ifelse(M == 0, "bl", paste0("m", M))) %>%
  #  select(-c(EXAMDATE.bl, M))
  
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
  #ab_pet <- adnimerge %>% select(RID, VISCODE, AV45, PIB, FBB) %>% mutate(RID = as.numeric(RID)) %>%
  #  filter(!if_all(c(AV45, PIB, FBB), is.na))
  
  ab_pet_ucb <- ADNIMERGE2::UCBERKELEY_AMY_6MM
  ab_pet_pib <- ADNIMERGE2::PIBPETSUVR %>% rowwise() %>% mutate(SUMMARY_SUVR = mean(ACG, FRC, PAR, PRC)) %>%
    mutate(AMYLOID_STATUS = SUMMARY_SUVR > 1.21)
  
  ab_pet <- ab_pet_ucb %>% select(RID, VISCODE2, TRACER, SUMMARY_SUVR, AMYLOID_STATUS)
  
  ab_pet <- ab_pet_pib %>% mutate(TRACER = "PIB") %>% rename(VISCODE2 = VISCODE) %>%
    select(RID, VISCODE2, TRACER, SUMMARY_SUVR, AMYLOID_STATUS) %>% rbind(ab_pet)
  
  # Diagnoses
  # CI = cognitively impaired, by diagnosis 
  # Not sure if we should only sue the cognitive tests for this one
  #dx.df <- adnimerge %>% select(RID, VISCODE, DX) %>% mutate(RID = as.numeric(RID)) %>%
  #  group_by(RID) %>% mutate(CI = any(DX %in% c("MCI", "Dementia"))) 
  
  dx.df <- ADNIMERGE2::DXSUM %>% distinct(RID, VISCODE2, DIAGNOSIS)
  dx.df <- ADNIMERGE2::CDR %>% select(RID, VISCODE2, CDGLOBAL) %>% right_join(dx.df, by =c("RID", "VISCODE2"))
  dx.df <- ADNIMERGE2::MMSE %>% distinct(RID, VISCODE2, .keep_all=TRUE) %>% select(RID, VISCODE2, MMSCORE) %>%
    right_join(dx.df, by =c("RID", "VISCODE2")) %>% rowwise() %>%
    mutate(CI = any(DIAGNOSIS %in% c("MCI", "Dementia"), CDGLOBAL > 0.0, MMSCORE < 26, na.rm=TRUE)) %>% arrange(RID, VISCODE2)
    
  # AB_pos
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
  
  ab_csf <- full_join(upennabeta, upennms, by = c("RID", "VISCODE2"), suffix=c(".re", ".ms")) %>%
    rowwise() %>% mutate(AB_pos.csf = any(AB_pos.re, AB_pos.ms, na.rm=TRUE))
  
  #ab_pet <- ab_pet %>% rowwise() %>% mutate(AB_pos.pet = any(AV45 > 1.11 | PIB > 1.22 | FBB > 1.08, na.rm=TRUE),
  #                                          VISCODE = as.character(VISCODE))
  
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