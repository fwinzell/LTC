fit_dpm2 <- function(scale_y=FALSE, scale_t=FALSE, include_hippo=FALSE) {
  source("~/R/LTC/utils/adni_data_loaders.r")
  
  get_ab_df <- source("~/R/LTC/get_ab_df.r")$value
  ab_df <- get_ab_df()
  
  # Update: subjects with PTID: 381_S_#### should be removed, ADNIMERGE2 is not updated ...
  remove_rids <- filter(ADNIMERGE2::REGISTRY, grepl("^381_S_.*", ADNIMERGE2::REGISTRY$PTID)) %>% select(RID, VISCODE2)
  
  dates.df <- ADNIMERGE2::REGISTRY %>% select(RID, VISCODE2, EXAMDATE) %>% distinct(RID, VISCODE2, .keep_all = TRUE)
  
  dx.df <- ADNIMERGE2::DXSUM %>% select(RID, VISCODE2, DIAGNOSIS, DXCONFID) %>% drop_na(DIAGNOSIS) %>%
    mutate(DX.bl = ifelse(VISCODE2 == "bl", DIAGNOSIS, NA)) %>%
    group_by(RID) %>% mutate(
      DX.bl = if (all(is.na(DX.bl))) NA else first(na.omit(DX.bl))
    ) %>% ungroup() %>% arrange(RID, VISCODE2)
  
  missing.bl <- filter(dx.df, is.na(DX.bl) & VISCODE2 == "sc") %>% mutate(DX.bl = DIAGNOSIS) %>%
    select(RID, DX.bl) %>% distinct() %>% rename(DX.sc = DX.bl)
  
  dx.df <- dx.df %>% left_join(missing.bl, by = "RID") %>% group_by(RID) %>% 
    mutate(DX.bl = ifelse(is.na(DX.bl), DX.sc, DX.bl)) %>% select(-DX.sc)
  
  dx.bl <- dx.df %>% select(RID, DX.bl) %>% distinct()
  dx.df <- select(dx.df, -DX.bl)
  
  #dx.df <- dx.df %>% filter(VISCODE2=="bl") %>% rename(DX.bl = DIAGNOSIS) %>% select(-VISCODE2) %>%
  #  right_join(dx.df, by="RID") %>% distinct(RID, VISCODE2, .keep_all = TRUE)
  
  cdr.df <- ADNIMERGE2::CDR %>% select(RID, VISCODE2, CDGLOBAL) %>% distinct(RID, VISCODE2, .keep_all = TRUE)
  
  #dx.df <- ADNIMERGE2::MMSE %>% select(RID, VISCODE2, MMSCORE) %>% distinct(RID, VISCODE2, .keep_all = TRUE) %>%
  #  right_join(dx.df, by =c("RID", "VISCODE2")) %>%
  #  rowwise() %>%
  #  mutate(CI = any(DIAGNOSIS %in% c("MCI", "Dementia"), CDGLOBAL > 0.0, MMSCORE < 26, na.rm=TRUE)) %>% arrange(RID, VISCODE2)
  
  adas.df <- ADNIMERGE2::ADAS %>% select(RID, VISCODE2, TOTAL13) %>% distinct(RID, VISCODE2, .keep_all = TRUE) %>% 
    rename(ADAS13 = TOTAL13) %>% drop_na(ADAS13)
  mmse.df <- ADNIMERGE2::MMSE %>% select(RID, VISCODE2, MMSCORE) %>% distinct(RID, VISCODE2, .keep_all = TRUE) %>% 
    rename(MMSE = MMSCORE) %>% drop_na(MMSE)
  amy.df <- ADNIMERGE2::UCBERKELEY_AMY_6MM %>% filter(TRACER == "FBP") %>% 
    select(RID, VISCODE2, SUMMARY_SUVR) %>% rename(FBP_SUVR = SUMMARY_SUVR) %>% distinct(RID, VISCODE2, .keep_all = TRUE)
  
  adni_dpm <- full_join(adas.df, mmse.df, by=c("RID", "VISCODE2")) %>% 
    full_join(amy.df, by=c("RID", "VISCODE2")) %>%
    left_join(dx.df, by = c("RID", "VISCODE2")) %>% 
    left_join(dates.df, by = c("RID", "VISCODE2")) %>%
    left_join(cdr.df, by = c("RID", "VISCODE2")) %>%
    left_join(dx.bl, by = "RID") %>%
    distinct(RID, VISCODE2, .keep_all = TRUE) %>% filter(! VISCODE2 %in% c("f", "uns1"))
  
  # Time by EXAMDATES
  adni_dpm <- ADNIMERGE2::REGISTRY %>% group_by(RID) %>% arrange(EXAMDATE) %>% 
    summarise(EXAMDATE.bl = EXAMDATE[1]) %>% right_join(adni_dpm, by="RID") %>% drop_na(EXAMDATE.bl)
  
  adni_dpm <- mutate(adni_dpm, Months = interval(EXAMDATE.bl, EXAMDATE) %/% months(1)) %>%
    mutate(VisM = as.numeric(ifelse(VISCODE2 %in% c("bl", "sc"), "0", gsub("m", "", VISCODE2)))) %>% 
    mutate(Months = ifelse(is.na(Months), VisM, Months)) %>% select(-VisM) %>% 
    mutate(Years = Months/12) 
  
  if (include_hippo) {
    hippo.df <- ucsf_longitudinal_all(only_vol=TRUE) %>%
      select(RID, Years, ST29SV, ST88SV, DX.bl) %>%
      mutate(Hippocampus = ST29SV+ST88SV) %>% select(RID, Years, Hippocampus) %>%
      distinct(RID, Years, .keep_all=TRUE)
    
    #hippo.df$Hippocampus <- round(hippo.df$Hippocampus*1000, 1)
    #ggplot(data = hippo.df, aes(x=Years, y=Hippocampus, group = RID, color = DX.bl)) +
    #  geom_line() +
    #  labs(x="") +
    #  scale_color_brewer(palette = "YlOrRd", name = "Baseline diagnosis") +
    #  theme_classic()
    
    adni_dpm <- left_join(adni_dpm, hippo.df, by= c("RID", "Years"))
  }
  
  adni_dpm <- anti_join(adni_dpm, remove_rids, by=c("RID", "VISCODE2"))
  
  adni_dpm %>% distinct(RID, DX.bl) %>% select(DX.bl) %>% table(useNA = "ifany")
  
  # Impute missing baseline diagnosis by CDGLOBAL
  adni_dpm <- adni_dpm %>% 
    mutate(NAs = if_all(c(ADAS13, MMSE, FBP_SUVR, DIAGNOSIS), is.na)) %>%
    filter(!NAs) %>%
    mutate(DX.bl = ifelse(is.na(DX.bl), 
                          ifelse(CDGLOBAL >= 1.0, "Dementia", 
                                 ifelse(CDGLOBAL >= 0.5, "MCI", "CN")), as.character(DX.bl))) %>%
    drop_na(DX.bl) %>%
    mutate(CN = ifelse(DX.bl == "CN", 1, 0),
           MCI = ifelse(DX.bl == "MCI", 1, 0),
           AD = ifelse(DX.bl == "Dementia", 1, 0)) %>%
    rowwise() %>%
    mutate(DXCI = DIAGNOSIS %in% c("MCI", "Dementia") & DXCONFID %in% c("Highly Confident", NA)) %>%
    mutate(CI = any(DXCI, CDGLOBAL > 0.0, MMSE < 26, na.rm=TRUE)) %>% arrange(RID, EXAMDATE)
  
  
  # NOTE: AB defined as AB- = 0, AB+ = 1
  adni_dpm <- ab_df %>% select(RID, VISCODE2, AB, AB_pos.pet) %>% right_join(adni_dpm, by=c("RID", "VISCODE2")) %>%
    group_by(RID) %>% mutate(valid_ab = any(!is.na(AB))) %>% filter(valid_ab) %>% select(-valid_ab) 
  
  # Remove any subjects with negative AB-PET status and cognitive impairment
  adni_dpm <- adni_dpm %>%  mutate(invalid = !AB_pos.pet & CI) %>% group_by(RID) %>% 
    mutate(invalid = any(invalid)) %>% ungroup() %>%
    filter(invalid == FALSE | is.na(invalid))
    
  #removals <- filter(adni_dpm, invalid)
  
  # Impute any missing AB status by carrying negatives backwards, and positive forwards
  # assuming monoticity
  missing_ab <- adni_dpm %>% filter(is.na(AB)) %>% select(RID) %>% unique() %>% unlist()
  table(adni_dpm$AB, useNA = "ifany")
  
  for (rid in missing_ab) {
    subj <- adni_dpm %>% filter(RID == rid) %>% arrange(EXAMDATE)
    idxs <- which(!is.na(subj$AB))
    if (length(idxs) > 0) {
      last_i = 1
      next_i = c(idxs, nrow(subj)+1)
      count = 1
      for (i in idxs) {
        count = count + 1
        if (subj$AB[i] == 0) {
          # Carry backwards
          subj[last_i:i, "AB"] <- subj$AB[i]
          last_i <- i
        } else {
          # Carry forwards
          subj[i:(next_i[count]-1), "AB"] <- subj$AB[i]
          last_i <- next_i[count]
        }
      }
      adni_dpm <- rows_update(adni_dpm, subj, by = c("RID", "VISCODE2"))
    }
  }
  

  # Remove any visits with CI and negative AB CSF status
  adni_dpm <- adni_dpm %>%  mutate(invalid = !AB & CI) %>% filter(invalid == FALSE | is.na(invalid)) %>% 
    group_by(RID) %>% mutate(hasbl = any(VISCODE2 %in% c("bl", "sc"))) %>% ungroup() %>% filter(hasbl) 
  
  ab.bl <- adni_dpm %>% select(RID, VISCODE2, AB) %>% drop_na(AB) %>% filter(VISCODE2 %in% c("bl", "sc")) %>% 
    mutate(negAB.bl = 1-AB) %>% select(RID, negAB.bl) %>% distinct()
  
  # Indicator for negative AB baseline status
  adni_dpm <- left_join(adni_dpm, ab.bl, by = "RID") %>%
    drop_na(negAB.bl)
    #mutate(negAB.bl = ifelse(is.na(negAB.bl), 0, negAB.bl))
  
  rids <- unique(adni_dpm$RID)
  
  if (scale_y) {
    adni_dpm <- mutate(adni_dpm, 
                       ADAS13 = ADAS13/85,
                       MMSE = MMSE/30,
                       FBP_SUVR = FBP_SUVR/max(FBP_SUVR, na.rm=TRUE))
  }
  
  # Convert to long format
  cols <- c("RID", "Months", "DX.bl", "CN", "MCI", "AD", "negAB.bl")
  adni_dpm_long <- rbind(adni_dpm[c("ADAS13", cols)] |> mutate(scale = "ADAS13") |> rename(value = ADAS13), 
                         adni_dpm[c("MMSE", cols)] |> mutate(scale = "MMSE") |> rename(value = MMSE),
                         adni_dpm[c("FBP_SUVR", cols)] |> mutate(scale = "FBP_SUVR") |> rename(value = FBP_SUVR))
  if (include_hippo) {
    adni_dpm_long <- rbind(adni_dpm_long, 
                           adni_dpm[c("Hippocampus", cols)] |> mutate(scale = "Hippocampus") |> rename(value = Hippocampus))
  }
  adni_dpm_long <- drop_na(adni_dpm_long, value)
  
  adni_dpm_long <- adni_dpm_long %>% mutate(
    scale = as.factor(scale),
    DX.bl = as.factor(DX.bl),
    RID = as.numeric(RID),
    Years = Months/12
  )
    
  adni_dpm_long$t <- adni_dpm_long$Years
  if (scale_t) {
    adni_dpm_long$t <- scale(adni_dpm_long$t)
  } 

  cols <- c("ADAS13", "MMSE", "FBP_SUVR")
  if (include_hippo) { 
    cols <- c(cols, "Hippocampus")
    }
  v.df <- adni_dpm %>% filter(CN == 1 & Months == 0.0) %>% ungroup() %>% select(all_of(cols)) %>% colMeans(na.rm=T)
  
  if (scale_t) {
    mu <- attr(adni_dpm_long$t, 'scaled:center')
    sg <- attr(adni_dpm_long$t, 'scaled:scale')
  } else {
    mu <- 0
    sg <- 1
  }
  
  
  fixed_start_coef_y <- c(l.scaleADAS13 = 0.5, #0.05,
                          l.scaleMMSE = -0.05, #-0.04,
                          l.scaleFBP_SUVR = 0.5, #0.05,
                          #l.scaleADAS13 = 0.05,
                          #l.scaleMMSE = -0.04,
                          #l.scaleFBP_SUVR = 0.05,
                          s.MCI = (6-mu)/sg,
                          s.AD = (12-mu)/sg,
                          s.negAB.bl = (-1-mu)/sg, 
                          g.scaleADAS13 = 1.5,
                          g.scaleMMSE = 1,
                          g.scaleFBP_SUVR = 2,
                          v.scale=v.df) 
                          #v.scaleADAS13 = v.df["ADAS13"],
                          #v.scaleMMSE = v.df["MMSE"],
                          #v.scaleFBP_SUVR = v.df["FBP_SUVR"])
  
  if (include_hippo) {
    fixed_start_coef_y <- c(fixed_start_coef_y,
                            l.scaleHippocampus = -0.01,
                            g.scaleHippocampus = 1)
  }
  
  ctrl <- nlmeControl(maxIter = 200, # 50, Can be increased
                      pnlsMaxIter = 20, 
                      msMaxIter = 500, # 50-500
                      minScale = 0.001,
                      tolerance = 1e-4,
                      niterEM = 100, # 25-100
                      pnlsTol = 0.001,
                      msTol = 1e-5,
                      msVerbose = FALSE,
                      apVar = TRUE,
                      minAbsParApVar = 0.05,
                      natural = TRUE)
  
   
  dpm_model <- progmod(value ~ progmod::exp_model(t, l, s, g, v),
                       data = adni_dpm_long,
                       fixed = list(l ~ scale + 0,
                                    s ~ MCI + AD + negAB.bl + 0,
                                    g ~ scale + 0,
                                    v ~ scale + 0),
                       random = list(s ~ 1,
                                     v ~ scale),
                       groups = ~ RID,
                       start = fixed_start_coef_y,
                       weights = varIdent(form = ~ 1 | scale),
                       method = "REML",
                       control = ctrl)
  
  
  adni_dpm$fixed_shift_multi <- with(adni_dpm,
                                     MCI * (fixed.effects(dpm_model)["s.MCI"]*sg)+mu +
                                       AD * (fixed.effects(dpm_model)["s.AD"]*sg)+mu +
                                       negAB.bl * (fixed.effects(dpm_model)["s.negAB.bl"])*sg+mu)
  
  pred_rand <- random.effects(dpm_model)
  adni_dpm$random_shift_multi <- (pred_rand[match(adni_dpm$RID, rownames(pred_rand)), 's.(Intercept)']*sg)+mu
  
  
  adni_dpm |> mutate(time_shift = fixed_shift_multi + random_shift_multi,
                     DX.bl = factor(DX.bl, levels = c("CN", "MCI", "Dementia"))) |> 
    select(RID, Years, time_shift, ADAS13, MMSE, FBP_SUVR, DX.bl) -> dpm_out
  
  gp_adas <- ggplot(data = dpm_out, aes(x=Years+time_shift, y=ADAS13, group = RID, color = DX.bl)) +
    geom_line() +
    labs(x="") +
    scale_color_brewer(palette = "YlOrRd", name = "Baseline diagnosis") + 
    theme_classic()
  
  gp_pet <- ggplot(data = dpm_out, aes(x=Years+time_shift, y=FBP_SUVR, group = RID, color = DX.bl)) +
    geom_line() +
    labs(x="") +
    scale_color_brewer(palette = "YlOrRd", name = "Baseline diagnosis") +
    theme_classic()
  
  gp_mmse <- ggplot(data = dpm_out, aes(x=Years+time_shift, y=MMSE, group = RID, color = DX.bl)) +
    geom_line() +
    labs(x="") +
    scale_color_brewer(palette = "YlOrRd", name = "Baseline diagnosis") +
    theme_classic()
  
  gp <- ggarrange(gp_adas, gp_mmse, gp_pet, ncol = 1, common.legend = TRUE, legend = "right")
  plot(gp)
  
  return(dpm_out)
}
