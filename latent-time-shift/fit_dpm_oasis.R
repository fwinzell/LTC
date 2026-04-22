library(dplyr)
library(progmod)
library(ggplot2)

source("~/R/LTC/utils/oasis_data_loaders.R")

plot_raw <- function(oasis_dpm) {
  gp_moca <- ggplot(data = oasis_dpm, aes(x=Years, y=oasisMOCA, group = oasisID, color = as.factor(DX.bl))) +
    geom_line() +
    labs(x="") +
    scale_color_brewer(palette = "YlOrRd", name = "Baseline diagnosis") + 
    theme_classic()
  
  gp_pet <- ggplot(data = oasis_dpm, aes(x=Years, y=PiB_SUVR, group = oasisID, color = as.factor(DX.bl))) +
    geom_point() +
    geom_line() +
    labs(x="") +
    scale_color_brewer(palette = "YlOrRd", name = "Baseline diagnosis") +
    theme_classic()
  
  gp_cdr <- ggplot(data = oasis_dpm, aes(x=Years, y=CDRSUM, group = oasisID, color = as.factor(DX.bl))) +
    geom_line() +
    labs(x="") +
    scale_color_brewer(palette = "YlOrRd", name = "Baseline diagnosis") +
    theme_classic()
  
  gp <- ggarrange(gp_moca, gp_cdr, gp_pet, ncol = 1, common.legend = TRUE, legend = "right")
  plot(gp)
}

fit_dpm <- function() {
  oasis_dpm <- get_dpm_data(use_centiloid = TRUE) 
  
  rids <- unique(oasis_dpm$OASISID)
  scale_t = FALSE
  
  # Convert to long format
  oasis_dpm <- as.data.frame(oasis_dpm) %>% 
    mutate(RID = as.numeric(gsub("OAS", "", OASISID)))
  
  cols <- c("OASISID", "RID", "Years", "DX.bl", "CN", "PMCI", "MCI", "AD")
  oasis_dpm_long <- rbind(oasis_dpm[c("CDRSUM", cols)] |> mutate(scale = "CDRSUM") |> rename(value = CDRSUM), 
                         oasis_dpm[c("MMSE", cols)] |> mutate(scale = "MMSE") |> rename(value = MMSE),
                         oasis_dpm[c("Centiloid_SUVR", cols)] |> mutate(scale = "Centiloid_SUVR") |> rename(value = Centiloid_SUVR)) %>% 
    #oasis_dpm[c("PiB_SUVR", cols)] |> mutate(scale = "PiB_SUVR") |> rename(value = PiB_SUVR)) %>%
    drop_na(value) %>% select(OASISID, RID, Years, scale, value, DX.bl, CN, PMCI, MCI, AD) %>% 
    mutate(scale = as.factor(scale), DX.bl = as.factor(DX.bl))

  
  oasis_dpm_long$t <- oasis_dpm_long$Years
  if (scale_t) {
    adni_dpm_long$t <- scale(adni_dpm_long$t)
  } 
  
  v.df <- oasis_dpm %>% filter(CN == 1 & Months == 0) %>% ungroup() %>% select(CDRSUM, MMSE, Centiloid_SUVR) %>% colMeans(na.rm=T)
  
  if (scale_t) {
    mu <- attr(adni_dpm_long$t, 'scaled:center')
    sg <- attr(adni_dpm_long$t, 'scaled:scale')
  } else {
    mu <- 0
    sg <- 1
  }
  
  
  fixed_start_coef_y <- c(l.scaleMMSE = -0.05, 
                          l.scaleCDRSUM = 0.25,
                          l.scaleCentiloid_SUVR = 10,
                          s.PMCI = (3-mu)/sg,
                          s.MCI = (6-mu)/sg,
                          s.AD = (12-mu)/sg,
                          g.scaleMMSE = 1,
                          g.scaleCDRSUM = 1,
                          g.scaleCentiloid_SUVR = 2,
                          v.df)
  
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
                       data = oasis_dpm_long,
                       fixed = list(l ~ scale + 0,
                                    s ~ PMCI + MCI + AD + 0,
                                    g ~ scale + 0,
                                    v ~ scale + 0),
                       random = list(s ~ 1,
                                     v ~ scale),
                       groups = ~ RID,
                       start = fixed_start_coef_y,
                       weights = varIdent(form = ~ 1 | scale),
                       method = "REML",
                       control = ctrl)
  
  
  oasis_dpm$fixed_shift_multi <- with(oasis_dpm,
                                     PMCI * (fixed.effects(dpm_model)["s.PMCI"]*sg)+mu +
                                       MCI * (fixed.effects(dpm_model)["s.MCI"]*sg)+mu +
                                       AD * (fixed.effects(dpm_model)["s.AD"]*sg)+mu)
  
  pred_rand <- random.effects(dpm_model)
  oasis_dpm$random_shift_multi <- (pred_rand[match(oasis_dpm$RID, rownames(pred_rand)), 's.(Intercept)']*sg)+mu
  
  
  oasis_dpm |> mutate(time_shift = fixed_shift_multi + random_shift_multi, 
                      DX.bl = factor(DX.bl, labels = c("CN", "Impaired", "MCI", "Dementia"))) |> 
    select(RID, Years, time_shift, CDRSUM, MMSE, Centiloid_SUVR, DX.bl) -> dpm_out
  
  gp_mmse <- ggplot(data = dpm_out, aes(x=Years+time_shift, y=MMSE, group = RID, color = DX.bl)) +
    geom_line() +
    labs(x="") +
    scale_color_brewer(palette = "YlOrRd", name = "Baseline diagnosis") + 
    theme_classic()
  
  gp_cdr <- ggplot(data = dpm_out, aes(x=Years+time_shift, y=CDRSUM, group = RID, color = DX.bl)) +
    geom_line() +
    labs(x="") +
    scale_color_brewer(palette = "YlOrRd", name = "Baseline diagnosis") +
    theme_classic()
  
  gp_pet <- ggplot(data = dpm_out, aes(x=Years+time_shift, y=Centiloid_SUVR, group = RID, color = DX.bl)) +
    geom_line() +
    geom_point() +
    labs(x="") +
    scale_color_brewer(palette = "YlOrRd", name = "Baseline diagnosis") +
    theme_classic()
  
  gp <- ggarrange(gp_mmse, gp_cdr, gp_pet, ncol = 1, common.legend = TRUE, legend = "right")
  plot(gp)
  
  return(dpm_out)
}
