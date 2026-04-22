library(dplyr)
library(progmod)
library(ggplot2)

source("~/R/LTC/utils/nacc_data_loaders.R")

plot_raw <- function(nacc_dpm) {
  gp_moca <- ggplot(data = nacc_dpm, aes(x=Years, y=NACCMOCA, group = NACCID, color = as.factor(NACCDXBL))) +
    geom_line() +
    labs(x="") +
    scale_color_brewer(palette = "YlOrRd", name = "Baseline diagnosis") + 
    theme_classic()
  
  gp_pet <- ggplot(data = nacc_dpm, aes(x=Years, y=PiB_SUVR, group = NACCID, color = as.factor(NACCDXBL))) +
    geom_point() +
    geom_line() +
    labs(x="") +
    scale_color_brewer(palette = "YlOrRd", name = "Baseline diagnosis") +
    theme_classic()
  
  gp_cdr <- ggplot(data = nacc_dpm, aes(x=Years, y=CDRSUM, group = NACCID, color = as.factor(NACCDXBL))) +
    geom_line() +
    labs(x="") +
    scale_color_brewer(palette = "YlOrRd", name = "Baseline diagnosis") +
    theme_classic()
  
  gp <- ggarrange(gp_moca, gp_cdr, gp_pet, ncol = 1, common.legend = TRUE, legend = "right")
  plot(gp)
}

fit_dpm <- function() {
  nacc_dpm <- get_dpm_data() 
  
  
  rids <- unique(nacc_dpm$NACCID)
  scale_t = FALSE
  
  # Convert to long format
  nacc_dpm <- as.data.frame(nacc_dpm) %>% 
    mutate(RID = as.numeric(gsub("NACC", "", NACCID)))
  
  cols <- c("NACCID", "Years", "NACCDXBL", "CN", "PMCI", "MCI", "AD", "negAB.bl")
  nacc_dpm_long <- rbind(nacc_dpm[c("CDRSUM", cols)] |> mutate(scale = "CDRSUM") |> rename(value = CDRSUM), 
                         #nacc_dpm[c("NACCMMSE", cols)] |> mutate(scale = "NACCMMSE") |> rename(value = NACCMMSE),
                         nacc_dpm[c("NACCMOCA", cols)] |> mutate(scale = "NACCMOCA") |> rename(value = NACCMOCA)) %>% 
    #nacc_dpm[c("PiB_SUVR", cols)] |> mutate(scale = "PiB_SUVR") |> rename(value = PiB_SUVR)) %>%
    drop_na(value) %>% select(NACCID, Years, scale, value, NACCDXBL, CN, PMCI, MCI, AD, negAB.bl)
  
  nacc_dpm_long <- nacc_dpm_long %>% mutate(
    scale = as.factor(scale),
    NACCDXBL = as.factor(NACCDXBL),
    negAB.bl = as.numeric(negAB.bl),
    RID = as.numeric(gsub("NACC", "", NACCID))
  )
  
  nacc_dpm_long$t <- nacc_dpm_long$Years
  if (scale_t) {
    adni_dpm_long$t <- scale(adni_dpm_long$t)
  } 
  
  
  v.df <- nacc_dpm %>% filter(CN == 1 & Months == 0) %>% ungroup() %>% select(CDRSUM, NACCMOCA) %>% colMeans(na.rm=T)
  
  if (scale_t) {
    mu <- attr(adni_dpm_long$t, 'scaled:center')
    sg <- attr(adni_dpm_long$t, 'scaled:scale')
  } else {
    mu <- 0
    sg <- 1
  }
  
  
  
  fixed_start_coef_y <- c(l.scaleNACCMOCA = -0.05, 
                          #l.scaleNACCMMSE = -0.05, 
                          #l.scalePiB_SUVR = 0.5, 
                          l.scaleCDRSUM = 0.25,
                          s.PMCI = (3-mu)/sg,
                          s.MCI = (6-mu)/sg,
                          s.AD = (12-mu)/sg,
                          s.negAB.bl = (-1-mu)/sg, 
                          g.scaleNACCMOCA = 1,
                          #g.scaleNACCMMSE = 1,
                          g.scaleCDRSUM = 1,
                          #g.scalePiB_SUVR = 2,
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
                       data = nacc_dpm_long,
                       fixed = list(l ~ scale + 0,
                                    s ~ PMCI + MCI + AD + negAB.bl + 0,
                                    g ~ scale + 0,
                                    v ~ scale + 0),
                       random = list(s ~ 1,
                                     v ~ scale),
                       groups = ~ RID,
                       start = fixed_start_coef_y,
                       weights = varIdent(form = ~ 1 | scale),
                       method = "REML",
                       control = ctrl)
  
  
  nacc_dpm$fixed_shift_multi <- with(nacc_dpm,
                                     PMCI * (fixed.effects(dpm_model)["s.PMCI"]*sg)+mu +
                                       MCI * (fixed.effects(dpm_model)["s.MCI"]*sg)+mu +
                                       AD * (fixed.effects(dpm_model)["s.AD"]*sg)+mu +
                                       negAB.bl * (fixed.effects(dpm_model)["s.negAB.bl"])*sg+mu)
  
  pred_rand <- random.effects(dpm_model)
  nacc_dpm$random_shift_multi <- (pred_rand[match(nacc_dpm$RID, rownames(pred_rand)), 's.(Intercept)']*sg)+mu
  
  
  nacc_dpm |> mutate(time_shift = fixed_shift_multi + random_shift_multi,
                     NACCDXBL = factor(NACCDXBL, labels = c("CN", "Impaired", "MCI", "Dementia"))) |> 
    select(RID, Years, time_shift, CDRSUM, NACCMOCA, NACCDXBL) -> dpm_out
  
  gp_moca <- ggplot(data = dpm_out, aes(x=Years+time_shift, y=NACCMOCA, group = RID, color = NACCDXBL)) +
    geom_line() +
    labs(x="") +
    scale_color_brewer(palette = "YlOrRd", name = "Baseline diagnosis") + 
    theme_classic()
  
  gp_cdr <- ggplot(data = dpm_out, aes(x=Years+time_shift, y=CDRSUM, group = RID, color = NACCDXBL)) +
    geom_line() +
    labs(x="") +
    scale_color_brewer(palette = "YlOrRd", name = "Baseline diagnosis") +
    theme_classic()
  
  gp <- ggarrange(gp_moca, gp_cdr, ncol = 1, common.legend = TRUE, legend = "right")
  plot(gp)
  
  return(dpm_out)
}
