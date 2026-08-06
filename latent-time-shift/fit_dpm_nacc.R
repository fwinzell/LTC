library(dplyr)
library(progmod)
library(ggplot2)
library(tidyr)
library(ggpubr)

source("~/R/LTC/utils/nacc_data_loaders.R")

plot_raw <- function(nacc_dpm) {
  gp_moca <- ggplot(data = nacc_dpm, aes(x=Years, y=NACCMOCA, group = NACCID, color = as.factor(NACCDXBL))) +
    geom_line() +
    labs(x="") +
    scale_color_brewer(palette = "YlOrRd", name = "Baseline diagnosis") + 
    theme_classic()
  
  gp_pet <- ggplot(data = nacc_dpm, aes(x=Years, y=GAAIN_SUVR, group = NACCID, color = as.factor(NACCDXBL))) +
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

fit_dpm <- function(remove_ceiling=FALSE, include_ab_pet=FALSE, scale_t = FALSE) {
  nacc_dpm <- get_dpm_data() 
  
  plot_raw(nacc_dpm)
  
  rids <- unique(nacc_dpm$NACCID)
  
  nacc_dpm <- as.data.frame(nacc_dpm) %>% 
    mutate(RID = as.numeric(gsub("NACC", "", NACCID)))
  
  if (remove_ceiling) {
    # A clear ceiling effect can be observed. This removes subsequent visits after
    # a max/min value have been observed
    nacc_dpm <- nacc_dpm %>% 
      group_by(NACCID) %>% arrange(Years) %>%
      filter(is.na(CDRSUM) | CDRSUM < 18 | cumsum(!is.na(CDRSUM) & CDRSUM == 18) == 1) %>%
      filter(is.na(NACCMOCA) | NACCMOCA > 0 | cumsum(!is.na(NACCMOCA) & NACCMOCA == 0) == 1) %>%
      ungroup()
  }
  
  #nacc_dpm <- nacc_dpm %>% filter(!is.na(CDRSUM) & !is.na(NACCMOCA))
  
  cols <- c("NACCID", "Years", "NACCDXBL", "CN", "PMCI", "MCI", "AD", "negAB.bl")
  nacc_dpm_long <- rbind(nacc_dpm[c("CDRSUM", cols)] |> mutate(scale = "CDRSUM") |> rename(value = CDRSUM), 
                         #nacc_dpm[c("NACCMMSE", cols)] |> mutate(scale = "NACCMMSE") |> rename(value = NACCMMSE),
                         nacc_dpm[c("NACCMOCA", cols)] |> mutate(scale = "NACCMOCA") |> rename(value = NACCMOCA))
  
  if (include_ab_pet) {
    nacc_dpm_long <- rbind(nacc_dpm_long, 
                           nacc_dpm[c("GAAIN_SUVR", cols)] |> mutate(scale = "GAAIN_SUVR") |> rename(value = GAAIN_SUVR))
  }
  
  nacc_dpm_long <- nacc_dpm_long %>% drop_na(value) %>% 
    select(NACCID, Years, scale, value, NACCDXBL, CN, PMCI, MCI, AD, negAB.bl)
  
  if (FALSE) {
    nacc_dpm_long <- rbind(
      nacc_dpm_long %>%
      group_by(NACCID) %>%
      filter(scale == "CDRSUM") %>%
        filter(value < 18 | cumsum(value == 18) == 1),
      nacc_dpm_long %>%
        group_by(NACCID) %>%
        filter(scale == "NACCMOCA") %>%
        filter(value > 0 | cumsum(value == 0) == 1)
    )
  }
  
  nacc_dpm_long <- nacc_dpm_long %>% mutate(
    scale = as.factor(scale),
    NACCDXBL = as.factor(NACCDXBL),
    negAB.bl = as.numeric(negAB.bl),
    RID = as.numeric(gsub("NACC", "", NACCID))
  )
  
  obs_by_id <- nacc_dpm_long %>%
    filter(!is.na(value)) %>%
    count(RID, scale) %>%
    tidyr::pivot_wider(
      names_from = scale,
      values_from = n,
      values_fill = 0
    )
  zero_ids <- obs_by_id %>% filter(NACCMOCA==0 | CDRSUM==0) %>% select(RID) %>% unlist()
  #nacc_dpm_long <- filter(nacc_dpm_long, !(RID %in% zero_ids))
  
  nacc_dpm_long$t <- nacc_dpm_long$Years
  if (scale_t) {
    nacc_dpm_long$t <- scale(nacc_dpm_long$t)
  } 
  
  if (scale_t) {
    mu <- attr(nacc_dpm_long$t, 'scaled:center')
    sg <- attr(nacc_dpm_long$t, 'scaled:scale')
  } else {
    mu <- 0
    sg <- 1
  }
  
  if (include_ab_pet) {
    v.df <- nacc_dpm %>% filter(CN == 1 & Months == 0) %>% ungroup() %>% 
      select(CDRSUM, NACCMOCA, GAAIN_SUVR) %>% colMeans(na.rm=T)
    
    
    fixed_start_coef_y <- c(l.scaleNACCMOCA = -0.05, 
                            l.scaleGAAIN_SUVR = 0.5, 
                            l.scaleCDRSUM = 0.25,
                            s.PMCI = (3-mu)/sg,
                            s.MCI = (6-mu)/sg,
                            s.AD = (12-mu)/sg,
                            s.negAB.bl = (-1-mu)/sg, 
                            g.scaleNACCMOCA = 1,
                            g.scaleCDRSUM = 1,
                            g.scaleGAAIN_SUVR = 2,
                            v.df)
  } else {
    v.df <- nacc_dpm %>% filter(CN == 1 & Months == 0) %>% ungroup() %>% 
      select(CDRSUM, NACCMOCA) %>% colMeans(na.rm=T)
    
    
    fixed_start_coef_y <- c(l.scaleNACCMOCA = -0.05, 
                            #l.scaleNACCMMSE = -0.05, 
                            l.scaleCDRSUM = 0.25,
                            s.PMCI = (3-mu)/sg,
                            s.MCI = (6-mu)/sg,
                            s.AD = (12-mu)/sg,
                            s.negAB.bl = (-1-mu)/sg, 
                            g.scaleNACCMOCA = 1,
                            #g.scaleNACCMMSE = 1,
                            g.scaleCDRSUM = 1,
                            v.df)
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
  
  if (include_ab_pet) {
    nacc_dpm |> mutate(time_shift = fixed_shift_multi + random_shift_multi,
                       NACCDXBL = factor(NACCDXBL, labels = c("CN", "Impaired", "MCI", "Dementia"))) |> 
      select(RID, Years, time_shift, CDRSUM, NACCMOCA, GAAIN_SUVR, NACCDXBL) -> dpm_out
    
    gp_ab_pet <- ggplot(data = dpm_out, aes(x=Years+time_shift, y=GAAIN_SUVR, group = RID, color = NACCDXBL)) +
      geom_line() +
      labs(x="") +
      scale_color_brewer(palette = "YlOrRd", name = "Baseline diagnosis") +
      theme_classic()
  } else {
    nacc_dpm |> mutate(time_shift = fixed_shift_multi + random_shift_multi,
                       NACCDXBL = factor(NACCDXBL, labels = c("CN", "Impaired", "MCI", "Dementia"))) |> 
      select(RID, Years, time_shift, CDRSUM, NACCMOCA, NACCDXBL) -> dpm_out
  }
  
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
  
  if (include_ab_pet) {
    gp <- ggarrange(gp_moca, gp_cdr, gp_ab_pet, ncol = 1, common.legend = TRUE, legend = "right")
  } else {
    gp <- ggarrange(gp_moca, gp_cdr, ncol = 1, common.legend = TRUE, legend = "right") 
  }
  plot(gp)
  
  return(dpm_out)
}

# Run NACC dpm like this:
# nacc_dpm_out <- fit_dpm(remove_ceiling = FALSE, include_ab_pet = FALSE)
