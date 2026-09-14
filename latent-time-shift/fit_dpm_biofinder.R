library(dplyr)
library(progmod)
library(ggplot2)

source("~/R/LTC/utils/biofinder_data_loaders.R")

plot_raw <- function(dpm_df) {
  library(ggpubr)
  gp_adas <- ggplot(data = dpm_df, aes(x=Years, y=adas, group = sid, color = diag)) +
    geom_line() +
    labs(x="") +
    scale_color_brewer(palette = "YlOrRd", name = "Baseline diagnosis") + 
    theme_classic()
  
  gp_pet <- ggplot(data = dpm_df, aes(x=Years, y=fnc_suvr, group = sid, color = diag)) +
    geom_point() +
    geom_line() +
    labs(x="") +
    scale_color_brewer(palette = "YlOrRd", name = "Baseline diagnosis") +
    theme_classic()
  
  gp_mmse <- ggplot(data = dpm_df, aes(x=Years, y=mmse, group = sid, color = diag)) +
    geom_line() +
    labs(x="") +
    scale_color_brewer(palette = "YlOrRd", name = "Baseline diagnosis") +
    theme_classic()
  
  gp <- ggarrange(gp_adas, gp_mmse, gp_pet, ncol = 1, common.legend = TRUE, legend = "right")
  plot(gp)
}

plot_long <- function(long_df, id_col="RID", diag_col="DX.bl") {
  library(ggpubr)
  scales <- unique(long_df$scale)
  
  plotList <- c()
  for (s in scales) {
    plot_df <- filter(long_df, scale == s)
    gp <- ggplot(data = plot_df, aes(x=t, y=value, group = .data[[id_col]], color = .data[[diag_col]])) +
      geom_point(size=1) +
      geom_line() +
      labs(x="", y=s) +
      scale_color_brewer(palette = "YlOrRd", name = "Baseline diagnosis") + 
      theme_classic()
    plotList <- c(plotList, gp)
  } 
  
  gp <- ggarrange(plotlist = plotList, ncol = 1, common.legend = TRUE, legend = "right")
  plot(gp)
}

dpm_df <- get_biofinder_dpm()


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

cols <- c("adas", "mmse", "fnc_suvr")
v.df <- dpm_df %>% filter(CU == 1 & Visit == 0) %>% ungroup() %>% select(all_of(cols)) %>% colMeans(na.rm=T)

scale_t=FALSE
if (scale_t) {
  mu <- attr(adni_dpm_long$t, 'scaled:center')
  sg <- attr(adni_dpm_long$t, 'scaled:scale')
} else {
  mu <- 0
  sg <- 1
}


fixed_start_coef_y <- c(l.scaleadas = 0.05, #0.05,
                        l.scalemmse = -0.1, #-0.04,
                        l.scalefnc_suvr = 0.5, #0.05,
                        s.SCD = (3-mu)/sg,
                        s.MCI = (6-mu)/sg,
                        s.AD = (12-mu)/sg,
                        s.negAB.bl = (-1-mu)/sg, 
                        g.scaleadas = 0.75,
                        g.scalemmse = 1,
                        g.scalefnc_suvr = 2,
                        v.scale=v.df) 


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
                     data = dpm_long,
                     fixed = list(l ~ scale + 0,
                                  s ~ SCD + MCI + AD + negAB.bl + 0,
                                  g ~ scale + 0,
                                  v ~ scale + 0),
                     random = list(s ~ 1,
                                   v ~ 1),
                     groups = ~ sid,
                     start = fixed_start_coef_y,
                     weights = varIdent(form = ~ 1 | scale),
                     method = "REML",
                     control = ctrl)

dpm_df$fixed_shift_multi <- with(dpm_df,
                                 SCD * (fixed.effects(dpm_model)["s.SCD"]*sg)+mu +
                                   MCI * (fixed.effects(dpm_model)["s.MCI"]*sg)+mu +
                                     AD * (fixed.effects(dpm_model)["s.AD"]*sg)+mu +
                                     negAB.bl * (fixed.effects(dpm_model)["s.negAB.bl"])*sg+mu)

pred_rand <- random.effects(dpm_model)
dpm_df$random_shift_multi <- (pred_rand[match(dpm_df$sid, rownames(pred_rand)), 's.(Intercept)']*sg)+mu


dpm_df |> mutate(time_shift = fixed_shift_multi + random_shift_multi) |> 
  select(sid, Years, time_shift, adas, mmse, fnc_suvr, diag) -> dpm_out

gp_adas <- ggplot(data = dpm_out, aes(x=Years+time_shift, y=adas, group = sid, color = diag)) +
  geom_line() +
  labs(x="") +
  scale_color_brewer(palette = "YlOrRd", name = "Baseline diagnosis") + 
  theme_classic()

gp_pet <- ggplot(data = dpm_out, aes(x=Years+time_shift, y=fnc_suvr, group = sid, color = diag)) +
  geom_line() +
  labs(x="") +
  scale_color_brewer(palette = "YlOrRd", name = "Baseline diagnosis") +
  theme_classic()

gp_mmse <- ggplot(data = dpm_out, aes(x=Years+time_shift, y=mmse, group = sid, color = diag)) +
  geom_line() +
  labs(x="") +
  scale_color_brewer(palette = "YlOrRd", name = "Baseline diagnosis") +
  theme_classic()

gp <- ggarrange(gp_adas, gp_mmse, gp_pet, ncol = 1, common.legend = TRUE, legend = "right")
plot(gp)

write.csv(dpm_out, "~/R/EDAP-data/BioFINDER/DPM_BioFINDER.csv")


