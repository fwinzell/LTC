library(ADNIMERGE)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(dplyr)
library(Rtsne)
library(paletteer)
library(tidyr)
library(broom)
library(nlme)

library(icenReg)
library(survival)

source("~/R/LTC/utils/oasis_data_loaders.R")
run <- "oasis_exp_km_ab_2"
load(paste("~/R/EDAP-data/LTC4/", run, ".Rdata", sep = ""))

mri_data <- load_oasis_mri_data(unified_norm=TRUE)
ab_pos_ids <- get_ab_pos_ids()
mri_controls <- mri_data %>% filter(DX.bl == "CN" & !(OASISID %in% ab_pos_ids))

dpm_df <- read.csv("~/R/EDAP-data/OASIS/DPM_OASIS.csv", header=TRUE)
mri_data <- dpm_df %>% select(RID, time_shift) %>% distinct() %>% right_join(mri_data, by="RID") %>%
  mutate(Time = Years + time_shift)

Clusters <- data.frame(
  Cluster = oasisLTC@Cluster,
  RID = oasisLTC@RID
)

Clusters %>% left_join(mri_data, by = "RID") %>% drop_na(Cluster) -> mri_data

source("~/R/LTC/utils/analysis_utils.R")

brain_regions <- get_brain_regions() 

limbic_vars <- filter(brain_regions, Region == "limbic")

##### make plots ######

mean_fn <- list(
  mu = ~ mean(.x, na.rm=TRUE),
  sigma = ~ sd(.x, na.rm=TRUE)
)

mean_and_sd <- function(vars) {
  vars <- unname(vars)
  mean_df <- mri_controls %>% select(RID, Months, all_of(vars)) %>% 
    mutate(v = rowSums(across(all_of(vars)))) %>%
    arrange(Months) %>% distinct(RID, .keep_all = TRUE) %>%
    summarise(across(v, mean_fn)) 
  return(mean_df)
}

# Hippocampus
vars <- c("Left.Hippocampus_volume", "Right.Hippocampus_volume") %>% setNames(c("Left", "Right"))
hippo <- mri_data %>% select(RID, DX.bl, Time, Cluster, all_of(vars)) %>% mutate(v = Left + Right)

mean_df <- mean_and_sd(vars)

hippo_plot <- survival_analysis(hippo, mu=mean_df$v_mu, sigma=mean_df$v_sigma, title_name="Hippocampus")
plot(hippo_plot$plot)


# Entorhinal
vars <- grepv("*entorhinal", colnames(mri_data)) %>% setNames(c("Left", "Right"))
entor <- mri_data %>% select(RID, DX.bl, Time, Cluster, all_of(vars)) %>% mutate(v = Left + Right)

mean_df <- mean_and_sd(vars)

ent_plot <- survival_analysis(entor, mu=mean_df$v_mu, sigma=mean_df$v_sigma, title_name="Entorhinal")
plot(ent_plot$plot)

# Accumbens Area
vars <- grepv("Accumbens", colnames(mri_data)) %>% setNames(c("Left", "Right"))
accum <- mri_data %>% select(RID, DX.bl, Time, Cluster, all_of(vars)) %>% mutate(v = Left + Right)

mean_df <- mean_and_sd(vars)

accum_plot <- survival_analysis(accum, mu=mean_df$v_mu, sigma=mean_df$v_sigma, title_name="Accumbens Area")

# Left-temporal
vars <- grepv("lh.*temporal_volume", colnames(mri_data))

lat_temp <- mri_data %>% select(RID, DX.bl, Time, Cluster, all_of(vars)) %>% 
  mutate(v = rowSums(across(all_of(vars))))

mean_df <- mean_and_sd(vars)

lt_p <- survival_analysis(lat_temp, mu=mean_df$v_mu, sigma=mean_df$v_sigma, title_name="Left Temporal", z_threshold = -1.28)


# Parietal
vars <- grepv("*parietal", colnames(mri_data))
parietal <- mri_data %>% select(RID, DX.bl, Time, Cluster, all_of(vars)) %>%
  mutate(v = rowSums(across(all_of(vars)))) %>% drop_na(v)

mean_df <- mean_and_sd(vars)

par_p <- survival_analysis(parietal, mu=mean_df$v_mu, sigma=mean_df$v_sigma, "Parietal", z_threshold = -1.28)


# Fusiform
vars <- grepv("*fusiform", colnames(mri_data)) %>% setNames(c("Left", "Right"))

fusi <- mri_data %>% select(RID, DX.bl, Time, Cluster, all_of(vars)) %>%
  mutate(v = Left + Right)
mean_df <- mean_and_sd(vars)

fusi_p <- survival_analysis(fusi, mu=mean_df$v_mu, sigma=mean_df$v_sigma, "Fusiform", z_threshold = -1.28)


# Insula
vars <- grepv("*insula", colnames(mri_data)) %>% setNames(c("Left", "Right"))

insul <- mri_data %>% select(RID, DX.bl, Time, Cluster, all_of(vars)) %>%
  mutate(v = Left + Right) %>% drop_na(v)
mean_df <- mean_and_sd(vars)

insul_p <- survival_analysis(insul, mu=mean_df$v_mu, sigma=mean_df$v_sigma, "Insula", z_threshold = -1)


plot <- ggarrange(
  hippo_plot$plot,
  ent_plot$plot,
  lt_p$plot,
  par_p$plot,
  fusi_p$plot,
  insul_p$plot,
  ncol = 3, nrow =2, common.legend = TRUE,
  legend = "bottom"
)

plot(plot)
ggsave("~/R/EDAP-data/plots/LTC4/oasis_mri_surv_plots_ab.png", plot, width = 10, height = 6, dpi =500)

#### make table ####

format_summary <- function(summary) {
  df <- summary %>% lapply(function(col) {
    if (is.numeric(col)) sprintf("$%.3f$", col) else as.character(col)
  }) %>% data.frame() %>%
    mutate(CI = sprintf("(%s, %s)", .[[2]], .[[3]])) %>%
    select(matches("HR.*"), CI, matches("^p.*")) 
  
  suffix <- sub("^HR\\.", "", names(df)[1])
  df <- rename(df, !!paste0("CI.", suffix) := CI)
  
  return(df)
}

df <- cbind( format_summary(hippo_plot$summary),
             format_summary(ent_plot$summary),
             format_summary(lt_p$summary),
             format_summary(par_p$summary),
             format_summary(fusi_p$summary),
             format_summary(insul_p$summary) )



## make latex table
library(knitr)

df1 <- df[, 1:(ncol(df)/3)]
df1 <- cbind(Cluster = gsub("Cluster", "", rownames(hippo_plot$summary)), df1)
df2 <- df[, (ncol(df)/3 + 1):(2*ncol(df)/3)]
df2 <- cbind(Cluster = gsub("Cluster", "", rownames(hippo_plot$summary)), df2)
df3 <- df[, (2*ncol(df)/3 + 1):ncol(df)]
df3 <- cbind(Cluster = gsub("Cluster", "", rownames(hippo_plot$summary)), df3)

kable(df1, format = "latex", escape=FALSE, booktabs=TRUE)

kable(df2, format = "latex", escape=FALSE, booktabs=TRUE)

kable(df3, format = "latex", escape=FALSE, booktabs=TRUE)




