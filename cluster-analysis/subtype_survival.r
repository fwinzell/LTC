library(ADNIMERGE2)
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

source("~/R/LTC/utils/adni_data_loaders.R")
run <- "exp_km_ab"
load(paste("~/R/EDAP-data/LTC4/", run, ".Rdata", sep = ""))

ucsf_data <- ucsf_longitudinal_all(only_vol = TRUE, filter_n = 0)

ab_df <- get_ab_df()
ab_df %>% select(RID, AB_any) %>% group_by(RID) %>% 
  mutate(AB_any = if_any(AB_any)) %>% distinct() %>% filter(AB_any) %>% 
  select(RID) %>% unlist() -> ab_pos_rids_any

mri_controls <- ucsf_data %>% filter(!(RID %in% ab_pos_rids_any)) %>% filter(DX.bl == "CN")

dpm_df <- read.csv("~/R/EDAP-data/ADNI/DPM_ADNI.csv", header=TRUE)
ucsf_data <- dpm_df %>% select(RID, time_shift) %>% distinct() %>% right_join(ucsf_data, by="RID") %>%
  mutate(Time = Years + time_shift)

Clusters <- data.frame(
  Cluster = adniLTC@Cluster,
  RID = adniLTC@RID
)

Clusters %>% left_join(ucsf_data, by = "RID") %>% drop_na(Cluster) -> ucsf_data

datadic <- datadict_as_tibble(ADNIMERGE2::DATADIC) %>% filter(TBLNAME == "UCSFFSL") %>%
  select(FLDNAME, TEXT) %>% filter(grepl("ST\\d+[A-Z]{2}", FLDNAME)) %>%
  mutate(Region = sub(".*\\s", "", TEXT))

get_mean_df <- function(vars) {
  
}


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

ucsf_data$Cluster <- relevel(ucsf_data$Cluster, ref="A")

# Hippocampus
vars <- filter(datadic, grepl("*Hippocampus", Region)) %>% select(FLDNAME) %>% unlist() %>% 
  setNames(c("Left", "Right"))
hippo <- ucsf_data %>% select(RID, DX.bl, Time, Cluster, all_of(vars)) %>% mutate(v = Left + Right)

mean_df <- mean_and_sd(vars)

hippo_plot <- survival_analysis(hippo, mu=mean_df$v_mu, sigma=mean_df$v_sigma, title_name="Hippocampus")
plot(hippo_plot$plot)


# Entorhinal
vars <- filter(datadic, grepl("*Entorhinal", Region)) %>% filter(grepl("*CV", FLDNAME)) %>%
  select(FLDNAME) %>% unlist() %>%setNames(c("Left", "Right"))

mean_df <- mean_and_sd(vars)

entor <- ucsf_data %>% select(RID, DX.bl, Time, Cluster, all_of(vars)) %>% mutate(v = Left + Right)
ent_plot <- survival_analysis(entor, mu=mean_df$v_mu, sigma=mean_df$v_sigma, "Entorhinal")
plot(ent_plot$plot)

# Accumbens Area
vars <- filter(datadic, grepl("*Accumbens", Region)) %>% select(FLDNAME) %>% unlist() %>% 
  setNames(c("Left", "Right"))

mean_df <- mean_and_sd(vars)

accum <- ucsf_data %>% select(RID, DX.bl, Time, Cluster, all_of(vars)) %>% mutate(v = Left + Right)
accum_plot <- survival_analysis(accum, "Accumbens Area", mu=mean_df$v_mu, sigma=mean_df$v_sigma)


# Left-temporal
#ucsf_data$Cluster <- relevel(ucsf_data$Cluster, ref="B")
vars <- filter(datadic, grepl("Left.*Temporal$", Region)) %>% filter(grepl("*CV", FLDNAME)) %>%
  select(FLDNAME) %>% unlist() %>% unname()

mean_df <- mean_and_sd(vars)

lat_temp <- ucsf_data %>% select(RID, DX.bl, Time, Cluster, all_of(vars)) %>% 
  mutate(v = rowSums(across(all_of(vars))))
lt_p <- survival_analysis(lat_temp, "Left Temporal", mu=mean_df$v_mu, sigma=mean_df$v_sigma, z_threshold = -1.28)
plot(lt_p$plot)

# Parietal
vars <- filter(datadic, grepl("*Parietal$", Region)) %>% filter(grepl("*CV", FLDNAME)) %>%
  select(FLDNAME) %>% unlist() %>% unname()

mean_df <- mean_and_sd(vars)

parietal <- ucsf_data %>% select(RID, DX.bl, Time, Cluster, all_of(vars)) %>%
  mutate(v = rowSums(across(all_of(vars)))) %>% drop_na(v)
par_p <- survival_analysis(parietal, mu=mean_df$v_mu, sigma=mean_df$v_sigma,  "Parietal")
plot(par_p$plot)

# Fusiform
vars <- filter(datadic, grepl("*Fusiform", Region)) %>% filter(grepl("*CV", FLDNAME)) %>%
  select(FLDNAME) %>% unlist() %>%setNames(c("Left", "Right"))

mean_df <- mean_and_sd(vars)

fusi <- ucsf_data %>% select(RID, DX.bl, Time, Cluster, all_of(vars)) %>%
  mutate(v = Left + Right)
fusi_p <- survival_analysis(fusi, mu=mean_df$v_mu, sigma=mean_df$v_sigma,  "Fusiform", z_threshold = -1.28)

# Occipital
#vars <- filter(datadic, grepl("*Occipital", Region)) %>% filter(grepl("*CV", FLDNAME)) %>%
#  select(FLDNAME) %>% unlist() %>%setNames(c("Left", "Right"))

#mean_df <- mean_and_sd(vars)

#occipt <- ucsf_data %>% select(RID, DX.bl, Time, Cluster, all_of(vars)) %>%
#  mutate(v = Left + Right) %>% drop_na(v)
#occip_p <- survival_analysis(occipt, mu=mean_df$v_mu, sigma=mean_df$v_sigma,  "Lateral Occipital")
#plot(occip_p$plot)

# Insula
vars <- filter(datadic, grepl("*Insula", Region)) %>% filter(grepl("*CV", FLDNAME)) %>%
  select(FLDNAME) %>% unlist() %>%setNames(c("Left", "Right"))

mean_df <- mean_and_sd(vars)

insul <- ucsf_data %>% select(RID, DX.bl, Time, Cluster, all_of(vars)) %>%
  mutate(v = Left + Right) %>% drop_na(v)

# Increase z-score threshold to get events for all clusters
insul_p <- survival_analysis(insul, "Insula", mu=mean_df$v_mu, sigma=mean_df$v_sigma, z_threshold = -1.15)
plot(insul_p$plot)
insul_p$summary

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
ggsave("~/R/EDAP-data/plots/LTC4/mri_surv_plots_ab.png", plot, width = 10, height = 6, dpi =500)

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

df1 <- df[, 1:(ncol(df)/2)]
df1 <- cbind(Cluster = gsub("Cluster", "", rownames(hippo_plot$summary)), df1)
df2 <- df[, (ncol(df)/2 + 1):ncol(df)]
df2 <- cbind(Cluster = gsub("Cluster", "", rownames(hippo_plot$summary)), df2)


kable(df1, format = "latex", escape=FALSE, booktabs=TRUE)


kable(df2, format = "latex", escape=FALSE, booktabs=TRUE)





