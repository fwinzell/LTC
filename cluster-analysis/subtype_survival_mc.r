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

multi_cohort_df <- read.csv("~/R/EDAP-data/MULTI_COHORT.csv", header = TRUE)
all.vars <- c(grepv("^(RH_|LH_|CC_)", colnames(multi_cohort_df)), "BRAINSTEM")

run <- "exp_km_ab_ao"
load(paste("~/R/EDAP-data/LTC_MC/", run, ".Rdata", sep = ""))

adni_dl <- new.env()
source("~/R/LTC/utils/adni_data_loaders.R", local=adni_dl)

oasis_dl <- new.env()
source("~/R/LTC/utils/oasis_data_loaders.R", local=oasis_dl)

nacc_dl <- new.env()
source("~/R/LTC/utils/nacc_data_loaders.R", local=nacc_dl)

mc_dl <- new.env()
source("~/R/LTC/utils/multi_cohort_loader.R", local=mc_dl)

source("~/R/LTC/utils/analysis_utils.R")

ab_pos_rids <- c(paste0("ADNI_", adni_dl$get_ab_pos_ids()),
                 paste0("NACC_", nacc_dl$get_ab_pos_ids()),
                 gsub("OAS", "OASIS_", oasis_dl$get_ab_pos_ids()))

mc_mri <- mc_dl$multi_cohort_mri()
mri_controls <- mc_mri %>% filter(DX.bl == "CN" & !(RID %in% ab_pos_rids))

# Might need to change this line - accidentially made RID numeric...
#crids <- levels(multiLTC@betas$RID)[multiLTC@RID]
Clusters <- data.frame(
  Cluster = multiLTC@Cluster,
  RID = multiLTC@RID
)

Clusters %>% left_join(multi_cohort_df, by = "RID") %>% drop_na(Cluster) -> multi_cohort_df

multi_cohort_df %>% distinct(RID, Cluster, Cohort) %>% select(Cluster, Cohort) %>% table()

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
vars <- c("LH_HIPPOCAMPUS", "RH_HIPPOCAMPUS") %>% setNames(c("Left", "Right"))
hippo <- multi_cohort_df %>% select(RID, DX.bl, Time, Cluster, all_of(vars)) %>% mutate(v = Left + Right)

mean_df <- mean_and_sd(vars)

hippo_plot <- survival_analysis(hippo, mu=mean_df$v_mu, sigma=mean_df$v_sigma, title_name="Hippocampus")
plot(hippo_plot$plot)


# Entorhinal
vars <- c("LH_ENTORHINAL", "RH_ENTORHINAL") %>% setNames(c("Left", "Right"))
entor <- multi_cohort_df %>% select(RID, DX.bl, Time, Cluster, all_of(vars)) %>% mutate(v = Left + Right)

mean_df <- mean_and_sd(vars)

ent_plot <- survival_analysis(entor, mu=mean_df$v_mu, sigma=mean_df$v_sigma, title_name="Entorhinal")
plot(ent_plot$plot)

# Accumbens Area
vars <- grepv("ACCUMBENS", colnames(multi_cohort_df)) %>% setNames(c("Left", "Right"))
accum <- multi_cohort_df %>% select(RID, DX.bl, Time, Cluster, all_of(vars)) %>% mutate(v = Left + Right)

mean_df <- mean_and_sd(vars)

accum_plot <- survival_analysis(accum, mu=mean_df$v_mu, sigma=mean_df$v_sigma, title_name="Accumbens Area", z_threshold = -1.28)

# Left-temporal
vars <- grepv("LH.*TEMPORAL", colnames(multi_cohort_df))

lat_temp <- multi_cohort_df %>% select(RID, DX.bl, Time, Cluster, all_of(vars)) %>% 
  mutate(v = rowSums(across(all_of(vars))))

mean_df <- mean_and_sd(vars)

lt_p <- survival_analysis(lat_temp, mu=mean_df$v_mu, sigma=mean_df$v_sigma, title_name="Left Temporal")


# Parietal
vars <- grepv("*PARIETAL", colnames(multi_cohort_df))
parietal <- multi_cohort_df %>% select(RID, DX.bl, Time, Cluster, all_of(vars)) %>%
  mutate(v = rowSums(across(all_of(vars)))) %>% drop_na(v)

mean_df <- mean_and_sd(vars)

par_p <- survival_analysis(parietal, mu=mean_df$v_mu, sigma=mean_df$v_sigma, "Parietal")


# Fusiform
vars <- c("LH_FUSIFORM", "RH_FUSIFORM") %>% setNames(c("Left", "Right"))

fusi <- multi_cohort_df %>% select(RID, DX.bl, Time, Cluster, all_of(vars)) %>%
  mutate(v = Left + Right)
mean_df <- mean_and_sd(vars)

fusi_p <- survival_analysis(fusi, mu=mean_df$v_mu, sigma=mean_df$v_sigma, "Fusiform")


# Insula
vars <- c("LH_INSULA", "RH_INSULA") %>% setNames(c("Left", "Right"))

insul <- multi_cohort_df %>% select(RID, DX.bl, Time, Cluster, all_of(vars)) %>%
  mutate(v = Left + Right) %>% drop_na(v)
mean_df <- mean_and_sd(vars)

insul_p <- survival_analysis(insul, mu=mean_df$v_mu, sigma=mean_df$v_sigma, "Insula", z_threshold = -1.28)


plot <- ggarrange(
  hippo_plot$plot,
  ent_plot$plot,
  accum_plot$plot,
  lt_p$plot,
  par_p$plot,
  fusi_p$plot,
  #insul_p$plot,
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




