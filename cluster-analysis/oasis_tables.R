library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(tidyverse)

library(aricode)
library(infotheo)
library(paletteer)

library(lme4)
library(stats)
library(reshape2)
library(tibble)
library(RColorBrewer)

library(lmerTest)
library(emmeans)

source("~/R/LTC/utils/oasis_data_loaders.R")

run <- "oasis_exp_km_ab"
load(paste("~/R/EDAP-data/LTC4/", run, ".Rdata", sep = ""))

mri_data <- load_oasis_mri_data(unified_norm=TRUE)
dpm_df <- read.csv("~/R/EDAP-data/OASIS/DPM_OASIS.csv", header=TRUE)
mri_data <- dpm_df %>% select(RID, time_shift) %>% distinct() %>% right_join(mri_data, by="RID") %>%
  mutate(Time = Years + time_shift)

Clusters <- data.frame(
  Cluster = oasisLTC@Cluster,
  RID = oasisLTC@RID
)

Clusters %>% left_join(mri_data, by = "RID") %>% drop_na(Cluster) -> mri_data


# Demographics
oasis_demog <- get_demographics()
oasis_dx <- get_diagnoses()

oasis_demog <- oasis_dx %>% distinct(RID, DX.bl) %>% right_join(oasis_demog, by="RID")

table_df <- select(oasis_demog, OASISID, RID, DX.bl, AgeatEntry, GENDER, EDUC, APOE4) %>%
  filter(RID %in% oasisLTC@RID)

# AB status
oasis_ab <- get_oasis_ab(FALSE) %>% group_by(OASISID, tracer) %>% summarise(
  PET_SUVR = max(PET_SUVR),
  .groups = "drop_last"
) %>% pivot_wider(names_from = tracer, values_from = PET_SUVR)

table_df <- left_join(table_df, oasis_ab, by="OASISID")


# Tau-PET
oasis_tau <- get_tau_pet(get_braak=TRUE) %>% group_by(OASISID) %>%
  summarise(
    across(Braak1_2:Tauopathy, max)
  ) %>% ungroup()

##### DATASET TABLE #####
table_df <- left_join(table_df, oasis_tau, by="OASISID")

mean_df <- table_df %>% group_by(DX.bl) %>% summarise(
  Age = round(mean(AgeatEntry),3),
  EDUC = round(mean(EDUC),2),
  FBP = round(mean(AV45, na.rm=T), 3),
  PIB = round(mean(PIB, na.rm=T), 3),
  Braak1_2 = round(mean(Braak1_2, na.rm=T), 3),
  Braak3_4 = round(mean(Braak3_4, na.rm=T), 3),
  Braak5_6 = round(mean(Braak5_6, na.rm=T), 3),
  Tauopathy = round(mean(Tauopathy, na.rm=T), 3)
) %>% ungroup()

mean_df <- table_df %>% summarise(
  Age = round(mean(AgeatEntry),3),
  EDUC = round(mean(EDUC),2),
  FBP = round(mean(AV45, na.rm=T), 3),
  PIB = round(mean(PIB, na.rm=T), 3),
  Braak1_2 = round(mean(Braak1_2, na.rm=T), 3),
  Braak3_4 = round(mean(Braak3_4, na.rm=T), 3),
  Braak5_6 = round(mean(Braak5_6, na.rm=T), 3),
  Tauopathy = round(mean(Tauopathy, na.rm=T), 3),
  DX.bl = "All"
) %>% rbind(mean_df) %>% tibble::column_to_rownames("DX.bl")

sd_df <- table_df %>% group_by(DX.bl) %>% summarise(
  Age = round(sd(AgeatEntry),3),
  EDUC = round(sd(EDUC),2),
  FBP = round(sd(AV45, na.rm=T), 3),
  PIB = round(sd(PIB, na.rm=T), 3),
  Braak1_2 = round(sd(Braak1_2, na.rm=T), 3),
  Braak3_4 = round(sd(Braak3_4, na.rm=T), 3),
  Braak5_6 = round(sd(Braak5_6, na.rm=T), 3),
  Tauopathy = round(sd(Tauopathy, na.rm=T), 3)
) %>% ungroup()

sd_df <- table_df %>% summarise(
  Age = round(sd(AgeatEntry),3),
  EDUC = round(sd(EDUC),2),
  FBP = round(sd(AV45, na.rm=T), 3),
  PIB = round(sd(PIB, na.rm=T), 3),
  Braak1_2 = round(sd(Braak1_2, na.rm=T), 3),
  Braak3_4 = round(sd(Braak3_4, na.rm=T), 3),
  Braak5_6 = round(sd(Braak5_6, na.rm=T), 3),
  Tauopathy = round(sd(Tauopathy, na.rm=T), 3),
  DX.bl = "All"
) %>% rbind(sd_df) %>% tibble::column_to_rownames("DX.bl")

perc_df <- table_df %>% group_by(DX.bl) %>% summarise(
  n = n(),
  Male = round(mean(GENDER == 1, na.rm=T)*100,1),
  Female = 100-Male,
  APOE4.0 = round(mean(APOE4 == 0, na.rm=T)*100,1),
  APOE4.1 = round(mean(APOE4 == 1, na.rm=T)*100,1),
  APOE4.2 = round(mean(APOE4 == 2, na.rm=T)*100,1),
) %>% ungroup()

perc_df <- table_df %>% summarise(
  n = n(),
  Male = round(mean(GENDER == 1, na.rm=T)*100,1),
  Female = 100-Male,
  APOE4.0 = round(mean(APOE4 == 0, na.rm=T)*100,1),
  APOE4.1 = round(mean(APOE4 == 1, na.rm=T)*100,1),
  APOE4.2 = round(mean(APOE4 == 2, na.rm=T)*100,1),
  DX.bl = "All"
) %>% rbind(perc_df) %>% tibble::column_to_rownames("DX.bl")


mean_df <- t(mean_df) %>% as.data.frame() %>%
  tibble::rownames_to_column("Variable") 
sd_df <- t(sd_df) %>% as.data.frame() %>%
  tibble::rownames_to_column("Variable") 
perc_df <- t(perc_df) %>% as.data.frame() %>%
  tibble::rownames_to_column("Variable") 

latex_strings <- mapply(
  function(m, s) sprintf("$%.3f \\pm %.3f$", m, s),
  mean_df[-1], sd_df[-1],
  SIMPLIFY = FALSE
) |> as.data.frame() |> mutate(Variable = mean_df$Variable)

latex_strings <- rbind(perc_df, latex_strings)

variable_names <- c("n", "Sex (M)", "Sex (F)", "\\textit{APOE} $\\epsilon$4 (0)", "\\textit{APOE} $\\epsilon$4 (1)", "\\textit{APOE} $\\epsilon$4 (2)",
                    "Age (Y)", "Educat. (Y)", 
                    "FBP", "PiB", 
                    "Braak I-II", "Braak III-IV", "Braak V-VI", "Cereb. Cortex Ref.")
latex_strings$Variable <- variable_names

library(knitr)
latex_strings <- select(latex_strings, Variable, CN, Impaired, MCI, Dementia, All)
kable(latex_strings, format = "latex", escape = FALSE, booktabs = TRUE)

##### CLUSTER ANALYSIS #####

source("~/R/LTC/utils/analysis_utils.R")

oasis_demog <- inner_join(oasis_demog, Clusters, by="RID")
oasis_dx <- left_join(oasis_dx, Clusters, by="RID") %>% drop_na(Cluster)
oasis_dpm <- read.csv("~/R/EDAP-data/OASIS/DPM_OASIS.csv", header=TRUE) %>% inner_join(Clusters, by="RID") %>%
  distinct(RID, time_shift)

##### AGE #####

age_df <- oasis_demog %>% select(RID, AgeatEntry, Cluster) %>% left_join(oasis_dpm, by="RID") %>%
  mutate(onset = AgeatEntry - time_shift)

res_age <- anova_tukey_test(age_df, "AgeatEntry")
res_onset <- anova_tukey_test(age_df, "onset")

##### Education ####

res_educ <- anova_tukey_test(oasis_demog, "EDUC")

##### APOE4 ######
apoe_df <- select(oasis_demog, RID, APOE4, Cluster)
apoe_summary <- apoe_df %>% summarise(.by=c(Cluster, APOE4), n = n()) %>% group_by(Cluster) %>% mutate(prop = n/sum(n))

chi_df_apoe <- cat_chi_test(apoe_df, varname = "APOE4")

##### Diagnosis #####

oasis_dx <- oasis_dx %>% mutate(DX = factor(DX, levels = c("CN", "Impaired", "MCI", "Dementia"), ordered = TRUE),
                                DX.bl = factor(DX.bl, levels = c("CN", "Impaired", "MCI", "Dementia"), ordered = TRUE)) %>%
  group_by(RID) %>% mutate(
    DX.highest = max(DX)
  ) %>% ungroup()


dx_df <- select(oasis_dx, RID, DX.bl, DX.highest, Cluster) %>% distinct() %>%
  group_by(RID) %>%
  mutate(
    CN2CI = ifelse(DX.bl == "CN", DX.highest != "CN", NA),
    MCI2AD = ifelse(DX.bl == "MCI", DX.highest == "Dementia", NA)
  ) %>% ungroup()

chi_df_dx <- cat_chi_test(dx_df, varname = "DX.bl")

dx_summary <- dx_df %>% summarise(.by=c(Cluster, DX.bl), n = n()) %>% group_by(Cluster) %>% mutate(prop = n/sum(n))

dx_transitions <- dx_df %>% select(Cluster, CN2CI, MCI2AD) %>% group_by(Cluster) %>% 
  mutate(nCN2CI = sum(CN2CI, na.rm=T),
         nMCI2AD = sum(MCI2AD, na.rm=T)) %>%
  mutate(CN2CI = sum(CN2CI, na.rm=T)/sum(!is.na(CN2CI)),
         MCI2AD = sum(MCI2AD, na.rm=T)/sum(!is.na(MCI2AD))) %>% distinct()


#### Co-pathologies ######

oasis_cp <- get_copath() %>% inner_join(Clusters, by="RID") %>% distinct()

lbd_summary <- oasis_cp %>% summarise(.by=c(Cluster, lbdis), n = n()) %>% group_by(Cluster) %>% mutate(prop = n/sum(n))
# Only 2 cases - do not include

cv_summary <- oasis_cp %>% summarise(.by=c(Cluster, cvd), n = n()) %>% group_by(Cluster) %>% mutate(prop = n/sum(n))

amndem_summary <- oasis_cp %>% summarise(.by=c(Cluster, amndem), n = n()) %>% group_by(Cluster) %>% mutate(prop = n/sum(n))

mri_data <- load_oasis_mri_data() %>% filter(RID %in% Clusters$RID) %>% left_join(Clusters, by="RID")

#### TDP-43 MRI signature ####
# (Left/Right Middle Temporal + Left/Right Inferior Temporal)/ Hippocampus

tdp.df <- mri_data %>% group_by(RID) %>% filter(Months == max(Months)) %>% ungroup() %>%
  mutate(TDP=(lh_middletemporal_volume+rh_middletemporal_volume+lh_inferiortemporal_volume+rh_inferiortemporal_volume)/(Left.Hippocampus_volume + Right.Hippocampus_volume)) %>% 
  select(RID, TDP, Cluster) %>% drop_na(TDP) %>% distinct(RID, .keep_all = TRUE)

res_tdp <- anova_tukey_test(tdp.df, colname="TDP")

##### AB-PET ######
oasis_ab <- get_oasis_ab() %>% mutate(RID = as.numeric(gsub("^OAS", "", OASISID))) %>% inner_join(Clusters, by="RID")
oasis_ab <- left_join(oasis_ab, oasis_dpm, by="RID") %>% mutate(Years = round(days_to_visit/30.5)/12,
                                                                  Time = Years + time_shift)
oasis_ab <- select(age_df, RID, AgeatEntry) %>% rename(AGE=AgeatEntry) %>% right_join(oasis_ab, by="RID")

ab_summary <- oasis_ab %>% summarise(.by=c(Cluster, tracer), mean = mean(Centiloid_SUVR), sd = sd(Centiloid_SUVR)) 

res_ab <- data.frame()
for (tr in unique(oasis_ab$tracer)) {
  p.ab <- filter(oasis_ab, tracer==tr) %>% get_pairwise_p_values(response = "Centiloid_SUVR", time_var = "Time", by_time = TRUE, correct_for_age = TRUE) 
  
  p.ab <- data.frame(p.ab) %>% separate(contrast, into = c("Cluster1", "Cluster2"), sep = " - ") %>% arrange(Cluster1, Cluster2) %>%
    mutate(p.abpet = round(p.value,3), tracer = tr)
  
  res_ab <- bind_rows(p.ab, res_ab)
}


##### Tau-PET ######
oasis_tau <- get_tau_pet(get_braak = TRUE) %>% mutate(RID = as.numeric(gsub("^OAS", "", OASISID))) %>% left_join(Clusters, by="RID") %>%
  drop_na(Cluster)
oasis_tau <- left_join(oasis_tau, oasis_dpm, by="RID") %>% mutate(Years = round(days_to_visit/30.5)/12,
                                                                  Time = Years + time_shift)

res_tau <- lapply(c("Braak1_2", "Braak3_4", "Braak5_6", "Tauopathy"), 
                  function(col) anova_tukey_test(oasis_tau, col))

##### num MRI ######
n_mri <- mri_data %>% group_by(RID) %>% summarise(n_mri = n()) %>% ungroup() %>% left_join(Clusters, by = "RID") %>%
  group_by(Cluster) %>% summarise(mean = mean(n_mri),
                                  sd = sd(n_mri))

#### GENDER ######
gender <- oasis_demog %>% select(RID, Cluster, GENDER) %>% group_by(Cluster) %>% summarise(Male = round(mean(GENDER == 1, na.rm=T)*100,1),
                                                                     Female = 100-Male)

table_df <- res_age$summary %>% select(Cluster, n) %>% pivot_wider(names_from=Cluster, values_from = n, names_prefix = "mean_") %>% mutate(Variable = "n")
table_df <- res_age$summary %>% select(-n) %>% pivot_wider(names_from=Cluster, values_from = c(mean,sd)) %>% mutate(Variable = "Age (Y)") %>% bind_rows(table_df)
table_df <- res_onset$summary %>% select(-n) %>% pivot_wider(names_from=Cluster, values_from = c(mean,sd)) %>% mutate(Variable = "Onset (Y)") %>% bind_rows(table_df)
table_df <- res_educ$summary %>% select(-n) %>% pivot_wider(names_from=Cluster, values_from = c(mean,sd)) %>% mutate(Variable = "Educat. (Y)") %>% bind_rows(table_df)
table_df <- n_mri %>% pivot_wider(names_from=Cluster, values_from = c(mean,sd)) %>% mutate(Variable = "Num. MRI") %>% bind_rows(table_df)
table_df <- gender %>% select(-Female) %>% pivot_wider(names_from=Cluster, values_from = Male, names_prefix = "mean_") %>% mutate(Variable = "Male (%)") %>% bind_rows(table_df)
table_df <- gender %>% select(-Male) %>% pivot_wider(names_from=Cluster, values_from = Female, names_prefix = "mean_") %>% mutate(Variable = "Female (%)") %>% bind_rows(table_df)

table_df_2 <- apoe_summary %>% drop_na(APOE4) %>% pivot_wider(names_from = Cluster, values_from=c(n, prop)) %>% arrange(APOE4) %>% 
  mutate(Variable = c("APOE e4 (0)", "APOE e4 (1)", "APOE e4 (2)")) %>% select(-APOE4)

table_df_2 <- dx_summary %>% pivot_wider(names_from = Cluster, values_from=c(n, prop)) %>% arrange(DX.bl) %>%
  mutate(Variable=as.character(DX.bl)) %>% select(-DX.bl) %>% bind_rows(table_df_2)

table_df_2 <- dx_transitions %>% select(Cluster, CN2CI, nCN2CI) %>% rename(n=nCN2CI, prop=CN2CI) %>%
  pivot_wider(names_from = Cluster, values_from=c(n, prop)) %>% mutate(Variable = "CN - CI") %>% bind_rows(table_df_2)

table_df_2 <- dx_transitions %>% select(Cluster, MCI2AD, nMCI2AD) %>% rename(n=nMCI2AD, prop=MCI2AD) %>%
  pivot_wider(names_from = Cluster, values_from=c(n, prop)) %>% mutate(Variable = "MCI - AD") %>% bind_rows(table_df_2)

table_df_2 <- filter(cv_summary, cvd==1) %>% pivot_wider(names_from=Cluster, values_from = c(n, prop)) %>% mutate(Variable = "CVD") %>% select(-cvd) %>%
  bind_rows(table_df_2)

table_df <- ab_summary %>% pivot_wider(names_from=Cluster, values_from = c(mean,sd)) %>% 
  mutate(Variable = ifelse(tracer == "AV45", "FBP", "PiB")) %>% select(-tracer) %>% bind_rows(table_df)

table_df <- res_tau[[1]]$summary %>% select(-n) %>% pivot_wider(names_from=Cluster, values_from = c(mean,sd)) %>% mutate(Variable = "Braak I-II") %>% bind_rows(table_df)
table_df <- res_tau[[2]]$summary %>% select(-n) %>% pivot_wider(names_from=Cluster, values_from = c(mean,sd)) %>% mutate(Variable = "Braak III-IV") %>% bind_rows(table_df)
table_df <- res_tau[[3]]$summary %>% select(-n) %>% pivot_wider(names_from=Cluster, values_from = c(mean,sd)) %>% mutate(Variable = "Braak V-VI") %>% bind_rows(table_df)
table_df <- res_tau[[4]]$summary %>% select(-n) %>% pivot_wider(names_from=Cluster, values_from = c(mean,sd)) %>% mutate(Variable = "Cereb. Cortex Ref.") %>% bind_rows(table_df)

res_onset$tukey


#### Make LaTeX tables ####
cs <- as.character(sort(unique(Clusters$Cluster)))
library(knitr)

format_pm <- function(mean, sd) {
  sprintf("$%.1f \\pm %.2f$", mean, sd)
}

for (c in cs) {
  table_df[[c]] <- sprintf("$%.1f \\pm %.2f$", table_df[[paste0("mean_", c)]], table_df[[paste0("sd_", c)]])
  table_df_2[[c]] <- sprintf("$%.1f$\\%% ($n=%.0f$)", 
                             as.numeric(table_df_2[[paste0("prop_", c)]])*100, 
                             as.numeric(table_df_2[[paste0("n_", c)]]))
}

order <- c("n", "Age (Y)", "Onset (Y)", "APOE e4 (0)", "APOE e4 (1)",
           "APOE e4 (2)", "Educat. (Y)", "Num. MRI", "Male (%)", "Female (%)",
           "CN", "Impaired", "MCI", "Dementia", "CN - CI", "MCI - AD",
           "FBP", "PiB",
           "Braak I-II", "Braak III-IV", "Braak V-VI", "Cereb. Cortex Ref.")

latex_table <- select(table_df, Variable, all_of(cs))
latex_table <- select(table_df_2, Variable, all_of(cs)) %>% bind_rows(latex_table) %>%
  arrange(factor(Variable, levels = order))

kable(latex_table, format = "latex", escape = FALSE, booktabs = TRUE)

stats_df <- res_age$tukey %>% select(Cluster1, Cluster2, lwr, upr, p.adj) %>% mutate(Variable = "Age")
stats_df <- res_onset$tukey %>% select(Cluster1, Cluster2, lwr, upr, p.adj) %>% mutate(Variable = "Onset") %>% bind_rows(stats_df)
stats_df <- res_educ$tukey %>% select(Cluster1, Cluster2, lwr, upr, p.adj) %>% mutate(Variable = "Educat.") %>% bind_rows(stats_df)

stats_df <- res_tau[[1]]$tukey %>% select(Cluster1, Cluster2, lwr, upr, p.adj) %>% mutate(Variable = "Braak I-II") %>% bind_rows(stats_df)
stats_df <- res_tau[[2]]$tukey %>% select(Cluster1, Cluster2, lwr, upr, p.adj) %>% mutate(Variable = "Braak III-IV") %>% bind_rows(stats_df)
stats_df <- res_tau[[3]]$tukey %>% select(Cluster1, Cluster2, lwr, upr, p.adj) %>% mutate(Variable = "Braak V-VI") %>% bind_rows(stats_df)
stats_df <- res_tau[[4]]$tukey %>% select(Cluster1, Cluster2, lwr, upr, p.adj) %>% mutate(Variable = "Cereb. Cortex Ref.") %>% bind_rows(stats_df)

stats_df <- res_tdp$tukey %>% select(Cluster1, Cluster2, lwr, upr, p.adj) %>% mutate(Variable = "TDP-34 MRI") %>% bind_rows(stats_df)

stats_df <- res_ab %>% select(Cluster1, Cluster2, p.abpet, tracer) %>% rename(p.adj = p.abpet) %>% mutate(Variable = ifelse(tracer == "AV45", "FBP", "PiB")) %>%
  select(-tracer) %>% bind_rows(stats_df)

stats_df$p.val = round(stats_df$p.adj, 3)

select(stats_df, Variable, p.val) %>% arrange(Variable)

chi_summary <- chi_df_dx %>% mutate(Variable = "Diagnosis")
chi_summary <- chi_df_apoe %>% mutate(Variable = "APOE4") %>% bind_rows(chi_summary)


#### Plots ####
library(ggalluvial)

p <- ggplot(age_df, aes(x = Cluster, y = onset, fill = Cluster)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.2) + 
  geom_point(size = 1, position = position_jitter(width = 0.2, height = 0), aes(color = Cluster)) +
  scale_color_paletteer_d("ggthemes::Tableau_10") +
  scale_fill_paletteer_d("ggthemes::Tableau_10") +
  ylim(40, 92) +
  theme_classic()
p <- p + labs(x = "Clusters", y ="Estimated age of onset")
plot(p)



#Alluvial DX plot
# count transitions
plot_df <- dx_df %>% mutate(DX.bl = recode(DX.bl, "Dementia" = "AD", "Impaired" = "MCI"),
                            DX.highest = recode(DX.highest, "Dementia" = "AD", "Impaired" = "MCI")) %>%
  group_by(Cluster) %>%
  count(DX.bl, DX.highest, name = "Freq")

gp_alluv_list <- lapply(cs, function(c) {
  gp_alluv <- plot_df %>% filter(Cluster == c) %>%
    ggplot(aes(axis1 = DX.bl, axis2 = DX.highest, y = Freq)) +
    geom_alluvium(aes(fill = DX.highest), width = 1/12) +
    geom_stratum(width = 1/8, fill = "grey90", color = "grey40") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("Baseline", "Final diagnosis"), expand = c(.1, .05)) +
    labs(x = NULL, y = "Count", title = paste0("Cluster ", c)) +
    scale_fill_brewer(palette = "YlOrRd", name = "Diagnosis") +
    theme_classic()
})

gp_alluv <- ggarrange(plotlist=gp_alluv_list, ncol=2, nrow=1, common.legend = TRUE, legend="right")
plot(gp_alluv)


tp <- ggplot(tdp.df, aes(x=Cluster, y = TDP, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  #geom_violin(alpha=0.5) +
  geom_point(size = 0.8, position = position_jitter(width = 0.2, height = 0), aes(color=Cluster)) +
  labs(x = "Cluster", y = "TDP-43 MRI ratio") +
  scale_color_paletteer_d("ggthemes::Tableau_10") +
  scale_fill_paletteer_d("ggthemes::Tableau_10") +
  theme_classic()
tp


abpet <- oasis_ab %>% select(RID, Cluster, tracer, Centiloid_SUVR) %>% group_by(RID, tracer) %>% 
  dplyr::summarize(Centiloid_SUVR = max(Centiloid_SUVR),
                   Cluster = Cluster[1], .groups = "drop_last") %>% ungroup() 

gp_fbp <- filter(abpet, tracer == "AV45") %>%
  ggplot(aes(x = Cluster, y = Centiloid_SUVR, fill = Cluster)) + 
  geom_violin(alpha = 0.5) + 
  geom_point(size = 0.8, position = position_jitter(width = 0.2, height = 0), color = "black") +
  scale_color_paletteer_d("ggthemes::Tableau_10") +
  scale_fill_paletteer_d("ggthemes::Tableau_10") +
  geom_hline(yintercept=20.6, linetype="dashed", color="black") +
  theme_classic() + 
  labs(x = "Clusters", y ="FBP SUVR")
plot(gp_fbp)

gp_pib <- filter(abpet, tracer == "PIB") %>%
  ggplot(aes(x = Cluster, y = Centiloid_SUVR, fill = Cluster)) + 
  geom_violin(alpha = 0.5) + 
  geom_point(size = 0.8, position = position_jitter(width = 0.2, height = 0), color = "black") +
  scale_color_paletteer_d("ggthemes::Tableau_10") +
  scale_fill_paletteer_d("ggthemes::Tableau_10") +
  geom_hline(yintercept=16.4, linetype="dashed", color="black") +
  theme_classic() + 
  labs(x = "Clusters", y ="PiB SUVR")
plot(gp_pib)

ab.plot <- ggarrange(gp_fbp, gp_pib, nrow=1, common.legend = TRUE, legend="right")


# Save plots
ggsave(p, filename = paste("~/R/EDAP-data/plots/LTC4/oasis_", run, "_onset_age.png", sep=""), width = 5, height = 4, dpi = 300)
ggsave(gp_alluv, filename = paste("~/R/EDAP-data/plots/LTC4/oasis_", run, "_dx.png", sep=""), width = 8, height = 6, dpi = 300)
ggsave(ab.plot, filename = paste("~/R/EDAP-data/plots/LTC4/oasis_", run, "_ab_violin.png"), width = 9, height =4, dpi=600)
ggsave(tp, filename = paste("~/R/EDAP-data/plots/LTC4/oasis_", run, "_tdp.png", sep=""), width = 5, height = 4, dpi = 300)






