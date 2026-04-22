library(ADNIMERGE2)
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

source("~/R/LTC/utils/adni_data_loaders.R")

source("~/R/LTC/utils/analysis_utils.R")

run <- "exp_km_ab"
load(paste("~/R/EDAP-data/LTC4/", run, ".Rdata", sep = ""))

ucsf_data <- ucsf_longitudinal_all(only_vol = TRUE, filter_n = 0)
dpm_df <- read.csv("~/R/EDAP-data/ADNI/DPM_ADNI.csv", header=TRUE)
ucsf_data <- dpm_df %>% select(RID, time_shift) %>% distinct() %>% right_join(ucsf_data, by="RID")

Clusters <- data.frame(
  Cluster = adniLTC@Cluster,
  RID = adniLTC@RID
)

ucsf_data <- filter(ucsf_data, RID %in% adniLTC@RID)


adni_demog <- ADNIMERGE2::PTDEMOG %>% filter(RID %in% adniLTC@RID) %>%
  mutate(VISDATE = parse_date_time(VISDATE, orders=c("Ymd")),
         PTDOB = parse_date_time(PTDOB, orders=c("m/y")), 
         AGE = interval(PTDOB, VISDATE) %/% years(1)) %>% group_by(RID) %>%
  mutate(AGE = min(AGE, na.rm=TRUE),
         PTEDUCAT = max(PTEDUCAT, na.rm=TRUE)) %>% ungroup() %>%
  select(RID, AGE, PTGENDER, PTEDUCAT) %>% na.omit() %>% distinct() %>% arrange(RID, AGE)

# ----- AGE -----
print("AGE ~ Cluster")
adni_demog |> group_by(RID) |> select(RID, AGE) |> 
  left_join(Clusters, by = "RID") |> drop_na(Cluster) -> age_df

age_summary <- age_df |> group_by(Cluster) |> summarise(mean = mean(AGE), sd = sd(AGE), n = n())

anova_res <- aov(AGE ~ Cluster, data=age_df)
tukey_res <- TukeyHSD(anova_res)
tukey_df <- data.frame(tukey_res$Cluster)
tukey_df <- tukey_df %>% rownames_to_column("Clusters") %>% 
  separate(Clusters, into = c("Cluster2", "Cluster1"), sep = "-") %>% arrange(Cluster1, Cluster2)

# Onset age (time shift)
ucsf_data |> group_by(RID) |> select(RID, time_shift) |> distinct() |> 
  left_join(Clusters, by = "RID") |> drop_na(Cluster) -> onset_df
age_df |> select(RID, AGE) |> inner_join(onset_df, by = "RID", keep = FALSE) -> onset_df
onset_df$AGE <- onset_df$AGE - onset_df$time_shift

onset_summary <- onset_df |> group_by(Cluster) |> summarise(mean = mean(AGE), sd = sd(AGE), n = n())

anova_res <- aov(AGE ~ Cluster, data=onset_df)
tukey_res <- TukeyHSD(anova_res)
tukey_df_2 <- data.frame(tukey_res$Cluster)
tukey_df_2 <- tukey_df_2 %>% rownames_to_column("Clusters") %>% 
  separate(Clusters, into = c("Cluster2", "Cluster1"), sep = "-") %>% arrange(Cluster1, Cluster2)
colnames(tukey_df_2)[3:6] <- paste(colnames(tukey_df)[3:6], ".onset", sep="")

# Demographics
adni_demog <- adni_demog %>% left_join(Clusters, by = "RID") 
demo_summary <- adni_demog %>% group_by(Cluster) %>% summarise(Edu.mean = mean(PTEDUCAT, na.rm=T),
                                                 Edu.sd = sd(PTEDUCAT, na.rm=T),
                                                 Ratio.male = mean(PTGENDER == "Male", na.rm=T),
                                                 .groups = "drop")

# Education stats
anova_res <- aov(PTEDUCAT ~ Cluster, data=adni_demog)
tukey_res <- TukeyHSD(anova_res)
tukey_df_3 <- data.frame(tukey_res$Cluster)
tukey_df_3 <- tukey_df_3 %>% rownames_to_column("Clusters") %>% 
  separate(Clusters, into = c("Cluster2", "Cluster1"), sep = "-") %>% arrange(Cluster1, Cluster2)
colnames(tukey_df_3)[3:6] <- paste(colnames(tukey_df)[3:6], ".educat", sep="")

##### APOE e4 ######
print("APOE ~ Cluster")

apoe_df <- ADNIMERGE2::APOERES %>% filter(RID %in% adniLTC@RID) %>% 
  select(RID, GENOTYPE) %>% 
  mutate(APOE4 = str_count(GENOTYPE, "4")) %>%
  left_join(Clusters, by="RID")
apoe_df$APOE4 <- as.factor(apoe_df$APOE4)

apoe_summary <- data.frame()
for(c in sort(unique(apoe_df$Cluster))) {
  apoe_df %>% filter(Cluster == c) %>% group_by(APOE4) %>% summarise(n = n()) %>% mutate(prop = n/sum(n)) -> summ_df
  summ_df$Cluster <- c
  apoe_summary <- rbind(apoe_summary, summ_df)
}

chi_df_apoe <- cat_chi_test(apoe_df, varname = "APOE4")

##### DX #####
dx.df <- ADNIMERGE2::DXSUM %>% select(RID, VISCODE2, DIAGNOSIS, DXCONFID) %>% drop_na(DIAGNOSIS) %>%
  mutate(DX.bl = ifelse(VISCODE2 == "bl", DIAGNOSIS, NA)) %>%
  group_by(RID) %>% mutate(
    DX.bl = first(na.omit(DX.bl))
  ) %>% ungroup() %>%
  mutate(DIAGNOSIS = factor(DIAGNOSIS, levels = c("CN", "MCI", "Dementia"), ordered = TRUE),
         DX.bl = factor(DX.bl, levels = c("CN", "MCI", "Dementia"), ordered = TRUE)) %>%
  group_by(RID) %>% mutate(
    DX.highest = max(DIAGNOSIS)
  ) %>% ungroup()
  
dx.bl <- select(dx.df, RID, DX.bl, DX.highest) %>% distinct() %>%
  left_join(Clusters, by = "RID") %>% drop_na(Cluster) %>%
  group_by(RID) %>%
  mutate(
    CN2CI = ifelse(DX.bl == "CN", DX.highest != "CN", NA),
    MCI2AD = ifelse(DX.bl == "MCI", DX.highest == "Dementia", NA)
  ) %>% ungroup()

chi_df_dx <- cat_chi_test(dx.bl, varname = "DX.bl")

dx_summary <- data.frame()
for(c in sort(unique(dx.bl$Cluster))) {
  dx.bl %>% filter(Cluster == c) %>% group_by(DX.bl) %>% summarise(n = n()) %>% mutate(prop = n/sum(n)) -> summ_df
  summ_df$Cluster <- c
  dx_summary <- rbind(dx_summary, summ_df)
}

dx_transitions <- dx.bl %>% select(Cluster, CN2CI, MCI2AD) %>% group_by(Cluster) %>% 
  mutate(nCN2CI = sum(CN2CI, na.rm=T),
         nMCI2AD = sum(MCI2AD, na.rm=T)) %>%
  mutate(CN2CI = sum(CN2CI, na.rm=T)/sum(!is.na(CN2CI)),
         MCI2AD = sum(MCI2AD, na.rm=T)/sum(!is.na(MCI2AD))) %>% distinct()


#### WMHs ####
# White matter leisons
wmh_df <- ADNIMERGE2::UCD_WMH
wmh_df$RID <- as.numeric(wmh_df$RID)
wmh_df <- wmh_df |> select(RID, VISCODE2, TOTAL_WMH) |> left_join(Clusters, by = "RID") |>
          drop_na(Cluster) 

visits <- wmh_df$VISCODE2
visits <- gsub("scmri", "0", visits)
visits <- gsub("bl", "0", visits)
visits <- gsub("sc", "0", visits)
visits <- sapply(visits, function(x) gsub("m", "", x))
wmh_df$M <- as.numeric(visits)
wmh_df |> distinct(RID, M, .keep_all = TRUE) |> arrange(RID, M) -> wmh_df

wmh_df <- adni_demog %>% select(RID, AGE) %>% unique() %>% right_join(wmh_df, by="RID")
wmh_df <- ucsf_data |> select(RID, time_shift) |> distinct(RID, time_shift) |> 
          inner_join(wmh_df, by = c("RID" = "RID")) |> mutate(
            Years = M/12,
            Time = Years + time_shift) 

library(lmerTest)
library(stringr)

combined_df <- data.frame()
for (c in sort(unique(wmh_df$Cluster))) {
  wmh_df <- mutate(wmh_df, Cluster = relevel(Cluster, c))
  fit_mm <- lmer(TOTAL_WMH ~ Time * Cluster + AGE + (Time|RID), data = wmh_df)
  summ_df <- summary(fit_mm) |> coef() |> as.data.frame()
  summ_df[grepl("^Cluster[A-Z]", rownames(summ_df)), ] |> select("Pr(>|t|)") -> cluster_df
  colnames(cluster_df) <- c("p (cluster)")
  cluster_df$Cluster1 <- c
  cluster_df$Cluster2 <- rownames(cluster_df) %>% str_remove("Cluster")
  summ_df[grepl("Time:Cluster[A-Z]", rownames(summ_df)), ] |> select("Pr(>|t|)") -> inter_df
  colnames(inter_df) <- c("p (interaction)")
  inter_df$Cluster1 <- c
  inter_df$Cluster2 <- rownames(inter_df) %>% str_remove("Time:Cluster")
  
  combined_df <- rbind(combined_df, left_join(cluster_df, inter_df, by = c("Cluster1", "Cluster2")))
}
# Bonferroni correction for multiple comparisons per cluster
n_comp <- length(unique(wmh_df$Cluster))-1
combined_df <- mutate(combined_df, across(where(is.numeric), ~ round(pmin(. * n_comp,1.0),3) ))


##### SAA+ #####
amp_syn <- read.csv("~/R/EDAP-data/ADNI/AMPRION_ASYN_SAA_04Mar2025.csv", header = TRUE)

amp_syn <- amp_syn %>% select(RID, Result) %>% right_join(Clusters, by = "RID") %>% na.omit() %>%
  filter(Result != "Indeterminate") %>%
  mutate(SAA = ifelse(Result %in% c("Detected-1", "Detected-2"), 1, 0))

saa_summary <- amp_syn %>% group_by(Cluster) %>% summarise(Ratio.detected = mean(SAA),
                                                           n = n(),
                                                          .groups = "drop")
print(saa_summary)


chi_df_saa <- cat_chi_test(amp_syn, varname="SAA")


#### ABETA/AV45 ####

# CSF biomarkers
get_ab_df <- source("~/R/EDAP-R/LTC/get_ab_df.r")$value
ab_df <- get_ab_df() %>% filter(RID %in% adniLTC@RID)

ab_df <- ab_df %>% select(RID, A4240.re, A4240.ms) %>% group_by(RID) %>% 
  summarise(A4240.ms=min(A4240.ms), A4240.re=min(A4240.re)) %>%
  left_join(Clusters, by="RID") %>% drop_na(Cluster)

ab_ia_summary <- ab_df %>% group_by(Cluster) %>% summarise(Mean = mean(A4240.re, na.rm=T),
                                                        SD = sd(A4240.re, na.rm=T),
                                                            .groups = "drop")
ab_ms_summary <- ab_df %>% group_by(Cluster) %>% summarise(Mean = mean(A4240.ms, na.rm=T),
                                                            SD = sd(A4240.ms, na.rm=T),
                                                            .groups = "drop")


#ucsf_data |> select(RID, time_shift) |> distinct(RID, time_shift) |> 
#  inner_join(upennabeta, by = c("RID" = "RID")) |> mutate(Time = M + time_shift) -> upennabeta

ab_df <- get_ab_df() %>% filter(RID %in% adniLTC@RID) %>% left_join(age_df, by="RID") %>%
  mutate(Months = as.numeric(gsub("[^0-9.-]", "", VISCODE2))) %>% 
  mutate(Months = ifelse(is.na(Months), 0, Months),
         Years = Months/12)
ab_df <- ucsf_data %>% select(RID, time_shift) %>% distinct() %>% right_join(ab_df, by="RID") %>% 
  mutate(Time = Years + time_shift)

p.ab.ms <- get_pairwise_p_values(ab_df, "A4240.ms", time_var = "Time", by_time = FALSE)
p.ab.ms <- data.frame(p.ab.ms) %>% separate(contrast, into = c("Cluster1", "Cluster2"), sep = " - ") %>% arrange(Cluster1, Cluster2) %>%
  mutate(p.ab = round(p.value,3))

p.ab.ia <- get_pairwise_p_values(ab_df, "A4240.re", time_var = "Time", by_time = FALSE)
p.ab.ia <- data.frame(p.ab.ia) %>% separate(contrast, into = c("Cluster1", "Cluster2"), sep = " - ") %>% arrange(Cluster1, Cluster2) %>%
  mutate(p.ab = round(p.value,3))

# PET
ab_pet <- ab_df %>% select(RID, time_shift, Years, Time, AGE, VISCODE2, Cluster, TRACER, SUMMARY_SUVR) %>% 
  pivot_wider(names_from = TRACER, values_from = SUMMARY_SUVR, id_cols= c(RID, time_shift, Years, Time, AGE, VISCODE2, Cluster)) %>% 
  select(-`NA`) %>% filter(!(is.na(FBP) & is.na(FBB) & is.na(PIB))) 

ab.pet.surv <- ab_pet %>% select(RID, Cluster, FBP, Time) %>% group_by(RID) %>% 
  mutate(AB = any(FBP > 1.11), 
         Time = max(Time), 
         FBP = max(FBP)) %>% unique() %>% na.omit()
ab_pet_summary <- ab.pet.surv %>% group_by(Cluster) %>% summarise(Mean = mean(FBP, na.rm=T),
                                                                  SD = sd(FBP, na.rm=T),
                                                                  .groups = "drop") 

p.amy.pet <- get_pairwise_p_values(ab_pet, "FBP", time_var = "Time", by_time = TRUE)  
p.amy.pet <- data.frame(p.amy.pet) %>% separate(contrast, into = c("Cluster1", "Cluster2"), sep = " - ") %>% arrange(Cluster1, Cluster2) %>%
  mutate(p.abpet = round(p.value,3))

# AB-status
#get_ab_status <- source("~/R/EDAP-R/LTC/get_ab_status.r")$value
#ab_status <- get_ab_status() %>% filter(RID %in% adniLTC@RID) %>% right_join(Clusters, by = "RID")

#ab_summary <- ab_status %>% group_by(Cluster) %>% 
#  dplyr::summarize(AB_pos=mean(AB, na.rm=T), 
#            N = sum(!is.na(AB))) 

#### TDP-43 MRI signature ####
#ST40CV, ST99CV = Left/Right Middle Temporal 
#ST32CV, ST91CV = Left/Right Inferior Temporal
# divide by hippocampus
tdp.df <- ucsf_data %>% group_by(RID) %>% filter(Months == max(Months)) %>% ungroup() %>%
  mutate(TDP=(ST40CV+ST99CV+ST32CV+ST91CV)/(ST29SV+ST88SV)) %>% select(RID, TDP) %>%
  left_join(Clusters, by = "RID") %>% drop_na(TDP) %>% distinct(RID, .keep_all = TRUE)

anova_res <- aov(TDP ~ Cluster, data=tdp.df)
tukey_res <- TukeyHSD(anova_res)
tukey_df_tdp <- data.frame(tukey_res$Cluster)
tukey_df_tdp <- tukey_df_tdp %>% rownames_to_column("Clusters") %>% 
  separate(Clusters, into = c("Cluster2", "Cluster1"), sep = "-") %>% arrange(Cluster1, Cluster2)

#### Tau-PET ####
get_tau_pet <- source("~/R/EDAP-R/LTC/get_tau_pet.r")$value

ucbtaupet <- get_tau_pet() %>% filter(RID %in% adniLTC@RID) %>%
  select(RID, VISCODE2, Years, CTX_ENTORHINAL_SUVR, META_TEMPORAL_SUVR, CTX_INSULA_SUVR, CTX_SUPERIORFRONTAL_SUVR, CTX_SUPERIORPARIETAL_SUVR) %>%
  left_join(Clusters, by = "RID")

ucbtaupet <- ucsf_data |> select(RID, time_shift) |> distinct(RID, time_shift) |> 
    right_join(ucbtaupet, by = "RID") |> mutate(Time = Years + time_shift) 

table(ucbtaupet$Cluster)

p.tau.pet <- lapply(c("CTX_ENTORHINAL_SUVR", "META_TEMPORAL_SUVR", "CTX_INSULA_SUVR", "CTX_SUPERIORFRONTAL_SUVR", "CTX_SUPERIORPARIETAL_SUVR"), 
                    function(col) get_pairwise_p_values(ucbtaupet, col, time_var = "Time", by_time = FALSE, correct_for_age = FALSE))
p.tau.pet <- data.frame(p.tau.pet) %>% separate(contrast, into = c("Cluster1", "Cluster2"), sep = " - ") %>% arrange(Cluster1, Cluster2) %>%
  rename(p.braak12 = p.value,
         p.braak34 = p.value.1,
         p.insu = p.value.2,
         p.front = p.value.3,
         p.pariet = p.value.4) %>% select(Cluster1, Cluster2, p.braak12, p.braak34, p.insu, p.front, p.pariet)

tau_summary <- ucbtaupet %>% select(RID, CTX_ENTORHINAL_SUVR, META_TEMPORAL_SUVR, CTX_INSULA_SUVR, CTX_SUPERIORFRONTAL_SUVR, CTX_SUPERIORPARIETAL_SUVR) %>%
  group_by(RID) %>% summarise_all(max) %>% left_join(Clusters, by = "RID") %>% select(-RID) %>% 
  group_by(Cluster) %>% summarise(across(everything(), list(Mean = mean, SD = sd)))

n_mri <- ucsf_data %>% group_by(RID) %>% summarise(n_mri = n()) %>% ungroup() %>% left_join(Clusters, by = "RID") %>%
  group_by(Cluster) %>% summarise(mean = mean(n_mri),
                                   sd = sd(n_mri))

#### Combine results ####
table_df <- data.frame("Variable" = c("n", "Age", "Onset", "APOE e4 (0)", "APOE e4 (1)",
                                      "APOE e4 (2)", "Educat.", "Num. MRI", "Gender",
                                      "CN", "MCI", "AD", "CN - MCI/AD", "MCI - AD",
                                      "SAA+", "AB42/AB40 (MS)", "AB42/AB40 (IA)", "AB-PET",
                                      "Braak I-II", "Braak III-IV", "Insula", "Sup. Frontal", "Sup. Parietal"))
for (clust in sort(unique(Clusters$Cluster))) {
  print(clust)
  age_c <- age_summary %>% filter(Cluster == clust)
  onset_c <- onset_summary %>% filter(Cluster == clust)
  apoe_c <- apoe_summary %>% filter(Cluster == clust)
  demo_c <- demo_summary %>% filter(Cluster == clust)
  dx_c <- filter(dx_summary, Cluster == clust)
  dx_t_c <- filter(dx_transitions, Cluster == clust)
  saa_c <- saa_summary %>% filter(Cluster == clust)
  n_mri_c <- n_mri %>% filter(Cluster == clust)
  ab4240_ms_c <- ab_ms_summary %>% filter(Cluster == clust)
  ab4240_ia_c <- ab_ia_summary %>% filter(Cluster == clust)
  abpet_c <- ab_pet_summary %>% filter(Cluster == clust)
  tau_c <- tau_summary %>% filter(Cluster == clust)
  col1 <- setNames(data.frame(c(round(age_c$n[1], 0),
                                round(age_c$mean[1], 2),
                                round(onset_c$mean[1], 2),
                                round(apoe_c$prop, 3)*100,
                                round(demo_c$Edu.mean[1], 2),
                                round(n_mri_c$mean[1], 2),
                                round(demo_c$Ratio.male[1]*100, 1),
                                round(dx_c$prop, 3)*100,
                                round(dx_t_c$CN2CI, 3)*100,
                                round(dx_t_c$MCI2AD, 3)*100,
                                round(saa_c$Ratio.detected*100, 1),
                                #round(ab_c$AB_pos*100, 1),
                                round(ab4240_ms_c$Mean, 3),
                                round(ab4240_ia_c$Mean, 3),
                                round(abpet_c$Mean, 2),
                                round(tau_c$CTX_ENTORHINAL_SUVR_Mean, 3),
                                round(tau_c$META_TEMPORAL_SUVR_Mean, 3),
                                round(tau_c$CTX_INSULA_SUVR_Mean, 3),
                                round(tau_c$CTX_SUPERIORFRONTAL_SUVR_Mean, 3),
                                round(tau_c$CTX_SUPERIORPARIETAL_SUVR_Mean, 3)
                                )), clust)
  col2 <- setNames(data.frame(c(round(age_c$n[1], 0), 
                                round(age_c$sd[1], 2), 
                                round(onset_c$sd[1], 2),
                                round(apoe_c$n, 0),
                                round(demo_c$Edu.sd[1], 2),
                                round(n_mri_c$sd[1], 2),
                                round((1-demo_c$Ratio.male)*100, 1),
                                round(dx_c$n, 0),
                                round(dx_t_c$nCN2CI, 0),
                                round(dx_t_c$nMCI2AD, 0),
                                round(saa_c$n),
                                #round(ab_c$N),
                                round(ab4240_ms_c$SD, 3),
                                round(ab4240_ia_c$SD, 3),
                                round(abpet_c$SD, 2),
                                round(tau_c$CTX_ENTORHINAL_SUVR_SD, 3),
                                round(tau_c$META_TEMPORAL_SUVR_SD, 3),
                                round(tau_c$CTX_INSULA_SUVR_SD, 3),
                                round(tau_c$CTX_SUPERIORFRONTAL_SUVR_SD, 3),
                                round(tau_c$CTX_SUPERIORPARIETAL_SUVR_SD, 3)
                                )), clust)
  table_df <- cbind(table_df, col1, col2)
}

stats_table <- left_join(tukey_df, tukey_df_2, by = c("Cluster1", "Cluster2")) %>% 
  left_join(tukey_df_3, by = c("Cluster1", "Cluster2"))%>% 
  left_join(., chi_df_apoe, by = c("Cluster1", "Cluster2")) %>% 
  left_join(., chi_df_dx, by = c("Cluster1", "Cluster2"), suffix = c(".apoe", ".dx")) %>% 
  select(Cluster1, Cluster2, everything())

biomark_table <- stats_table %>% select(Cluster1, Cluster2) %>% 
  left_join(combined_df, by = c("Cluster1", "Cluster2"))  %>%
  left_join(., chi_df_saa, by = c("Cluster1", "Cluster2")) %>% 
  left_join(tukey_df_tdp, by = c("Cluster1", "Cluster2")) %>%
  #left_join(res.df.r, by = c("Cluster1", "Cluster2")) %>%
  select(Cluster1, Cluster2, everything()) %>%
  mutate_if(is.numeric, round, 3)
biomark_table <- p.ab.ia %>% select(Cluster1, Cluster2, p.ab) %>% right_join(biomark_table, by = c("Cluster1", "Cluster2"))
biomark_table <- p.ab.ms %>% select(Cluster1, Cluster2, p.ab) %>% right_join(biomark_table, by = c("Cluster1", "Cluster2"), suffix = c(".ms", ".ia"))
biomark_table <- p.amy.pet %>% select(Cluster1, Cluster2, p.abpet) %>% right_join(biomark_table, by = c("Cluster1", "Cluster2"))


print(table_df)
print(stats_table)
print(biomark_table)

#### Make LaTeX tables ####
cs <- as.character(sort(unique(Clusters$Cluster)))
library(knitr)

format_pm <- function(mean, sd) {
  sprintf("$%.1f \\pm %.2f$", mean, sd)
}

latex_df <- data.frame(table_df) %>% filter(Variable %in% c("Age", "Onset", "Educat.", "Num. MRI"))
for (col in cs) {
  sd_col <- paste0(col, ".1")
  latex_df[[col]] <- sprintf("$%.1f \\pm %.2f$", latex_df[[col]], latex_df[[sd_col]])
}

# Percentage variables APOE4, SAA+, DX
latex_df_2 <- data.frame(table_df) %>% 
  filter(Variable %in% c("APOE e4 (0)", "APOE e4 (1)", "APOE e4 (2)", 
                         "CN", "MCI", "AD", "CN - MCI/AD", "MCI - AD", "SAA+"))
for (col in cs) {
  sd_col <- paste0(col, ".1")
  latex_df_2[[col]] <- sprintf("$%.1f$\\%% ($n=%.0f$)", as.numeric(latex_df_2[[col]]), as.numeric(latex_df_2[[sd_col]]))
}

latex_df_3 <- data.frame(table_df) %>% filter(Variable == "Gender") 
male <- latex_df_3 %>% select(all_of(cs)) %>% mutate(Variable = "Sex (M)")
latex_df_3 <- latex_df_3 %>% select(all_of(paste0(cs, ".1"))) %>% setNames(cs) %>%
  mutate(Variable = "Sex (F)") %>% rbind(male)
latex_df_3[cs] <- sapply(latex_df_3[cs], function(x) sprintf("$%.1f$\\%%", as.numeric(x)))
  
#for (col in cs) {
#  fem_col <- paste0(col, ".1")
  #latex_df_3[[col]] <- sprintf("($%.1f$\\%%, $%.1f$\\%%)", as.numeric(latex_df_3[[col]]), as.numeric(latex_df_3[[fem_col]]))
#}

latex_df_4 <- data.frame(table_df) %>% 
  filter(Variable %in% c("AB42/AB40 (MS)", "AB42/AB40 (IA)", "AB-PET", "Braak I-II", 
                         "Braak III-IV", "Insula", "Sup. Frontal", "Sup. Parietal"))
for (col in cs) {
  sd_col <- paste0(col, ".1")
  latex_df_4[[col]] <- sprintf("$%.3f \\pm %.3f$", latex_df_4[[col]], latex_df_4[[sd_col]])
}

latex_df <- latex_df[, !grepl("\\.1$", names(latex_df))]
latex_df_2 <- latex_df_2[, !grepl("\\.1$", names(latex_df_2))]
latex_df_4 <- latex_df_4[, !grepl("\\.1$", names(latex_df_4))]

latex_df <- rbind(latex_df, latex_df_3, latex_df_2, latex_df_4) 
colnames(latex_df) <- c("", sprintf("\\textbf{%s} ($n=%.0f$)", cs, table_df[1, cs])) 

kable(latex_df, format = "latex", escape = FALSE, booktabs = TRUE)

stats_table %>% mutate(Cluster = sprintf("%s-%s", Cluster1, Cluster2)) %>%
  select(c(Cluster, lwr, upr, p.adj, lwr.onset, upr.onset, p.adj.onset, p.adj.educat, adj_p.apoe, adj_p.dx)) -> latex_stats

for (col in colnames(latex_stats)[2:ncol(latex_stats)]) {
  if(grepl("upr|lwr", col)) {
    latex_stats[[col]] <- sprintf("$%.2f$", latex_stats[[col]])
  } else {
    latex_stats[[col]] <- sprintf("$%.3f$", latex_stats[[col]])
    }
}

latex_stats <- latex_stats %>% mutate(
  CI_age = sprintf("(%s, %s)", lwr, upr),
  CI_onset = sprintf("(%s, %s)", lwr.onset, upr.onset)
) %>% select(c(Cluster, CI_age, p.adj, CI_onset, p.adj.onset, p.adj.educat, adj_p.apoe, adj_p.dx))

kable(latex_stats, format = "latex", escape = FALSE, booktabs = TRUE)

latex_biomk <- biomark_table %>% mutate(Cluster = sprintf("%s-%s", Cluster1, Cluster2)) %>%
  select(c(Cluster, p.ab.ia, p.ab.ms, p.abpet, adj_p, `p (cluster)`, `p (interaction)`, lwr, upr, p.adj)) 

for (col in colnames(latex_biomk)[2:ncol(latex_biomk)]) {
  if(grepl("upr|lwr", col)) {
    latex_biomk[[col]] <- sprintf("$%.2f$", latex_biomk[[col]])
  } else {
    latex_biomk[[col]] <- sprintf("$%.3f$", latex_biomk[[col]])
  }
}

latex_biomk <- latex_biomk  %>% mutate(CI_tdp = sprintf("(%s, %s)", lwr, upr)) %>% 
  select(c(Cluster, p.ab.ia, p.ab.ms,  p.abpet, adj_p, `p (cluster)`, `p (interaction)`, CI_tdp, p.adj)) 

kable(latex_biomk, format = "latex", escape = FALSE, booktabs = TRUE)

latex_tau_pet <- p.tau.pet %>% mutate(Cluster = sprintf("%s-%s", Cluster1, Cluster2)) %>% 
  select(Cluster, p.braak12, p.braak34, p.insu, p.front, p.pariet)
for (col in colnames(latex_tau_pet)[2:ncol(latex_tau_pet)]) {
  latex_tau_pet[[col]] <- sprintf("$%.3f$", latex_tau_pet[[col]])
}
kable(latex_tau_pet, format = "latex", escape = FALSE, booktabs = TRUE)


#library(writexl)
#write_xlsx(table_df, paste("~/R/EDAP-data/", kml3d.run, "_table.xlsx"))
#write_xlsx(stats_table, paste("~/R/EDAP-data/", kml3d.run, "_stats.xlsx"))

#### Plots ####
library(ggalluvial)

p <- ggplot(age_df, aes(x = Cluster, y = AGE, fill = Cluster)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.2) + 
  geom_point(size = 1, position = position_jitter(width = 0.2, height = 0), aes(color = Cluster)) +
  scale_color_paletteer_d("ggthemes::Tableau_10") +
  scale_fill_paletteer_d("ggthemes::Tableau_10") +
  ylim(40, 92) +
  theme_classic()
p <- p + labs(x = "Clusters", y ="Age")
plot(p)

po <- ggplot(onset_df, aes(x = Cluster, y = AGE, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  geom_point(size=1, position = position_jitter(width = 0.2, height = 0), aes(color = Cluster)) +
  scale_color_paletteer_d("ggthemes::Tableau_10") +
  scale_fill_paletteer_d("ggthemes::Tableau_10") +
  theme_classic() +
  ylim(40, 92)
po <- po + labs(x = "Clusters", y ="Estimated age of onset")
plot(po)

#Alluvial DX plot
# count transitions
plot_df <- dx.bl %>% mutate(DX.bl = recode(DX.bl, "Dementia" = "AD"),
                            DX.highest = recode(DX.highest, "Dementia" = "AD")) %>%
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

gp_alluv <- ggarrange(plotlist=gp_alluv_list, ncol=2, nrow=2, common.legend = TRUE, legend="right")
plot(gp_alluv)

wmh_df$Cluster <- factor(wmh_df$Cluster, levels = sort(levels(wmh_df$Cluster)))
gp_lm <- ggplot(wmh_df, aes(x = Time, y = TOTAL_WMH, color = Cluster, fill = Cluster)) +
  geom_smooth(method = "lm", se = TRUE) +  # se = TRUE adds confidence intervals
  labs(x = "Years (shifted)", y = "White Matter Hyperintensity", color = "Cluster", fill = "Cluster") +
  theme(legend.position = "right") +
  scale_color_paletteer_d("ggthemes::Tableau_10") +
  scale_fill_paletteer_d("ggthemes::Tableau_10") +
  theme_classic() +
  scale_x_continuous(limits = c(0, NA), expand = c(0,0))
plot(gp_lm)

gp_wmh <- ggplot(wmh_df, aes(x=Time, y=TOTAL_WMH, colour = Cluster, fill=Cluster, group=RID)) +
  geom_point(shape=21, color="black") +
  geom_line() +
  facet_wrap(~Cluster) +
  labs(x = "Years (shifted)", y = "White Matter Hyperintensity", color = "Cluster", fill = "Cluster") +
  scale_color_paletteer_d("ggthemes::Tableau_10") +
  scale_fill_paletteer_d("ggthemes::Tableau_10") +
  theme_classic() 
plot(gp_wmh)


tp <- ggplot(tdp.df, aes(x=Cluster, y = TDP, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  #geom_violin(alpha=0.5) +
  geom_point(size = 0.8, position = position_jitter(width = 0.2, height = 0), aes(color=Cluster)) +
  labs(x = "Cluster", y = "TDP-43 MRI ratio") +
  scale_color_paletteer_d("ggthemes::Tableau_10") +
  scale_fill_paletteer_d("ggthemes::Tableau_10") +
  theme_classic()
tp

sap <- ggplot(saa_summary, aes(x=Cluster, y=Ratio.detected, fill=Cluster)) +
  geom_col() +
  labs(x = "Cluster", y = "SAA+ (%)") + 
  scale_color_paletteer_d("ggthemes::Tableau_10") +
  scale_fill_paletteer_d("ggthemes::Tableau_10") +
  theme_classic()
sap
  

gp.ab.ms <- ggplot(ab_df, aes(x = Cluster, y = A4240.ms, fill = Cluster)) + 
  geom_violin(alpha = 0.5) + 
  geom_point(size = 0.8, position = position_jitter(width = 0.2, height = 0), color = "black") +
  scale_color_paletteer_d("ggsci::lanonc_lancet") +
  scale_fill_paletteer_d("ggsci::lanonc_lancet") +
  geom_hline(yintercept=0.138, linetype="dashed", color="black") +
  theme_classic() + 
  labs(x = "Clusters", y ="AB42/AB40 (MS)")
plot(gp.ab.ms)

gp.ab.ia <- ggplot(ab_df, aes(x = Cluster, y = A4240.re, fill = Cluster)) + 
  geom_violin(alpha = 0.5) + 
  geom_point(size = 0.8, position = position_jitter(width = 0.2, height = 0), color = "black") +
  scale_color_paletteer_d("ggthemes::Tableau_10") +
  scale_fill_paletteer_d("ggthemes::Tableau_10") +
  geom_hline(yintercept=0.066, linetype="dashed", color="black") +
  theme_classic() + 
  labs(x = "Clusters", y ="AB42/AB40 (IA)")
plot(gp.ab.ia)

ab_pet %>% select(RID, Cluster, FBP) %>% group_by(RID) %>% dplyr::summarize(FBP = max(FBP))
gp.abpet <- ab_pet %>% select(RID, Cluster, FBP) %>% group_by(RID) %>% dplyr::summarize(FBP = max(FBP),
                                                                                  Cluster = Cluster[1]) %>%
  ggplot(aes(x = Cluster, y = FBP, fill = Cluster)) + 
  geom_violin(alpha = 0.5) + 
  geom_point(size = 0.8, position = position_jitter(width = 0.2, height = 0), color = "black") +
  scale_color_paletteer_d("ggthemes::Tableau_10") +
  scale_fill_paletteer_d("ggthemes::Tableau_10") +
  geom_hline(yintercept=1.11, linetype="dashed", color="black") +
  theme_classic() + 
  labs(x = "Clusters", y ="AB-PET")
plot(gp.abpet)

ab.plot <- ggarrange(gp.ab.ia, gp.abpet, nrow=1, common.legend = TRUE, legend="right")


# Save plots
ggsave(p, filename = paste("~/R/EDAP-data/plots/LTC4/", run, "_age.png", sep=""), width = 5, height = 4, dpi = 300)
ggsave(po, filename = paste("~/R/EDAP-data/plots/LTC4/", run, "_onset.png", sep=""), width = 5, height = 4, dpi = 300)
ggsave(gp_alluv, filename = paste("~/R/EDAP-data/plots/LTC4/", run, "_dx.png", sep=""), width = 8, height = 6, dpi = 300)
ggsave(gp_wmh, filename = paste("~/R/EDAP-data/plots/LTC4/", run, "_wmh.png", sep=""), width = 7, height = 5, dpi = 300)
ggsave(ab.plot, filename = paste("~/R/EDAP-data/plots/LTC4/", run, "_ab_violin.png"), width = 9, height =4, dpi=600)

plot_of_ages <- ggarrange(p, po, nrow = 1, common.legend = TRUE, legend = "right")
plot(plot_of_ages)
ggsave(plot_of_ages, 
       filename = paste("~/R/EDAP-data/plots/", run, "_p_ages_row.png", sep=""), 
                        width = 9, height = 4, dpi = 500)

ggsave(tp, filename = paste("~/R/EDAP-data/plots/LTC4/", run, "_tdp.png", sep=""), width = 7, height = 5, dpi = 300)
ggsave(sap, filename = paste("~/R/EDAP-data/plots/LTC4/", run, "_saa.png", sep=""), width = 6.5, height = 5, dpi = 300)

