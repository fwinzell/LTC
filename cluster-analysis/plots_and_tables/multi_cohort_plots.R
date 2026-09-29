library(ADNIMERGE2)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(tidyverse)

#library(aricode)
#library(infotheo)
library(paletteer)

library(lme4)
library(stats)
library(reshape2)
library(tibble)
library(RColorBrewer)

library(lmerTest)
library(emmeans)

#### Set up #####
adni_dl <- new.env()
source("~/R/LTC/utils/adni_data_loaders.R", local=adni_dl)

oasis_dl <- new.env()
source("~/R/LTC/utils/oasis_data_loaders.R", local=oasis_dl)

source("~/R/LTC/utils/analysis_utils.R")

multi_cohort_df <- read.csv("~/R/EDAP-data/MULTI_COHORT.csv", header = TRUE) %>%
  filter_out(Cohort == "NACC")
all.vars <- c(grepv("^(RH_|LH_|CC_)", colnames(multi_cohort_df)), "BRAINSTEM")

run <- "exp_km_ab_ao"
load(paste("~/R/EDAP-data/LTC_MC/", run, ".Rdata", sep = ""))

Clusters <- data.frame(
  Cluster = multiLTC@Cluster,
  RID = multiLTC@RID
)

multi_cohort_df <- Clusters %>% left_join(multi_cohort_df, by = "RID") %>% drop_na(Cluster) 


#### Demographics ####
adni_demog <- adni_dl$get_demographics() %>%
  select(-GENOTYPE) %>% rename(GENDER = PTGENDER,
                               EDUCATION = PTEDUCAT)
adni_dx <- adni_dl$get_diagnoses()

adni_demog <- adni_dx %>% distinct(RID, DX.bl) %>% right_join(adni_demog, by="RID") %>%
  mutate(RID = paste0('ADNI_', RID))

oasis_demog <- oasis_dl$get_demographics()
oasis_dx <- oasis_dl$get_diagnoses()

oasis_demog <- oasis_dx %>% distinct(RID, DX.bl) %>% right_join(oasis_demog, by="RID") %>%
  select(RID, DX.bl, AgeatEntry, GENDER, EDUC, APOE4) %>% mutate(RID = paste0('OASIS_', RID),
                                                                 GENDER = recode(GENDER, `1` = 'Male', `2` = 'Female')) %>%
  rename(AGE = AgeatEntry,
         EDUCATION = EDUC) %>% drop_na(DX.bl)

mc_demog <- rbind(adni_demog, oasis_demog)

# Age, Age of onset
age_df <- multi_cohort_df %>% distinct(RID, time_shift, Cluster) %>% right_join(mc_demog, by='RID') %>%
  drop_na(Cluster) %>%
  distinct(RID, AGE, time_shift, Cluster) %>% mutate(ONSET = AGE - time_shift)

p <- ggplot(age_df, aes(x = Cluster, y = AGE, fill = Cluster)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.2) + 
  geom_point(size = 1, position = position_jitter(width = 0.2, height = 0), aes(color = Cluster)) +
  scale_color_paletteer_d("ggthemes::Tableau_10") +
  scale_fill_paletteer_d("ggthemes::Tableau_10") +
  ylim(40, 92) +
  theme_classic()
p <- p + labs(x = "Clusters", y ="Age")
plot(p)

po <- ggplot(age_df, aes(x = Cluster, y = ONSET, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  geom_point(size=1, position = position_jitter(width = 0.2, height = 0), aes(color = Cluster)) +
  scale_color_paletteer_d("ggthemes::Tableau_10") +
  scale_fill_paletteer_d("ggthemes::Tableau_10") +
  theme_classic() +
  ylim(40, 92)
po <- po + labs(x = "Clusters", y ="Estimated age of onset")
plot(po)


# Amyloid-beta
adni_ab_df <- adni_dl$get_ab_df() %>% 
  mutate(RID = paste0('ADNI_', RID))
oasis_ab_df <- oasis_dl$get_oasis_ab(use_centiloid = FALSE) %>%
  mutate(RID = gsub('OAS', 'OASIS_', OASISID))

FBP_df <- adni_ab_df %>% filter(TRACER == 'FBP') %>% group_by(RID) %>%
  mutate(FBP = max(SUMMARY_SUVR)) %>% ungroup() %>% distinct(RID, FBP)

FBP_df <- oasis_ab_df %>% filter(tracer == 'AV45') %>% group_by(RID) %>%
  mutate(FBP = max(PET_SUVR)) %>% ungroup() %>% distinct(RID, FBP) %>% rbind(FBP_df)

FBP_df <- left_join(FBP_df, Clusters, by='RID') %>% drop_na(Cluster)

gp.abpet <- FBP_df %>%
  ggplot(aes(x = Cluster, y = FBP, fill = Cluster)) + 
  geom_violin(alpha = 0.5) + 
  geom_point(size = 0.8, position = position_jitter(width = 0.2, height = 0), color = "black") +
  scale_color_paletteer_d("ggthemes::Tableau_10") +
  scale_fill_paletteer_d("ggthemes::Tableau_10") +
  geom_hline(yintercept=1.11, linetype="dashed", color="black") +
  theme_classic() + 
  labs(x = "Clusters", y ="AB-PET")
plot(gp.abpet)

PIB_df <- adni_ab_df %>% filter(TRACER == 'PIB') %>% group_by(RID) %>%
  mutate(PIB = max(SUMMARY_SUVR)) %>% ungroup() %>% distinct(RID, PIB)

PIB_df <- oasis_ab_df %>% filter(tracer == 'PIB') %>% group_by(RID) %>%
  mutate(PIB = max(PET_SUVR)) %>% ungroup() %>% distinct(RID, PIB) %>% rbind(PIB_df)

PIB_df <- left_join(PIB_df, Clusters, by='RID') %>% drop_na(Cluster)

gp.abpet <- PIB_df %>%
  ggplot(aes(x = Cluster, y = PIB, fill = Cluster)) + 
  geom_violin(alpha = 0.5) + 
  geom_point(size = 0.8, position = position_jitter(width = 0.2, height = 0), color = "black") +
  scale_color_paletteer_d("ggthemes::Tableau_10") +
  scale_fill_paletteer_d("ggthemes::Tableau_10") +
  geom_hline(yintercept=1.11, linetype="dashed", color="black") +
  theme_classic() + 
  labs(x = "Clusters", y ="AB-PET")
plot(gp.abpet)


