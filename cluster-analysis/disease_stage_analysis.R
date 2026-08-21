library(ADNIMERGE2)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(dplyr)
library(Rtsne)
library(paletteer)
library(tidyr)
library(broom)
library(lme4)
library(lmerTest)

multi_cohort_df <- read.csv("~/R/EDAP-data/MULTI_COHORT.csv", header = TRUE)

run <- "exp_km_ab_ao"
load(paste("~/R/EDAP-data/LTC_MC/", run, ".Rdata", sep = ""))

adni_dl <- new.env()
source("~/R/LTC/utils/adni_data_loaders.R", local=adni_dl)

oasis_dl <- new.env()
source("~/R/LTC/utils/oasis_data_loaders.R", local=oasis_dl)

mc_dl <- new.env()
source("~/R/LTC/utils/multi_cohort_loader.R", local=mc_dl)

source("~/R/LTC/utils/analysis_utils.R")

ab_pos_rids <- c(paste0("ADNI_", adni_dl$get_ab_pos_ids()),
                 gsub("OAS", "OASIS_", oasis_dl$get_ab_pos_ids()))

mc_mri <- mc_dl$multi_cohort_mri()
mri_controls <- mc_mri %>% filter(DX.bl == "CN" & !(RID %in% ab_pos_rids))

Clusters <- data.frame(
  Cluster = multiLTC@Cluster,
  RID = multiLTC@RID
)

Clusters %>% left_join(multi_cohort_df, by = "RID") %>% drop_na(Cluster) -> multi_cohort_df

tran <- lapply(unique(Clusters$Cluster), function(c) {
  multi_cohort_df %>% filter(Cluster == c) %>% select(Time) %>% unlist() %>% quantile(c(0.05, 0.95))
})

tran

multi_cohort_df$Stage <- cut(multi_cohort_df$Time, 
                      breaks = c(-Inf, 0, 5, 10, 12.5, 15, Inf),
                      include.lowest = TRUE)
mri_cols <- c(grepv("^(RH_|LH_|CC_)", colnames(multi_cohort_df)), "BRAINSTEM")

multi_cohort_df %>%
  count(RID, Time) %>%
  filter(n > 1) %>% nrow()

ucsf_data %>%
  count(RID, Months) %>%
  filter(n > 1) %>% nrow()

mc_df_z <- data.frame()

for (varname in mri_cols) {
  mu <- mri_controls %>% select(RID, Months, all_of(varname)) %>% na.omit() %>%
    arrange(Months) %>% distinct(RID, .keep_all = TRUE) %>%
    summarise(Mean = mean(.data[[varname]], na.rm=TRUE),
              SD = sd(.data[[varname]], na.rm=TRUE)) 
  
  #z_col = paste0(varname, "_z")
  df <- multi_cohort_df %>% select(RID, Time, DX.bl, Cluster, Stage, all_of(varname)) %>%
    drop_na(Cluster) %>% mutate(!!varname := (.data[[varname]]-mu$Mean)/mu$SD) 
  
  if (ncol(mc_df_z) == 0) {
    mc_df_z <- df
  } else {
    mc_df_z <- left_join(mc_df_z, df, by=c('RID', 'Time', 'Cluster', 'DX.bl', 'Stage'))
  }
}

stages <- levels(multi_cohort_df$Stage)
df_stage <- multi_cohort_df %>% filter(Stage == stages[4]) %>% 
  mutate(hippocampus = LH_HIPPOCAMPUS + RH_HIPPOCAMPUS)

library(lme4)
library(lmerTest)
                
m0 <- lmer(
  hippocampus ~ Time + (1 | RID),
  data = df_stage,
  REML = FALSE
)

m1 <- lmer(
  hippocampus ~ Time + Cluster + (1 | RID),
  data = df_stage,
  REML = FALSE
)

anova(m0, m1)
               