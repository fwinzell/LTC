


multi_cohort_df <- read.csv("~/R/EDAP-data/MULTI_COHORT_4.csv", header = TRUE)
all.vars <- c(grepv("^(RH_|LH_|CC_)", colnames(multi_cohort_df)), "BRAINSTEM")

adni_dl <- new.env()
source("~/R/LTC/utils/adni_data_loaders.R", local=adni_dl)

oasis_dl <- new.env()
source("~/R/LTC/utils/oasis_data_loaders.R", local=oasis_dl)

adni_copath <- adni_loader$get_copath() %>% 
  mutate(RID = paste0("ADNI_", RID)) %>% 
  mutate(
    WMH_ICV = TOTAL_WMH / CEREBRUM_TCV,
    log_WMH = log1p(TOTAL_WMH),
    log_WMH_ICV = log1p(WMH_ICV)
  ) %>%
  rename(WMH = TOTAL_WMH)

wmh_transformed <- adni_copath %>%
  select(RID, TOTAL_WMH, log_WMH, WMH_ICV, log_WMH_ICV) %>%
  pivot_longer(
    cols = -RID,
    names_to = "Transformation",
    values_to = "WMH"
  )

ggplot(wmh_transformed, aes(x = WMH)) +
  geom_histogram(bins = 40, na.rm = TRUE) +
  facet_wrap(~ Transformation, scales = "free") +
  theme_classic() +
  labs(
    x = "WMH",
    y = "Count"
  )

ggplot(wmh_transformed, aes(sample = WMH)) +
  stat_qq(na.rm = TRUE) +
  stat_qq_line(na.rm = TRUE) +
  facet_wrap(~ Transformation, scales = "free") +
  theme_classic()
 

## OASIS 
oasis_copath <- oasis_loader$get_copath() %>%
  rename(WMH = WM.hypointensities_volume) %>%
  mutate(WMH = WMH / 1e3, # convert to cm^3
         ICV = IntraCranialVol / 1e3,
          WMH_ICV = WMH / ICV,
         log_WMH = log1p(WMH), 
         log_WMH_ICV = log1p(WMH_ICV))

wmh_transformed <- oasis_copath %>%
  select(RID, WMH, log_WMH, WMH_ICV, log_WMH_ICV) %>%
  pivot_longer(
    cols = -RID,
    names_to = "Transformation",
    values_to = "WMH"
  )


oasis_wmh <- select(oasis_copath, RID, Years, WMH, log_WMH, WMH_ICV, log_WMH_ICV) %>%
  mutate(Cohort = "OASIS")
adni_wmh <- select(adni_copath, RID, Years, WMH, log_WMH, WMH_ICV, log_WMH_ICV) %>%
  mutate(Cohort = "ADNI")

wmh_df <- rbind(oasis_wmh, adni_wmh)

ggplot(wmh_df, aes(x = WMH, fill = Cohort)) +
  geom_density(na.rm = TRUE, alpha = 0.5) +
  theme_classic() +
  labs(
    x = "WMH",
    y = "Density"
  )

ggplot(wmh_df, aes(x = WMH_ICV, fill = Cohort)) +
  geom_density(na.rm = TRUE, alpha = 0.5) +
  theme_classic() +
  labs(
    x = "WMH_ICV",
    y = "Density"
  )

ggplot(wmh_df, aes(x = log_WMH, fill = Cohort)) +
  geom_density(na.rm = TRUE, alpha = 0.5) +
  theme_classic() +
  labs(
    x = "log_WMH",
    y = "Density"
  )

ggplot(wmh_df, aes(x = log_WMH_ICV, fill = Cohort)) +
  geom_density(na.rm = TRUE, alpha = 0.5) +
  theme_classic() +
  labs(
    x = "log_WMH_ICV",
    y = "Density"
  )


