##### BioFINDER

biofinder <- read.csv("~/R/EDAP-data/BioFINDER/data_for_Filip/data_filip.csv")

cortical_vols <- readxl::read_excel("~/R/EDAP-data/BioFINDER/data_for_Filip/cortical_vols_filip.xlsx") %>%
  filter(!is.na(sid))
cortical_vols <- filter(cortical_vols, !if_all(colnames(cortical_vols), is.na)) %>%
  mutate(mri_date = str_extract(csv_icv__index, "\\d{8}"),
         mri_date = ymd(mri_date)) 

counts <- cortical_vols %>% group_by(sid) %>% summarise(n_obs = n()) 
mean(counts$n_obs)
sd(counts$n_obs)

library(ggplot2)

gp1 <- ggplot(counts, aes(x = n_obs)) +
  geom_bar(aes(y = after_stat(count) / sum(after_stat(count))), fill = "steelblue") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = seq(min(counts$n_obs), max(counts$n_obs), by = 1)) +
  labs(x = "n_obs", y = "Percent of total", title = "BioFINDER coritcal vols")

subcortical_vols <- grepv("samseg_vols", colnames(biofinder))
subcortical_df <- biofinder %>% select(sid, Visit, cognitive_status_baseline_variable, diagnosis_baseline_variable,
                                       all_of(subcortical_vols), icv_mm3) %>%
  mutate(mri_date = str_extract(csv_samseg_vols__index, "\\d{8}"),
         mri_date = ymd(mri_date)) 

subcortical_df <- filter(subcortical_df, !if_all(all_of(subcortical_vols[subcortical_vols != "csv_samseg_vols__index"]), is.na))

counts <- subcortical_df %>% group_by(sid) %>% summarise(n_obs = n()) 
mean(counts$n_obs)
sd(counts$n_obs)

gp2 <- ggplot(counts, aes(x = n_obs)) +
  geom_bar(aes(y = after_stat(count) / sum(after_stat(count))), fill = "darkblue") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = seq(min(counts$n_obs), max(counts$n_obs), by = 1)) +
  labs(x = "n_obs", y = "Percent of total", title = "BioFINDER subcoritcal vols")

#### After pre-processing

# Data loading
source("~/R/LTC/utils/biofinder_data_loaders.R")

mri_df <- get_mri_data(normalize = TRUE) %>% select(-OPTICCHIASM)

ab_df <- get_ab_df()

ab_pos_ids <- ab_df %>% group_by(sid) %>% mutate(AB_any = any(AB)) %>% ungroup() %>%
  filter(AB_any) %>% select(sid) %>% unlist()

mri_df <- filter(mri_df, sid %in% ab_pos_ids)

dpm_res <- read.csv("~/R/EDAP-data/BioFINDER/DPM_BioFINDER.csv") %>% distinct(sid, time_shift)

mri_df <- left_join(mri_df, dpm_res, by="sid") %>% filter(!is.na(time_shift)) %>%
  mutate(Time = Years + time_shift) %>%
  rename(RID = sid)

# Count observations
counts <- mri_df %>%
  group_by(RID) %>%
  summarise(n_obs = n())

mean(counts$n_obs)
sd(counts$n_obs)

gp3 <- ggplot(counts, aes(x = n_obs)) +
  geom_bar(aes(y = after_stat(count) / sum(after_stat(count))), fill = "purple") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = seq(min(counts$n_obs), max(counts$n_obs), by = 1)) +
  labs(x = "n_obs", y = "Percent of total", title = "BioFINDER after pre-p.")

#### ADNI+OASIS
multi_cohort_df <- read.csv("~/R/EDAP-data/MULTI_COHORT_4.csv", header = TRUE)

# Filter out NACC
multi_cohort_df <- filter_out(multi_cohort_df, Cohort == "NACC")

all.vars <- c(grepv("^(RH_|LH_|CC_)", colnames(multi_cohort_df)), "BRAINSTEM")

# Count observations
counts <- multi_cohort_df %>%
  group_by(RID) %>%
  summarise(n_obs = n())

mean(counts$n_obs)
sd(counts$n_obs)

gp4 <- ggplot(counts, aes(x = n_obs)) +
  geom_bar(aes(y = after_stat(count) / sum(after_stat(count))), fill = "darkred") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = seq(min(counts$n_obs), max(counts$n_obs), by = 1)) +
  labs(x = "n_obs", y = "Percent of total", title = "ADNI+OASIS")

library(ggpubr)

gp <- ggarrange(gp1, gp2, gp3, gp4, ncol=1)
plot(gp)


