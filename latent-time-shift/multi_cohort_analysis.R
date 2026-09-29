library(dplyr)
library(ggplot2)
library(Rtsne)
library(cluster)
library(mclust)

multi_cohort_df <- read.csv("~/R/EDAP-data/MULTI_COHORT_4.csv", header = TRUE)

# Filter out NACC
multi_cohort_df <- filter_out(multi_cohort_df, Cohort == "NACC")

all.vars <- c(grepv("^(RH_|LH_|CC_)", colnames(multi_cohort_df)), "BRAINSTEM")


# Data loading
source("~/R/LTC/utils/biofinder_data_loaders.R")

biofinder_df <- get_mri_data_updated(normalize = TRUE) #%>% select(-OPTICCHIASM)

ab_df <- get_ab_df()

ab_pos_ids <- ab_df %>% group_by(sid) %>% mutate(AB_any = any(AB)) %>% ungroup() %>%
  filter(AB_any) %>% select(sid) %>% unlist()

biofinder_df <- filter(biofinder_df, sid %in% ab_pos_ids)

dpm_res <- read.csv("~/R/EDAP-data/BioFINDER/DPM_BioFINDER.csv") %>% distinct(sid, time_shift)

biofinder_df <- left_join(biofinder_df, dpm_res, by="sid") %>% filter(!is.na(time_shift)) %>%
  mutate(Time = Years + time_shift,
         Cohort = "BF2") %>%
  rename(RID = sid)

bf.vars <- c(grepv("^(RH_|LH_|CC_)", colnames(biofinder_df)), "BRAINSTEM", "OPTICCHIASM")
mri.vars <- intersect(all.vars, bf.vars)

common_cols <- intersect(names(multi_cohort_df), names(biofinder_df))

multi_cohort_df <- rbind(multi_cohort_df[, common_cols], biofinder_df[, common_cols])


# Have a look at time-shift
# Hmm this is a problem - need to control for diagnosis?
# Move this to a new script
ts <- distinct(multi_cohort_df, RID, time_shift, Cohort, DX.bl)

ggplot(ts, aes(x=time_shift, fill = Cohort)) +
  geom_histogram(alpha = 0.5) +
  theme_classic()

ts %>% group_by(Cohort) %>% summarise(mean = mean(time_shift),
                                      median = median(time_shift))

aov_result <- aov(time_shift ~ Cohort, data = ts)
summary(aov_result)

kruskal.test(time_shift ~ Cohort, data = ts)

pairwise.wilcox.test(ts$time_shift, ts$Cohort, p.adjust.method = "BH")



