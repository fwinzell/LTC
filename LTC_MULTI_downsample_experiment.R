library(ADNIMERGE2)
library(progmod)
library(tidyr)
library(dplyr)
library(tibble)
library(ClusterR)

library(progress)
library(purrr)

library(caret)
library(stringr)

library(ggplot2)
library(ggpubr)
library(lubridate)

library(parallel)
library(foreach)
library(doParallel)
library(fst)

# Useful functions for loading data and plotting trajectories
# Clustering functions
source("~/R/LTC/utils/clustering.R")

# Model fitting
source("~/R/LTC/utils/model_utils.R")

# Extra utils for clustering and visualization
source("~/R/LTC/utils/cluster_utils.R")

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

downsample_visits <- function(df, target_mean, id = "RID", time = "Years",
                              min_visits = 1, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  out <- df
  
  n_subjects <- length(unique(out[[id]]))
  target_n_visits <- round(target_mean * n_subjects)
  n_remove <- nrow(out) - target_n_visits
  
  if (n_remove <= 0)
    return(out)
  
  for (i in seq_len(n_remove)) {
    
    n_visits <- table(out[[id]])
    eligible <- names(n_visits[n_visits > min_visits])
    
    if (length(eligible) == 0) {
      warning("Minimum number of visits reached.")
      break
    }
    
    selected_id <- sample(eligible, 1)
    idx <- which(as.character(out[[id]]) == selected_id)
    
    last_visit <- idx[which.max(out[[time]][idx])]
    out <- out[-last_visit, ]
  }
  
  out
}

mri_downsampled <- downsample_visits(
  multi_cohort_df,
  target_mean = 2.3,
  id = "RID",
  time = "Years",
  min_visits = 2,
  seed = 123
)

table(table(multi_cohort_df$RID))
table(table(mri_downsampled$RID))

## BioFINDER
source("~/R/LTC/utils/biofinder_data_loaders.R")
mri_df <- get_mri_data_updated(normalize = TRUE) #%>% select(-OPTICCHIASM)
ab_df <- get_ab_df()
ab_pos_ids <- ab_df %>% group_by(sid) %>% mutate(AB_any = any(AB)) %>% ungroup() %>%
  filter(AB_any) %>% select(sid) %>% unlist()
mri_df <- filter(mri_df, sid %in% ab_pos_ids)
dpm_res <- read.csv("~/R/EDAP-data/BioFINDER/DPM_BioFINDER.csv") %>% distinct(sid, time_shift)
mri_df <- left_join(mri_df, dpm_res, by="sid") %>% filter(!is.na(time_shift)) %>%
  mutate(Time = Years + time_shift) %>%
  rename(RID = sid)

round(table(table(mri_df$RID))/length(unique(mri_df$RID)), 2)
round(table(table(mri_downsampled$RID))/length(unique(mri_downsampled$RID)), 2)


# Count observations
counts <- mri_downsampled %>%
  group_by(RID) %>%
  summarise(n_obs = n())

mean(counts$n_obs)
sd(counts$n_obs)


match_visit_distribution <- function(df, reference_df,
                                     id = "RID",
                                     time = "Years",
                                     n_subjects = NULL,
                                     seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # ---------------------------------------------------------
  # 1. Visit counts in reference cohort
  # ---------------------------------------------------------
  
  ref_counts <- table(table(reference_df[[id]]))
  ref_prop <- ref_counts / sum(ref_counts)
  
  # ---------------------------------------------------------
  # 2. Select subjects from the target cohort
  # ---------------------------------------------------------
  
  subjects <- unique(df[[id]])
  
  if (is.null(n_subjects)) {
    n_subjects <- length(unique(reference_df[[id]]))
  }
  
  selected_subjects <- sample(subjects, n_subjects)
  
  out <- df[df[[id]] %in% selected_subjects, ]
  
  # ---------------------------------------------------------
  # 3. Current visit counts
  # ---------------------------------------------------------
  
  current_counts <- table(table(out[[id]]))
  
  # Visit-count categories present in reference
  target_counts <- round(ref_prop * n_subjects)
  
  # Correct rounding so total = n_subjects
  target_counts[names(target_counts)[1]] <-
    target_counts[names(target_counts)[1]] +
    (n_subjects - sum(target_counts))
  
  # ---------------------------------------------------------
  # 4. Adjust each participant to the desired distribution
  # ---------------------------------------------------------
  
  # Participants grouped by current number of visits
  visit_n <- table(out[[id]])
  
  for (target_n in names(target_counts)) {
    
    target_n <- as.integer(target_n)
    n_target <- target_counts[as.character(target_n)]
    
    current_subjects <- names(visit_n[visit_n == target_n])
    
    # Already enough
    if (length(current_subjects) >= n_target)
      next
    
    # Need additional subjects with MORE visits
    need <- n_target - length(current_subjects)
    
    donors <- names(visit_n[visit_n > target_n])
    
    if (length(donors) == 0)
      stop("Not enough subjects with more visits to create target distribution.")
    
    donors <- sample(donors, min(length(donors), need))
    
    # Reduce each donor to target_n visits
    for (sid in donors) {
      
      idx <- which(as.character(out[[id]]) == sid)
      
      # Order chronologically
      idx <- idx[order(out[[time]][idx])]
      
      # Keep earliest target_n visits
      keep <- idx[seq_len(target_n)]
      
      out <- out[-setdiff(idx, keep), ]
    }
    
    # Recalculate visit counts
    visit_n <- table(out[[id]])
  }
  
  out
}

mc_matched <- match_visit_distribution(multi_cohort_df, mri_df)

round(table(table(mc_matched$RID))/length(unique(mc_matched$RID)), 2)

