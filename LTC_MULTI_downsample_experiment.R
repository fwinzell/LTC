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

fit_inital = TRUE # set to FALSE to load previous initial model fitting
# 1. Load dataset
#multi_cohort_df_ <- read.csv("~/R/EDAP-data/MULTI_COHORT.csv", header = TRUE)
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


# Count observations
counts <- mri_downsampled %>%
  group_by(RID) %>%
  summarise(n_obs = n())

mean(counts$n_obs)
sd(counts$n_obs)

