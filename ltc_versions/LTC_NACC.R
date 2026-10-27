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
source("~/R/LTC/utils/nacc_data_loaders.R")

# Clustering functions
source("~/R/LTC/utils/clustering.R")

# Model fitting
source("~/R/LTC/utils/model_utils.R")

# Extra utils for clustering and visualization
source("~/R/LTC/utils/cluster_utils.R")

# MCDP model fitting
source("~/R/LTC/latent-time-shift/fit_dpm_nacc.R")

plot_raw_mri <- function(varname, nacc_mri_data) {
  df <- nacc_mri_data %>% select(RID, Time, all_of(varname), NACCDXBL) %>%
    rename(t = Time, y = varname)
  
  gp <- ggplot(df, aes(x = t, y = y, color = NACCDXBL)) +
    geom_line(aes(group = RID),  alpha=1) +
    geom_point() +
    labs(title = "", x = "Years (time shifted)", y = varname, color = NULL) +
    scale_color_brewer(palette = "YlOrRd") +
    theme_minimal() +
    ylim(0, NA)
  
  plot(gp)
  
}

#### Initial setup ####
fit_mcdp = FALSE # set to FALSE to load previous MCDP run
fit_inital = FALSE # set to FALSE to load previous initial model fitting

# 1. Load datasets
nacc_mri_scan <- load_nacc_mri()
nacc_mri_scan$RID <- as.numeric(gsub("^NACC", "", nacc_mri_scan$NACCID))

length(unique(nacc_mri_scan$NACCID))

# 2. Fit MCDP model to estimate time-shift
if (fit_mcdp) {
  dpm_df <- fit_dpm()
  write.csv(dpm_df, "~/R/EDAP-data/NACC/DPM_NACC.csv")
} else {
  dpm_df <- read.csv("~/R/EDAP-data/NACC/DPM_NACC.csv", header=TRUE)
}

# Transform MRI data to latent time scale
dpm_df |> select(RID, time_shift) |>
  unique() |> inner_join(nacc_mri_scan, join_by(RID)) -> nacc_mri_scan

rids <- unique(nacc_mri_scan$RID)

nacc_mri_scan$Time = nacc_mri_scan$Years + nacc_mri_scan$time_shift

mri_cols <- c("LEFT_HIPPO", "RIGHT_HIPPO", colnames(nacc_mri_scan)[grep("^[LR]H_[^_]*_GVOL", colnames(nacc_mri_scan))])

# 3. Fit inital NLMMs
beta_df <- data.frame(RID = rids)
for (varname in mri_cols) {
  dsubset <- select(nacc_mri_scan, RID, Time, all_of(varname))
  dsubset <- dsubset %>% rename(y = varname, t = Time) %>% drop_na(y)
  result <- exp_nlmms_fn(varname = varname, dsubset = dsubset)
  if (!is.null(result)) {
    beta_df <- full_join(result, beta_df, by = c("RID"))
  }
}

for (varname in mri_cols) {
  plot_raw_mri(varname, nacc_mri_scan) 
}

numCores <- detectCores()
registerDoParallel(numCores)



if (fit_inital) {
  write_fst(nacc_mri_scan, "~/R/LTC/tmp/ucsf_ab.fst", compress=0)
  system.time(
    foreach(i = seq_along(all.vars)) %dopar% {
      varname = all.vars[i]
      
      dsubset <- read_fst("~/R/LTC/tmp/ucsf_ab.fst", columns = c("RID", "Time", "DX.bl", varname))
      dsubset <- dsubset %>% rename(y = varname, t = Time) %>% drop_na(y)
      result <- exp_nlmms_sample_fn(varname, dsubset, n_samples = 25, verbose=FALSE)
      
      saveRDS(result, file = paste0("~/R/LTC/tmp/out/", varname, ".rds"))
      NULL
    }
  )
  
  results <- lapply(all.vars, function(x) readRDS(paste0("~/R/LTC/tmp/out/", x, ".rds")))
  
  results <- results[!sapply(results, is.null)]
  
  betaList <- lapply(results, `[[`, "betas") 
  beta_df <- purrr::reduce(betaList, full_join, by = "RID")
  
  nlmmBasic <- list(
    betas = beta_df,
    bic = sapply(results, `[[`, "bic"),
    aic = sapply(results, `[[`, "aic"),
    logLikes = sapply(results, `[[`, "logLike")
  )
  
  save(nlmmBasic, file = "~/R/EDAP-data/LTC4/nlmmBasic.Rdata")
} else {
  load("~/R/EDAP-data/LTC4/nlmmBasic.Rdata")
}




