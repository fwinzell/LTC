library(ADNIMERGE2)
library(progmod)
library(tidyr)
library(dplyr)
library(tibble)

library(progress)
library(purrr)

library(caret)
library(stringr)

library(ggplot2)
library(ggpubr)
library(lubridate)

multi_cohort <- function() {
# This merges ADNI, NACC and OASIS cohorts for LTC modeling and clustering
# ADNI: n=561 
# NACC: n=325 
# OASIS: n=320 -> 1206 subjects 

### ADNI ####

# Useful functions for loading data and plotting trajectories
source("~/R/LTC/utils/adni_data_loaders.R")

# Clustering functions
source("~/R/LTC/utils/clustering.R")

# Model fitting
source("~/R/LTC/utils/model_utils.R")

# Extra utils for clustering and visualization
source("~/R/LTC/utils/cluster_utils.R")

# 1. Load datasets
ucsf_all <- ucsf_longitudinal_all(only_vol=TRUE, filter_n=0)

# Find AB positives
ab_df <- get_ab_df()

# Any positive Abeta biomarker at any time point
ab_df %>% select(RID, AB_any) %>% group_by(RID) %>% 
  mutate(AB_any = if_any(AB_any)) %>% distinct() %>% filter(AB_any) %>% 
  select(RID) %>% unlist() -> ab_pos_rids_any

ucsf_ab <- ucsf_all %>% filter(RID %in% ab_pos_rids_any)

length(unique(ucsf_ab$RID))

# 2. Load MCDP model time-shift estimates
adni_dpm_df <- read.csv("~/R/EDAP-data/ADNI/DPM_ADNI.csv", header=TRUE)

# Transform MRI data to latent time scale
adni_dpm_df |> select(RID, time_shift) |>
  unique() |> inner_join(ucsf_ab, join_by(RID)) -> ucsf_ab

ucsf_ab$Time = ucsf_ab$Years + ucsf_ab$time_shift

adni_idxs <- grepl("^ST\\d+(CV|SV|SA|TA)$", colnames(ucsf_ab)) 
adni_mri_vars <- colnames(ucsf_ab)[adni_idxs == 1]

#### NACC ####
# Useful functions for loading data and plotting trajectories
source("~/R/LTC/utils/nacc_data_loaders.R")

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

# 1. Load datasets
nacc_mri_scan <- load_nacc_mri()
nacc_mri_scan$RID <- as.numeric(gsub("^NACC", "", nacc_mri_scan$NACCID))

length(unique(nacc_mri_scan$NACCID))

# 2. Load MCDP model time-shift estimates
nacc_dpm_df <- read.csv("~/R/EDAP-data/NACC/DPM_NACC.csv", header=TRUE)

# Transform MRI data to latent time scale
nacc_dpm_df |> select(RID, time_shift) |>
  unique() |> inner_join(nacc_mri_scan, join_by(RID)) -> nacc_mri_scan

rids <- unique(nacc_mri_scan$RID)

nacc_mri_scan$Time = nacc_mri_scan$Years + nacc_mri_scan$time_shift

nacc_mri_vars <- c(colnames(nacc_mri_scan)[grep("^(LEFT|RIGHT|CC)_", colnames(nacc_mri_scan))], 
                   colnames(nacc_mri_scan)[grep("^[LR]H_[^_]*_GVOL", colnames(nacc_mri_scan))],
                   "BRAIN_STEM")

#### OASIS ####
# Useful functions for loading data and plotting trajectories
source("~/R/LTC/utils/oasis_data_loaders.R")

# MCDP model fitting
source("~/R/LTC/latent-time-shift/fit_dpm_oasis.R")

# 1. Load datasets
oasis_mri_data <- load_oasis_mri_data(unified_norm=FALSE)
#mri_data$RID <- as.numeric(gsub("^OAS", "", mri_data$OASISID))
oasis_mri_data$DX.bl <- factor(oasis_mri_data$DX.bl, labels = c("CN", "Impaired", "MCI", "Dementia"))

ab_pos_ids <- get_ab_pos_ids()
oasis_mri_data <- filter(oasis_mri_data, OASISID %in% ab_pos_ids)

length(unique(oasis_mri_data$OASISID))

# 2. Load MCDP model time-shift estimates
oasis_dpm_df <- read.csv("~/R/EDAP-data/OASIS/DPM_OASIS.csv", header=TRUE)


# Transform MRI data to latent time scale
oasis_dpm_df |> select(RID, time_shift) |>
  unique() |> inner_join(oasis_mri_data, join_by(RID)) -> oasis_mri_data

rids <- unique(oasis_mri_data$RID)

oasis_mri_data$Time = oasis_mri_data$Years + oasis_mri_data$time_shift

#oasis_mri_vars <- select(oasis_mri_data, "lh_bankssts_volume":"Right.Thalamus.Proper_volume") %>% colnames()
oasis_mri_vars <- c(grepv("^[lr]h_[a-z]+_volume$", colnames(oasis_mri_data)), 
              grepv("^(Left|Right)\\.[a-zA-Z._]+_volume$", colnames(oasis_mri_data)),
              grepv("^CC_[a-zA-Z._]+_volume$", colnames(oasis_mri_data)),
              "Brain.Stem_volume")


#### Merging ####

nacc_mri_scan <- nacc_mri_scan %>% rename(DX.bl = NACCDXBL) %>%
  select(c("RID", "time_shift", "DX.bl", "Months", "Years", "Time", all_of(nacc_mri_vars)))

oasis_mri_data <- oasis_mri_data %>% 
  select(c("RID", "time_shift", "DX.bl", "Months", "Years", "Time", all_of(oasis_mri_vars)))

# RID
ucsf_ab$RID <- paste0("ADNI_", ucsf_ab$RID)
nacc_mri_scan$RID <- paste0("NACC_", nacc_mri_scan$RID)
oasis_mri_data$RID <- paste0("OASIS_", oasis_mri_data$RID)

# MRI variables
datadic <- datadict_as_tibble(ADNIMERGE2::DATADIC) %>% filter(TBLNAME == "UCSFFSL") %>%
  select(FLDNAME, TEXT) %>% filter(grepl("ST\\d+[A-Z]{2}", FLDNAME)) %>%
  mutate(Region = sub(".*\\s", "", TEXT))

convert_name <- function(x) {
  words <- regmatches(x, gregexpr("[A-Z][a-z]*", x))[[1]]
  
  if (length(words) == 0) {
    return(toupper(x))
  }
  
  if (words[1] == "Right") {
    region <- toupper(paste0(words[-1], collapse = ""))
    paste0("RH_", region)
  } else if (words[1] == "Left") {
    region <- toupper(paste0(words[-1], collapse = ""))
    paste0("LH_", region)
  } else {
    toupper(paste0(words, collapse = ""))
  }
}

datadic$Formatted <- sapply(datadic$Region, convert_name)

lookup <- setNames(datadic$Formatted, datadic$FLDNAME)
names(ucsf_ab) <- ifelse(names(ucsf_ab) %in% names(lookup), 
                         lookup[names(ucsf_ab)], 
                         names(ucsf_ab))
names(ucsf_ab) <- gsub("CORPUSCALLOSUM", "CC_", names(ucsf_ab))

names(nacc_mri_scan) <- gsub("_GVOL", "", names(nacc_mri_scan))
names(nacc_mri_scan) <- gsub("LEFT_", "LH_", names(nacc_mri_scan))
names(nacc_mri_scan) <- gsub("RIGHT_", "RH_", names(nacc_mri_scan))
names(nacc_mri_scan) <- gsub("_PROPER", "", names(nacc_mri_scan))

nacc_mri_scan <- rename(nacc_mri_scan,
                        RH_ACCUMBENSAREA=RH_ACCUMBENS_AREA,
                        LH_ACCUMBENSAREA=LH_ACCUMBENS_AREA,
                        RH_CEREBELLUMCORTEX=RH_CEREBELLUM_CORTEX,
                        LH_CEREBELLUMCORTEX=LH_CEREBELLUM_CORTEX,
                        RH_CEREBELLUMWM=RH_CEREBELLUM_WHITE_MATTER,
                        LH_CEREBELLUMWM=LH_CEREBELLUM_WHITE_MATTER,
                        BRAINSTEM=BRAIN_STEM)
cc_name <- function(x) {
  if (grepl("^CC", x)) {
    x <- gsub("_", "", x)
    gsub("^CC", "CC_", x)
  } else {
    x
  }
}

names(nacc_mri_scan) <- sapply(names(nacc_mri_scan), cc_name)

#oasis_conv <- lapply(names(oasis_mri_data), convert_name)
names(oasis_mri_data) <- sapply(names(oasis_mri_data), function(name) {
  if (name %in% oasis_mri_vars) {
    name <- convert_name(name)
    name <- cc_name(name)
  }
  name
})
names(oasis_mri_data) <- gsub("_VOLUME", "", names(oasis_mri_data))
oasis_mri_data <- rename(oasis_mri_data,
                         LH_ACCUMBENSAREA=LH_ACCUMBENS,
                         RH_ACCUMBENSAREA=RH_ACCUMBENS,
                         LH_THALAMUS=LH_THALAMUSPROPER,
                         RH_THALAMUS=RH_THALAMUSPROPER,
                         RH_CEREBELLUMWM=RH_CEREBELLUMWHITEMATTER,
                         LH_CEREBELLUMWM=LH_CEREBELLUMWHITEMATTER)

ucsf_ab$Cohort <- "ADNI"
nacc_mri_scan$Cohort <- "NACC"
oasis_mri_data$Cohort <- "OASIS"

all_columns <- base::intersect(base::intersect(colnames(ucsf_ab), colnames(nacc_mri_scan)), 
                               colnames(oasis_mri_data))

missing_adni <- setdiff(colnames(ucsf_ab), all_columns)
missing_nacc <- setdiff(colnames(nacc_mri_scan), all_columns)
missing_oasis <- setdiff(colnames(oasis_mri_data), all_columns)

# Missing Ok

multi_cohort_df <- rbind(
  select(ucsf_ab, all_of(all_columns)),
  select(nacc_mri_scan, all_of(all_columns)),
  select(oasis_mri_data, all_of(all_columns))
)


# Have a look

#ggplot(multi_cohort_df, aes(x = RH_HIPPOCAMPUS, fill = Cohort)) +
#  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
#  theme_minimal()

# Save as .csv
write.csv(multi_cohort_df, file = "~/R/EDAP-data/MULTI_COHORT_2.csv", row.names = FALSE)
}


multi_cohort_mri <- function() {
  # This merges ADNI, NACC and OASIS MRI data
  ### ADNI ####
  
  # Useful functions for loading data and plotting trajectories
  source("~/R/LTC/utils/adni_data_loaders.R")
  
  # Clustering functions
  source("~/R/LTC/utils/clustering.R")
  
  # Model fitting
  source("~/R/LTC/utils/model_utils.R")
  
  # Extra utils for clustering and visualization
  source("~/R/LTC/utils/cluster_utils.R")
  
  # 1. Load datasets
  ucsf_all <- ucsf_longitudinal_all(only_vol=TRUE, filter_n=0)
  
  adni_idxs <- grepl("^ST\\d+(CV|SV|SA|TA)$", colnames(ucsf_all)) 
  adni_mri_vars <- colnames(ucsf_all)[adni_idxs == 1]
  
  #### NACC ####
  # Useful functions for loading data and plotting trajectories
  source("~/R/LTC/utils/nacc_data_loaders.R")
  
  # 1. Load datasets
  nacc_mri_scan <- load_nacc_mri()
  nacc_mri_scan$RID <- as.numeric(gsub("^NACC", "", nacc_mri_scan$NACCID))
  
  nacc_mri_vars <- c(colnames(nacc_mri_scan)[grep("^(LEFT|RIGHT|CC)_", colnames(nacc_mri_scan))], 
                     colnames(nacc_mri_scan)[grep("^[LR]H_[^_]*_GVOL", colnames(nacc_mri_scan))],
                     "BRAIN_STEM")
  
  #### OASIS ####
  # Useful functions for loading data and plotting trajectories
  source("~/R/LTC/utils/oasis_data_loaders.R")
  
  # 1. Load datasets
  oasis_mri_data <- load_oasis_mri_data(unified_norm=FALSE)
  #mri_data$RID <- as.numeric(gsub("^OAS", "", mri_data$OASISID))
  oasis_mri_data$DX.bl <- factor(oasis_mri_data$DX.bl, labels = c("CN", "Impaired", "MCI", "Dementia"))

  
  #oasis_mri_vars <- select(oasis_mri_data, "lh_bankssts_volume":"Right.Thalamus.Proper_volume") %>% colnames()
  oasis_mri_vars <- c(grepv("^[lr]h_[a-z]+_volume$", colnames(oasis_mri_data)), 
                      grepv("^(Left|Right)\\.[a-zA-Z._]+_volume$", colnames(oasis_mri_data)),
                      grepv("^CC_[a-zA-Z._]+_volume$", colnames(oasis_mri_data)),
                      "Brain.Stem_volume")
  
  
  #### Merging ####
  
  nacc_mri_scan <- nacc_mri_scan %>% rename(DX.bl = NACCDXBL) %>%
    select(c("RID", "DX.bl", "Months", "Years", all_of(nacc_mri_vars)))
  
  oasis_mri_data <- oasis_mri_data %>% 
    select(c("RID", "DX.bl", "Months", "Years", all_of(oasis_mri_vars)))
  
  # RID
  ucsf_all$RID <- paste0("ADNI_", ucsf_all$RID)
  nacc_mri_scan$RID <- paste0("NACC_", nacc_mri_scan$RID)
  oasis_mri_data$RID <- paste0("OASIS_", oasis_mri_data$RID)
  
  # MRI variables
  datadic <- datadict_as_tibble(ADNIMERGE2::DATADIC) %>% filter(TBLNAME == "UCSFFSL") %>%
    select(FLDNAME, TEXT) %>% filter(grepl("ST\\d+[A-Z]{2}", FLDNAME)) %>%
    mutate(Region = sub(".*\\s", "", TEXT))
  
  convert_name <- function(x) {
    words <- regmatches(x, gregexpr("[A-Z][a-z]*", x))[[1]]
    
    if (length(words) == 0) {
      return(toupper(x))
    }
    
    if (words[1] == "Right") {
      region <- toupper(paste0(words[-1], collapse = ""))
      paste0("RH_", region)
    } else if (words[1] == "Left") {
      region <- toupper(paste0(words[-1], collapse = ""))
      paste0("LH_", region)
    } else {
      toupper(paste0(words, collapse = ""))
    }
  }
  
  datadic$Formatted <- sapply(datadic$Region, convert_name)
  
  lookup <- setNames(datadic$Formatted, datadic$FLDNAME)
  names(ucsf_all) <- ifelse(names(ucsf_all) %in% names(lookup), 
                           lookup[names(ucsf_all)], 
                           names(ucsf_all))
  names(ucsf_all) <- gsub("CORPUSCALLOSUM", "CC_", names(ucsf_all))
  
  names(nacc_mri_scan) <- gsub("_GVOL", "", names(nacc_mri_scan))
  names(nacc_mri_scan) <- gsub("LEFT_", "LH_", names(nacc_mri_scan))
  names(nacc_mri_scan) <- gsub("RIGHT_", "RH_", names(nacc_mri_scan))
  names(nacc_mri_scan) <- gsub("_PROPER", "", names(nacc_mri_scan))
  
  nacc_mri_scan <- rename(nacc_mri_scan,
                          RH_ACCUMBENSAREA=RH_ACCUMBENS_AREA,
                          LH_ACCUMBENSAREA=LH_ACCUMBENS_AREA,
                          RH_CEREBELLUMCORTEX=RH_CEREBELLUM_CORTEX,
                          LH_CEREBELLUMCORTEX=LH_CEREBELLUM_CORTEX,
                          RH_CEREBELLUMWM=RH_CEREBELLUM_WHITE_MATTER,
                          LH_CEREBELLUMWM=LH_CEREBELLUM_WHITE_MATTER,
                          BRAINSTEM=BRAIN_STEM)
  cc_name <- function(x) {
    if (grepl("^CC", x)) {
      x <- gsub("_", "", x)
      gsub("^CC", "CC_", x)
    } else {
      x
    }
  }
  
  names(nacc_mri_scan) <- sapply(names(nacc_mri_scan), cc_name)
  
  #oasis_conv <- lapply(names(oasis_mri_data), convert_name)
  names(oasis_mri_data) <- sapply(names(oasis_mri_data), function(name) {
    if (name %in% oasis_mri_vars) {
      name <- convert_name(name)
      name <- cc_name(name)
    }
    name
  })
  names(oasis_mri_data) <- gsub("_VOLUME", "", names(oasis_mri_data))
  oasis_mri_data <- rename(oasis_mri_data,
                           LH_ACCUMBENSAREA=LH_ACCUMBENS,
                           RH_ACCUMBENSAREA=RH_ACCUMBENS,
                           LH_THALAMUS=LH_THALAMUSPROPER,
                           RH_THALAMUS=RH_THALAMUSPROPER,
                           RH_CEREBELLUMWM=RH_CEREBELLUMWHITEMATTER,
                           LH_CEREBELLUMWM=LH_CEREBELLUMWHITEMATTER)
  
  ucsf_all$Cohort <- "ADNI"
  nacc_mri_scan$Cohort <- "NACC"
  oasis_mri_data$Cohort <- "OASIS"
  
  all_columns <- base::intersect(base::intersect(colnames(ucsf_all), colnames(nacc_mri_scan)), 
                                 colnames(oasis_mri_data))
  
  multi_cohort_mri <- rbind(
    select(ucsf_all, all_of(all_columns)),
    select(nacc_mri_scan, all_of(all_columns)),
    select(oasis_mri_data, all_of(all_columns))
  )
  
  return(multi_cohort_mri)
}



