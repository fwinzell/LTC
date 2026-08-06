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

dataset <- "adni"
source("~/R/LTC/utils/analysis_utils.R")

if (dataset == "adni") {
  source("~/R/LTC/utils/adni_data_loaders.R")
  run <- "exp_km_ab"
  load(paste("~/R/EDAP-data/LTC4/", run, ".Rdata", sep = ""))
  
  ab_df <- get_ab_df()
  
  ab_df %>% select(RID, AB_any) %>% group_by(RID) %>% 
    mutate(AB_any = if_any(AB_any)) %>% distinct() %>% filter(AB_any) %>% 
    select(RID) %>% unlist() -> ab_pos_rids_any
  
  mri_data <- ucsf_longitudinal_all(only_vol = TRUE, filter_n = 0)
  mri_controls <- mri_data %>% filter(!(RID %in% ab_pos_rids_any)) %>% filter(DX.bl == "CN")
  
  dpm_df <- read.csv("~/R/EDAP-data/ADNI/DPM_ADNI.csv", header=TRUE)
  mri_data <- dpm_df %>% select(RID, time_shift) %>% distinct() %>% right_join(mri_data, by="RID") %>%
    mutate(Time = Years + time_shift) %>% drop_na(Time)
  
  Clusters <- data.frame(
    Cluster = adniLTC@Cluster,
    RID = adniLTC@RID
  )
  
  mri_data <- left_join(mri_data, Clusters, by = "RID")
  
  datadic <- datadict_as_tibble(ADNIMERGE2::DATADIC) %>% filter(TBLNAME == "UCSFFSL") %>%
    select(FLDNAME, TEXT) %>% filter(grepl("ST\\d+[A-Z]{2}", FLDNAME)) %>%
    mutate(Region = sub(".*\\s", "", TEXT))
  
  tran <- lapply(unique(Clusters$Cluster), function(c) {
    mri_data %>% filter(Cluster == c) %>% select(Time) %>% unlist() %>% quantile(c(0.05, 0.95))
  })
  
  tran
  
  mri_data$Stage <- cut(mri_data$Time, 
                         breaks = c(-Inf, 0, 5, 10, 12.5, 15, Inf),
                         include.lowest = TRUE)
  
  mri_cols <- grepv("^ST\\d+[A-Z]*", colnames(mri_data)) 
  
} else if(dataset == "oasis") {
  source("~/R/LTC/utils/oasis_data_loaders.R")
  run <- "oasis_exp_km_ab_2"
  load(paste("~/R/EDAP-data/LTC4/", run, ".Rdata", sep = ""))
  
  mri_data <- load_oasis_mri_data(unified_norm = TRUE)
  ab_pos_ids <- get_ab_pos_ids()
  #mri_data <- filter(mri_data, OASISID %in% ab_pos_ids)
  
  mri_controls <- mri_data %>% filter(DX.bl == "CN" & !(OASISID %in% ab_pos_ids))
  mri_data <- filter(mri_data, OASISID %in% ab_pos_ids)
  
  dpm_df <- read.csv("~/R/EDAP-data/OASIS/DPM_OASIS.csv", header=TRUE)
  mri_data <- dpm_df %>% select(RID, time_shift) %>% distinct() %>% right_join(mri_data, by="RID") %>%
    mutate(Time = Years + time_shift) %>% drop_na(Time)
  
  Clusters <- data.frame(
    Cluster = oasisLTC@Cluster,
    RID = oasisLTC@RID
  )
  
  mri_data <- left_join(mri_data, Clusters, by = "RID")
  
  tran <- lapply(unique(Clusters$Cluster), function(c) {
    mri_data %>% filter(Cluster == c) %>% select(Time) %>% unlist() %>% quantile(c(0.05, 0.95))
  })
  
  mri_data$Stage <- cut(mri_data$Time, 
                         breaks = c(-Inf, 0, 5, 10, 12.5, 15, Inf),
                         include.lowest = TRUE)
  
  mri_cols <- c(grepv("^[lr]h_[a-z]+_volume$", colnames(mri_data)), 
                grepv("^(Left|Right)\\.[a-zA-Z._]+_volume$", colnames(mri_data)))
  
}

#region_names <- data.frame()
#for (varname in adniLTC@varNames) {
#  descr <- unlist(datadic[which(datadic$FLDNAME == varname), "Region"])
#  region_names <- rbind(region_names, data.frame(Variable=varname,
#                                                 Region=descr,
#                                                 Mean=mean(mri_data[[varname]], na.rm = TRUE), 
#                                                 SD=sd(mri_data[[varname]], na.rm = TRUE)))
#}
#region_names |> arrange(desc(Mean)) -> region_names

brain_regions <- get_brain_regions()

all_means <- data.frame()
all_collect <- data.frame()


for (varname in mri_cols) {
  #mu <- mri_data %>% filter(Time < quantile(Time, 0.025)) %>% select(all_of(varname)) %>% 
  #  summarise(Mean = mean(.data[[varname]], na.rm=TRUE),
  #            SD = sd(.data[[varname]], na.rm=TRUE)) 
  
  mu <- mri_controls %>% select(RID, Months, all_of(varname)) %>% na.omit() %>%
    arrange(Months) %>% distinct(RID, .keep_all = TRUE) %>%
    summarise(Mean = mean(.data[[varname]], na.rm=TRUE),
              SD = sd(.data[[varname]], na.rm=TRUE)) 
  
  df <- mri_data %>% select(RID, Time, DX.bl, Cluster, Stage, all_of(varname)) %>%
    drop_na(Cluster) %>% mutate(z = (.data[[varname]]-mu$Mean)/mu$SD) %>%
    group_by(RID, Stage) %>% 
    slice_min(abs(Time - median(Time)), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  collect <- select(df, RID, Time, Cluster, Stage, z) %>% 
    mutate(Variable = varname)
  
  all_collect <- rbind(all_collect, collect)
  
  means <- df %>%
    group_by(Cluster, Stage) %>%
    summarise(Mean = mean(z, na.rm = TRUE),
              SD = sd(z, na.rm=TRUE),
              n = n(), .groups = "drop") %>%
    mutate(Variable = varname)
  
  if (ncol(all_means) == 0) {
    all_means <- means
  } else {
    all_means <- rbind(all_means, means) 
  }
}


collect.wide <- pivot_wider(all_collect, names_from = Variable, values_from = z, values_fn = mean)

stages = levels(mri_data$Stage)
clusters = levels(mri_data$Cluster)

data <- filter(all_collect, Stage ==stages[1]) 

library(fpc)

bhatta_fun <- function(stage_1, cluster_a, stage_2, cluster_b) {
  data <- filter(collect.wide, Stage %in% c(stage_1, stage_2), Cluster %in% c(cluster_a, cluster_b)) %>% na.omit() 
  
  pca <- data %>% select(all_of(mri_cols)) %>% prcomp(rank. = 10)
  
  if (cluster_a != cluster_b) {
    res <- data.frame(
      Cluster = data$Cluster,
      pca$x
    )
    
    Sigma_1 <- res %>% filter(Cluster == cluster_a) %>% select(-Cluster) %>% cov()
    Mu_1 <- res %>% filter(Cluster == cluster_a) %>% select(-Cluster) %>% colMeans(na.rm=TRUE)
    
    Sigma_2 <- res %>% filter(Cluster == cluster_b) %>% select(-Cluster) %>% cov()
    Mu_2 <- res %>% filter(Cluster == cluster_b) %>% select(-Cluster) %>% colMeans(na.rm=TRUE)  
  } else {
    res <- data.frame(
      Stage = data$Stage,
      pca$x
    )
    
    Sigma_1 <- res %>% filter(Stage == stage_1) %>% select(-Stage) %>% cov()
    Mu_1 <- res %>% filter(Stage == stage_1) %>% select(-Stage) %>% colMeans(na.rm=TRUE)
    
    Sigma_2 <- res %>% filter(Stage == stage_2) %>% select(-Stage) %>% cov()
    Mu_2 <- res %>% filter(Stage == stage_2) %>% select(-Stage) %>% colMeans(na.rm=TRUE)
  }
  
  
  D_b <- fpc::bhattacharyya.dist(mu1 = Mu_1, mu2 = Mu_2, Sigma1 = Sigma_1, Sigma2 = Sigma_2)
  return(exp(-D_b))
}


bdists <- lapply(stages,  bhatta_fun, cluster_a="A", cluster_b="B", stage_2=stages[6])
bdists



#### Bad idea #####

library(MANOVA.RM)

eeg <- spread(EEG, feature, resp)
head(eeg)
fit <- multRM(cbind(brainrate, complexity) ~ sex * region, data = eeg, subject = "id", within = "region", iter = 1000)
summary(fit)

fit <- multRM(cbind(ST104CV, ST105CV) ~ Cluster * Stage, data = collect.wide, subject = "RID", iter=1000)

stage_1 <- filter(all_means, Stage == "[-Inf,0]")

                               