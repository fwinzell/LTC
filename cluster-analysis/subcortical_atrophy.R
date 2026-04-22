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
  
  mri_cols <- grepv("^ST\\d+SV*", colnames(mri_data)) 
  
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

brain_regions <- get_brain_regions()

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
  
}


collect.wide <- pivot_wider(all_collect, names_from = Variable, values_from = z, values_fn = mean)

if (dataset == "adni") {
  bap <- brain_regions |> select(Variables, Lobe, BraakStage)
  bap <- datadic %>% select(FLDNAME, Region) %>% right_join(bap, by = c("Region" = "Variables"))
  
  all_collect <- left_join(all_collect, bap, by = c("Variable" = "FLDNAME"))
} 

for (col in mri_cols) {
  df <- filter(all_collect, Variable == col)
  varname = df$Region[1]
  gp <- ggplot(df, aes(x=Time, y=z, color=Cluster, fill=Cluster, group=RID)) +
    geom_line() +
    geom_point() +
    scale_color_paletteer_d("ggthemes::Tableau_10") +
    scale_fill_paletteer_d("ggthemes::Tableau_10") +
    #geom_smooth(method="loess") +
    geom_hline(yintercept = 0, linetype=2) +
    geom_ribbon(
      inherit.aes = FALSE,
      aes(x = Time, ymin = -1.96, ymax = 1.96),
      fill = "grey70",
      alpha = 0.25
    ) +
    labs(title=varname) +
    theme_classic()
  plot(gp)
}

