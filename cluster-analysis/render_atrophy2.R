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

dataset <- "oasis"
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
    slice_min(abs(Time - median(Time)), n = 1, with_ties=FALSE) %>%
    ungroup()
  
  means <- df %>% 
    group_by(Cluster, Stage) %>%
    summarise(Mean = mean(z, na.rm = TRUE),
              SD = sd(z, na.rm=TRUE),
              n = n(), .groups = "drop") %>%
    mutate(Variable = varname)
  
  all_means <- rbind(all_means, means) 
}

#all_means <- pivot_longer(all_means, cols = matches("^ST\\d+[A-Z]*"), names_to = "Variable", values_to = "Mean")

if (dataset == "adni") {
  bap <- brain_regions |> select(Variables, Lobe, BraakStage)
  bap <- datadic %>% select(FLDNAME, Region) %>% right_join(bap, by = c("Region" = "Variables"))
  
  all_means <- left_join(all_means, bap, by = c("Variable" = "FLDNAME"))
  
  all_means <- all_means %>% mutate(Type = sub(".*([A-Z]{2})$", "\\1", Variable)) %>% mutate(Type = as.factor(Type)) 
  all_means <- all_means %>% filter(Type == "CV") %>% mutate(Region = tolower(Region))
  
  workbench_names <- read.csv("/Users/filipwinzell/Workbench/atlas/Desikan_connectome_workbench.csv") %>% 
    mutate(Region = sub("ctx-lh-", "left", X)) %>% mutate(Region = sub("ctx-rh-", "right", Region))
} else if (dataset == "oasis") {
  all_means <- mutate(all_means, Region = gsub("_", "-", str_remove(Variable, "_volume"))) %>%
    select(-Variable)
  workbench_names <- read.csv("/Users/filipwinzell/Workbench/atlas/Desikan_connectome_workbench.csv") %>% 
    mutate(Region = sub("ctx-", "", X)) 
}


distinct(all_means, Cluster, n)
############## WORKBENCH PROJECTIONS #############################

# directory in which the workbench software is stored and in which all surface renderings should be stored
dir.workbench.software = "/Users/filipwinzell/Workbench/software/workbench"#paste0(dir.root.olink, "AHBA_correlations/WorkBench_projections/software/workbench/")

# directory in which you want all surface renderings to be stored
dir.workbench = paste0('/Users/filipwinzell/Workbench/Surface_renderings/LTC4/', dataset, '/') 

# directory in which the atlases are stored
dir.atlas = '/Users/filipwinzell/Workbench/atlas'#paste0(dir.root.olink, "AHBA_correlations/WorkBench_projections/atlas")

# render to Desikan atlas function
render_to_fs <- function(path_to_vector_txt_file, output_folder, out_file){
  
  fs_dlabel=paste0(dir.atlas, "/Desikan.dlabel.nii")
  #fs_dscalar=paste0(dir.atlas, "/cifti/Schaefer2018_200Parcels_7Networks_order.dscalar.nii")
  fs_pscalar=paste0(dir.atlas, "/Desikan.pscalar.nii")
  
  # render pet mean change
  command1="#!/bin/sh"
  command2=paste0("export PATH=$PATH:", dir.workbench.software, "/bin_macosx64")
  command3=paste0("wb_command -cifti-convert -from-text ", path_to_vector_txt_file, " ",fs_pscalar, " ", out_file)
  
  writeLines(c(command1, command2, command3), paste0(output_folder, "/tmp_render_to_workbench.sh"))
  bash_command=paste0("bash ", paste0(output_folder, "/tmp_render_to_workbench.sh"))
  system(bash_command)
  
}


# render to Schaefer atlas function
render_to_schaefer <- function(path_to_vector_txt_file, output_folder, out_file){
  
  Atlas.dlabel=paste0(dir.atlas, "/Schaefer2018_200Parcels_7Networks_order.dlabel.nii")
  Atlas.dscalar=paste0(dir.atlas, "/Schaefer2018_200Parcels_7Networks_order.dscalar.nii")
  Atlas.pscalar=paste0(dir.atlas, "/Schaefer2018_200Parcels_7Networks_order.pscalar.nii")
  
  # render pet mean change
  command1="#!/bin/sh"
  command2= paste0("export PATH=$PATH:", dir.workbench.software, "/bin_macosx64")
  command3=paste0("wb_command -cifti-convert -from-text ", path_to_vector_txt_file, " ",Atlas.pscalar, " ", out_file)
  
  writeLines(c(command1, command2, command3), paste0(output_folder, "/bash_render_to_workbench.sh"))
  bash.command=paste0("bash ", paste0(output_folder, "/bash_render_to_workbench.sh"))
  system(bash.command)
  
}


clusters <- unique(all_means$Cluster)
stages <- unique(all_means$Stage)

for (c in clusters) {
  stage_idx = 1
  for (s in stages) {
    wb_df <- all_means %>% filter(Cluster == c & Stage == s) 
    workbench_names %>% select(ROI, Label, Region) %>% inner_join(wb_df, by = "Region") -> wb_df
    
    missing <- setdiff(workbench_names$ROI, wb_df$ROI)
    if (length(missing) > 0) {
      wb_df <- workbench_names %>% filter(ROI %in% missing) %>% select(ROI, Label, Region) %>%
        mutate(Cluster = c, Stage = s, Mean = NA, SD = NA, n = NA) %>% rbind(wb_df)
    }
    wb_df <- arrange(wb_df, ROI)
    
    atr.vec <- wb_df$Mean
    atr.vec[is.na(atr.vec)] <- 0
    
    outfile <- paste("exp_km_ab", c, stage_idx, sep="_")
    stage_idx = stage_idx + 1
    
    outfile_txt=paste0(dir.workbench, "tab_", outfile,".txt")
    write.table(atr.vec, file = outfile_txt, row.names = F, col.names = F)
    outfile_cifti=paste0(dir.workbench, outfile,".pscalar.nii")
    render_to_fs(outfile_txt, dir.workbench, outfile_cifti)
    
  }
}


#matches <- cluster_a[cluster_a$Region %in% workbench_names$Region, ]
                               