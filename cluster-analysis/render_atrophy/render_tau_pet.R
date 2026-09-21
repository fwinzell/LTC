library(dplyr)
library(tidyr)
library(ADNIMERGE2)
library(lme4)
library(stats)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(tibble)
library(kml3d)
library(RColorBrewer)
library(paletteer)


source("~/R/LTC/utils/adni_data_loaders.R")
run <- "exp_km_ab"
load(paste("~/R/EDAP-data/LTC4/", run, ".Rdata", sep = ""))

Clusters <- data.frame(
  Cluster = adniLTC@Cluster,
  RID = adniLTC@RID
)

dpm_df <- read.csv("~/R/EDAP-data/ADNI/DPM_ADNI.csv", header=TRUE)

# TAU PET
## Partial Volume Correction (PVC)
# Tau data corrected for partial volume effects using the Geometric Transfer Matrix (GTM) approach. 
# Using the MRI closest in time to tau scan, the GTM approach models all FreeSurfer-defined ROIs as well as regions in which offtarget
# binding is common (e.g., choroid plexus in FTP; meninges in MK6240) to reduce
# contamination from these regions into neighboring regions of interest. 

#tau.pet.pvc <- ADNIMERGE::ucberkeleyav1451_pvc
#tau.pet.pvc <- tau.pet.pvc |> left_join(Clusters, by = c("RID" = "RID")) |> filter(!is.na(Cluster))

#visits <- tau.pet.pvc$VISCODE
#visits <- gsub("bl", "0", visits)
#visits <- sapply(visits, function(x) gsub("m", "", x))
#tau.pet.pvc$M <- as.numeric(visits)
#tau.pet.pvc |> distinct(RID, M, .keep_all = TRUE) |> arrange(RID, EXAMDATE) -> tau.pet.pvc

# Some subjects have missing M values. Calculate the replacement using the EXAMDATE 
# and other M values, drop any single cases with missing M

#tau.pet.pvc$M.bl <- epitools::as.month(as.character(tau.pet.pvc$EXAMDATE))$stratum
#tau.pet.pvc <- tau.pet.pvc %>% group_by(RID) %>% mutate(M.bl = (M.bl - min(M.bl))/30.5) %>%
#  mutate(M.start = round(mean(na.omit(M-M.bl)))) %>% mutate(M = round(if_else(is.na(M), M.bl+M.start, M))) %>%
#  drop_na(M)

tau.pet.pvc <- get_tau_pet()
ab.df <- get_ab_df() 
dx.df <- get_diagnoses()

ab_neg_rids <- filter(ab.df, VISCODE2=="bl" & !AB_any) %>% distinct(RID) %>% unlist()
controls <- filter(dx.df, RID %in% ab_neg_rids & DIAGNOSIS == "CN" & !CI) %>% distinct(RID) %>% unlist()
tau.controls <- filter(tau.pet.pvc, RID %in% controls)

tau.pet.pvc <- tau.pet.pvc |> left_join(Clusters, by = c("RID" = "RID")) |> filter(!is.na(Cluster))

tau.pet.pvc <- distinct(dpm_df, RID, time_shift) %>% inner_join(tau.pet.pvc, by="RID") %>%
  mutate(Time = Years + time_shift)

tau.pet.pvc %>% distinct(RID, Cluster) %>% count(Cluster)
length(unique(tau.pet.pvc$RID))

tau.pet.pvc$Stage <- cut(tau.pet.pvc$Time, 
                      breaks = c(-Inf, 0, 5, 10, 12.5, Inf),
                      include.lowest = TRUE)
table(tau.pet.pvc[c("Cluster", "Stage")])

vars <- colnames(tau.pet.pvc)[grep("^CTX_[LR]H_.*_SUVR$", names(tau.pet.pvc))]

tau.pet.pvc <- mutate(tau.pet.pvc, Cluster = factor(as.character(Cluster)))

all_means <- data.frame()
for (varname in vars) {
  #mu <- mri_data %>% filter(Time < quantile(Time, 0.025)) %>% select(all_of(varname)) %>% 
  #  summarise(Mean = mean(.data[[varname]], na.rm=TRUE),
  #            SD = sd(.data[[varname]], na.rm=TRUE)) 
  
  mu <- tau.controls %>% select(RID, Months, all_of(varname)) %>% na.omit() %>%
    arrange(Months) %>% distinct(RID, .keep_all = TRUE) %>%
    summarise(Mean = mean(.data[[varname]], na.rm=TRUE),
              SD = sd(.data[[varname]], na.rm=TRUE)) 
  
  df <- tau.pet.pvc %>% select(RID, Time, Cluster, Stage, all_of(varname)) %>%
    drop_na(Cluster) %>% mutate(z = (.data[[varname]]-mu$Mean)/mu$SD) %>%
    group_by(RID, Stage) %>% 
    slice_min(abs(Time - median(Time)), n = 1) %>%
    ungroup()
  
  means <- df %>%
    group_by(Cluster, Stage) %>%
    summarise(Mean = mean(z, na.rm = TRUE),
              SD = sd(z, na.rm=TRUE),
              n = n(), .groups = "drop") %>%
    mutate(Region = varname)
  
  all_means <- rbind(all_means, means) 
}

############## WORKBENCH PROJECTIONS #############################

# directory in which the workbench software is stored and in which all surface renderings should be stored
dir.workbench.software = "/Users/filipwinzell/Workbench/software/workbench"#paste0(dir.root.olink, "AHBA_correlations/WorkBench_projections/software/workbench/")

# directory in which you want all surface renderings to be stored
dir.workbench = '/Users/filipwinzell/Workbench/Surface_renderings/'#paste0(dir.root.olink, "AHBA_correlations/WorkBench_projections/surface_renderings/")

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

workbench_names <- read.csv("/Users/filipwinzell/Workbench/atlas/Desikan_connectome_workbench.csv") %>% 
  mutate(Region = gsub("-", "_", toupper(X))) 

clusters <- unique(all_means$Cluster)
stages <- unique(all_means$Stage)

for (c in clusters) {
  stage_idx = 1
  for (s in stages) {
    wb_df <- all_means %>% filter(Cluster == c & Stage == s) %>% mutate(Region = gsub("_SUVR", "", Region))
    workbench_names %>% select(ROI, Label, Region) %>% left_join(wb_df, by = "Region") -> wb_df
    
    atr.vec <- wb_df$Mean
    atr.vec[is.na(atr.vec)] <- 0
    
    outfile <- paste("Tau_suvr_exp_km_ab", c, stage_idx,sep="_")
    stage_idx = stage_idx + 1
    
    outfile_txt=paste0(dir.workbench, "tab_", outfile,".txt")
    write.table(atr.vec, file = outfile_txt, row.names = F, col.names = F)
    outfile_cifti=paste0(dir.workbench, outfile,".pscalar.nii")
    render_to_fs(outfile_txt, dir.workbench, outfile_cifti)
    
  }
}




