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


multi_cohort_df <- read.csv("~/R/EDAP-data/MULTI_COHORT_4.csv", header = TRUE)

run <- "exp_km_ab_ao_2"
load(paste("~/R/EDAP-data/LTC_MC/new/", run, ".Rdata", sep = ""))

adni_dl <- new.env()
source("~/R/LTC/utils/adni_data_loaders.R", local=adni_dl)

oasis_dl <- new.env()
source("~/R/LTC/utils/oasis_data_loaders.R", local=oasis_dl)

mc_dl <- new.env()
source("~/R/LTC/utils/multi_cohort_loader.R", local=mc_dl)

ab_pos_rids <- c(paste0("ADNI_", adni_dl$get_ab_pos_ids()),
                 gsub("OAS", "OASIS_", oasis_dl$get_ab_pos_ids()))


Clusters <- data.frame(
  Cluster = multiLTC@Cluster,
  RID = multiLTC@RID
)

Clusters %>% left_join(multi_cohort_df, by = "RID") %>% drop_na(Cluster) -> multi_cohort_df

# TAU PET
## Partial Volume Correction (PVC)
# Tau data corrected for partial volume effects using the Geometric Transfer Matrix (GTM) approach. 
# Using the MRI closest in time to tau scan, the GTM approach models all FreeSurfer-defined ROIs as well as regions in which offtarget
# binding is common (e.g., choroid plexus in FTP; meninges in MK6240) to reduce
# contamination from these regions into neighboring regions of interest. 

# Some subjects have missing M values. Calculate the replacement using the EXAMDATE 
# and other M values, drop any single cases with missing M


adni_tau <- adni_dl$get_tau_pet() %>%
  mutate(RID = paste0("ADNI_", RID),
         Cohort = "ADNI")
oasis_tau <- oasis_dl$get_tau_pet(get_braak = FALSE)

vars <- colnames(adni_tau)[grep("^CTX_[LR]H_.*_SUVR$", names(adni_tau))]

source("~/R/LTC/cluster-analysis/render_atrophy/convert_tau_names.R")

oasis_cols <- grepv("^PET_fSUVR_rsf_[LR]_CTX", colnames(oasis_tau))
oasis_tau <- select(oasis_tau, OASISID, tracer, days_to_visit, all_of(oasis_cols)) %>%
  mutate(Years = days_to_visit/365.25,
         RID = gsub("^OAS", "OASIS_", OASISID),
         Cohort = "OASIS")
names(oasis_tau) <- unname(sapply(colnames(oasis_tau), convert_oasis_tau, simplify = "array"))

common <- intersect(names(oasis_tau), names(adni_tau))
tau_df <- bind_rows(select(oasis_tau, all_of(common)), select(adni_tau, all_of(common))) %>%
  select(RID, Cohort, Years, everything())


adni_dx <- adni_dl$get_diagnoses()

adni_pos_ids <- adni_dl$get_ab_pos_ids()
adni_controls <- filter(adni_dx, DIAGNOSIS == "CN" & !CI) %>% 
  filter(!(RID %in% adni_pos_ids)) %>% distinct(RID) %>% 
  mutate(RID = paste0("ADNI_", RID)) %>% unlist() 

oasis_dx <- oasis_dl$get_diagnoses()
oasis_pos_ids <- oasis_dl$get_ab_pos_ids()
oasis_controls <- filter(oasis_dx, DX == "CN") %>%
  filter(!(OASISID %in% oasis_pos_ids)) %>% distinct(OASISID) %>% 
  mutate(OASISID = gsub("^OAS", "OASIS_", OASISID)) %>% unlist()

tau.controls <- filter(tau_df, RID %in% c(adni_controls, oasis_controls))

ggplot(tau_df, aes(x=Years, y=CTX_LH_INFERIORPARIETAL_SUVR, fill = Cohort, color = Cohort)) +
  geom_point(aes(group = RID)) +
  facet_wrap(~Cohort) +
  theme_classic()


tau_df <- multi_cohort_df %>% distinct(RID, time_shift) %>% right_join(tau_df, by = "RID") %>%
  mutate(Time = Years + time_shift) %>% left_join(Clusters, by = "RID") %>% drop_na(Cluster)

tau_df %>% distinct(RID, Cluster) %>% count(Cluster)
length(unique(tau_df$RID))

tau_df$Stage <- cut(tau_df$Time, 
                      breaks = c(-Inf, 3, 6, 9, 12, Inf),
                      include.lowest = TRUE)
table(tau_df[c("Cluster", "Stage")], useNA = "ifany")

tau_df <- mutate(tau_df, Cluster = factor(as.character(Cluster)))

all_means <- data.frame()
for (varname in vars) {
  #mu <- mri_data %>% filter(Time < quantile(Time, 0.025)) %>% select(all_of(varname)) %>% 
  #  summarise(Mean = mean(.data[[varname]], na.rm=TRUE),
  #            SD = sd(.data[[varname]], na.rm=TRUE)) 
  
  mu <- tau.controls %>% select(RID, Years, all_of(varname)) %>% na.omit() %>%
    arrange(Years) %>% distinct(RID, .keep_all = TRUE) %>%
    summarise(Mean = mean(.data[[varname]], na.rm=TRUE),
              SD = sd(.data[[varname]], na.rm=TRUE)) 
  
  df <- tau_df %>% select(RID, Time, Cluster, Stage, all_of(varname)) %>%
    drop_na(Cluster) %>% mutate(z = (.data[[varname]]-mu$Mean)/mu$SD) %>%
    group_by(RID, Stage) %>% 
    slice_min(abs(Time - median(Time)), n = 1) %>%
    ungroup()
  
  means <- df %>%
    group_by(Cluster, Stage) %>%
    summarise(Mean = mean(z, na.rm = TRUE),
              Median = median(z, na.rm=TRUE),
              SD = sd(z, na.rm=TRUE),
              n = n(), 
              t = Mean / (SD / sqrt(n)),
              .groups = "drop") %>%
    mutate(Region = varname)
  
  all_means <- rbind(all_means, means) 
}

############## WORKBENCH PROJECTIONS #############################

# directory in which the workbench software is stored and in which all surface renderings should be stored
dir.workbench.software = "/Users/filipwinzell/Workbench/software/workbench"#paste0(dir.root.olink, "AHBA_correlations/WorkBench_projections/software/workbench/")

# directory in which you want all surface renderings to be stored
dir.workbench = '/Users/filipwinzell/Workbench/Surface_renderings/TAU/'#paste0(dir.root.olink, "AHBA_correlations/WorkBench_projections/surface_renderings/")

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
    
    atr.vec <- wb_df$Median
    atr.vec[is.na(atr.vec)] <- 0
    
    outfile <- paste("Tau_med_suvr_exp_", c, stage_idx,sep="_")
    stage_idx = stage_idx + 1
    
    outfile_txt=paste0(dir.workbench, "tab_", outfile,".txt")
    write.table(atr.vec, file = outfile_txt, row.names = F, col.names = F)
    outfile_cifti=paste0(dir.workbench, outfile,".pscalar.nii")
    render_to_fs(outfile_txt, dir.workbench, outfile_cifti)
    
  }
}




