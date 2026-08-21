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

multi_cohort_df <- read.csv("~/R/EDAP-data/MULTI_COHORT.csv", header = TRUE)

run <- "exp_km_ab_ao"
load(paste("~/R/EDAP-data/LTC_MC/", run, ".Rdata", sep = ""))

adni_dl <- new.env()
source("~/R/LTC/utils/adni_data_loaders.R", local=adni_dl)

oasis_dl <- new.env()
source("~/R/LTC/utils/oasis_data_loaders.R", local=oasis_dl)

mc_dl <- new.env()
source("~/R/LTC/utils/multi_cohort_loader.R", local=mc_dl)

source("~/R/LTC/utils/analysis_utils.R")

ab_pos_rids <- c(paste0("ADNI_", adni_dl$get_ab_pos_ids()),
                 gsub("OAS", "OASIS_", oasis_dl$get_ab_pos_ids()))

mc_mri <- mc_dl$multi_cohort_mri()
mri_controls <- mc_mri %>% filter(DX.bl == "CN" & !(RID %in% ab_pos_rids))

Clusters <- data.frame(
  Cluster = multiLTC@Cluster,
  RID = multiLTC@RID
)

Clusters %>% left_join(multi_cohort_df, by = "RID") %>% drop_na(Cluster) -> multi_cohort_df

tran <- lapply(unique(Clusters$Cluster), function(c) {
  multi_cohort_df %>% filter(Cluster == c) %>% select(Time) %>% unlist() %>% quantile(c(0.05, 0.95))
})

tran

multi_cohort_df$Stage <- cut(multi_cohort_df$Time, 
                      breaks = c(-Inf, 0, 5, 10, 12.5, 15, Inf),
                      include.lowest = TRUE)
mri_cols <- c(grepv("^(RH_|LH_|CC_)", colnames(multi_cohort_df)), "BRAINSTEM")

all_means <- data.frame()

for (varname in mri_cols) {
  mu <- mri_controls %>% select(RID, Months, all_of(varname)) %>% na.omit() %>%
    arrange(Months) %>% distinct(RID, .keep_all = TRUE) %>%
    summarise(Mean = mean(.data[[varname]], na.rm=TRUE),
              SD = sd(.data[[varname]], na.rm=TRUE)) 
  
  df <- multi_cohort_df %>% select(RID, Time, DX.bl, Cluster, Stage, all_of(varname)) %>%
    drop_na(Cluster) %>% mutate(z = (.data[[varname]]-mu$Mean)/mu$SD) %>%
    group_by(RID, Stage) %>% 
    slice_min(abs(Time - median(Time)), n = 1, with_ties=FALSE) %>%
    ungroup()
  
  means <- df %>% 
    group_by(Cluster, Stage) %>%
    summarise(Mean = mean(z, na.rm = TRUE),
              SD = sd(z, na.rm=TRUE),
              n = n(), .groups = "drop") %>%
    mutate(Region = varname)
  
  all_means <- rbind(all_means, means) 
}

#all_means <- pivot_longer(all_means, cols = matches("^ST\\d+[A-Z]*"), names_to = "Variable", values_to = "Mean")
workbench_names <- read.csv("/Users/filipwinzell/Workbench/atlas/Desikan_connectome_workbench.csv") %>%
  mutate(Region = sub("-", "_", toupper(sub("ctx-", "", X)))) 

distinct(all_means, Cluster, n)
############## WORKBENCH PROJECTIONS #############################

# directory in which the workbench software is stored and in which all surface renderings should be stored
dir.workbench.software = "/Users/filipwinzell/Workbench/software/workbench"#paste0(dir.root.olink, "AHBA_correlations/WorkBench_projections/software/workbench/")

# directory in which you want all surface renderings to be stored
dir.workbench = paste0('/Users/filipwinzell/Workbench/Surface_renderings/LTC_MC/') 

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
                               