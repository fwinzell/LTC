library(ADNIMERGE2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(stats)
#library(GauPro)

source("~/R/LTC/utils/adni_data_loaders.R")

ucsf_data <- ucsf_longitudinal_all(only_vol=TRUE, filter_n=0)

dpm_df <- read.csv("~/R/EDAP-data/ADNI/DPM_ADNI.csv", header=TRUE)

dpm_df |> select(RID, time_shift) |>
  unique() |> inner_join(ucsf_data, join_by(RID)) -> ucsf_data

ucsf_data$Time = ucsf_data$Years + ucsf_data$time_shift
ucsf_data <- mutate(ucsf_data, DX.bl = factor(DX.bl, levels = c("CN", "MCI", "Dementia")))

get_time_shift_gp <- function(varname, title, df = ucsf_data) {
  suff <- sub(".*([A-Z]{2})$", "\\1", varname)
  print(suff)
  if (suff == "SV") {
    ynm <- "Subcortical Vol."
  } else if (suff == "CV") {
    ynm <- "Cortical Vol."
  } else if (suff == "HS") {
    ynm <- "Hippocampal Vol."
  } else if (suff == "TA") {
    ynm <- "Thickness"
  } else if (suff == "SA") {
    ynm <- "Surface Area"
  }
  
  traj <- ggplot(df, aes(x = Years, y = !!sym(varname), group = RID, color = DX.bl)) +
    geom_line(alpha=0.85) +
    labs(title = title, x = "Years", y = ynm, color = NULL) +
    scale_color_brewer(palette = "YlOrRd") +
    guides(color = guide_legend(override.aes = list(size = 7.5))) +
    theme_classic() +
    scale_x_continuous(limits = c(0, NA), expand = c(0,0)) 
  
  traj_ts <- ggplot(df) +
    geom_line(aes(x = Time, y = !!sym(varname), group = RID, color = DX.bl),  alpha=0.85) +
    labs(title = " ", x = "Years (time shifted)", y = ynm, color = NULL) +
    scale_color_brewer(palette = "YlOrRd") + 
    guides(color = guide_legend(override.aes = list(size = 7.5))) +
    theme_classic() +
    scale_x_continuous(limits = c(0, NA), expand = c(0,0)) 
  
  return(list(traj, traj_ts))
}


gp1 <- get_time_shift_gp("ST24CV", "Left Entorhinal")
gp2 <- get_time_shift_gp("ST40CV", "Left Middle Temporal")
gp3 <- get_time_shift_gp("ST88SV", "Right Hippocampus")

merged_list <- list()
for (i in 1:2) {
  m = 3*(i-1)+1
  merged_list[[m]] <- gp1[[i]]
  merged_list[[m+1]] <- gp2[[i]]
  merged_list[[m+2]] <- gp3[[i]]
}

comb <- ggarrange(plotlist = merged_list, ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
plot(comb)

ggsave("~/R/EDAP-data/plots/LTC4/mri-time-shift.png", comb, width = 10, height = 7, dpi=1000)


