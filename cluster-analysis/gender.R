library(ADNIMERGE2)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(dplyr)
library(Rtsne)
library(paletteer)
library(tidyr)
library(broom)
library(nlme)

library(icenReg)
library(survival)

source("~/R/LTC/utils/adni_data_loaders.R")
run <- "exp_km_ab"
load(paste("~/R/EDAP-data/LTC4/", run, ".Rdata", sep = ""))

ucsf_data_norm <- ucsf_longitudinal_all(only_vol = TRUE, filter_n = 0, normalize = TRUE)
ucsf_data_raw <- ucsf_longitudinal_all(only_vol = TRUE, filter_n = 0, normalize = FALSE)

dpm_df <- read.csv("~/R/EDAP-data/ADNI/DPM_ADNI.csv", header=TRUE)
ucsf_data <- dpm_df %>% select(RID, time_shift) %>% distinct() %>% right_join(ucsf_data, by="RID") %>%
  mutate(Time = Years + time_shift)

Clusters <- data.frame(
  Cluster = adniLTC@Cluster,
  RID = adniLTC@RID
)

Clusters %>% left_join(ucsf_data, by = "RID") %>% drop_na(Cluster) -> ucsf_data

adni_demog <- get_demographics()

ucsf_data_norm <- adni_demog %>% select(RID, PTGENDER) %>% right_join(ucsf_data_norm, by="RID")
ucsf_data_raw <- adni_demog %>% select(RID, PTGENDER) %>% right_join(ucsf_data_raw, by="RID")

mri_cols <- grepv("^ST\\d+(CV|SV|SA|TA)$", colnames(ucsf_data))

col <- mri_cols[1]
for (col in mri_cols) {
  gp <- ggplot(ucsf_data_raw, aes(y=!!sym(col), x=PTGENDER, fill=PTGENDER)) +
    geom_boxplot()
  plot(gp)
}



