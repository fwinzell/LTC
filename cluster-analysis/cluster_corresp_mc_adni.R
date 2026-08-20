library(clue)
library(dplyr)
library(mclust)

source("~/R/LTC/utils/cluster_utils.R")

load("~/R/EDAP-data/LTC4/exp_km_ab.Rdata")

Clusters_adni <- data.frame(
  Cluster_adni = adniLTC@Cluster,
  RID = adniLTC@RID
) 


load("~/R/EDAP-data/LTC4/oasis_exp_km_ab.Rdata")

Clusters_oasis <- data.frame(
  Cluster_oasis = oasisLTC@Cluster,
  RID = oasisLTC@RID
)


run <- "exp_km_ab_ao"
load(paste("~/R/EDAP-data/LTC_MC/", run, ".Rdata", sep = ""))

Clusters <- data.frame(
  Cluster_mc = multiLTC@Cluster,
  RID = multiLTC@RID
) %>% mutate(Cohort = gsub("_.*", "", RID))


# ADNI comparison

comb <- Clusters %>% filter(Cohort == 'ADNI') %>% mutate(RID = as.numeric(gsub("ADNI_", "", RID))) %>%
  left_join(Clusters_adni, by='RID')

table(comb$Cluster_mc, comb$Cluster_adni)

# OASIS comparison

comb <- Clusters %>% filter(Cohort == 'OASIS') %>% mutate(RID = as.numeric(gsub("OASIS_", "", RID))) %>%
  left_join(Clusters_oasis, by='RID')

table(comb$Cluster_mc, comb$Cluster_oasis)




