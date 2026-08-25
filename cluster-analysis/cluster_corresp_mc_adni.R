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

comb <- Clusters_new %>% filter(Cohort == 'ADNI') %>% mutate(RID = as.numeric(gsub("ADNI_", "", RID))) %>%
  left_join(Clusters_adni, by='RID')

table(comb$Cluster_new, comb$Cluster_adni)

cm_df <- as.data.frame(table(comb$Cluster_new, comb$Cluster_adni)) %>%
  rename(`ADNI+OASIS` = Var1,
         ADNI = Var2)

ggplot(cm_df, aes(x = ADNI, y = `ADNI+OASIS`, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), color = "black", size = 5) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Confusion Matrix", fill = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 1))

# OASIS comparison

comb <- Clusters_new %>% filter(Cohort == 'OASIS') %>% mutate(RID = as.numeric(gsub("OASIS_", "", RID))) %>%
  left_join(Clusters_oasis, by='RID')

table(comb$Cluster_new, comb$Cluster_oasis)


### New ADNI ####

load(paste("~/R/EDAP-data/LTC_MC/new/", run, ".Rdata", sep = ""))

Clusters_new <- data.frame(
  Cluster_new = multiLTC@Cluster,
  RID = multiLTC@RID
) %>% mutate(Cohort = gsub("_.*", "", RID))

comb <- left_join(Clusters_new, Clusters, by='RID')

table(comb$Cluster_new, comb$Cluster_mc)

cm_df <- as.data.frame(table(comb$Cluster_mc, comb$Cluster_adni)) %>%
  rename(`ADNI+OASIS` = Var1,
         ADNI = Var2)

ggplot(cm_df, aes(x = ADNI, y = `ADNI+OASIS`, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), color = "black", size = 5) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Confusion Matrix", fill = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 1))


