library(clue)
library(dplyr)
library(mclust)

source("~/R/LTC/utils/cluster_utils.R")

load("~/R/EDAP-data/LTC_MC/new/exp_km_ab_ao.RData")

Clusters <- data.frame(
  Cluster_0 = multiLTC@Cluster,
  RID = multiLTC@RID
)
n <- length(unique(Clusters$Cluster_0))

for (ii in 1:5) {
  run <- sprintf("exp_km_ab_ao_%d", ii)
  load(paste("~/R/EDAP-data/LTC_MC/cross_validation/", run, ".Rdata", sep = ""))
  
  thisC <- data.frame(Cluster = multiLTC@Cluster, RID = multiLTC@RID)
  
  colnames(thisC) <- c(sprintf("Cluster_%d", ii), "RID")
  Clusters <- Clusters %>% left_join(thisC, by = "RID")
  Clusters[, sprintf("Cluster_%d", ii)] <- align_clusters(Clusters[, sprintf("Cluster_%d", ii)], Clusters$Cluster_0)
  
}


rm(multiLTC)

table(Clusters$Cluster_0, Clusters$Cluster_1)
table(Clusters$Cluster_0, Clusters$Cluster_2)
table(Clusters$Cluster_0, Clusters$Cluster_3)
table(Clusters$Cluster_0, Clusters$Cluster_4)
table(Clusters$Cluster_0, Clusters$Cluster_5)


adjustedRandIndex(Clusters$Cluster_0, Clusters$Cluster_1)
adjustedRandIndex(Clusters$Cluster_0, Clusters$Cluster_2)
adjustedRandIndex(Clusters$Cluster_0, Clusters$Cluster_3)
adjustedRandIndex(Clusters$Cluster_0, Clusters$Cluster_4)
adjustedRandIndex(Clusters$Cluster_0, Clusters$Cluster_5)

#Clusters %>% select(RID, Cluster_0, Cluster_4) %>% na.omit() %>% summarise(ari = adjustedRandIndex(Cluster_0, Cluster_4))
#adjustedRandIndex(Clusters$Cluster_0, Clusters$Cluster_4)

aris <- Clusters %>% select(-RID, -Cluster_0) %>% sapply(function(col) {
  adjustedRandIndex(Clusters$Cluster_0, col)
})
print(mean(aris))
print(sd(aris))

sim_df <- data.frame()
for(c in unique(Clusters$Cluster_0)) {
  
  df <- Clusters %>% filter(Cluster_0 == c)
  ratio <- df %>% summarise(across(-1, ~ {
    valid <- !is.na(.x)
    mean(.x[valid] == c)
  }))
  ratio$Cluster <- c
  sim_df <- rbind(sim_df, ratio)
}


