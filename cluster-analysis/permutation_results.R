
source("~/R/LTC/utils/cluster_utils.R")

run <- "exp_km_ab_ao_2"
load(paste("~/R/EDAP-data/LTC_MC/new/", run, ".Rdata", sep = ""))

Clusters <- data.frame(
  Cluster = multiLTC@Cluster,
  RID = multiLTC@RID
)

for (cl in LETTERS[1:4]) {
  load(paste0("~/R/EDAP-data/LTC_MC/permuted/exp_km_ab_ao_", cl, ".Rdata"))
  
  thisC <- data.frame(Cluster = multiLTC@Cluster, RID = multiLTC@RID)
  
  colnames(thisC) <- c(paste0("Cluster_", cl), "RID")
  Clusters <- Clusters %>% left_join(thisC, by = "RID")
  Clusters[, paste0("Cluster_", cl)] <- align_clusters_new(Clusters[, paste0("Cluster_", cl)], Clusters$Cluster)
  
}

align = Clusters[, paste0("Cluster_", cl)]
to = Clusters$Cluster


table(Clusters$Cluster, Clusters$Cluster_A)
table(Clusters$Cluster, Clusters$Cluster_B)
table(Clusters$Cluster, Clusters$Cluster_C)
table(Clusters$Cluster, Clusters$Cluster_D)

calculate_accuracy <- function(ref, test) {
  tab <- table(ref, test)
  sum(diag(tab))/sum(tab)
}

calculate_accuracy(Clusters$Cluster, Clusters$Cluster_A)
calculate_accuracy(Clusters$Cluster, Clusters$Cluster_B)
calculate_accuracy(Clusters$Cluster, Clusters$Cluster_C)
calculate_accuracy(Clusters$Cluster, Clusters$Cluster_D)

Clusters$RID <- as.character(Clusters$RID)
missc_A <- filter(Clusters, as.character(Cluster) != as.character(Cluster_A)) %>% select(RID)
missc_B <- filter(Clusters, as.character(Cluster) != as.character(Cluster_B)) %>% select(RID)
missc_C <- filter(Clusters, as.character(Cluster) != as.character(Cluster_C)) %>% select(RID)
missc_D <- filter(Clusters, as.character(Cluster) != as.character(Cluster_D)) %>% select(RID)

ffs <- rbind(missc_A, missc_B, missc_C, missc_D) %>% distinct(RID)
true <- filter_out(Clusters, RID %in% ffs$RID)



