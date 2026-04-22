run <- "exp_km_ab"
load(paste("~/R/EDAP-data/LTC4/", run, ".Rdata", sep = ""))

# Extra utils for clustering and visualization
source("~/R/LTC/utils/cluster_utils.R")

dendro <- plot_dendrogram(adniLTC, save=TRUE)


