
load("~/R/EDAP-data/LTC_MC/cross_validation/exp_km_ab_ao_1.Rdata")

Clusters_1 <- data.frame(
  Cluster = multiLTC@Cluster,
  RID = multiLTC@RID
)

load("~/R/EDAP-data/LTC_MC/cross_validation/exp_km_ab_ao_2.Rdata")

Clusters_2 <- data.frame(
  Cluster = multiLTC@Cluster,
  RID = multiLTC@RID
)

for (fold in folds) {
  rids <- multi_cohort_df$RID[fold]
  overlap <- length(intersect(rids, Clusters_1$RID))
  print(paste0("Ratioed: ", overlap/nrow(Clusters_1)))
}

for (fold in folds) {
  rids <- multi_cohort_df$RID[fold]
  overlap <- length(intersect(rids, Clusters_2$RID))
  print(paste0("Ratioed: ", overlap/nrow(Clusters_2)))
}

fold_1 <- setdiff(multi_cohort_df$RID, Clusters_1$RID)
fold_2 <- setdiff(multi_cohort_df$RID, Clusters_2$RID)

fold_rids <- unique(c(fold_1, fold_2))


#### Cross-validation ####
k=3
set.seed(123)
dffs <- multi_cohort_df %>% filter_out(RID %in% fold_rids) %>% distinct(RID, DX.bl)
new_folds <- caret::createFolds(dffs$DX.bl, k = 3)

idx_folds <- lapply(new_folds, function(rids) {
  which(multi_cohort_df$RID %in% dffs[rids, "RID"])
})


names(idx_folds) <- c("Fold3", "Fold4", "Fold5")

idxs_1 <- which(multi_cohort_df$RID %in% fold_1)
idxs_2 <- which(multi_cohort_df$RID %in% fold_2)

idx_folds$Fold1 <- idxs_1
idx_folds$Fold2 <- idxs_2

folds <- idx_folds[sort(names(idx_folds))]

sum(unlist(lapply(folds, length)))

rm(Clusters_1)
rm(Clusters_2)
rm(dffs)
rm(multiLTC)
rm(new_folds)


