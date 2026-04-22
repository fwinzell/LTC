library(ADNIMERGE2)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(dplyr)
library(Rtsne)
library(paletteer)
library(ggsci)
library(tidyr)

library(gbm)
library(xgboost)
library(caret)
library(pROC)

dataset = "oasis"

source(paste0("~/R/LTC/utils/", dataset, "_data_loaders.R"))

if (dataset == "adni") {
  run <- "exp_km_ab"
  load(paste("~/R/EDAP-data/LTC4/", run, ".Rdata", sep = ""))
  
  mri_data <- ucsf_longitudinal_all(only_vol = TRUE, filter_n = 0)
  dpm_df <- read.csv("~/R/EDAP-data/ADNI/DPM_ADNI.csv", header=TRUE)
  mri_data <- dpm_df %>% select(RID, time_shift) %>% distinct() %>% right_join(mri_data, by="RID") %>%
    mutate(Time = Years + time_shift)
  
  Clusters <- data.frame(
    Cluster = adniLTC@Cluster,
    RID = adniLTC@RID
  )
  
  mri_data <- Clusters %>% left_join(mri_data, by = "RID") %>% drop_na(Cluster) %>%
    select(RID, Cluster, Time, all_of(adniLTC@varNames))
  
  mri_cols <- adniLTC@varNames
} else if(dataset == "oasis") {
  run <- "oasis_exp_km_ab_2"
  load(paste("~/R/EDAP-data/LTC4/", run, ".Rdata", sep = ""))
  
  mri_data <- load_oasis_mri_data(unified_norm = TRUE)
  
  dpm_df <- read.csv("~/R/EDAP-data/OASIS/DPM_OASIS.csv", header=TRUE)
  mri_data <- dpm_df %>% select(RID, time_shift) %>% distinct() %>% right_join(mri_data, by="RID") %>%
    mutate(Time = Years + time_shift) %>% drop_na(Time)
  
  Clusters <- data.frame(
    Cluster = oasisLTC@Cluster,
    RID = oasisLTC@RID
  )
  
  mri_data <- left_join(mri_data, Clusters, by = "RID") %>% drop_na(Cluster)
  
  mri_cols <- oasisLTC@varNames
  #mri_cols <- c(grepv("^[lr]h_[a-z]+_volume$", colnames(mri_data)), 
  #              grepv("^(Left|Right)\\.[a-zA-Z._]+_volume$", colnames(mri_data)))
}



create_strat_folds <- function(df, k = 5) {
  df %>% select(RID, Cluster) %>% unique() -> subjects  
  
  folds <- caret::createFolds(subjects$Cluster, k = 5)
  
  idx_folds <- lapply(folds, function(rids) {
    which(df$RID %in% subjects[rids, "RID"])
  })
  
  return(idx_folds)
}


set.seed(123)
folds <- create_strat_folds(mri_data, k=5) #caret::createFolds(cldp$Cluster, k = 5) #groupKFold(cldp$RID, k = 5)

i = 1
n_clusters <- length(unique(mri_data$Cluster))
cv_conf_df <- data.frame()
acc = 0
rocAll <- data.frame()
varAll <- data.frame()
aucAll <- data.frame()
for (f_idxs in folds) {
  
  train_data <- mri_data[-f_idxs, ]
  test_data <- mri_data[f_idxs, ]
  
  train_labels <- as.numeric(train_data$Cluster) - 1
  train_matrix <- xgb.DMatrix(as.matrix(train_data[,mri_cols]), label = train_labels)
  test_labels <- as.numeric(test_data$Cluster) - 1
  test_matrix <- xgb.DMatrix(as.matrix(test_data[,mri_cols]), label = test_labels)
  
  xgb_model <- xgboost(
    data = train_matrix,
    objective = "multi:softprob",
    num_class = n_clusters,
    nrounds = 200,
    max_depth = 10,
    eta = 0.1
  )
  
  pred_probs <- predict(xgb_model, test_matrix) 
  pred_matrix <- matrix(pred_probs, ncol=n_clusters, byrow=T, dimnames = list(c(), levels(test_data$Cluster)))
  pred_labels <- max.col(pred_matrix)
  
  res_df <- data.frame(Cluster = test_data$Cluster, RID = test_data$RID, Param = test_data$Time,
                       Prediciton = factor(pred_labels, levels = 1:n_clusters, labels = levels(test_data$Cluster)),
                       pred_matrix)
                       
  
  conf_mat <- confusionMatrix(res_df$Prediciton, res_df$Cluster)
  overall_acc <- conf_mat$overall['Accuracy']
  print(paste("Fold", i, "Accuracy", overall_acc))
  acc = acc + overall_acc
  
  conf_df <- as.data.frame(conf_mat$table)
  #conf_df$Freq <- conf_df$Freq / rowSums(conf_mat$table)[conf_df$Reference]
  #conf_df$Freq <- round(conf_df$Freq, 3)
  conf_df |> group_by(Reference) |> reframe(Prediction = Prediction, Freq = Freq) -> conf_df #/ sum(Freq)) -> conf_df
  
  conf_df$Prediction <- factor(conf_df$Prediction, levels = rev(levels(conf_df$Prediction)))
  
  cv_conf_df <- rbind(cv_conf_df, cbind(conf_df, Fold = i))
  
  conf_plot <- ggplot(conf_df, aes(x = Reference, y = Prediction, fill = round(Freq, 3))) +
    geom_tile(color = "white") +
    geom_text(aes(label = Freq), size = 5) +  # Show numbers inside the tiles
    scale_fill_gradient(low = "white", high = "blue") +
    labs(title = paste("Confusion Matrix Fold: ", i, sep=""), x = "Actual Class", y = "Predicted Class") +
    theme_minimal() +
    theme(axis.text.x = element_text(hjust = 1))
  #plot(conf_plot)
  
  # AUC
  #multi_class_auc <- multiclass.roc(response = res_df$Cluster, predictor = as.matrix(pred$predictions))
  #aucs = c(aucs, multi_class_auc$auc)
  
  rocs <- lapply(levels(res_df$Cluster), function(cls) {
    one_vs_rest <- factor(ifelse(res_df$Cluster == cls, 1, 0), levels=c(0,1))
    rocky <- roc(response = one_vs_rest, predictor = res_df[, cls], levels = c(0, 1), direction = "<")
    return(rocky)
  })
  
  roc_df <- do.call(rbind, lapply(seq_along(rocs), function(j) {
    data.frame(
      Sensitivity = rev(rocs[[j]]$sensitivities), 
      Specificity = rev(1 - rocs[[j]]$specificities), 
      Class = levels(res_df$Cluster)[j]
    )
  }))
  
  #roc_df$Specificity <- round(roc_df$Specificity, 2)
  #roc_df <- roc_df %>% group_by(Class, Specificity) %>% summarise(Sensitivity = mean(Sensitivity))
  
  spec_grid <- seq(0, 1, length.out = 101)
  roc_df %>% group_by(Class) %>% do({
    interp <- approx(
      x = .$Specificity,
      y = .$Sensitivity,
      xout = spec_grid,
      rule = 2
    )
    data.frame(
      Specificity = interp$x,
      Sensitivity = interp$y,
      Class = unique(.$Class)
    )
  }) -> roc_df
  
  
  rocAll <- rbind(rocAll, cbind(roc_df, Fold = i))
  
  auc_df <- data.frame(Class = levels(res_df$Cluster), 
                       AUC = sapply(rocs, function(roc) roc$auc))
  aucAll <- rbind(aucAll, cbind(auc_df, Fold = i))
  
  i = i + 1
}

cv_conf_df %>% pivot_wider(names_from = Fold, values_from = Freq) %>% 
  mutate(Mean = round(rowMeans(select(., where(is.numeric))), 3), 
         SD = round(apply(select(., where(is.numeric)), 1 , sd), 3)) -> cv_conf_df

rocAll %>% pivot_wider(names_from = Fold, values_from = Sensitivity, names_prefix = "Fold ") %>%
  ungroup() %>%
  mutate(TPR = rowMeans(select(., contains("Fold"))),
         Lower_CI = apply(select(., contains("Fold")), 1, function(x) quantile(x, 0.025)),
         Upper_CI = apply(select(., contains("Fold")), 1, function(x) quantile(x, 0.975))) -> rocAll

aucAll %>% group_by(Class) %>% summarise(AUC = mean(AUC)) -> aucs
aucs$Specificity <- 0.25 #seq(from = 0.25, to = 0.75, length.out = opt)
aucs$TPR <- rocAll %>% filter(Specificity == 0.25) %>% pull(TPR)


conf_p <- ggplot(cv_conf_df, aes(x = Reference, y = Prediction, fill = Mean)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste(Mean, SD, sep="\n ±")), size = 5) +  # Show numbers inside the tiles
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "", x = "Actual Class", y = "Predicted Class") +
  theme_minimal() +
  theme(legend.position = "none")
plot(conf_p)

mauc = mean(aucAll$AUC)
sdauc = sd(aucAll$AUC)
print(paste("Average Multi-label AUC", mauc, sep =": "))
print(paste("SD AUC", sdauc, sep =": "))

# For exp_km: AUC=0.7485 (0.0961)

acc = acc / length(folds)
print(paste("Average Accuracy", acc, sep =": "))

aucs$Class <- as.factor(aucs$Class)
rocp <- ggplot(rocAll, aes(x = Specificity, y = TPR, color = Class)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI, fill = Class), alpha = 0.2, color = NA) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(title = "", 
       x = "False Positive Rate", 
       y = "True Positive Rate") +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_color_paletteer_d("ggthemes::Tableau_10") +
  scale_fill_paletteer_d("ggthemes::Tableau_10") +
  scale_x_continuous(limits = c(0, 1.01), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1.01), expand = c(0,0)) +
  coord_fixed() +
  geom_text_repel(data = aucs, 
                  aes(label = paste0(Class, " (AUC = ", round(AUC, 2), ")"),
                      color = Class),
                  size = 4, 
                  nudge_x = 0.2,   # Adjust text position
                  nudge_y = -0.05, 
                  segment.size = 0.3) +
    guides(color = guide_legend(override.aes = list(label = ""))) # Remove text from legend
plot(rocp)


get_brain_regions_fn <- source("~/R/EDAP-R/get_brain_regions.R")$value
brain_regions <- get_brain_regions_fn()

varAll %>% group_by(Cluster) %>% summarise_all(mean) %>% 
  select(-Fold) %>% pivot_longer(cols = -Cluster, names_to = "MR.var", values_to = "Importance") -> summary_vars

sapply(summary_vars$MR.var, function(varname) {
  tail(strsplit(attr(mri_data[[varname]], 'label'), " ")[[1]],n=1)
}) -> summary_vars$Variables

brain_regions |> select(Variables, Lobe) |> left_join(summary_vars, by = "Variables") |> drop_na() -> summary_vars
summary_vars |> group_by(Cluster, Lobe) %>% summarise(Importance = mean(Importance)) -> summary_vars
  
varp <- summary_vars |> ggplot(aes(x = Lobe, y = Importance, fill = Cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Mean Importance of MR Variables by Lobe", x = "Lobe", y = "Mean Importance") +
  theme_minimal() +
  theme(legend.position = "right") +
  coord_flip() 

plot(varp)

ggsave(conf_p, filename = "~/R/EDAP-data/plots/LTC4/ltc_xgb_confmat.png", width = 5, height = 5, dpi = 300)
ggsave(rocp, filename = "~/R/EDAP-data/plots/LTC4/oasis_ltc_xgb_roc.png", width = 5, height = 4, dpi = 300)
#ggsave(varp, filename = "~/R/EDAP-plots/rf_varimp.png", width = 7, height = 10, dpi = 300)


