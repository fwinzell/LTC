library(dplyr)
library(ggplot2)
library(Rtsne)
library(cluster)
library(mclust)

save_dir <- "~/R/EDAP-data/LTC_MC/new/"

run <- "exp_km_ab_ao_2"
load(paste("~/R/EDAP-data/LTC_MC/new/", run, ".Rdata", sep = ""))

Clusters <- data.frame(
  Cluster = multiLTC@Cluster,
  RID = multiLTC@RID
)


visualize <- function(df, param_cols, group_col="Cluster") {
  df <- filter_out(df, is.na(.data[[group_col]]))
  
  # Keep complete cases for this analysis
  df_complete <- df %>%
    filter(complete.cases(across(all_of(param_cols))))
  
  #### T-SNE analysis ####
  
  X <- scale(df_complete[, param_cols])
  
  set.seed(123)
  
  tsne <- Rtsne(
    X,
    dims = 2,
    perplexity = 30,
    pca = TRUE,
    check_duplicates = FALSE
  )
  
  tsne_df <- df_complete %>%
    select(RID, all_of(group_col)) %>%
    mutate(
      TSNE1 = tsne$Y[, 1],
      TSNE2 = tsne$Y[, 2]
    )
  
  p_tsne <- ggplot(tsne_df, aes(TSNE1, TSNE2, colour = .data[[group_col]])) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_paletteer_d("ggthemes::Tableau_10") +
    scale_fill_paletteer_d("ggthemes::Tableau_10") +
    theme_classic() +
    labs(
      x = "t-SNE 1",
      y = "t-SNE 2",
      colour = group_col
    )
  
  p_tsne
}


#### Try fitting new Naive cluster models #####
source("~/R/LTC/utils/model_utils.R")
library(purrr)

multi_cohort_df <- read.csv("~/R/EDAP-data/MULTI_COHORT_4.csv", header = TRUE)

# Filter out NACC
multi_cohort_df <- filter_out(multi_cohort_df, Cohort == "NACC")

all.vars <- c(grepv("^(RH_|LH_|CC_)", colnames(multi_cohort_df)), "BRAINSTEM")

# Load BioFINDER data

# Data loading
source("~/R/LTC/utils/biofinder_data_loaders.R")

biofinder_df <- get_mri_data_updated(normalize = TRUE) #%>% select(-OPTICCHIASM)

ab_df <- get_ab_df()

ab_pos_ids <- ab_df %>% group_by(sid) %>% mutate(AB_any = any(AB)) %>% ungroup() %>%
  filter(AB_any) %>% select(sid) %>% unlist()

biofinder_df <- filter(biofinder_df, sid %in% ab_pos_ids)

dpm_res <- read.csv("~/R/EDAP-data/BioFINDER/DPM_BioFINDER.csv") %>% distinct(sid, time_shift)

biofinder_df <- left_join(biofinder_df, dpm_res, by="sid") %>% filter(!is.na(time_shift)) %>%
  mutate(Time = Years + time_shift,
         Cohort = "BF2") %>%
  rename(RID = sid)

bf.vars <- c(grepv("^(RH_|LH_|CC_)", colnames(biofinder_df)), "BRAINSTEM", "OPTICCHIASM")
mri.vars <- intersect(all.vars, bf.vars)

common_cols <- intersect(names(multi_cohort_df), names(biofinder_df))

multi_cohort_df <- rbind(multi_cohort_df[, common_cols], biofinder_df[, common_cols])

multi_cohort_df <- mutate(multi_cohort_df, 
                          Cohort = gsub("[0-9_]", "", RID)) %>%
  mutate(Cohort = factor(ifelse(Cohort == "BF", "BF2", "ADNI/OASIS")))

# Fit/load models
if (TRUE) {
  results <- list()
  for (i in seq_along(mri.vars)) {
    varname = mri.vars[i]
    dsubset <- select(multi_cohort_df, RID, Cohort, Time, all_of(varname)) %>%
      rename(t=Time, y=varname) %>% drop_na(y)
    results[[i]] <- exp_nlmms_cohort_params(varname, dsubset, n_samples=10, verbose=FALSE)
  }

  save(results, file = "~/R/EDAP-data/LTC_MC/cluster_naive_models_2.Rdata") 
} else {
  load("~/R/EDAP-data/LTC_MC/cluster_naive_models.Rdata")  
}


results <- results[!sapply(results, is.null)]

randList <- lapply(results, `[[`, "random")
fixList <- lapply(results, `[[`, "fixed")

rand_df <- purrr::reduce(randList, full_join, by = "RID")

# Do not do this...
func_slope <- data.frame(RID=rand_df$RID)
for (i in seq_along(fixList)) {
  fix <- fixList[[i]]
  rnd <- randList[[i]]
  slop <- data.frame(RID = rnd$RID,
                     -exp(fix[grepl("l", names(fix))])*exp(exp(fix[grepl("g", names(fix))])*(1+rnd[grepl("gi", names(rnd))]))
                     )
  func_slope <- left_join(func_slope, slop, by="RID")
}

func_slope

# Run assingment test

df <- rand_df %>%
  full_join(Clusters, by = "RID") %>%
  mutate(Cohort = gsub("[0-9_]", "", RID)) %>%
  mutate(Cluster = factor(Cluster),
         Cohort = factor(Cohort)) 

gi_cols <- grepv(".gi", colnames(df))
df <- select(df, RID, Cluster, Cohort, all_of(gi_cols))

visualize(df, gi_cols, "Cluster")
visualize(df, gi_cols, "Cohort")

# Parameter columns
param_cols <- gi_cols #setdiff(names(func_slope), "RID")

# Keep complete cases for this analysis
df_complete <- df %>%
  filter(complete.cases(across(all_of(param_cols))))

set.seed(123)

run_assignment_test <- function(df_complete, 
                                classifier="centroid", k=5, 
                                scale_params=TRUE) {
  # 80/20 split
  # Stratified split is preferable
  train_idx <- unlist(
    lapply(
      split(seq_len(nrow(df_complete)), df_complete$Cluster),
      function(idx) {
        sample(idx, floor(0.8 * length(idx)))
      }
    )
  )
  
  train <- df_complete[train_idx, ]
  test  <- df_complete[-train_idx, ]
  
  # Scale using training set only
  train_means <- sapply(train[, param_cols], mean)
  train_sds   <- sapply(train[, param_cols], sd)
  
  # Avoid division by zero
  train_sds[train_sds == 0] <- 1
  
  if (scale_params) {
    X_train <- scale(
      train[, param_cols],
      center = train_means,
      scale = train_sds
    )
    
    X_test <- scale(
      test[, param_cols],
      center = train_means,
      scale = train_sds
    ) 
  } else {
    X_train <- train[, param_cols]
    X_test <- test[, param_cols]
  }
  
  if (classifier == "centroid") {
    # Centroids
    centroids <- sapply(
      levels(df_complete$Cluster),
      function(cl) {
        colMeans(
          X_train[train$Cluster == cl, , drop = FALSE]
        )
      }
    )
    
    centroids <- t(centroids)
    
    # Distances
    distance_matrix <- sapply(
      seq_len(nrow(centroids)),
      function(k) {
        rowSums(
          (sweep(X_test, 2, centroids[k, ], "-"))^2
        )^0.5
      }
    )
    
    colnames(distance_matrix) <- rownames(centroids)
    
    predicted <- colnames(distance_matrix)[
      max.col(-distance_matrix)
    ]
    
    test_results <- test %>%
      select(RID, Cluster) %>%
      mutate(
        Predicted = factor(predicted, levels = levels(Cluster)),
        Distance = apply(distance_matrix, 1, min),
        Margin = apply(distance_matrix, 1, function(x) {
          sx <- sort(x)
          sx[2] - sx[1]
        }
        )
      )
  } else if (classifier == "knn") {
    # kNN
    predicted <- knn(
      train = X_train,
      test = X_test,
      cl = train$Cluster,
      k = k
    )
    
    test_results <- test %>%
      select(RID, Cluster) %>%
      mutate(
        Original = test$Cluster,
        Predicted = factor(predicted, levels = levels(Cluster))
      )
  } else {
    test_results = NULL
  }
  
  return(test_results)
}

run_cohort_test <- function(df_complete) {
  # 80/20 split
  # Stratified split is preferable
  train_idx <- unlist(
    lapply(
      split(seq_len(nrow(df_complete)), df_complete$Cohort),
      function(idx) {
        sample(idx, floor(0.8 * length(idx)))
      }
    )
  )
  
  train <- df_complete[train_idx, ]
  test  <- df_complete[-train_idx, ]
  
  # Scale using training set only
  train_means <- sapply(train[, param_cols], mean)
  train_sds   <- sapply(train[, param_cols], sd)
  
  # Avoid division by zero
  train_sds[train_sds == 0] <- 1
  
  X_train <- scale(
    train[, param_cols],
    center = train_means,
    scale = train_sds
  )
  
  X_test <- scale(
    test[, param_cols],
    center = train_means,
    scale = train_sds
  )
  
  # Centroids
  centroids <- sapply(
    levels(df_complete$Cohort),
    function(co) {
      colMeans(
        X_train[train$Cohort == co, , drop = FALSE]
      )
    }
  )
  
  centroids <- t(centroids)
  
  # Distances
  distance_matrix <- sapply(
    seq_len(nrow(centroids)),
    function(k) {
      rowSums(
        (sweep(X_test, 2, centroids[k, ], "-"))^2
      )^0.5
    }
  )
  
  colnames(distance_matrix) <- rownames(centroids)
  
  predicted <- colnames(distance_matrix)[
    max.col(-distance_matrix)
  ]
  
  test_results <- test %>%
    select(RID, Cohort) %>%
    mutate(
      Predicted = factor(predicted, levels = levels(Cohort)),
      Distance = apply(distance_matrix, 1, min),
      Margin = apply(distance_matrix, 1, function(x) {
        sx <- sort(x)
        sx[2] - sx[1]
      }
      )
    )
  return(test_results)
}

get_assignment_metrics <- function(distance_matrix) {
  
  sorted <- t(apply(distance_matrix, 1, sort))
  
  data.frame(
    min_distance = sorted[, 1],
    second_distance = sorted[, 2],
    margin = sorted[, 2] - sorted[, 1],
    relative_margin =
      (sorted[, 2] - sorted[, 1]) / sorted[, 1]
  )
}

run_centroid_cv <- function(df_complete, k_folds = 5) {
  df_complete <- df_complete %>% filter_out(is.na(Cluster)) 
  subjects <- df_complete %>% select(RID, Cluster) %>% unique() 
  
  folds <- caret::createFolds(subjects$Cluster, k = k_folds)
  
  resultList <- list()
  for (ii in 1:k_folds) {
    train = df_complete[-folds[[ii]], ]
    val = df_complete[folds[[ii]], ]
    
    # Scale using training set only
    train_means <- sapply(train[, param_cols], mean)
    train_sds   <- sapply(train[, param_cols], sd)
    
    # Avoid division by zero
    train_sds[train_sds == 0] <- 1
    
    X_train <- scale(train[, param_cols], center = train_means, scale = train_sds)
      
    X_val <- scale(val[, param_cols], center = train_means, scale = train_sds)
    
    # Centroids
    centroids <- sapply(
      levels(df_complete$Cluster),
      function(cl) {
        colMeans(X_train[train$Cluster == cl, , drop = FALSE]) 
        })
      
    centroids <- t(centroids)
      
    # Distances
    distance_matrix <- sapply(
      seq_len(nrow(centroids)),
      function(k) {
        rowSums((sweep(X_val, 2, centroids[k, ], "-"))^2)^0.5
      }
    )
      
    colnames(distance_matrix) <- rownames(centroids)
    
    metrics <- get_assignment_metrics(distance_matrix)
    
    predicted <- colnames(distance_matrix)[max.col(-distance_matrix)]
      
    val_results <- val %>%
      select(RID, Cluster) %>%
      mutate(
        Predicted = factor(predicted, levels = levels(Cluster)),
        Correct = Predicted == Cluster
      ) %>% cbind(metrics)
    
    resultList[[ii]] <- val_results
  }
  
  results_df <- do.call("rbind", resultList)
  
  gp1 <- ggplot(results_df,
         aes(relative_margin, fill = Correct)) +
    geom_density(alpha = 0.4) +
    theme_classic() +
    labs(
      x = "Relative distance margin",
      y = "Density"
    )
  
  gp2 <- ggplot(results_df,
         aes(min_distance, fill = Correct)) +
    geom_density(alpha = 0.4) +
    theme_classic() +
    labs(
      x = "Distance to nearest centroid",
      y = "Density"
    )
  
  plot(gp1)
  plot(gp2)
  
  return(results_df)
}

df_cohort <- df_complete %>% mutate(Cohort = factor(ifelse(Cohort == "BF", "BF2", "ADNI/OASIS")))

bap <- run_cohort_test(df_cohort)
table(Cohort = bap$Cohort, Predicted = bap$Predicted)

library(randomForest)

rf <- randomForest(
  x = df_cohort[, param_cols],
  y = df_cohort$Cohort,
  ntree = 500
)

rf

test_results <- run_assignment_test(df_complete)
table(
  Original = test_results$Cluster,
  Predicted = test_results$Predicted
)

# Sensitivity / recall for each original cluster
test_results %>%
  group_by(Cluster) %>%
  summarise(
    n = n(),
    correct = sum(Cluster == Predicted),
    recall = mean(Cluster == Predicted)
  )

n_reps <- 100

results <- vector("list", n_reps)
df_complete <- filter_out(df_complete, is.na(Cluster))

for (r in seq_len(n_reps)) {
  
  test_results <- run_assignment_test(df_complete, classifier = "centroid")
  
  recall <- test_results %>%
    group_by(Cluster) %>%
    summarise(
      recall = mean(Cluster == Predicted)
    ) %>% pivot_wider(names_from = Cluster, values_from = recall, names_prefix = "recall.")
  
  results[[r]] <- data.frame(
    Rep = r,
    Accuracy = mean(test_results$Cluster == test_results$Predicted),
    ARI = adjustedRandIndex(test_results$Cluster, test_results$Predicted),
    recall
  )
}

results_df <- bind_rows(results)
results_df %>%
  summarise(
    mean_accuracy = mean(Accuracy),
    sd_accuracy = sd(Accuracy),
    mean_ARI = mean(ARI),
    sd_ARI = sd(ARI),
    median_ARI = median(ARI),
    lower_ARI = quantile(ARI, 0.025),
    upper_ARI = quantile(ARI, 0.975),
    
    mean_recall.A = mean(recall.A),
    mean_recall.B = mean(recall.B),
    mean_recall.C = mean(recall.C),
    mean_recall.D = mean(recall.D)
  )

visualize(params_df = select(df, -Cluster), Clusters = Clusters)


# Try kNN classifier instead

library(class)
library(mclust)

set.seed(123)

n_reps <- 100

results <- vector("list", n_reps)

for (r in seq_len(n_reps)) {
  
 test_results <- run_assignment_test(df_complete, classifier = "knn", k=1)
  
 results[[r]] <- data.frame(
   Rep = r,
   Accuracy = mean(test_results$Cluster == test_results$Predicted),
   ARI = adjustedRandIndex(test_results$Cluster, test_results$Predicted)
 )
}

results_knn <- bind_rows(results)

results_knn %>%
  summarise(
    mean_accuracy = mean(Accuracy),
    sd_accuracy = sd(Accuracy),
    mean_ari = mean(ARI),
    sd_ari = sd(ARI)
  )



# calculate distance to centroids
X <- scale(df_complete[, param_cols])

centroids_all <- sapply(
  levels(df_complete$Cluster),
  function(cl) {
    colMeans(
      X[df_complete$Cluster == cl, , drop = FALSE]
    )
  }
) %>%
  t()

centroid_dist <- sapply(
  seq_len(nrow(X)),
  function(i) {
    cl <- as.character(df_complete$Cluster[i])
    
    sqrt(sum(
      (X[i, ] - centroids_all[cl, ])^2
    ))
  }
)

distance_df <- data.frame(
  Cluster = df_complete$Cluster,
  Distance = centroid_dist
)

ggplot(distance_df, aes(Cluster, Distance)) +
  geom_boxplot() +
  theme_classic() +
  labs(
    x = "Cluster",
    y = "Distance to cluster centroid"
  )


# Centroid distance analysis
df_complete <- df %>%
  filter(complete.cases(across(all_of(param_cols))))


cv_predictions <- run_centroid_cv(df_complete, k_folds = 6) 


thresholds <- seq(0, 1, by = 0.02)

threshold_results <- lapply(thresholds, function(th) {
  
  retained <- cv_predictions$relative_margin >= th
  
  data.frame(
    threshold = th,
    coverage = mean(retained),
    accuracy = mean(
      cv_predictions$Correct[retained]
    )
  )
}) %>%
  bind_rows()

ggplot(threshold_results,
       aes(coverage, accuracy)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  labs(
    x = "Proportion assigned",
    y = "Accuracy among assigned participants"
  )




