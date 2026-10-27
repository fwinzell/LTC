library(ADNIMERGE2)
library(progmod)
library(tidyr)
library(dplyr)
library(tibble)
library(ClusterR)
library(kml3d)
library(progress)
library(purrr)
library(mclust)
library(caret)
library(stringr)

library(ggplot2)
library(ggpubr)
library(lubridate)

library(parallel)
library(foreach)
library(doParallel)
library(fst)

# Useful functions for loading data and plotting trajectories
source("~/R/LTC/utils/ucsf_data_loaders.R")

# Clustering functions
source("~/R/LTC/utils/clustering.R")

# Model fitting
source("~/R/LTC/utils/model_utils.R")

# Extra utils for clustering and visualization
source("~/R/LTC/utils/cluster_utils.R")

create_strat_folds <- function(df, k = 5) {
  df %>% select(RID, DX.bl) %>% unique() -> subjects  
  
  folds <- caret::createFolds(subjects$DX.bl, k = 5)
  
  idx_folds <- lapply(folds, function(rids) {
    which(df$RID %in% subjects[rids, "RID"])
  })
  
  return(idx_folds)
}

#### Initial steps ####
fit_mcdp = FALSE # set to FALSE to load previous MCDP run
fit_inital = FALSE # set to FALSE to load previous initial model fitting
parallell = TRUE # set to TRUE for parallell execution with foreach
# 1. Load datasets
ucsf_all <- ucsf_longitudinal_all(only_vol=TRUE, filter_n=0)

# Find AB positives
get_ab_df <- source("~/R/LTC/get_ab_df.r")$value
ab_df <- get_ab_df()

# Any positive Abeta biomarker at any time point
ab_df %>% select(RID, AB_any) %>% group_by(RID) %>% 
  mutate(AB_any = if_any(AB_any)) %>% distinct() %>% filter(AB_any) %>% 
  select(RID) %>% unlist() -> ab_pos_rids_any

ucsf_ab <- ucsf_all %>% filter(RID %in% ab_pos_rids_any)

length(unique(ucsf_ab$RID))

# 2. Fit MCDP model to estimate time-shift
if (fit_mcdp) {
  fit_dpm2 <- source("~/R/LTC/latent-time-shift/fit_dpm.R")$value
  dpm_df <- fit_dpm2(scale_t=FALSE, scale_y=FALSE)
  write.csv("~/R/EDAP-data/ADNI/DPM_ADNI.csv")
} else {
  dpm_df <- read.csv("~/R/EDAP-data/ADNI/DPM_ADNI.csv", header=TRUE)
}

# Transform MRI data to latent time scale
dpm_df |> select(RID, time_shift) |>
  unique() |> inner_join(ucsf_ab, join_by(RID)) -> ucsf_ab

rids <- unique(ucsf_ab$RID)

ucsf_ab$Time = ucsf_ab$Years + ucsf_ab$time_shift

all.idxs <- grepl("^ST\\d+(CV|SV|SA|TA)$", colnames(ucsf_ab)) 
all.vars <- colnames(ucsf_ab)[all.idxs == 1]

# 3. Fit inital NLMMs
numCores <- detectCores()
registerDoParallel(numCores)

if (fit_inital) {
  write_fst(ucsf_ab, "~/R/LTC/tmp/ucsf_ab.fst", compress=0)
  system.time(
    foreach(i = seq_along(all.vars)) %dopar% {
      varname = all.vars[i]
      dsubset <- read_fst("~/R/LTC/tmp/ucsf_ab.fst", columns = c("RID", "Time", "DX.bl", varname))
      dsubset <- dsubset %>% rename(y = varname, t = Time) %>% drop_na(y)
      result <- exp_nlmms_sample_fn(varname, dsubset, n_samples = 25, verbose=FALSE)
      
      saveRDS(result, file = paste0("~/R/LTC/tmp/out/", varname, ".rds"))
      NULL
    }
  )
  
  results <- lapply(all.vars, function(x) readRDS(paste0("~/R/tmp/out/", x, ".rds")))
  
  results <- results[!sapply(results, is.null)]
  
  betaList <- lapply(results, `[[`, "betas") 
  beta_df <- purrr::reduce(betaList, full_join, by = "RID")
  
  nlmmBasic <- list(
    betas = beta_df,
    bic = sapply(results, `[[`, "bic"),
    aic = sapply(results, `[[`, "aic"),
    logLikes = sapply(results, `[[`, "logLike")
  )
  
  save(nlmmBasic, file = "~/R/EDAP-data/LTC4/nlmmBasic.Rdata")
} else {
  load("~/R/EDAP-data/LTC4/nlmmBasic.Rdata")
}

cat("Inital models fitted to", round(dim(nlmmBasic$betas)[2]/length(all.vars)*100), "% of variables") 
cat("Mean BIC: ", mean(nlmmBasic$bic))



#### KM, cross-validation ####
folds <- create_strat_folds(ucsf_ab, k=5)

for(ii in 1:1) {

  ucsf_clust <- ucsf_ab[-folds[[ii]], ]
  rids <- unique(ucsf_clust$RID)
  
  max_clusters <- 8
  
  clusterList <- list()
  clusterList[[1]] <- data.frame(RID = rids, Cluster = 1, probs=NA) # One cluster
  beta_0 <- nlmmBasic$betas %>% filter(RID %in% rids) %>%
    mutate_all(~replace(., is.na(.), 0))

  clusterList[[2]] <- betaKM(beta_0) # First clustering
  table(clusterList[[2]]$Cluster)
  clusterPairs <- list(c(1,2)) # Keep track of active cluster pairs that can be split
  curr_c = 2
  
  treeIdx <- c(1,2)
  
  #for(c in unique(cluster_df$Cluster)) {
  #  clusterList[[c]] <- cluster_df$RID[which(cluster_df$Cluster == c)]
  #}
  ucsf_clust <- ucsf_clust %>% left_join(clusterList[[curr_c]], by = "RID") 
  nlmmBest <- fit_cluster_nlmms_foreach(ucsf_clust, parallell=parallell)
  
  cat("2 Cluster models fitted to", round(length(nlmmBest$bic)/length(all.vars)*100), "% of variables \n") 
  
  while(length(clusterPairs)>0) {
    ucsf_clust <- ucsf_ab %>% left_join(clusterList[[curr_c]], by = "RID") %>%
      drop_na(Cluster)
    #%>% select(RID, Month.bl, DX.bl, time_shift, Cluster, all_of(vars[1:25]))
    
    clabels <- unique(ucsf_clust$Cluster)
    next_label <- max(clabels)+1
    if (next_label > max_clusters) {
      print("Maximum number of clusters reached without convergence")
      break
    }
    
    nlmmCandidates <- list()
    nlmmCandidates[[1]] <- nlmmBest
    next_c = 2
    
    # Collect cluster candidates
    newPairs = c()
    # Pop the first cluster pair
    pair = clusterPairs[[1]]
    for(c in pair){
      cat("Splitting cluster ", c, " of ", pair, "\n")
      ucsf_subset <- ucsf_clust %>% filter(Cluster == c)
      crids <- unique(ucsf_subset$RID)
      beta_subset <- nlmmBest$betas %>% filter(RID %in% crids) %>%
        mutate_all(~replace(., is.na(.), 0)) %>%
        select(where(~ sd(.x, na.rm = TRUE) != 0))
      
      c_df <- betaKM(beta_subset)
      print(table(c_df$Cluster))
      
      # Only one cluster
      if (any(c_df$Cluster == 0)) {
        print("only one cluster - continues")
        #nlmmCandidates[[next_c]] <- nlmmCandidates[[1]]
        #next_c = next_c + 1
      } else {
        c_df <- c_df %>% mutate(Cluster = ifelse(Cluster == 1, c, next_label))
        
        newPairs <- c(newPairs, list(c(c, next_label)))
        
        next_label = next_label+1 # Update next label
        # Add "Old" branch from other cluster to the clustering
        prev_c_df <- clusterList[[curr_c]] %>% filter(!(RID %in% crids))
        #c_df <- c_df %>% select(-P) %>% rbind(prev_c_df)
        c_df <- rbind(c_df, prev_c_df)
        
        nlmmCandidates[[next_c]] <- ucsf_clust %>% select(-Cluster) %>% left_join(c_df, by = "RID") %>%
          fit_cluster_nlmms_foreach(n_samples=5, parallell=parallell)
        
        append_c <- length(clusterList)+1
        clusterList[[append_c]] <- c_df
        next_c = next_c + 1
      }
    }
    # All combined if we have done two splits
    if(length(newPairs) == 2) {
      # Merge previous two clustering
      print(paste("Evaluating 4 cluster solution with of pairs", newPairs[1], newPairs[2]))
      prev_c_df <- clusterList[[append_c-1]] %>% filter(!(RID %in% crids))
      curr_c_df <- clusterList[[append_c]] %>% filter(RID %in% crids)
      c_df <- rbind(prev_c_df, curr_c_df)
      
      nlmmCandidates[[next_c]] <- ucsf_clust %>% select(-Cluster) %>% left_join(c_df, by = "RID") %>%
        fit_cluster_nlmms_foreach(n_samples=5, restart = FALSE, parallell=parallell)
      clusterList[[append_c+1]] <- c_df
      next_c = next_c + 1
    }
    
    # Find the best clustering
    #bics <- sapply(nlmmCandidates, function(x) mean(x$bic))
    
    best_idx <- evaluate_bics(nlmmCandidates)
    
    nlmmBest <- nlmmCandidates[[best_idx]]
    
    # Remove this pair
    clusterPairs <- clusterPairs[-1]
    
    if(best_idx == 1) {
      # do nothing, we are done with this pair
      new_best_c <- curr_c
    } else if(best_idx == 4) {
      # add both clusters to the list of cluster pairs
      new_best_c <- length(clusterList)
      clusterPairs <- c(clusterPairs, newPairs)
      cat("Added new pairs: ", newPairs[[1]], " and ", newPairs[[2]], "\n")
      treeIdx <- c(treeIdx, new_best_c)
    } else {
      # add only the new cluster pair
      new_best_c <- length(clusterList)+(best_idx-length(nlmmCandidates))
      clusterPairs <- c(clusterPairs, newPairs[best_idx-1])
      cat("Added new pair: ", newPairs[[best_idx-1]], "\n")
      treeIdx <- c(treeIdx, new_best_c)
    }
    
    # Make a better tree solution?
    curr_c <- new_best_c
    print("Current best:")
    print(table( clusterList[[curr_c]]$Cluster ))
  }
  
  cluster_df <- clusterList[[new_best_c]]
  
  # Assign new letter labels 
  ctab <- table(cluster_df$Cluster)
  mapping <- setNames(LETTERS[seq_along(ctab)], names(ctab))
  cluster_df$Cluster <- as.factor(mapping[as.character(cluster_df$Cluster)])
  
  # Make tree
  #treeIdx <- treeIdx[1:length(treeIdx)-1]
  adjMat <- matrix(data = NA, nrow = length(ctab), ncol = length(treeIdx),
                   dimnames = list(levels(cluster_df$Cluster), treeIdx))
  
  for (j in treeIdx) {
    tab <- table(clusterList[[j]]$Cluster)
    adjMat[1:length(tab), as.character(j)] <- tab
  }
  
  # Assign new letter labels - with decreasing frequency
  #freq <- sort(table(cluster_df$Cluster), decreasing = TRUE)
  #mapping <- setNames(LETTERS[seq_along(freq)], names(freq))
  #cluster_df$Cluster <- as.factor(mapping[as.character(cluster_df$Cluster)])
  
  table(cluster_df$Cluster)
  
  setClass("clusterObject",
           slots = c(
             RID = "numeric",
             Cluster = "factor",
             varNames = "character",
             probs = "numeric",
             betas = "data.frame",
             BIC = "numeric",
             AIC = "numeric",
             ll = "numeric",
             tree = "matrix"
           ))
  
  adniLTC <- new("clusterObject",
                 RID = cluster_df$RID,
                 Cluster = cluster_df$Cluster,
                 varNames = colnames(nlmmBest$betas)[-1],
                 probs = NA_real_,
                 betas = data.frame(nlmmBest$betas),
                 BIC = nlmmBest$bic,
                 AIC = nlmmBest$aic,
                 ll = nlmmBest$logLikes,
                 tree = adjMat)
  
  save(adniLTC, file = paste0("~/R/EDAP-data/LTC4/exp_km_ab_fold_", ii, ".Rdata"))
}



