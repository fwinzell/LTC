library(ADNIMERGE2)
library(progmod)
library(tidyr)
library(dplyr)
library(tibble)
library(ClusterR)

library(progress)
library(purrr)

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
# Clustering functions
source("~/R/LTC/utils/clustering.R")

# Model fitting
source("~/R/LTC/utils/model_utils.R")

# Extra utils for clustering and visualization
source("~/R/LTC/utils/cluster_utils.R")

fit_inital = FALSE # set to FALSE to load previous initial model fitting
# 1. Load dataset
#multi_cohort_df_ <- read.csv("~/R/EDAP-data/MULTI_COHORT.csv", header = TRUE)
multi_cohort_df <- read.csv("~/R/EDAP-data/MULTI_COHORT_3.csv", header = TRUE)

# Filter out NACC
multi_cohort_df <- filter_out(multi_cohort_df, Cohort == "NACC")

all.vars <- c(grepv("^(RH_|LH_|CC_)", colnames(multi_cohort_df)), "BRAINSTEM")

# 3. Fit inital NLMMs
numCores <- detectCores()
num_workers <- max(1L, min(8L, numCores-2L))
cl <- parallel::makeCluster(num_workers)
doParallel::registerDoParallel(num_workers)

dont.run <- c()

if (fit_inital) {
  write_fst(multi_cohort_df, "~/R/LTC/tmp/multi_cohort_df.fst", compress=0)
  system.time(
    foreach(i = seq_along(all.vars)) %dopar% {
      varname = all.vars[i]
      if (varname %in% dont.run) next
      
      dsubset <- read_fst("~/R/LTC/tmp/multi_cohort_df.fst", columns = c("RID", "Time", "DX.bl", varname))
      dsubset <- dsubset %>% rename(y = varname, t = Time) %>% drop_na(y)
      result <- exp_nlmms_sample_fn(varname, dsubset, n_samples = 25, verbose=FALSE)
      
      saveRDS(result, file = paste0("~/R/LTC/tmp/out/", varname, ".rds"))
      NULL
    }
  )
  
  results <- lapply(all.vars, function(x) readRDS(paste0("~/R/LTC/tmp/out/", x, ".rds")))
  
  results <- results[!sapply(results, is.null)]
  
  betaList <- lapply(results, `[[`, "betas") 
  beta_df <- purrr::reduce(betaList, full_join, by = "RID")
  
  nlmmBasic <- list(
    betas = beta_df,
    bic = sapply(results, `[[`, "bic"),
    aic = sapply(results, `[[`, "aic"),
    logLikes = sapply(results, `[[`, "logLike")
  )
  
  save(nlmmBasic, file = "~/R/EDAP-data/LTC_MC/new/nlmmBasic_AO.Rdata")
} else {
  load("~/R/EDAP-data/LTC_MC/new/nlmmBasic_AO.Rdata")
}

cat("Inital models fitted to", round((dim(nlmmBasic$betas)[2]-1)/length(all.vars)*100), "% of variables") 
cat("Mean BIC: ", mean(nlmmBasic$bic))

#### KM, no cross-validation ####
for(ii in 1:1) {
  set.seed(ii)
  mc_clust <- multi_cohort_df #%>% select(RID, time_shift, M, DX.bl, Months, Years, all_of(vars))
  rids <- unique(mc_clust$RID)
  
  max_clusters <- 8
  mcr <- 0.5 # minimal convergence rate
  
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
  mc_clust <- mc_clust %>% left_join(clusterList[[curr_c]], by = "RID") 
  valid_vars <- colnames(nlmmBasic$betas)[-1]
  nlmmBest <- fit_cluster_nlmms_foreach(data=mc_clust, run.vars=valid_vars, n_samples = 25, parallell = TRUE)
  #nlmmBest <- fit_cluster_nlmms_foreach_adni(mc_clust, verbose = FALSE, restart = TRUE, n_samples = 25, parallell = TRUE)
  
  cat("2 Cluster models fitted to", round(length(nlmmBest$bic)/length(all.vars)*100), "% of variables \n") 
  #save(nlmmBest, file = "~/R/EDAP-data/LTC_MC/nlmmBest0.Rdata")
  
  while(length(clusterPairs)>0) {
    mc_clust <- multi_cohort_df %>% left_join(clusterList[[curr_c]], by = "RID") %>%
      drop_na(Cluster)
    #%>% select(RID, Month.bl, DX.bl, time_shift, Cluster, all_of(vars[1:25]))
    
    clabels <- unique(mc_clust$Cluster)
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
      ucsf_subset <- mc_clust %>% filter(Cluster == c)
      crids <- unique(ucsf_subset$RID)
      beta_subset <- nlmmBest$betas %>% filter(RID %in% crids) %>%
        mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
        select(RID, where(~ is.numeric(.x) && sd(.x, na.rm = TRUE) != 0))
        #mutate_all(~replace(., is.na(.), 0)) %>%
        #select(where(~ sd(.x, na.rm = TRUE) != 0))
      
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
        
        nlmmCandidates[[next_c]] <- mc_clust %>% select(-Cluster) %>% left_join(c_df, by = "RID") %>%
          fit_cluster_nlmms_foreach(run.vars=valid_vars, n_samples = 25, parallell = TRUE)
        
        cat("New model fitted to", round(length(nlmmCandidates[[next_c]]$bic)/length(valid_vars)*100), "% of variables \n") 
        
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
      
      nlmmCandidates[[next_c]] <- mc_clust %>% select(-Cluster) %>% left_join(c_df, by = "RID") %>%
        fit_cluster_nlmms_foreach(run.vars=valid_vars, n_samples = 25, parallell = TRUE)
      
      cat("New model fitted to", round(length(nlmmCandidates[[next_c]]$bic)/length(valid_vars)*100), "% of variables \n") 
      
      clusterList[[append_c+1]] <- c_df
      next_c = next_c + 1
    }
    
    # Find the best clustering
    bics <- sapply(nlmmCandidates, function(x) mean(x$bic))
    
    if (length(nlmmCandidates) > 1) {
      best_idx <- evaluate_bics_2(nlmmCandidates, valid_vars, min_conv_rate=mcr)
    } else {
      best_idx = 1
    }
  
    
    nlmmBest <- nlmmCandidates[[best_idx]]
    save(nlmmBest, file = "~/R/EDAP-data/LTC_MC/new/nlmmBest_AO.Rdata")
    
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
      new_best_c <- length(clusterList)+(best_idx-ncol(all.bics))
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
  
  table(cluster_df$Cluster)
  
  setClass("clusterObject",
           slots = c(
             RID = "factor",
             Cluster = "factor",
             varNames = "character",
             probs = "numeric",
             betas = "data.frame",
             BIC = "numeric",
             AIC = "numeric",
             ll = "numeric",
             tree = "matrix"
           ))
  
  multiLTC <- new("clusterObject",
                 RID = cluster_df$RID,
                 Cluster = cluster_df$Cluster,
                 varNames = colnames(nlmmBest$betas)[-1],
                 probs = NA_real_,
                 betas = data.frame(nlmmBest$betas),
                 BIC = nlmmBest$bic,
                 AIC = nlmmBest$aic,
                 ll = nlmmBest$logLikes,
                 tree = adjMat)
  
  save(multiLTC, file = "~/R/EDAP-data/LTC_MC/new/exp_km_ab_ao.Rdata")
}

for(i in 1:length(clusterList)){
  print(i)
  print(table(clusterList[[i]]$Cluster))
}

for(i in 1:length(treeIdx)){
  print(table(clusterList[[treeIdx[i]]]$Cluster))
}
treeIdx


plot_dendrogram(multiLTC, save=TRUE)


cluster_df <- mutate(cluster_df, Cohort = gsub("_.*", "", RID))

