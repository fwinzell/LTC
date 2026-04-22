library(dplyr)
library(ClusterR)

betaKM <- function(beta_df, k=2, minClusterSize=10) {
  #beta.idxs <- grepl("^ST\\d+(CV|SV|HS|SA|TA)$", colnames(beta_df))
  beta.vars <- colnames(beta_df)[-1]
  rids <- unique(beta_df$RID)
  
  z_df <- beta_df
  z_df <- scale(beta_df[, beta.vars]) 
  # If a column is constant, scale generates NaNs, remove any of those columns
  z_df <- z_df[, !apply(is.na(z_df), 2, any)]
  
  clusters <- kmeans(z_df, centers=k, iter.max=25, nstart=100)
  
  # Find which cluster is larger
  largest_c <- as.numeric(names(which.max(table(clusters$cluster))))
  
  # Create a mapping so that the largest cluster becomes 1
  new_cluster_labels <- ifelse(clusters$cluster == largest_c, 1, 2)
  clusters$cluster <- new_cluster_labels
  
  cluster_df <- data.frame(RID = beta_df$RID, Cluster = clusters$cluster)
  
  # Set too small clusters to 0 - noise
  for(ki in 1:k) {
    if(clusters$size[ki] < minClusterSize) {
      cluster_df <- cluster_df %>% mutate(
        Cluster = ifelse(Cluster == ki, 0, Cluster))
    }
  }
  
  return(cluster_df)
}

random_clusters <- function(beta_df, k=2, minClusterSize=10) {
  cluster_labels <- LETTERS[1:k]
  rids <- beta_df$RID
  cluster_df <- data.frame(RID = rids, Cluster = sample(cluster_labels, length(rids), replace=TRUE))
  min_size = min(table(cluster_df$Cluster))
  if (min_size < minClusterSize) {
    cluster_df <- random_clusters(beta_df, k=k, minClusterSize = minClusterSize)
  } else {
    return(cluster_df)
  }
}

betaGMM <- function(beta_df, k=2, minClusterSize=10) {
  beta.idxs <- grepl("^ST\\d+(CV|SV|HS|SA|TA)$", colnames(beta_df))
  beta.vars <- colnames(beta_df)[beta.idxs == 1]
  rids <- unique(beta_df$RID)
  
  z_df <- beta_df
  z_df <- scale(beta_df[, beta.vars])
  # If a column is constant, scale generates NaNs, remove any of those columns
  z_df <- z_df[, !apply(is.na(z_df), 2, any)]
  
  gmm <- GMM(z_df, k, "maha_dist", "random_subset", km_iter=50, em_iter=100, verbose=FALSE)
  
  clusters <- predict_GMM(z_df, gmm$centroids, gmm$covariance_matrices, gmm$weights)
  
  probs <- apply(clusters$cluster_proba, 1, max)
  
  cluster_df <- data.frame(RID = beta_df$RID, Cluster = clusters$cluster_labels, 
                           P = probs)
  
  # Set too small clusters to 0 - noise
  for(ki in 1:k) {
    if(table(cluster_df$Cluster)[ki] < minClusterSize) {
      cluster_df <- cluster_df %>% mutate(
        Cluster = ifelse(Cluster == ki, 0, Cluster))
    }
  }
  
  return(cluster_df) 
}