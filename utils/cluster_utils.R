library(dplyr)
library(ggplot2)

evaluate_bics <- function(nlmmCandidates, all.vars) {
  all.bics <- sapply(nlmmCandidates, function(x) {
    idx <- which(all.vars %in% names(x$betas))
    bics <- rep(NA, length(all.vars))
    bics[idx] <- x$bic
    return(bics)
  } ) 
  
  
  # Check if any is significantly better than parent 
  p_values <- lapply(2:ncol(all.bics), function(col) {
    df <- na.omit(cbind(all.bics[, 1], all.bics[, col]))
    diff <- df[, 2] - df[, 1]
    wilcox.test(diff, alternative = "less")$p.value
  })
  
  better_idx <- which(p_values < 0.05)
  
  # If any is significantly better, select the best
  if (length(better_idx) == 0) {
    best_idx <- 1
  } else if (length(better_idx) > 1) {
    pairs <- combn(better_idx+1, 2, simplify = FALSE)
    best_idx <- lapply(pairs, function(cols) {
      idx = which.min(colMeans(na.omit(all.bics[, cols])))
      cols[idx]
    })
    best_idx <- as.numeric(names(which.max(table(unlist(best_idx)))))
  } else {
    best_idx <- better_idx
  }
  
  return(best_idx)
}

evaluate_bics_2 <- function(nlmmCandidates, all.vars, min_conv_rate=0.90) {
  if (length(nlmmCandidates) < 2) {
    return(1)
  }
  
  all.bics <- sapply(nlmmCandidates, function(x) {
    idx <- which(all.vars %in% names(x$betas))
    bics <- rep(NA, length(all.vars))
    bics[idx] <- x$bic
    return(bics)
  }, simplify = "matrix" ) 
  
  # Remove any clusterings that did not reach sufficient convergence rate
  conv_rates <- sapply(2:ncol(all.bics), function(col) {
    x <- all.bics[,col]
    1-length(which(is.na(x)))/length(all.vars)
  })
  idxs <- c(1, which(conv_rates>min_conv_rate)+1)
  all.bics <- all.bics[,idxs, drop = FALSE]
  
  if (ncol(all.bics) == 1) {
    return(1)
  }
  
  # Check if any is significantly better than parent 
  p_values <- lapply(2:ncol(all.bics), function(col) {
    df <- na.omit(cbind(all.bics[, 1], all.bics[, col]))
    diff <- df[, 2] - df[, 1]
    wilcox.test(diff, alternative = "less")$p.value
  })
  
  better_idx <- which(p_values < 0.05)
  
  # If any is significantly better, select the best
  if (length(better_idx) == 0) {
    best_idx <- 1
  } else if (length(better_idx) > 1) {
    pairs <- combn(better_idx+1, 2, simplify = FALSE)
    best_idx <- lapply(pairs, function(cols) {
      idx = which.min(colMeans(na.omit(all.bics[, cols])))
      cols[idx]
    })
    best_idx <- as.numeric(names(which.max(table(unlist(best_idx)))))
  } else {
    best_idx <- better_idx
  }
  
  return(best_idx)
}

plot_dendrogram <- function(LTCobj, save=TRUE, file.name="dendro") {
  # Generate dendrogram
  # Can move this to separate script?
  
  segments <- data.frame(x=0, y=1, xend=0, yend=2)
  segments <- rbind(segments,
                    data.frame(x=-LTCobj@tree[1,2], y=2, xend=LTCobj@tree[2,2], yend=2),
                    data.frame(x=-LTCobj@tree[1,2], y=2, xend=-LTCobj@tree[1,2], yend=3),
                    data.frame(x=LTCobj@tree[2,2], y=2, xend=LTCobj@tree[2,2], yend=3))
  
  prev_level <- na.omit(LTCobj@tree[,2])
  nextIdx = 3
  activeClusters = names(prev_level)
  xValues = setNames(c(-LTCobj@tree[1,2], LTCobj@tree[2,2]), activeClusters)
  inactiveClusters = list()
  negativeBranch = c("A")
  
  for (level in 3:ncol(LTCobj@tree)) {
    new_level = LTCobj@tree[,level]
    for (c in activeClusters) {
      if (new_level[c] != prev_level[c]) {
        # branch has branched
        new_pair = c(LTCobj@tree[c,level],LTCobj@tree[nextIdx,level])
        new_label = names(LTCobj@tree[,level])[nextIdx]
        prev_x = xValues[c]
        if (c %in% negativeBranch) {
          negativeBranch <- c(negativeBranch, new_label)
        } 
        
        segments <- rbind(segments,
                          data.frame(x=prev_x-new_pair[1], y=level, xend=prev_x+new_pair[2], yend=level),
                          data.frame(x=prev_x-new_pair[1], y=level, xend=prev_x-new_pair[1], yend=level+1),
                          data.frame(x=prev_x+new_pair[2], y=level, xend=prev_x+new_pair[2], yend=level+1))
        nextIdx = nextIdx+1
        activeClusters <- c(activeClusters, new_label)
        xValues[[c]] = prev_x-new_pair[1]
        xValues[[new_label]] = prev_x+new_pair[2]
      } else {
        # No branch, add to inactive clusters, remove from active
        inactiveClusters[[c]] <- c(xValues[c], level)
        activeClusters <- setdiff(activeClusters, c)
      }
    }
    prev_level = new_level
  }
  
  max_level = max(segments$yend)
  # add remaining segments for inactive clusters
  for (cxy in inactiveClusters) {
    segments <- rbind(segments, 
                      data.frame(x=cxy[1], y=cxy[2], xend=cxy[1], yend=max_level))
  }
  
  x_ticks <- segments %>% filter(yend == max_level) %>% select(x)
  x_labels <- unique(c(negativeBranch, rownames(LTCobj@tree)))
  x_labels <- sprintf("%s (n=%d)", x_labels, LTCobj@tree[x_labels,ncol(LTCobj@tree)])
  p <- ggplot(segments) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    scale_y_reverse(limits=c(max_level+0.3, 1)) +
    theme_classic() +
    #scale_x_discrete("Clusters", breaks = xbr$x, labels=c("A", "D", "C", "E", "B"), limits = c("0", "1500"), position = "bottom") 
    # Disable default x-axis elements
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    annotate("text",
             x = unlist(x_ticks), 
             y = max_level+0.2,        # slightly above the line (since y is reversed)
             label = unlist(x_labels),
             vjust = 0, size = 5) +
    labs(x = "Clusters", y = NULL) + coord_cartesian(clip = "off") +
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.12)))
  plot(p)
  
  w = length(x_labels)*1.5
  h = max_level
  if (save) {
    ggsave(p, filename = paste0("~/R/EDAP-data/plots/LTC/", file.name, ".png"), width = w, height = h, dpi = 300)
  }
  
  return(p)
}

align_clusters <- function(align, to) {
  tab <- table(to, align)
  
  c <- ncol(tab)
  n <- nrow(tab)
  M <- matrix(0, nrow=n, ncol=n)
  M[1:n, 1:c] <- tab
  
  assign <- clue::solve_LSAP(M, maximum = TRUE)
  
  newlvls <- c()
  for(idx in 1:c) {
    map <- which(assign == idx)
    newlvls <- c(newlvls, levels(to)[map])
  }
  levels(align) <- newlvls
  align <- factor(align, levels = LETTERS[1:c])
  return(align)
}
