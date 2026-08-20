library(dplyr)


get_brain_regions <- function(region.fname="~/R/EDAP-data/ADNI/reggs.xlsx", 
                              braak.fname="~/R/EDAP-data/ADNI/cho_stages.csv") {
  # load xlsx file
  library(readxl)
  brain_regions <- readxl::read_xlsx(region.fname)
  brain_regions$Lobe <- as.factor(tolower(brain_regions$Lobe))
  brain_regions$Region <- as.factor(tolower(brain_regions$Region))
  
  braak_stages <- read.csv(braak.fname)
  
  braak_stages$Stage <- as.factor(braak_stages$name)
  braak_stages$regions <- sub(".*-", "", braak_stages$id)
  
  brain_regions$BraakStage <- NA
  for (r in unique(braak_stages$regions)) {
    braak <- braak_stages[braak_stages$regions == r,]
    idxs <- grep(r, brain_regions$Variables, ignore.case=TRUE, fixed=FALSE)
    brain_regions$BraakStage[idxs] <- braak$name[1]
  }
  table(brain_regions$BraakStage)
  
  idxs <- grep("hippocampus", brain_regions$Variables, ignore.case=TRUE, fixed=FALSE)
  brain_regions$BraakStage[idxs] <- "III-IV"
  
  return(brain_regions)
}

survival_analysis <- function(df, title_name, mu, sigma, z_threshold = -1.6449) {
  library(icenReg)
  library(survival)
  
  #mu <- df %>% filter(Time < quantile(Time, 0.025)) %>% summarise(mean = mean(v),
  #                                                                sd = sd(v)) 
  df <- df %>% mutate(z = (v-mu)/sigma) %>% mutate(Event = ifelse(z < z_threshold, 1, 0))
  df$surv_time <- df$Time - min(df$Time)
  
  table(df[c("Cluster", "Event")])
  df %>% filter(Cluster == "F") %>% group_by(RID) %>% mutate(last_follow_up = max(Time)) %>%
    ungroup() %>% distinct(RID, last_follow_up) %>% select(last_follow_up) %>% hist()
  
  surv_data <- data.frame()
  for (rid in unique(df$RID)) {
    sub.df <- filter(df, RID == rid) %>% arrange(surv_time)
    cross_index = which(sub.df$z < z_threshold)[1]
    if (!is.na(cross_index)) {
      if (cross_index == 1) {
        time1 = 0
        time2 = min(sub.df$surv_time)
      } else {
        time1 = sub.df$surv_time[cross_index - 1]
        time2 = sub.df$surv_time[cross_index]
      } 
    } else {
      time1 = max(sub.df$surv_time)
      time2 = Inf
    }
    
    res = data.frame(
      RID = rid,
      cross_index = cross_index,
      time1 = time1,
      time2 = time2,
      event = !is.na(cross_index)
    )
    surv_data <- rbind(surv_data, res)
  }
  
  Clusters <- df %>% select(RID, Cluster) %>% distinct() 
  surv_data <- left_join(surv_data, Clusters, by="RID")
  
  counts <- table(surv_data[c("event", "Cluster")])
  if (any(counts == 0)) {
    idx = which(counts == 0, arr.ind=T)
    zeros <- data.frame(
      event = rownames(counts)[idx[, 'event']],
      Cluster = colnames(counts)[idx[, 'Cluster']]
    )
  for (i in seq_len(nrow(zeros))) {
    cat('Zeros encoutered in cluster ', zeros$Cluster[i], ' for event indicator: ', zeros$event[i])
  }
    print('Test with lower z-score threshold?')
    return(list('plot' = NA, "summary" = NA))
  }
  
  # 1. fit semi-parametric model
  fit <- ic_sp(Surv(time1, time2, type = "interval2") ~ Cluster,
               data = surv_data,
               bs_samples = 100)
  print(summary(fit))
  
  summ <- summary(fit)
  res.df <- data.frame(summ$summaryParameters) 
  res.df <- cbind(res.df, exp(confint(fit))) %>% select(`Exp.Est.`, `2.5 %`, `97.5 %`, p)
  colnames(res.df) <- paste(c("HR", "lwr", "upr", "p"), title_name, sep=".")
  
  
  # 2. test significance of clusters with a LRT 
  fit_null <- ic_sp(Surv(time1, time2, type = "interval2") ~ 1,
                    data = surv_data)
  
  # Extract log-likelihoods
  ll_null <- fit_null$llk
  ll <- fit$llk
  
  # Compute LRT statistic
  LRT <- 2 * (ll - ll_null)
  coefs <- length(coef(fit)) - length(coef(fit_null))
  p_value <- pchisq(LRT, df = coefs, lower.tail = FALSE)
  
  cat("Likelihood Ratio Test:\n",
      "  LRT =", round(LRT, 3), "\n",
      "  df =", coefs, "\n",
      "  p-value =", signif(p_value, 4), "\n")
  
  # 3. Plotting (use non-parametric model, just as a Kaplan Meier plot)
  fit_np <- ic_np(Surv(time1, time2, type = "interval2") ~ Cluster,
                  data = surv_data)
  
  surv_df <- lapply(names(fit_np$scurves), function(c) {
    data.frame(
      fit_np$scurves[[c]]$Tbull_ints,
      S = fit_np$scurves[[c]]$S_curves$baseline,
      Cluster = c
    )
  }) %>% dplyr::bind_rows()
  
  surv_df$surv_time <- surv_df %>% select(lower, upper) %>% apply(1, mean)
  
  tmax <- max(surv_data$time2[is.finite(surv_data$time2)])
  surv_df <- rbind(surv_df, data.frame(
    Cluster = levels(surv_data$Cluster),
    S = 1.0,
    surv_time = 0.0,
    upper = 0.0,
    lower = 0.0
  ), 
  data.frame(
    Cluster = levels(surv_data$Cluster),
    S = 0.0,
    surv_time = tmax,
    upper = tmax,
    lower = tmax
  )
  )
  surv_df$Time <- surv_df$surv_time + min(hippo$Time)
  surv_df$lower <- surv_df$lower + min(hippo$Time)
  surv_df$upper <- surv_df$upper + min(hippo$Time)
  
  gp <- ggplot(surv_df, aes(x=Time, y=S, color=Cluster)) +
    geom_step() +
    geom_step(aes(x=lower), linetype = 3) +
    geom_step(aes(x=upper), linetype=3) +
    scale_color_paletteer_d("ggthemes::Tableau_10") +
    scale_fill_paletteer_d("ggthemes::Tableau_10") +
    labs(title = title_name, y="") +
    theme_classic()
  
  return(list('plot' = gp, "summary" = res.df))
}


get_brain_regions <- function() {
  # load xlsx file
  library(readxl)
  brain_regions <- readxl::read_xlsx("~/R/EDAP-data/ADNI/reggs.xlsx")
  brain_regions$Lobe <- as.factor(tolower(brain_regions$Lobe))
  brain_regions$Region <- as.factor(tolower(brain_regions$Region))
  
  braak_stages <- read.csv("~/R/EDAP-data/ADNI/cho_stages.csv")
  
  braak_stages$Stage <- as.factor(braak_stages$name)
  braak_stages$regions <- sub(".*-", "", braak_stages$id)
  
  brain_regions$BraakStage <- NA
  for (r in unique(braak_stages$regions)) {
    braak <- braak_stages[braak_stages$regions == r,]
    idxs <- grep(r, brain_regions$Variables, ignore.case=TRUE, fixed=FALSE)
    brain_regions$BraakStage[idxs] <- braak$name[1]
  }
  table(brain_regions$BraakStage)
  
  idxs <- grep("hippocampus", brain_regions$Variables, ignore.case=TRUE, fixed=FALSE)
  brain_regions$BraakStage[idxs] <- "III-IV"
  
  return(brain_regions)
}

cat_chi_test <- function(df, varname) {
  cs <- as.character(sort(unique(df$Cluster)))
  cluster_pairs <- combn(cs, 2, simply = FALSE)
  df <- rename(df, x=all_of(varname))
  
  chi_res <- apply(cluster_pairs, 2, function(pair) {
    print(pair)
    subset <- df %>% filter(Cluster %in% pair)
    subset$Cluster <- factor(subset$Cluster)
    subset$x <- factor(subset$x)
    cont_tab <- table(subset$Cluster, subset$x)
    test_res <- chisq.test(cont_tab)
    data.frame(Cluster1 = pair[1], Cluster2 = pair[2], X.squared = test_res$statistic, p.value = test_res$p.value)
  })
  
  chi_df <- do.call(rbind, chi_res)
  chi_df$adj_p <- p.adjust(chi_df$p.value, method = "bonferroni")
  rownames(chi_df) <- NULL
  chi_df <- chi_df %>% arrange(Cluster1, Cluster2)
}


get_pairwise_p_values <- function(df, response, time_var, correct_for_age=TRUE, by_time = FALSE) {
  df <- df %>% rename(y = all_of(response),
                      t = all_of(time_var))
  
  if (correct_for_age) {
    fit <- lmer(y ~ t * Cluster + AGE + (1|RID), data = df)
  } else {
    fit <- lmer(y ~ t * Cluster + (1|RID), data = df)
  }
  
  if (by_time) {
    emm_time <- emmeans(fit, ~ Cluster|t)
    pair_res <- pairs(emm_time, adjust="tukey")
  } else {
    emm <- emmeans(fit, ~ Cluster)
    pair_res <- contrast(emm, method="pairwise", adjust="tukey")
  }
  
  return(pair_res)
}

anova_tukey_test <- function(df, colname) {
  df <- rename(df, y=all_of(colname))
  summary_df <- df |> group_by(Cluster) |> summarise(mean = mean(y), sd = sd(y), n = n())
  
  anova_res <- aov(y ~ Cluster, data=df)
  tukey_res <- TukeyHSD(anova_res)
  tukey_df <- data.frame(tukey_res$Cluster)
  tukey_df <- tukey_df %>% rownames_to_column("Clusters") %>% 
    separate(Clusters, into = c("Cluster2", "Cluster1"), sep = "-") %>% arrange(Cluster1, Cluster2)
  
  return(list("summary"=summary_df, "tukey"=tukey_df))
}

