#library(ADNIMERGE)
library(ADNIMERGE2)
library(progmod)
library(tidyr)
library(dplyr)
library(tibble)
library(ClusterR)
#library(kml3d)
library(progress)
library(purrr)
#library(mclust)
library(caret)
library(stringr)

library(ggplot2)
library(ggpubr)
library(lubridate)

library(parallel)
library(foreach)
library(doParallel)
library(fst)

ucsf_longitudinal_all <- function(only_vol=FALSE, filter_n=1) {
  # Update: subjects with PTID: 381_S_#### should be removed, ADNIMERGE2 is not updated ...
  remove_rids <- filter(ADNIMERGE2::REGISTRY, grepl("^381_S_.*", ADNIMERGE2::REGISTRY$PTID)) %>% select(RID, VISCODE2)
  
  # Load all Longitudinal UCSF datasets
  ucsf_data1 <- ADNIMERGE2::UCSFFSL
  ucsf_data2 <- ADNIMERGE2::UCSFFSL51
  
  cols <- intersect(colnames(ucsf_data1), colnames(ucsf_data2))
  ucsf_data1 <- ucsf_data1 %>% select(all_of(cols))
  ucsf_data2 <- ucsf_data2 %>% select(all_of(cols))
  
  ucsf_data <- rbind(ucsf_data1, ucsf_data2)
  
  ucsf_data3 <- ADNIMERGE2::UCSFFSL51ALL
  cols <- intersect(colnames(ucsf_data), colnames(ucsf_data3))
  ucsf_data <- ucsf_data3 %>% select(all_of(cols)) %>% rbind(ucsf_data)
  
  ucsf_data4 <- ADNIMERGE2::UCSFFSL51Y1
  cols <- intersect(colnames(ucsf_data), colnames(ucsf_data4))
  ucsf_data <- ucsf_data4 %>% select(all_of(cols)) %>% rbind(ucsf_data)
  
  ucsf_data$RID <- as.numeric(ucsf_data$RID)
  
  ucsf_data <- anti_join(ucsf_data, remove_rids, by=c("RID", "VISCODE2"))
  
  ucsf_data <- ADNIMERGE2::REGISTRY %>% filter(VISCODE2 == "bl") %>% select(RID, EXAMDATE) %>%
    rename(EXAMDATE.bl = EXAMDATE) %>% right_join(ucsf_data, by="RID") %>% drop_na(EXAMDATE.bl)
  
  ucsf_data <- mutate(ucsf_data, Months = interval(EXAMDATE.bl, EXAMDATE) %/% months(1))
  
  dx.df <- ADNIMERGE2::DXSUM %>% filter(VISCODE2 %in% c("bl", "sc")) %>% select(RID, VISCODE2, DIAGNOSIS) %>% 
    mutate(priority = case_when(
      VISCODE2 == "bl" ~ 1,
      VISCODE2 == "sc" ~ 2,
      TRUE ~ 3
    )) %>%
    arrange(RID, priority) %>%
    group_by(RID) %>%
    slice(1) %>%
    ungroup() %>%
    select(RID, DIAGNOSIS, VISCODE2) %>% rename(DX.bl = DIAGNOSIS)
  
  ucsf_data <- dx.df %>% select(-VISCODE2) %>% right_join(ucsf_data, by="RID") %>% drop_na(DX.bl)
  
  
  # Drop the following variables
  #ST128SV - WMHypoIntensities
  #ST125SV - RightVessel
  #ST68SV - NonWMHypoIntensities
  #ST66SV - LeftVessel
  #ST7SV - CSF
  #ST80SV - RightChoroidPlexus
  #ST21SV - LeftChoroidPlexus
  #Ventricles
  #ST8SV
  #ST96SV
  #ST37SV
  #ST9SV
  #ST127SV
  #ST30SV
  
  ucsf_data <- ucsf_data |> select(-ST8SV, -ST128SV, -ST125SV, -ST68SV, 
                                     -ST66SV, -ST7SV, -ST80SV, -ST21SV, 
                                     -ST96SV, -ST37SV, -ST9SV, -ST127SV, 
                                     -ST30SV, -ST89SV)
  
  # Filter Bad Quality MRI
  # Pass: All regions OK
  # Patial: Some region have failed, remove failed regions
  ucsf_partial <- filter(ucsf_data, OVERALLQC == "Partial")
  
  # Remove failed regions from partial
  datadic <- datadict_as_tibble(ADNIMERGE2::DATADIC) %>% filter(TBLNAME == "UCSFFSL") 
  qctable <- datadic %>% filter(str_detect(FLDNAME, "QC")) %>% select(FLDNAME, TEXT) %>%
    mutate(ST_codes = str_extract_all(TEXT, "ST\\d+")) %>%
    tidyr::unnest(ST_codes)
  qcList <- c("BGQC", "TEMPQC", "CWMQC", "FRONTQC", "INSULAQC", "OCCQC", "PARQC", "VENTQC")
  for (qc_var in qcList) {
    codes <- qctable %>% filter(FLDNAME == qc_var) %>% select(ST_codes) %>% unlist()
    ucsf_partial <- ucsf_partial %>% mutate(across(starts_with(codes),
                                                   ~ if_else(.data[[qc_var]] == "Fail", NA, .)))
  }
  
  ucsf_data <- filter(ucsf_data, OVERALLQC == "Pass") %>% rbind(ucsf_partial) %>%
    arrange(RID, Months)
  
  ucsf_data |> filter(!is.na(Months)) -> ucsf_data
  
  icv <- ucsf_data %>% select(RID, VISCODE2, ST10CV)
  ucsf_data <- select(ucsf_data, -ST10CV)
  
  vol.idxs <- grepl("^ST\\d+(CV|SV)$", colnames(ucsf_data))
  vol.vars <- colnames(ucsf_data)[vol.idxs == 1]
  all.idxs <- grepl("^ST\\d+(CV|SV|SA|TA)$", colnames(ucsf_data))
  all.vars <- colnames(ucsf_data)[all.idxs == 1]
  
  # Normalize by ICV
  for (varname in vol.vars) {
    ucsf_data[[varname]] <- ucsf_data[[varname]] / icv$ST10CV
  }
  
  if (only_vol) {
    ucsf_data <- ucsf_data %>% select(RID, DX.bl, Months, all_of(vol.vars)) %>%
      filter(!if_all(all_of(vol.vars), is.na))
  } else {
    ucsf_data <- ucsf_data %>% select(RID, DX.bl, Months, all_of(all.vars)) %>%
      filter(!if_all(all_of(all.vars), is.na))
  }
  
  filtered_ids <- ucsf_data |>
    group_by(RID) |>
    tally() |>
    filter(n >= filter_n) |>
    pull(RID)
  
  ucsf_data |> filter(RID %in% filtered_ids) -> ucsf_data
  
  ucsf_data$Years <- ucsf_data$Months/12
  return(ucsf_data)
}

ucsf_xs_all <- function(only_vol=FALSE, filter_n=0) {
  ucsf_data <- ADNIMERGE2::UCSFFSX7
  
  ucsf_data <- ucsf_data |> select(-ST8SV, -ST128SV, -ST125SV, -ST68SV, 
                                   -ST66SV, -ST7SV, -ST80SV, -ST21SV, 
                                   -ST96SV, -ST37SV, -ST9SV, -ST127SV, 
                                   -ST30SV, -ST89SV)
  
  ucsf_data <- ucsf_data %>%
    mutate(
      M = VISCODE2 %>%
        gsub("scmri|bl|sc|m", "0", .) %>%
        as.numeric()
    ) 
  
  # Calculate time since baseline
  # Remove participants without baseline visit (only screening not sufficient)
  
  #missing_m <- ucsf_data %>% filter(is.na(Months)) %>% select(RID, EXAMDATE)
  ucsf_data <- ADNIMERGE2::REGISTRY %>% filter(VISCODE2 == "bl") %>% select(RID, EXAMDATE) %>%
    rename(EXAMDATE.bl = EXAMDATE) %>% right_join(ucsf_data, by="RID") %>% drop_na(EXAMDATE.bl)
  
  ucsf_data <- mutate(ucsf_data, Months = interval(EXAMDATE.bl, EXAMDATE) %/% months(1))
  
  dx.df <- ADNIMERGE2::DXSUM %>% filter(VISCODE2 %in% c("bl", "sc")) %>% select(RID, VISCODE2, DIAGNOSIS) %>% 
    mutate(priority = case_when(
      VISCODE2 == "bl" ~ 1,
      VISCODE2 == "sc" ~ 2,
      TRUE ~ 3
    )) %>%
    arrange(RID, priority) %>%
    group_by(RID) %>%
    slice(1) %>%
    ungroup() %>%
    select(RID, DIAGNOSIS, VISCODE2) %>% rename(DX.bl = DIAGNOSIS)
  
  ucsf_data <- dx.df %>% select(-VISCODE2) %>% right_join(ucsf_data, by="RID") %>% drop_na(DX.bl)
  
  ucsf_data |> filter(OVERALLQC != "Fail" | is.na(OVERALLQC)) -> ucsf_data
  
  icv <- ucsf_data %>% select(RID, VISCODE2, ST10CV)
  ucsf_data <- select(ucsf_data, -ST10CV)

  vol.idxs <- grepl("^ST\\d+(CV|SV)$", colnames(ucsf_data))
  vol.vars <- colnames(ucsf_data)[vol.idxs == 1]
  all.idxs <- grepl("^ST\\d+(CV|SV|SA|TA)$", colnames(ucsf_data))
  all.vars <- colnames(ucsf_data)[all.idxs == 1]
  
  # Normalize by ICV
  for (varname in vol.vars) {
    ucsf_data[[varname]] <- ucsf_data[[varname]] / icv$ST10CV
  }
  
  if (only_vol) {
    ucsf_data <- ucsf_data %>% select(RID, M, DX.bl, Months, all_of(vol.vars)) %>%
      filter(!if_all(all_of(vol.vars), is.na))
  } else {
    ucsf_data <- ucsf_data %>% select(RID, M, DX.bl, Months, all_of(all.vars)) %>%
      filter(!if_all(all_of(all.vars), is.na))
  }
  
  filtered_ids <- ucsf_data |>
    group_by(RID) |>
    tally() |>
    filter(n >= filter_n) |>
    pull(RID)
  
  ucsf_data |> filter(RID %in% filtered_ids) -> ucsf_data
  
  ucsf_data$Years <- ucsf_data$Months/12
  return(ucsf_data)
  
}

betaKM <- function(beta_df, k=2, minClusterSize=10) {
  beta.idxs <- grepl("^ST\\d+(CV|SV|HS|SA|TA)$", colnames(beta_df))
  beta.vars <- colnames(beta_df)[beta.idxs == 1]
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

get_gp <- function(varname, df, model_fit) {
  
  tvec <- seq(min(df$t), max(df$t), length.out = 100)
  rids <- unique(df$RID)
  new.df <- data.frame(RID = rep(rids, each = 100), t = rep(tvec, times = length(rids)))
  
  new.df$y <- predict(model_fit, new.df)
  dx.df <- df %>% select(RID, DX.bl) %>% distinct()
  new.df <- left_join(new.df, dx.df, by="RID")
  
  traj_ts <- ggplot(df, aes(x = t, y = y, color = DX.bl)) +
    geom_line(aes(group = RID),  alpha=1) +
    labs(title = "", x = "Years (time shifted)", y = varname, color = NULL) +
    scale_color_brewer(palette = "YlOrRd") +
    geom_line(data=new.df, aes(group = RID), alpha=0.1) +
    theme_minimal() +
    ylim(0, NA)
  
  
  return(traj_ts)
}

get_start_estimates <- function(dsubset) {
  # Time should be in years
  dsubset$t_abs <- abs(dsubset$t)
  
  # 1) asymptote guess = early median
  mu_lower <- with(dsubset, tapply(y, cut(t, quantile(t, probs = c(0,0.05), na.rm=TRUE), include.lowest=TRUE), median, na.rm=TRUE))
  mu_lower <- as.numeric(mu_lower)
  
  # Minimal difference between asymptote and intercept to avoid scale 0
  epsilon <- mu_lower*0.001
  
  # 2) intercept = median at t=0
  y0 <- with(dsubset, tapply(y, cut(t_abs, quantile(t_abs, probs = c(0,0.05), na.rm=TRUE), include.lowest=TRUE), median, na.rm=TRUE))
  y0 <- min(as.numeric(y0), mu_lower-epsilon)
  
  # 3) late time median
  mu_upper <- with(dsubset, tapply(y, cut(t, quantile(t, probs = c(0.95,1), na.rm=TRUE), include.lowest=TRUE), median, na.rm=TRUE))
  mu_upper <- min(as.numeric(mu_upper), y0-epsilon)
  
  # Late time
  t_late = unname(quantile(dsubset$t, probs=0.99))
  
  # 4) Assign values
  v = mu_lower
  l = y0-mu_lower
  
  g = 1/t_late*log((mu_upper-mu_lower)/(y0-mu_lower))
  #a = l * exp(g)
  
  start_vals <- c(l=l, g=g, v=v)
  
  return(start_vals)
}

exp_model_expr <- deriv(
  expression(-exp(l) * exp(exp(g)*(1+gi)*t) + v),
  namevec = c('l', 'g', 'gi', 'v'),
  function.arg = c('t', 'l', 'g', 'gi', 'v'))

exp_model_expr_ts <- deriv(
  expression(-exp(l) * exp(exp(g)*(1+gi)*(t+s)) + v),
  namevec = c('l', 'g', 'gi', 's', 'v'),
  function.arg = c('t', 'l', 'g', 'gi', 's', 'v'))

exp_model_expr_u <- deriv(
  expression(-exp(l) * exp(exp(g)*t) + v),
  namevec = c('l', 'g', 'v'),
  function.arg = c('t', 'l', 'g', 'v'))

exp_nlmms_fn <- function(varname, ucsf_data) {
  ucsf_data |> select(RID, Time, DX.bl, all_of(varname)) |> 
    rename(y = varname, t = Time) %>% drop_na(y) -> dsubset
  
  mu_t <- mean(dsubset$t)
  dsubset$t <- (dsubset$t - mu_t)/10
  dsubset$y <- dsubset$y/quantile(dsubset$y, 0.99)
  
  start_vals <- get_start_estimates(dsubset)
  
  ctrl <- nlme::nlmeControl(maxIter=100, 
                            niterEM = 100, 
                            msVerbose = FALSE, 
                            pnlsMaxIter = 20, 
                            msMaxIter = 200, 
                            returnObject = FALSE,
                            tolerance = 1e-6,     # overall convergence (default ~1e-6)
                            pnlsTol = 1e-4,       # PNLS step tolerance (loosen if PNLS fails)
                            msTol = 1e-6,
                            minScale = 0.001,
                            opt = "nlminb")
  
  exp_model <- tryCatch( { nlme(
    y ~ exp_model_expr(t, l, g, gi, v),
    #y ~ exp_model_expr_u(t, l, g, v),
    data = dsubset,
    fixed = list(l ~ 1,
                 g ~ 1,
                 v ~ 1),                 # fixed effects for scale l and g and vertical additive effect (v)
    random = list(gi ~ 1,
                  v ~ 1),  # random vertical additive effect (v_i) and decline (g_i)
    groups = ~ RID,
    start = c(l = log(-start_vals['l']), v=start_vals['v'], g = log(start_vals['g'])),
    #start = fixed.effects(fix_model),
    control = ctrl,
    method = "REML"
  )
  }, error = function(e) {
    message("Fitting of model to ", varname, " failed: ", e$message)
    return(NULL)
  })
  
  if (!is.null(exp_model)) {
    #new_gp <- get_gp(varname, dsubset, exp_model)
    #plot(new_gp)
    
    rand <- ranef(exp_model)
    
    rand <- rownames_to_column(rand, var = "RID")
    
    decline <- rand %>% select(RID, gi) # or gi!!
    colnames(decline) <- c("RID", varname)
    
    decline$RID <- as.numeric(decline$RID)
    
    results <- list(
      betas = decline,
      bic = BIC(exp_model),
      aic = AIC(exp_model),
      logLike = exp_model$logLik
    )
    
    return(results)
  } else {
    return(NULL)
  }
  
}

exp_nlmms_sample_fn <- function(varname, dsubset, n_samples=100, verbose=FALSE) {
  mu_t <- mean(dsubset$t)
  dsubset$t <- (dsubset$t - mu_t)/10
  dsubset$y <- dsubset$y/quantile(dsubset$y, 0.99)
  rids <- unique(dsubset$RID)
  
  ctrl <- nlme::nlmeControl(maxIter=100, 
                            niterEM = 100, 
                            msVerbose = FALSE, 
                            pnlsMaxIter = 20, 
                            msMaxIter = 200, 
                            returnObject = FALSE,
                            tolerance = 1e-2,     # overall convergence (default ~1e-6)
                            pnlsTol = 1e-1,       # PNLS step tolerance (loosen if PNLS fails)
                            msTol = 1e-2,
                            minScale = 0.0001,
                            opt = "nlminb")
  
  # Bootstrap to get different starting values
  #start_vals <- matrix(NA, nrow = n_samples, ncol = 3)
  #colnames(start_vals) = c("l", "g", "v")
  best_ll = 0
  llrs = c()
  exp_model=NULL
  for (b in 1:n_samples) {
    # Resample rows with replacement
    start_vals <- dsubset %>% filter(RID %in% sample(RID, size = length(RID), replace = TRUE)) %>%
      do(as.data.frame(t(get_start_estimates(.)))) %>% unlist()
    
    test_model <- tryCatch( { nlme(
      #y ~ exp_model_expr_ts(t, l, g, gi, s, v),
      y ~ exp_model_expr(t, l, g, gi, v),
      #y ~ exp_model_expr_u(t, l, g, v),
      data = dsubset,
      fixed = list(l ~ 1,
                   g ~ 1,
                   v ~ 1),                 # fixed effects for scale l and g and vertical additive effect (v)
      random = list(gi ~ 1,
                    #s ~ 1,
                    v ~ 1),  # random vertical additive effect (v_i) and decline (g_i)
      groups = ~ RID,
      start = c(l = log(-start_vals['l']), v=start_vals['v'], g = log(start_vals['g'])),
      #start = fixed.effects(fix_model),
      control = ctrl,
      method = "REML"
    )
    }, error = function(e) {
      message("Fitting of model to ", varname, " failed: ", e$message)
      return(NULL)
    })
    if (!is.null(test_model)) {
      ll_test = logLik(test_model)
      # Likelihood ratios - if are less than 2 for 3 in a row, model is converged and we break the loop
      llrs = tail(c(llrs, abs(2*(ll_test-best_ll))), 3)
      if (length(llrs) > 2 & max(llrs) < 2) {
        if (verbose) print("Model converged")
        break
      }
      if (ll_test > best_ll) {
        exp_model <- test_model
        best_ll = ll_test
        if (verbose) {
          cat("New highest likelihood: ", best_ll, "\n")
          cat("BIC: ", BIC(exp_model), "\n")
        }
      }
      rm(test_model)
    }
  }
  
  if (!is.null(exp_model)) {
    ctrl$minScale <- 0.001
    ctrl$tolerance <- 1e-6
    ctrl$pnlsTol <- 1e-4
    ctrl$msTol <- 1e-6
    #ctrl$returnObject <- FALSE
    exp_model <- tryCatch( { 
      update(exp_model, control=ctrl)
    }, error = function(e) {
      message("Fitting of model to ", varname, " failed: ", e$message)
      return(exp_model)
    })
    
    #new_gp <- get_gp(varname, dsubset, exp_model)
    #plot(new_gp)
    
    rand <- ranef(exp_model)
    
    rand <- rownames_to_column(rand, var = "RID")
    
    decline <- rand %>% select(RID, gi) # or gi!!
    colnames(decline) <- c("RID", varname)
    
    decline$RID <- as.numeric(decline$RID)
    
    results <- list(
      betas = decline,
      bic = BIC(exp_model),
      aic = AIC(exp_model),
      logLike = exp_model$logLik
    )
    rm(exp_model)
    
    return(results)
  } else {
    return(NULL)
  }
  
}

fit_one_cluster_nlmm_fn <- function(varname, data) {
  data |> select(RID, Time, DX.bl, Cluster, all_of(varname)) |> 
    rename(y = varname, t = Time, c = Cluster) %>% mutate(c = as.factor(c)) %>% 
    na.omit() -> dsubset
  
  mu_t <- mean(dsubset$t)
  dsubset$t <- (dsubset$t - mu_t)/10
  dsubset$y <- dsubset$y/quantile(dsubset$y, 0.99)

  start_vals <- dsubset %>% group_by(c) %>% do(as.data.frame(t(get_start_estimates(.))))
  
  l0 <- setNames(start_vals$l, paste0("l.c", start_vals$c))
  g0 <- setNames(start_vals$g, paste0("g.c", start_vals$c))
  v0 <- setNames(start_vals$v, paste0("v.c", start_vals$c))
  #v0 <- c(v=mean(start_vals$v))
  
  #ctrl = nlmeControl(maxIter=100, pnlsMaxIter = 10, tolerance = 1e-5, pnlsTol = 1e-4, niterEM = 100)
  ctrl <- nlme::nlmeControl(niterEM = 100, msVerbose = FALSE, pnlsMaxIter = 20, msMaxIter = 200,
                            returnObject = FALSE, tolerance = 1e-3, pnlsTol = 1e-2)
  
  exp_model <- tryCatch( { nlme(
    y ~ exp_model_expr(t, l, g, gi, v),
    #y ~ exp_model_expr_u(t, l, g, v),
    data = dsubset,
    fixed = list(l ~ c - 1,
                 g ~ c - 1,
                 v ~ c - 1),                 # fixed effects for scale (a) and vertical additive effect (v)
    random = list(gi ~ 1,
                  v ~ 1),  # random vertical additive effect (v_i) and decline (g_i)
    groups = ~ RID,
    start = c(log(-l0), log(g0), v0),
    weights = varIdent(form = ~ 1 | c),      # allow different variances within different clusters
    control = ctrl,
    method = "REML"
  )
  }, error = function(e) {
    message("Fitting of model to ", varname, " failed: ", e$message)
  })
  
  if (!is.null(exp_model)) {
    
    rand <- ranef(exp_model)
    
    rand <- rownames_to_column(rand, var = "RID")
    
    decline <- rand %>% select(RID, gi)
    #decline <- rand %>% select(RID, `g.(Intercept)`)
    colnames(decline) <- c("RID", varname)
    
    decline$RID <- as.numeric(decline$RID)
    
    results <- list(
      betas = decline,
      bic = BIC(exp_model),
      aic = AIC(exp_model),
      logLike = exp_model$logLik
    )
    
    return(results)
  } 
}

fit_one_cluster_nlmm_sample_fn <- function(varname, dsubset, n_samples=25, verbose=FALSE) {
  mu_t <- mean(dsubset$t)
  dsubset$t <- (dsubset$t - mu_t)/10
  dsubset$y <- dsubset$y/quantile(dsubset$y, 0.99)
  
  ctrl <- nlme::nlmeControl(maxIter=100, 
                            niterEM = 100, 
                            msVerbose = FALSE, 
                            pnlsMaxIter = 20, 
                            msMaxIter = 200, 
                            returnObject = FALSE,
                            tolerance = 1e-2,    
                            pnlsTol = 1e-1,       
                            msTol = 1e-2,
                            minScale = 0.0001,
                            opt = "nlminb")
  
  # Bootstrap to get different starting values
  best_ll = 0
  llrs = c()
  exp_model=NULL
  for (b in 1:n_samples) {
    # Resample rows with replacement
    # Stratify by cluster
    start_vals <- dsubset %>% group_by(c) %>% filter(RID %in% sample(RID, size = length(RID), replace = TRUE)) %>%
                  do(as.data.frame(t(get_start_estimates(.)))) 
    
  
    l0 <- setNames(start_vals$l, paste0("l.c", start_vals$c))
    g0 <- setNames(start_vals$g, paste0("g.c", start_vals$c))
    v0 <- setNames(start_vals$v, paste0("v.c", start_vals$c))

    test_model <- tryCatch( { nlme(
      y ~ exp_model_expr(t, l, g, gi, v),
      data = dsubset,
      fixed = list(l ~ c - 1,
                   g ~ c - 1,
                   v ~ c - 1),                
      random = list(gi ~ 1,
                    v ~ 1),  
      groups = ~ RID,
      start = c(log(-l0), log(g0), v0),
      weights = varIdent(form = ~ 1 | c), # allow different variances within different clusters
      control = ctrl,
      method = "REML"
    )
    }, error = function(e) {
      message("Fitting of model to ", varname, " failed: ", e$message)
    })
    if (!is.null(test_model)) {
      ll_test = logLik(test_model)
      # Likelihood ratios - if are less than 2 for 3 in a row, model is converged and we break the loop
      llrs = tail(c(llrs, abs(2*(ll_test-best_ll))), 3)
      if (length(llrs) > 2 & max(llrs) < 2) {
        if (verbose) print("Model converged")
        break
      }
      if (ll_test > best_ll) {
        exp_model <- test_model
        best_ll = ll_test
        if (verbose) { 
          cat("New highest likelihood: ", best_ll, "\n")
          cat("BIC: ", BIC(exp_model), "\n")
        }
      }
      rm(test_model)
    }
  }
  
  if (!is.null(exp_model)) {
    ctrl$minScale <- 0.001
    ctrl$tolerance <- 1e-6
    ctrl$pnlsTol <- 1e-4
    ctrl$msTol <- 1e-6
    exp_model <- tryCatch( { 
      update(exp_model, control=ctrl)
    }, error = function(e) {
      message("Fitting of model to ", varname, " failed: ", e$message)
      return(exp_model)
    })
    
    rand <- ranef(exp_model)
    
    rand <- rownames_to_column(rand, var = "RID")
    
    decline <- rand %>% select(RID, gi)
    colnames(decline) <- c("RID", varname)
    
    decline$RID <- as.numeric(decline$RID)
    
    results <- list(
      betas = decline,
      bic = BIC(exp_model),
      aic = AIC(exp_model),
      logLike = exp_model$logLik
    )
    
    rm(exp_model)
    return(results)
  } 
}

fit_cluster_nlmms_mc <- function(ucsf_data) {
  all.idxs <- grepl("^ST\\d+(CV|SV|SA|TA)$", colnames(ucsf_data)) 
  all.vars <- colnames(ucsf_data)[all.idxs == 1]
  
  rids <- unique(ucsf_data$RID)
  
  n_clusters <- length(unique(ucsf_data$Cluster))
  
  resultsList <- mclapply(all.vars, fit_one_cluster_nlmm_fn, data=ucsf_data)
  
  return(list(betas = beta_df, 
              logLikes = logLikes,
              aic = AICs, 
              bic = BICs))
}

fit_cluster_nlmms_foreach <- function(ucsf_data, silent=FALSE, n_samples=25) {
  all.idxs <- grepl("^ST\\d+(CV|SV|SA|TA)$", colnames(ucsf_data)) 
  all.vars <- colnames(ucsf_data)[all.idxs == 1]
  
  rids <- unique(ucsf_data$RID)
  
  n_clusters <- length(unique(ucsf_data$Cluster))
  
  beta_df <- data.frame(RID = rids)
  logLikes <- c()
  BICs <- c()
  AICs <- c()
  
  write_fst(ucsf_data, "~/R/LTC/tmp/ucsf_data.fst", compress=0)
  foreach(i = seq_along(all.vars)) %dopar% {
    varname = all.vars[i]
    dsubset <- read_fst("~/R/LTC/tmp/ucsf_data.fst", 
                        columns = c("RID", "Time", "DX.bl", "Cluster", varname))
    dsubset <- dsubset %>% rename(y = varname, t = Time, c = Cluster) %>% mutate(c = as.factor(c)) %>% 
               na.omit() -> dsubset
    result <- fit_one_cluster_nlmm_sample_fn(varname, dsubset, n_samples = n_samples, verbose=FALSE)
    saveRDS(result, file = paste0("~/R/LTC/tmp/out/", varname, "_c.rds"))
    NULL
  }
  
  results <- lapply(all.vars, function(x) readRDS(paste0("~/R/LTC/tmp/out/", x, "_c.rds")))
  
  results <- results[!sapply(results, is.null)]
  
  betaList <- lapply(results, `[[`, "betas") 
  beta_df <- purrr::reduce(betaList, full_join, by = "RID")
  
  return(list(
    betas = beta_df,
    bic = sapply(results, `[[`, "bic"),
    aic = sapply(results, `[[`, "aic"),
    logLikes = sapply(results, `[[`, "logLike")
  ))
}

#### Initial steps ####
fit_mcdp = FALSE # set to FALSE to load previous MCDP run
fit_inital = FALSE # set to FALSE to load previous initial model fitting
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
  
  save(nlmmBasic, file = "~/R/EDAP-data/LTC4/nlmmBasic.Rdata")
} else {
  load("~/R/EDAP-data/LTC4/nlmmBasic.Rdata")
}

cat("Inital models fitted to", round(dim(nlmmBasic$betas)[2]/length(all.vars)*100), "% of variables") 
cat("Mean BIC: ", mean(nlmmBasic$bic))

#### KM, no cross-validation ####
for(ii in 1:1) {
  set.seed(ii)
  ucsf_clust <- ucsf_ab #%>% select(RID, time_shift, M, DX.bl, Months, Years, all_of(vars))
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
  nlmmBest <- fit_cluster_nlmms_foreach(ucsf_clust)
  
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
          fit_cluster_nlmms_foreach(n_samples=25)
        
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
        fit_cluster_nlmms_foreach(n_samples=25)
      clusterList[[append_c+1]] <- c_df
      next_c = next_c + 1
    }
    
    # Find the best clustering
    bics <- sapply(nlmmCandidates, function(x) mean(x$bic))
    aics <- sapply(nlmmCandidates, function(x) mean(x$aic))
    
    all.bics <- sapply(nlmmCandidates, function(x) {
      idx <- which(all.vars %in% names(x$betas))
      bics <- rep(NA, length(all.vars))
      bics[idx] <- x$bic
      return(bics)
    } ) 
    
    #c.idxs <- unlist(apply(all.bics, 1, which.min))
    
    #all.bics[which(is.na(all.bics))] <- max(all.bics, na.rm=T)
    
    #total_bics <- apply(all.bics, 2, sum)
    
    
    # Check if any is significantly better than parent 
    p_values <- lapply(2:4, function(col) {
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
    
    
    nlmmBest <- nlmmCandidates[[best_idx]]
    
    # Remove this pair
    #clusterPairs <- clusterPairs[!sapply(clusterPairs, function(x) all(x == pair))]
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
  
  save(adniLTC, file = "~/R/EDAP-data/LTC4/exp_km_ab.Rdata")
}

for(i in 1:length(clusterList)){
  print(i)
  print(table(clusterList[[i]]$Cluster))
}

for(i in 1:length(treeIdx)){
  print(table(clusterList[[treeIdx[i]]]$Cluster))
}
treeIdx


# Generate dendrogram
# Can move this to separate script?

segments <- data.frame(x=0, y=1, xend=0, yend=2)
segments <- rbind(segments,
                  data.frame(x=-adniLTC@tree[1,2], y=2, xend=adniLTC@tree[2,2], yend=2),
                  data.frame(x=-adniLTC@tree[1,2], y=2, xend=-adniLTC@tree[1,2], yend=3),
                  data.frame(x=adniLTC@tree[2,2], y=2, xend=adniLTC@tree[2,2], yend=3))

prev_level <- na.omit(adniLTC@tree[,2])
nextIdx = 3
activeClusters = names(prev_level)
negativeBranch = c("A")

for (level in 3:ncol(adniLTC@tree)) {
  new_level = adniLTC@tree[,level]
  for (c in activeClusters) {
    if (new_level[c] != prev_level[c]) {
      # branch has branched
      new_pair = c(adniLTC@tree[c,level],adniLTC@tree[nextIdx,level])
      if (c %in% negativeBranch) {
        prev_x = -prev_level[c]
        negativeBranch <- c(negativeBranch, names(adniLTC@tree[,level])[nextIdx])
      } else {
        prev_x = prev_level[c] 
      }
      segments <- rbind(segments,
                        data.frame(x=prev_x-new_pair[1], y=level, xend=prev_x+new_pair[2], yend=level),
                        data.frame(x=prev_x-new_pair[1], y=level, xend=prev_x-new_pair[1], yend=level+1),
                        data.frame(x=prev_x+new_pair[2], y=level, xend=prev_x+new_pair[2], yend=level+1))
    }
  }
}

max_level = max(segments$yend)
x_ticks <- segments %>% filter(yend == max_level) %>% select(x)
x_labels <- unique(c(negativeBranch, rownames(adniLTC@tree)))
x_labels <- sprintf("%s (n=%d)", x_labels, adniLTC@tree[x_labels,ncol(adniLTC@tree)])
p <- ggplot(segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  scale_y_reverse(limits=c(max_level+0.3, 1)) +
  theme_classic() +
  #scale_x_discrete("Clusters", breaks = xbr$x, labels=c("A", "D", "C", "E", "B"), limits = c("0", "1500"), position = "bottom") 
  # Disable default x-axis elements
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  ) +
  annotate("text",
           x = unlist(x_ticks), 
           y = max_level+0.2,        # slightly above the line (since y is reversed)
           label = unlist(x_labels),
           vjust = 0, size = 3) +
  labs(x = "Clusters", y = NULL)
p

ggsave(p, filename = "~/R/EDAP-data/plots/LTC4/dendro.png", width = 6.5, height = 4, dpi = 300)

