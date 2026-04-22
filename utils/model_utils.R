library(dplyr)
library(ggplot2)
library(progmod)
library(nlme)
library(purrr)
library(R.utils)

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

exp_nlmms_fn <- function(varname, dsubset) {
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
    
    test_model <- tryCatch( { 
      withTimeout(
        nlme(
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
          verbose = verbose,
          method = "REML"
        ), timeout=600
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

fit_cluster_nlmms_foreach_adni <- function(ucsf_data, verbose=FALSE, n_samples=25, restart=FALSE, parallell=TRUE) {
  all.idxs <- grepl("^ST\\d+(CV|SV|SA|TA)$", colnames(ucsf_data)) 
  all.vars <- colnames(ucsf_data)[all.idxs == 1]
  
  if (restart) {
    exist_vars <- grep('^ST\\d+[A-Z]{2}_c.rds', list.files("~/R/LTC/tmp/out/"), value=TRUE)
    exist_vars <- gsub('_c.rds', '', exist_vars)
    run.vars <- setdiff(all.vars, exist_vars)
  } else {
    run.vars <- all.vars
  }
  
  result <- fit_cluster_nlmms_foreach(data=ucsf_data, 
                                      run.vars=run.vars, 
                                      verbose=verbsoe, 
                                      n_samples=n_samples, 
                                      parallell=parallell,
                                      all.vars=all.vars)
  
  return(result)
}

get_remaining_vars <- function(vars) {
  exist_vars <- gsub('_c.rds', '', list.files("~/R/LTC/tmp/out/"))
  return(setdiff(vars, exist_vars))
}

fit_cluster_nlmms_foreach <- function(data, run.vars, verbose=FALSE, n_samples=25, parallell=TRUE, all.vars=c()) {
  if (length(all.vars) == 0) {
    all.vars=run.vars
  }
  
  rids <- unique(data$RID)
  n_clusters <- length(unique(data$Cluster))
  
  beta_df <- data.frame(RID = rids)
  logLikes <- c()
  BICs <- c()
  AICs <- c()
  
  write_fst(data, "~/R/LTC/tmp/cluster_data.fst", compress=0)
  if (parallell) {
    foreach(i = seq_along(run.vars)) %dopar% {
      varname = run.vars[i]
      dsubset <- read_fst("~/R/LTC/tmp/cluster_data.fst", 
                          columns = c("RID", "Time", "DX.bl", "Cluster", varname))
      dsubset <- dsubset %>% rename(y = varname, t = Time, c = Cluster) %>% mutate(c = as.factor(c)) %>% 
        na.omit() -> dsubset
      result <- fit_one_cluster_nlmm_sample_fn(varname, dsubset, n_samples = n_samples, verbose=verbose)
      saveRDS(result, file = paste0("~/R/LTC/tmp/out/", varname, "_c.rds"))
      NULL
    }
  } else {
    foreach(i = seq_along(run.vars)) %do% {
      varname = run.vars[i]
      dsubset <- select(data, RID, Time, DX.bl, Cluster, all_of(varname))
      dsubset <- dsubset %>% rename(y = varname, t = Time, c = Cluster) %>% mutate(c = as.factor(c)) %>% 
        na.omit() -> dsubset
      result <- fit_one_cluster_nlmm_sample_fn(varname, dsubset, n_samples = n_samples, verbose=verbose)
      saveRDS(result, file = paste0("~/R/LTC/tmp/out/", varname, "_c.rds"))
      NULL
    }
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

fit_cluster_nlmms_slow <- function(ucsf_data, silent=FALSE, n_samples=25) {
  all.idxs <- grepl("^ST\\d+(CV|SV|SA|TA)$", colnames(ucsf_data)) 
  all.vars <- colnames(ucsf_data)[all.idxs == 1]
  
  rids <- unique(ucsf_data$RID)
  
  n_clusters <- length(unique(ucsf_data$Cluster))
  
  beta_df <- data.frame(RID = rids)
  logLikes <- c()
  BICs <- c()
  AICs <- c()
  
  results <- list()
  for (i in seq_along(all.vars))  {
    varname = run.vars[i]
    dsubset <- select(ucsf_data, RID, Time, DX.bl, Cluster, all_of(varname))
    dsubset <- dsubset %>% rename(y = varname, t = Time, c = Cluster) %>% mutate(c = as.factor(c)) %>% 
      na.omit() -> dsubset
    results[[i]] <- fit_one_cluster_nlmm_sample_fn(varname, dsubset, n_samples = n_samples, verbose=TRUE)
  }
  
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



