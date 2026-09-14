library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(123)

# ============================================================
# Settings
# ============================================================

source("~/R/LTC/utils/model_utils.R")

multi_cohort_df <- read.csv("~/R/EDAP-data/MULTI_COHORT_4.csv", header = TRUE)
multi_cohort_df <- filter_out(multi_cohort_df, Cohort == "NACC")
all.vars <- c(grepv("^(RH_|LH_|CC_)", colnames(multi_cohort_df)), "BRAINSTEM")

start_vals <- sapply(all.vars, function(x) {
  dsubset <- multi_cohort_df %>% dplyr::select(RID, Time, all_of(x)) %>%
    rename(y = x, t = Time)
  get_start_estimates(dsubset)
})

n_patients <- length(unique(multi_cohort_df$RID))
n_visits_min <- min(table(multi_cohort_df$RID))
n_visits_max <- max(table(multi_cohort_df$RID))
t_min <- min(multi_cohort_df$Time)
t_max <- max(multi_cohort_df$Time)

n_variables = 6

# Set to TRUE for known clusters
# Set to FALSE for a continuous/no-cluster population
simulate_clusters <- TRUE

n_clusters <- 5

sigma_error <- 0

# SDs of patient-level random intercepts/slopes
sd_intercept <- sd(start_vals['l', ] + start_vals['v', ])
sd_slope <- sd(start_vals['g', ])

# Correlation between variables in patient-level effects
rho_intercept <- 0
rho_slope <- 0

# ============================================================
# 1. Disease trajectories
# ============================================================

# Each variable has its own baseline trajectory
# These can later be replaced by functions resembling
# your actual NLMM fixed effects.

trajectory_function <- function(time, l, g, v, cluster=NULL) {
  #mu <- l*exp(g*t) + v
  
  # --------------------------------------------------------
  # Cluster-specific deviations
  # --------------------------------------------------------
  
  if (!is.null(cluster)) {
    cluster_effects <- matrix(c(
      # cluster 1
      0.00, 0.00, 0.00,
      # cluster 2
      l*0.5, 0.00, 0.00,
      # cluster 3
      l*0.5, g*0.5, 0.00,
      # cluster 4
      l*0.5, g*0.5, v*0.5,
      # cluster 5
      l*0.25, g, 0.00
    ),
    nrow = 5,
    byrow = TRUE)
    
    mu <- (l+cluster_effects[cluster, 1]) * exp((g+cluster_effects[cluster, 2])*time) + v + cluster_effects[cluster, 3]
    
  } else {
    mu <- l*exp(g*time) + v
  }
  
  return(mu)
}

trajectory_function_ <- function(time, variable, cluster = NULL) {
  
  # Variable-specific baseline
  baseline_intercept <- c(
    1.0, 1.2, 0.9, 1.1, 0.8, 1.3
  )[variable]
  
  baseline_slope <- c(
    -0.03, -0.04, -0.025,
    -0.05, -0.035, -0.045
  )[variable]
  
  baseline_curvature <- c(
    -0.002, -0.003, -0.001,
    -0.002, -0.002, -0.003
  )[variable]
  
  mu <- baseline_intercept +
    baseline_slope * time +
    baseline_curvature * time^2
  
  # --------------------------------------------------------
  # Cluster-specific deviations
  # --------------------------------------------------------
  
  if (!is.null(cluster)) {
    
    # Different patterns for the five clusters.
    #
    # Each variable gets a somewhat different cluster
    # effect, so clustering is genuinely multivariate.
    
    cluster_effects <- matrix(
      c(
        # cluster 1
        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
        
        # cluster 2
        -0.05, -0.03, -0.02, -0.06, -0.04, -0.03,
        
        # cluster 3
        -0.10, -0.08, -0.05, -0.12, -0.07, -0.10,
        
        # cluster 4
        0.04, -0.05,  0.03, -0.08,  0.06, -0.04,
        
        # cluster 5
        -0.03,  0.06, -0.08,  0.05, -0.05,  0.08
      ),
      nrow = 5,
      byrow = TRUE
    )
    
    mu <- mu +
      cluster_effects[cluster, variable] * time
  }
  
  mu
}


# ============================================================
# 2. Patient-level multivariate random effects
# ============================================================

# Correlation matrices

R_intercept <- matrix(
  rho_intercept,
  n_variables,
  n_variables
)

diag(R_intercept) <- 1

R_slope <- matrix(
  rho_slope,
  n_variables,
  n_variables
)

diag(R_slope) <- 1


# Convert correlation matrices to covariance matrices

Sigma_intercept <-
  diag(rep(sd_intercept, n_variables)) %*%
  R_intercept %*%
  diag(rep(sd_intercept, n_variables))

Sigma_slope <-
  diag(rep(sd_slope, n_variables)) %*%
  R_slope %*%
  diag(rep(sd_slope, n_variables))


# Function for multivariate normal simulation

library(MASS)

random_intercepts <- MASS::mvrnorm(
  n = n_patients,
  mu = rep(0, n_variables),
  Sigma = Sigma_intercept
)

random_slopes <- MASS::mvrnorm(
  n = n_patients,
  mu = rep(0, n_variables),
  Sigma = Sigma_slope
)


# ============================================================
# 3. Assign true clusters
# ============================================================

if (simulate_clusters) {
  
  true_cluster <- sample(
    1:n_clusters,
    n_patients,
    replace = TRUE
  )
  
} else {
  
  # No true clusters
  true_cluster <- rep(NA_integer_, n_patients)
  
}


# ============================================================
# 4. Generate longitudinal data
# ============================================================

dat <- lapply(1:n_patients, function(i) {
  
  # Patient-specific number of observations
  n_visits <- sample(
    n_visits_min:n_visits_max,
    size = 1
  )
  
  # Disease times
  disease_time <- sort(
    runif(n_visits, 0, 10)
  )
  
  patient_data <- lapply(
    1:n_variables,
    function(var) {
      
      cluster_i <-
        if (simulate_clusters)
          true_cluster[i]
      else
        NULL
      
      mu <- trajectory_function(
        time = disease_time,
        l = start_vals['l', var],
        g = start_vals['g', var],
        v = start_vals['v', var],
        cluster = cluster_i
      )
      
      y <- mu +
        start_vals['l', var] * 
        exp(start_vals['g', var] * random_slopes[i, var] * disease_time) +
        random_intercepts[i, var] +
        rnorm(
          n_visits,
          mean = 0,
          sd = sigma_error
        )
      
      data.frame(
        patient = i,
        disease_time = disease_time,
        variable = paste0("V", var),
        outcome = y
      )
    }
  )
  
  bind_rows(patient_data)
  
}) %>%
  bind_rows()


# Add cluster information
dat <- dat %>%
  left_join(
    data.frame(
      patient = 1:n_patients,
      true_cluster = true_cluster
    ),
    by = "patient"
  )

# Wide format
dat_wide <- pivot_wider(dat, 
                        id_cols=c(patient, disease_time, true_cluster), 
                        names_from = variable, 
                        values_from = outcome)


# ============================================================
# 5. Look at the data
# ============================================================

head(dat_wide)

ggplot(
  dat_wide,
  aes(
    x = disease_time,
    y = V5,
    group = patient,
    colour = factor(true_cluster)
  )
) +
  geom_line(alpha = 0.2) +
  geom_smooth(
    aes(group = true_cluster),
    method = "loess",
    se = FALSE,
    linewidth = 1.2
  ) +
  theme_classic() +
  labs(
    colour = "True cluster",
    x = "Disease time",
    y = "Outcome"
  )

