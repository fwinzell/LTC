library(dplyr)
library(ggplot2)
library(tidyr)

set.seed(123)

# ============================================================
# Settings
# ============================================================

multi_cohort_df <- read.csv("~/R/EDAP-data/MULTI_COHORT_4.csv", header = TRUE)
multi_cohort_df <- filter_out(multi_cohort_df, Cohort == "NACC")
all.vars <- c(grepv("^(RH_|LH_|CC_)", colnames(multi_cohort_df)), "BRAINSTEM")

n_patients <- length(unique(multi_cohort_df$RID))
n_visits_min <- min(table(multi_cohort_df$RID))
n_visits_max <- max(table(multi_cohort_df$RID))
t_min <- min(multi_cohort_df$Time)
t_max <- max(multi_cohort_df$Time)

# Number of true clusters
n_clusters <- 5

dsubset <- multi_cohort_df %>% select(RID, Time, all_of(all.vars[1])) %>%
  filter(abs(Time) < 2) %>% rename(y = all.vars[1])

# Measurement noise
sigma_error <- 0.20

# Patient-level random effects
sigma_intercept <- mean(dsubset$y)
sigma_slope <- 0.10

# ============================================================
# 1. Define the underlying trajectory for each true cluster
# ============================================================

# Disease time ranges from 0 to 10
time_grid <- seq(t_min, t_max, length.out = 100)

# Define cluster-specific trajectories.
# Modify these to resemble your actual trajectories.

trajectory_function <- function(time, cluster) {
  
  # Baseline trajectory
  baseline <- 1 - 0.04 * time - 0.003 * time^2
  
  # Cluster-specific deviations
  if (cluster == 1) {
    # Slow progression
    effect <- -0.01 * time
    
  } else if (cluster == 2) {
    # Moderate progression
    effect <- -0.04 * time
    
  } else if (cluster == 3) {
    # Fast progression
    effect <- -0.08 * time
    
  } else if (cluster == 4) {
    # Early acceleration
    effect <- -0.08 * (time / (1 + exp(-(time - 4))))
    
  } else if (cluster == 5) {
    # Delayed progression
    effect <- -0.08 * (time / (1 + exp(-(time - 7))))
    
  }
  
  baseline + effect
}

# ============================================================
# 2. Assign patients to true clusters
# ============================================================

patient_info <- data.frame(
  patient = 1:n_patients,
  true_cluster = sample(
    1:n_clusters,
    n_patients,
    replace = TRUE
  )
)

# ============================================================
# 3. Create patient-specific random effects
# ============================================================

patient_info <- patient_info %>%
  mutate(
    random_intercept = rnorm(
      n_patients,
      mean = 0,
      sd = sigma_intercept
    ),
    
    random_slope = rnorm(
      n_patients,
      mean = 0,
      sd = sigma_slope
    )
  )

# ============================================================
# 4. Generate longitudinal observations
# ============================================================

dat <- lapply(1:n_patients, function(i) {
  
  # Different number of visits for each patient
  n_visits <- sample(
    n_visits_min:n_visits_max,
    1
  )
  
  # Patient-specific observation times
  disease_time <- sort(
    runif(n_visits, 0, 10)
  )
  
  cluster <- patient_info$true_cluster[i]
  
  # Population-level trajectory
  mu <- trajectory_function(
    disease_time,
    cluster
  )
  
  # Patient-specific deviation
  y <- mu +
    patient_info$random_intercept[i] +
    patient_info$random_slope[i] * disease_time +
    rnorm(
      n_visits,
      mean = 0,
      sd = sigma_error
    )
  
  data.frame(
    patient = i,
    disease_time = disease_time,
    outcome = y,
    true_cluster = cluster
  )
}) %>%
  bind_rows()

# ============================================================
# 5. Plot the simulated data
# ============================================================

ggplot(
  dat,
  aes(
    x = disease_time,
    y = outcome,
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
