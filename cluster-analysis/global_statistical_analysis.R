library(ADNIMERGE2)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(dplyr)
library(Rtsne)
library(paletteer)
library(tidyr)
library(broom)
library(lme4)
library(lmerTest)

multi_cohort_df <- read.csv("~/R/EDAP-data/MULTI_COHORT_4.csv", header = TRUE)

run <- "exp_km_ab_ao"
load(paste("~/R/EDAP-data/LTC_MC/new/", run, ".Rdata", sep = ""))

adni_dl <- new.env()
source("~/R/LTC/utils/adni_data_loaders.R", local=adni_dl)

oasis_dl <- new.env()
source("~/R/LTC/utils/oasis_data_loaders.R", local=oasis_dl)

mc_dl <- new.env()
source("~/R/LTC/utils/multi_cohort_loader.R", local=mc_dl)

source("~/R/LTC/utils/analysis_utils.R")

ab_pos_rids <- c(paste0("ADNI_", adni_dl$get_ab_pos_ids()),
                 gsub("OAS", "OASIS_", oasis_dl$get_ab_pos_ids()))

mc_mri <- mc_dl$multi_cohort_mri()
mri_controls <- mc_mri %>% filter(DX.bl == "CN" & !(RID %in% ab_pos_rids))

Clusters <- data.frame(
  Cluster = multiLTC@Cluster,
  RID = multiLTC@RID
)

Clusters %>% left_join(multi_cohort_df, by = "RID") %>% drop_na(Cluster) -> multi_cohort_df

mri_cols <- c(grepv("^(RH_|LH_|CC_)", colnames(multi_cohort_df)), "BRAINSTEM")

mc_df_z <- data.frame()

for (varname in mri_cols) {
  mu <- mri_controls %>% select(RID, Months, all_of(varname)) %>% na.omit() %>%
    arrange(Months) %>% distinct(RID, .keep_all = TRUE) %>%
    summarise(Mean = mean(.data[[varname]], na.rm=TRUE),
              SD = sd(.data[[varname]], na.rm=TRUE)) 
  
  #z_col = paste0(varname, "_z")
  df <- multi_cohort_df %>% select(RID, Time, DX.bl, Cluster, all_of(varname)) %>%
    drop_na(Cluster) %>% mutate(!!varname := (.data[[varname]]-mu$Mean)/mu$SD) 
  
  if (ncol(mc_df_z) == 0) {
    mc_df_z <- df
  } else {
    mc_df_z <- left_join(mc_df_z, df, by=c('RID', 'Time', 'Cluster', 'DX.bl'))
  }
}


library(lme4)
library(lmerTest)
library(emmeans)
library(splines)

var = mri_cols[1]

# Statistical test:
# H0: no cluster effect after accounting for disease stage and global severity
# H1: regional atrophy additionally differs between the clusters

mc_df_z$global_z <- rowMeans(
  mc_df_z[, setdiff(mri_cols, var)],
  na.rm = TRUE
)

preds_0 <- "ns(Time, df=3) + global_z + (1 | RID)"
preds_1 <- "ns(Time, df=3) + global_z + Cluster + (1 | RID)"

m0 <- lmer(
  reformulate(preds_0, response = var),
  data = mc_df_z,
  REML = FALSE
)

m1 <- lmer(
  reformulate(preds_1, response = var),
  data = mc_df_z,
  REML = FALSE
)

p_val = anova(m0, m1)$`Pr(>Chisq)`[2]

results <- lapply(mri_cols, function(var) {
  
  mc_df_z$global_z <- rowMeans(
    mc_df_z[, setdiff(mri_cols, var)],
    na.rm = TRUE
  )
  
  preds_0 <- "ns(Time, df=3) + global_z + (1 | RID)"
  preds_1 <- "ns(Time, df=3) + global_z + Cluster + (1 | RID)"
  
  m0 <- lmer(
    reformulate(preds_0, response = var),
    data = mc_df_z,
    REML = FALSE
  )
  
  m1 <- lmer(
    reformulate(preds_1, response = var),
    data = mc_df_z,
    REML = FALSE
  )
  
  
  cmp <- anova(m0, m1)
  
  data.frame(
    region = var,
    chisq = cmp$Chisq[2],
    df = cmp$Df[2],
    p = cmp$`Pr(>Chisq)`[2]
  )
})

results <- do.call(rbind, results)

results$p_fdr <- p.adjust(results$p, method = "BH")


results_2 <- lapply(mri_cols, function(var) {
  
  mc_df_z$global_z <- rowMeans(
    mc_df_z[, setdiff(mri_cols, var)],
    na.rm = TRUE
  )
  
  preds_1 <- "ns(Time, df=3) + global_z + Cluster + (1 | RID)"
  preds_2 <- "ns(Time, df=3) * Cluster + global_z + (1 | RID)"
  
  m1 <- lmer(
    reformulate(preds_1, response = var),
    data = mc_df_z,
    REML = FALSE
  )
  
  m2 <- lmer(
    reformulate(preds_2, response = var),
    data = mc_df_z,
    REML = FALSE
  )
  
  cmp <- anova(m1, m2)
  
  data.frame(
    region = var,
    chisq = cmp$Chisq[2],
    df = cmp$Df[2],
    p = cmp$`Pr(>Chisq)`[2]
  )
})

results_2 <- do.call(rbind, results_2)

results_2$p_fdr <- p.adjust(results_2$p, method = "BH")

table((results$p_fdr < 0.05))
table((results_2$p_fdr < 0.05))

make_dk_brain_plot <- function(result_df) {
  library(ggseg)
  library(ggplot2)
  
  desikan_atlas <- dk()
  all_labels <- desikan_atlas$core %>%
    select(hemi, label)
  
  result_df <- left_join(all_labels, result_df, by="label")
  
  gp_left <- ggplot() +
    geom_brain(atlas = dk(),
               data = result_df, 
               mapping = aes(fill=p_plot),
               position = position_brain(nrow=2, ncol =2),
               hemi = "left",
               view = c("lateral", "medial", "superior", "inferior")) +
    scale_fill_viridis_c(
      option = "magma",
      na.value = "grey90",
      name = expression(-log[10](FDR~p))
    ) +
    theme_void() +
    theme(legend.position = "bottom")
  
  gp_right <- ggplot() +
    geom_brain(atlas = dk(),
               data = result_df, 
               mapping = aes(fill=p_plot),
               position = position_brain(nrow=2, ncol =2),
               hemi = "right",
               view = c("lateral", "medial", "superior", "inferior")) +
    scale_fill_viridis_c(
      option = "magma",
      na.value = "grey90",
      name = expression(-log[10](FDR~p))
    ) +
    theme_void() +
    theme(legend.position = "bottom")
  
  
  gp <- ggarrange(gp_left, gp_right, ncol=2, common.legend = TRUE, legend = "bottom")
  plot(gp)
  return(gp)
}


library(stringr)

results <- mutate(results, 
                  label = str_to_lower(region),
                  p_plot = ifelse(p_fdr < 0.05, -log(p_fdr), NA),
                  p_inter = cut(p_fdr, breaks = c(0, 0.001, 0.005, 0.01, 0.05))
                  ) %>% select(-region)

gp1 <- make_dk_brain_plot(results)

results_2 <- mutate(results_2, 
                  label = str_to_lower(region),
                  p_plot = ifelse(p_fdr < 0.05, -log(p_fdr), NA),
                  p_inter = cut(p_fdr, breaks = c(0, 0.001, 0.005, 0.01, 0.05))
              ) %>% select(-region)

gp2 <- make_dk_brain_plot(results_2)

make_aseg_brain_plot <- function(result_df) {
  aseg_atlas <- aseg()
  all_labels <- aseg_atlas$core %>%
    select(hemi, label) %>%
    rename(label_ = label) %>%
    mutate(label = gsub("left-", "lh_", str_to_lower(label_))) %>%
    mutate(label = gsub("right-", "rh_", label)) %>%
    mutate(label = sub("_", "\u2021PH\u2021", label),
           label = gsub("[-_]", "", label),
           label = sub("\u2021PH\u2021", "_", label))
  
  result_df <- left_join(all_labels, result_df, by="label") %>%
    rename(dk_label = label, 
           label = label_)
  
  gp <- ggplot() +
    geom_brain(atlas = aseg(),
               data = result_df, 
               mapping = aes(fill=p_plot),
               position = position_brain(nrow=2),
               view = c("axial_3", "axial_5", "coronal_1", "sagittal")
    ) +
    scale_fill_viridis_c(
      option = "magma",
      na.value = "grey90",
      name = expression(-log[10](FDR~p))
    ) +
    theme_void() +
    theme(legend.position = "bottom")
  
  plot(gp)
  return(gp)
}

gp_aseg_1 <- make_aseg_brain_plot(results)
gp_aseg_2 <- make_aseg_brain_plot(results_2)

