
# Extra utils for clustering and visualization
source("~/R/LTC/utils/cluster_utils.R")

multi_cohort_df <- read.csv("~/R/EDAP-data/MULTI_COHORT_4.csv", header = TRUE)

# Filter out NACC
multi_cohort_df <- filter_out(multi_cohort_df, Cohort == "NACC")

all.vars <- c(grepv("^(RH_|LH_|CC_)", colnames(multi_cohort_df)), "BRAINSTEM")

load("~/R/EDAP-data/LTC_MC/new/nlmmBasic_AO_fp.Rdata")

run <- "exp_km_ab_ao"
load(paste("~/R/EDAP-data/LTC_MC/new/", run, ".Rdata", sep = ""))
Clusters <- data.frame(
  Cluster = multiLTC@Cluster,
  RID = multiLTC@RID
)

func_params <- nlmmBasic$func_params

silly <- silhouette_score(Clusters, func_params, plot=TRUE)

sil_df <- as.data.frame(silly)


