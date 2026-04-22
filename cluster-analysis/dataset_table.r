library(ADNIMERGE2)
library(dplyr)
library(lubridate)
library(stringr)

source("~/R/LTC/utils/ucsf_data_loaders.R")

run <- "exp_km_ab"
load(paste("~/R/EDAP-data/LTC4/", run, ".Rdata", sep = ""))

ucsf_data <- ucsf_longitudinal_all(only_vol = TRUE, filter_n = 0)

Clusters <- data.frame(
  Cluster = adniLTC@Cluster,
  RID = adniLTC@RID
)


adni_demog <- ADNIMERGE2::PTDEMOG %>% filter(RID %in% adniLTC@RID) %>%
  mutate(VISDATE = parse_date_time(VISDATE, orders=c("Ymd")),
         PTDOB = parse_date_time(PTDOB, orders=c("m/y")), 
         AGE = interval(PTDOB, VISDATE) %/% years(1)) %>% group_by(RID) %>%
  mutate(AGE = min(AGE, na.rm=TRUE),
         PTEDUCAT = max(PTEDUCAT, na.rm=TRUE)) %>% ungroup() %>%
  select(RID, AGE, PTGENDER, PTEDUCAT) %>% na.omit() %>% distinct() %>% arrange(RID, AGE)
 
dx.df <- ADNIMERGE2::DXSUM %>% select(RID, VISCODE2, DIAGNOSIS, DXCONFID) %>% drop_na(DIAGNOSIS) %>%
  mutate(DX.bl = ifelse(VISCODE2 == "bl", DIAGNOSIS, NA)) %>%
  group_by(RID) %>% mutate(
    DX.bl = first(na.omit(DX.bl))
  ) %>% ungroup() %>%
  mutate(DIAGNOSIS = factor(DIAGNOSIS, levels = c("CN", "MCI", "Dementia"), ordered = TRUE),
         DX.bl = factor(DX.bl, levels = c("CN", "MCI", "Dementia"), ordered = TRUE)) %>%
  group_by(RID) %>% mutate(
    DX.highest = max(DIAGNOSIS)
  ) %>% ungroup()

dx.bl <- select(dx.df, RID, DX.bl, DX.highest) %>% distinct() %>%
  left_join(Clusters, by = "RID") %>% drop_na(Cluster) %>%
  group_by(RID) %>%
  mutate(
    CN2CI = ifelse(DX.bl == "CN", DX.highest != "CN", NA),
    MCI2AD = ifelse(DX.bl == "MCI", DX.highest == "Dementia", NA)
  ) %>% ungroup()


apoe.df <- ADNIMERGE2::APOERES %>% filter(RID %in% adniLTC@RID) %>% 
  select(RID, GENOTYPE) %>% 
  mutate(APOE4 = str_count(GENOTYPE, "4"))


table_df <- adni_demog %>% left_join(apoe.df, by="RID") 
table_df <- dx.bl %>% select(RID, DX.bl, CN2CI, MCI2AD) %>% distinct() %>% right_join(table_df, by="RID")

dpm_df <- read.csv("~/R/EDAP-data/ADNI/DPM_ADNI.csv", header=TRUE)
table_df <- dpm_df %>% select(RID, time_shift) %>% distinct() %>% right_join(table_df, by="RID")

table_df <- ucsf_data %>% group_by(RID) %>% summarise(n_mri = n()) %>% right_join(table_df, by="RID")

table_df$Onset <- table_df$AGE - table_df$time_shift

# CSF AB

get_ab_df <- source("~/R/LTC/get_ab_df.r")$value
ab_df <- get_ab_df() %>% filter(RID %in% adniLTC@RID)

table_df <- ab_df %>% select(RID, A4240.re, A4240.ms) %>% group_by(RID) %>% 
  summarise(A4240.ms=min(A4240.ms), A4240.re=min(A4240.re)) %>%
  right_join(table_df, by="RID")

ab_pet <- ab_df %>% select(RID, TRACER, SUMMARY_SUVR) %>% drop_na(TRACER) %>%
  pivot_wider(id_cols=RID, names_from = TRACER, values_from = SUMMARY_SUVR, values_fn = max)
table_df <- left_join(table_df, ab_pet, by="RID")

# tau-PET

get_tau_pet <- source("~/R/LTC/cluster-analysis/get_tau_pet.r")$value

ucbtaupet <- get_tau_pet() %>% select(RID, VISCODE2, Months, Years, 
                                      CTX_ENTORHINAL_SUVR, META_TEMPORAL_SUVR, CTX_INSULA_SUVR, CTX_SUPERIORFRONTAL_SUVR, CTX_SUPERIORPARIETAL_SUVR)

ucbtaupet <- ucbtaupet %>% filter(RID %in% adniLTC@RID) %>% group_by(RID) %>%
  summarise(CTX_ENTORHINAL_SUVR = max(CTX_ENTORHINAL_SUVR), 
            META_TEMPORAL_SUVR = max(META_TEMPORAL_SUVR), 
            CTX_INSULA_SUVR = max(CTX_INSULA_SUVR), 
            CTX_SUPERIORFRONTAL_SUVR = max(CTX_SUPERIORFRONTAL_SUVR),
            CTX_SUPERIORPARIETAL_SUVR = max(CTX_SUPERIORPARIETAL_SUVR)) 

table_df <- left_join(table_df, ucbtaupet, by = "RID")

table_df %>% group_by(DX.bl) %>% summarise(
  n = n(),
  AGE.mean = mean(AGE),
  AGE.sd = sd(AGE),
  Male = mean(PTGENDER == "Male", na.rm=T)*100,
  Female = 100-Male,
  EDUCAT.mean = mean(PTEDUCAT),
  EDUCAT.sd = sd(PTEDUCAT),
  APOE4.0 = mean(APOE4 == 0, na.rm=T)*100,
  APOE4.1 = mean(APOE4 == 1, na.rm=T)*100,
  APOE4.2 = mean(APOE4 == 2, na.rm=T)*100,
  n_mri.mean = mean(n_mri),
  n_mri.sd = sd(n_mri),
  A4240.ms.mean = mean(A4240.ms, na.rm=T),
  A4240.ms.sd = sd(A4240.ms, na.rm=T),
  A4240.re.mean = mean(A4240.re, na.rm=T),
  A4240.re.sd = sd(A4240.re, na.rm=T),
  FBP.mean = mean(FBP, na.rm=T),
  FBP.sd = sd(FBP, na.rm=T),
  FBB.mean = mean(FBB, na.rm=T),
  FBB.sd = sd(FBB, na.rm=T),
  PIB.mean = mean(PIB, na.rm=T),
  PIB.sd = sd(PIB, na.rm=T),
  BRAAK1.mean = mean(CTX_ENTORHINAL_SUVR, na.rm=T),
  BRAAK1.sd = sd(CTX_ENTORHINAL_SUVR, na.rm=T),
  BRAAK2.mean = mean(META_TEMPORAL_SUVR, na.rm=T),
  BRAAK2.sd = sd(META_TEMPORAL_SUVR, na.rm=T),
  INSULA.mean = mean(CTX_INSULA_SUVR, na.rm=T),
  INSULA.sd = sd(CTX_INSULA_SUVR, na.rm=T),
  SUPFRONT.mean = mean(CTX_SUPERIORFRONTAL_SUVR, na.rm=TRUE),
  SUPFRONT.sd = sd(CTX_SUPERIORFRONTAL_SUVR, na.rm=TRUE),
  SUPPAR.mean = mean(CTX_SUPERIORPARIETAL_SUVR, na.rm=T),
  SUPPAR.sd = sd(CTX_SUPERIORPARIETAL_SUVR, na.rm=T),
) -> summary_all

table_df %>% ungroup() %>% summarise(
  n = n(),
  AGE.mean = mean(AGE),
  AGE.sd = sd(AGE),
  Male = mean(PTGENDER == "Male", na.rm=T)*100,
  Female = 100-Male,
  EDUCAT.mean = mean(PTEDUCAT),
  EDUCAT.sd = sd(PTEDUCAT),
  APOE4.0 = mean(APOE4 == 0, na.rm=T)*100,
  APOE4.1 = mean(APOE4 == 1, na.rm=T)*100,
  APOE4.2 = mean(APOE4 == 2, na.rm=T)*100,
  n_mri.mean = mean(n_mri),
  n_mri.sd = sd(n_mri),
  A4240.ms.mean = mean(A4240.ms, na.rm=T),
  A4240.ms.sd = sd(A4240.ms, na.rm=T),
  A4240.re.mean = mean(A4240.re, na.rm=T),
  A4240.re.sd = sd(A4240.re, na.rm=T),
  FBP.mean = mean(FBP, na.rm=T),
  FBP.sd = sd(FBP, na.rm=T),
  FBB.mean = mean(FBB, na.rm=T),  
  FBB.sd = sd(FBB, na.rm=T),
  PIB.mean = mean(PIB, na.rm=T),
  PIB.sd = sd(PIB, na.rm=T),
  BRAAK1.mean = mean(CTX_ENTORHINAL_SUVR, na.rm=T),
  BRAAK1.sd = sd(CTX_ENTORHINAL_SUVR, na.rm=T),
  BRAAK2.mean = mean(META_TEMPORAL_SUVR, na.rm=T),
  BRAAK2.sd = sd(META_TEMPORAL_SUVR, na.rm=T),
  INSULA.mean = mean(CTX_INSULA_SUVR, na.rm=T),
  INSULA.sd = sd(CTX_INSULA_SUVR, na.rm=T),
  SUPFRONT.mean = mean(CTX_SUPERIORFRONTAL_SUVR, na.rm=TRUE),
  SUPFRONT.sd = sd(CTX_SUPERIORFRONTAL_SUVR, na.rm=TRUE),
  SUPPAR.mean = mean(CTX_SUPERIORPARIETAL_SUVR, na.rm=T),
  SUPPAR.sd = sd(CTX_SUPERIORPARIETAL_SUVR, na.rm=T),
  DX.bl = "All"
) %>% rbind(summary_all) %>% tibble::column_to_rownames("DX.bl") -> summary_all

mean_cols <- colnames(summary_all)[grepl(".mean", colnames(summary_all))]
sd_cols <- gsub(".mean", ".sd", mean_cols)



mean_df <- summary_all %>% select(all_of(mean_cols)) %>%
  t() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("Variable") %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  mutate(Variable = gsub(".mean", "", Variable)) 
  
sd_df <- summary_all %>% select(all_of(sd_cols)) %>%
  t() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("Variable") %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  mutate(Variable = gsub(".sd", "", Variable))

latex_strings <- mapply(
  function(m, s) sprintf("$%.3f \\pm %.3f$", m, s),
  mean_df[-1], sd_df[-1],
  SIMPLIFY = FALSE
) |> as.data.frame() |> mutate(Variable = mean_df$Variable)

perc_df <- summary_all %>% select(-all_of(c(mean_cols, sd_cols, "n"))) %>%
  t() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("Variable") %>%
  mutate(across(where(is.numeric), ~ sprintf("$%.1f\\%%$", .)))

rest_df <- summary_all %>% select(all_of(c("n"))) %>%
  t() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("Variable") %>%
  mutate(across(where(is.numeric), ~ sprintf("$%d$", .)))
  
latex_strings <- rbind(rest_df, perc_df, latex_strings)

variable_names <- c("n", "Sex (M)", "Sex (F)", "\\textit{APOE} $\\epsilon$4 (0)", "\\textit{APOE} $\\epsilon$4 (1)", "\\textit{APOE} $\\epsilon$4 (2)",
                    "Age (Y)", "Educat. (Y)", "Num. MRI", 
                    "A$\\beta$42/A$\\beta$40 (MS)", "A$\\beta$42/A$\\beta$40 (IA)", 
                    "FBP", "FBB", "PiB", 
                    "Braak I-II", "Braak III-IV", "Insula", "Sup. Frontal", "Sup. Parietal")
latex_strings$Variable <- variable_names

library(knitr)
latex_strings <- select(latex_strings, Variable, CN, MCI, Dementia, All)
kable(latex_strings, format = "latex", escape = FALSE, booktabs = TRUE)





