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


source("~/R/LTC/utils/ucsf_data_loaders.R")
run <- "exp_km_ab"
load(paste("~/R/EDAP-data/LTC4/", run, ".Rdata", sep = ""))

get_ab_df <- source("~/R/LTC/get_ab_df.r")$value
ab_df <- get_ab_df()

ab_df %>% select(RID, AB_any) %>% group_by(RID) %>% 
  mutate(AB_any = if_any(AB_any)) %>% distinct() %>% filter(AB_any) %>% 
  select(RID) %>% unlist() -> ab_pos_rids_any

ucsf_data <- ucsf_longitudinal_all(only_vol = TRUE, filter_n = 0)
ucsf_controls <- ucsf_data %>% filter(!(RID %in% ab_pos_rids_any)) %>% filter(DX.bl == "CN")

dpm_df <- read.csv("~/R/EDAP-data/ADNI/DPM_ADNI.csv", header=TRUE)
ucsf_data <- dpm_df %>% select(RID, time_shift) %>% distinct() %>% right_join(ucsf_data, by="RID") %>%
  mutate(Time = Years + time_shift) %>% drop_na(Time)

Clusters <- data.frame(
  Cluster = adniLTC@Cluster,
  RID = adniLTC@RID
)

ucsf_data <- left_join(ucsf_data, Clusters, by = "RID")

datadic <- datadict_as_tibble(ADNIMERGE2::DATADIC) %>% filter(TBLNAME == "UCSFFSL") %>%
  select(FLDNAME, TEXT) %>% filter(grepl("ST\\d+[A-Z]{2}", FLDNAME)) %>%
  mutate(Region = sub(".*\\s", "", TEXT))


region_names <- data.frame()
for (varname in adniLTC@varNames) {
  descr <- unlist(datadic[which(datadic$FLDNAME == varname), "Region"])
  region_names <- rbind(region_names, data.frame(Variable=varname,
                                                 Region=descr,
                                                 Mean=mean(ucsf_data[[varname]], na.rm = TRUE), 
                                                 SD=sd(ucsf_data[[varname]], na.rm = TRUE)))
}

region_names |> arrange(desc(Mean)) -> region_names

get_brain_regions_fn <- source("~/R/EDAP-R/cluster-analysis/get_brain_regions.R")$value
brain_regions <- get_brain_regions_fn()

tran <- lapply(unique(Clusters$Cluster), function(c) {
  ucsf_data %>% filter(Cluster == c) %>% select(Time) %>% unlist() %>% quantile(c(0.05, 0.95))
})

tran

ucsf_data$Stage <- cut(ucsf_data$Time, 
                       breaks = c(-Inf, -1.5, 0, 5, 10, 15, Inf),
                       include.lowest = TRUE)


all_means <- data.frame()
all_collect <- data.frame()

all.idxs <- grepl("^ST\\d+[A-Z]*", colnames(ucsf_data)) 
all.vars <- colnames(ucsf_data)[all.idxs == 1]

for (varname in all.vars) {
  #mu <- ucsf_data %>% filter(Time < quantile(Time, 0.025)) %>% select(all_of(varname)) %>% 
  #  summarise(Mean = mean(.data[[varname]], na.rm=TRUE),
  #            SD = sd(.data[[varname]], na.rm=TRUE)) 
  
  mu <- ucsf_controls %>% select(RID, Months, all_of(varname)) %>% na.omit() %>%
    arrange(Months) %>% distinct(RID, .keep_all = TRUE) %>%
    summarise(Mean = mean(.data[[varname]], na.rm=TRUE),
              SD = sd(.data[[varname]], na.rm=TRUE)) 
  
  df <- ucsf_data %>% select(RID, Time, DX.bl, Cluster, Stage, all_of(varname)) %>%
    drop_na(Cluster) %>% mutate(z = (.data[[varname]]-mu$Mean)/mu$SD)
  
  collect <- select(df, RID, Time, Cluster, Stage, z) %>% 
    mutate(Variable = varname)
  
  all_collect <- rbind(all_collect, collect)
  
  means <- df %>%
    group_by(Cluster, Stage) %>%
    summarise(mean_values = mean(z, na.rm = TRUE), .groups = "drop") 
  colnames(means) <- c("Cluster", "Stage", varname)
  
  if (ncol(all_means) == 0) {
    all_means <- means
  } else {
    all_means <- left_join(all_means, means, by=c("Cluster", "Stage")) 
  }
}

all_means <- pivot_longer(all_means, cols = matches("^ST\\d+[A-Z]*"), names_to = "Variable", values_to = "Mean")

bap <- brain_regions |> select(Variables, Lobe, BraakStage)
bap <- datadic %>% select(FLDNAME, Region) %>% right_join(bap, by = c("Region" = "Variables"))

all_means <- left_join(all_means, bap, by = c("Variable" = "FLDNAME"))
all_collect <- left_join(all_collect, bap, by = c("Variable" = "FLDNAME"))

df <- all_means %>% select(Cluster, Stage, Region, Mean) %>% 
  mutate(Side = str_extract(Region, "(Left|Right)")) %>% drop_na(Side) %>%
  mutate(Region = str_remove(Region, "(Left|Right)")) %>%
  pivot_wider(names_from = Side, values_from = Mean)

df_all <- all_collect %>% select(RID, Time, Cluster, Stage, Region, z) %>%
  mutate(Side = str_extract(Region, "(Left|Right)")) %>% drop_na(Side) %>%
  mutate(Region = str_remove(Region, "(Left|Right)")) %>%
  pivot_wider(names_from = Side, values_from = z, values_fn = mean)


lat_plot <- ggplot(df, aes(x=Right, y=Left, color=Cluster)) +
  geom_point() +
  scale_color_paletteer_d("ggthemes::Tableau_10") +
  geom_abline(linetype=2) +
  facet_wrap(~Stage)
plot(lat_plot)


df_ <- filter(df_all, Stage=="(15, Inf]")
lat_plot <- ggplot(df_all, aes(x=Right, y=Left, color=Stage)) +
  geom_point() +
  #scale_color_paletteer_d("ggthemes::Tableau_10") +
  scale_color_brewer(palette = "YlOrRd") +
  geom_abline(linetype=2) +
  facet_wrap(~Cluster) +
  theme_minimal()
plot(lat_plot)


corr.results <- df_all %>%
  group_by(Cluster, Stage) %>%
  summarise(
    r     = cor.test(Left, Right, method = "pearson")$estimate,
    p_val = cor.test(Left, Right, method = "pearson")$p.value,
    n     = sum(!is.na(Left) & !is.na(Right)),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_val, method = "BH"))

wilcox.results <- df_all %>%
  group_by(Cluster, Stage) %>%
  summarise(
    p_val  = wilcox.test(Left, Right, paired = TRUE)$p.value,
    median_diff = median(Left - Right, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_val, method = "BH")) %>%
  mutate(p_val = round(p_val, digits=4),
         p_adj = round(p_adj, digits=4)) %>%
  filter(p_adj < 0.05)
