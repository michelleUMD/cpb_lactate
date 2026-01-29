source("/studies/cardiac/support/cpb/adequate/analyses/fangm/utils.R")
source("/studies/cardiac/support/cpb/adequate/analyses/fangm/xgboost/utils_xgboost.R")

# Define start and end lactate  ----
cldta_first_last <- get_first_peak_lactate_df(iv_evtof1_upper = 6)
get_first_last_lactate_plots(cldta_first_last)
plot_lactate_over_time_quartiles(cldta_first_last, n_bins = 20)
plot_time_hist(cldta_first_last, bins = 30)


# Feature engineering ---- 
cldta_feature_df <- cldta_first_last
cldta_feature_df <- cldta_feature_df %>%
  get_mean_summary(df_onpump_list[["ci"]], "ci") %>%
  get_mean_summary(df_onpump_list[["do2_interp"]], "do2_interp_m2") %>%
  get_mean_summary(df_onpump_list[["map"]], "map") %>%
  get_mean_summary(df_onpump_list[["svo2"]], "svo2") %>%
  get_mean_summary(df_onpump_list[["glucose"]], "gluc_io") %>%
  get_mean_summary(df_onpump_list[["hgb"]], "hgb_io") %>%
  get_mean_summary(df_onpump_list[["pao2"]], "pao2") %>% 
  get_mean_summary(df_onpump_list[["temp"]], "temp") %>% 
  get_mean_summary(df_onpump_list[["ohgb"]], "ohgb_io") %>% 
  get_mean_summary(df_onpump_list[["paco2"]], "paco2") %>% 
  get_min_summary(df_onpump_list[["map"]], "map") %>% 
  get_max_summary(df_onpump_list[["glucose"]], "gluc_io") %>%
  get_sd_summary(df_onpump_list[["map"]], "map") %>% 
  get_sd_summary(df_onpump_list[["ci"]], "ci")

atemp_omit_first_five <- df_onpump_list[["art_temp"]] %>% 
  filter(iv_evth > 5/60)
vtemp_omit_first_five <- df_onpump_list[["ven_temp"]] %>% 
  filter(iv_evth > 5/60)

cldta_feature_df <- cldta_feature_df %>%
  get_mean_summary(atemp_omit_first_five, "atemp") %>%
  get_mean_summary(vtemp_omit_first_five, "vtemp") 


cldta_feature_df <- cldta_feature_df %>%
  get_time_summary(df_onpump_list[["do2_interp"]], "do2_interp", cutoff = 500, direction = "below") %>% 
  get_time_summary(df_onpump_list[["do2_interp"]], "do2_interp", cutoff = 800, direction = "above") %>%
  get_time_summary(df_onpump_list[["map"]], "map", cutoff = 70, direction = "below") %>% 
  get_time_summary(df_onpump_list[["map"]], "map", cutoff = 85, direction = "above") %>% 
  get_time_summary(df_onpump_list[["svo2"]], "svo2", cutoff = 75, direction = "below") %>% 
  get_time_summary(df_onpump_list[["svo2"]], "svo2", cutoff = 85, direction = "above")

cldta_feature_df <- convert_to_cpb_frac(cldta_feature_df)

# Pre-CPB means 
cldta_feature_df <- cldta_feature_df %>% 
  get_precpb_mean_summary(df_list[["map"]], "map", within_x_hours = 2) %>% 
  get_precpb_mean_summary(df_list[["hgb"]], "hgb_io", within_x_hours = 24) %>% 
  get_precpb_mean_summary(df_list[["pao2"]], "pao2", within_x_hours = 6) %>% 
  get_precpb_mean_summary(df_list[["glucose"]], "gluc_io", within_x_hours = 24) 

# Post-CPB means 
cldta_feature_df <- cldta_feature_df %>%
  get_postcpb_mean_summary(df_list[["map"]], "map") %>%
  get_postcpb_mean_summary(df_list[["hgb"]], "hgb_io") %>%
  get_postcpb_mean_summary(df_list[["pao2"]], "pao2")


## Pressors ---- 
cldta_feature_df <- get_pressor_df(cldta_feature_df)
cldta_feature_df <- get_pressor_df_mcg(cldta_feature_df)
cldta_feature_df <- get_precpb_pressors_df(cldta_feature_df)
cldta_feature_df <- get_postcpb_pressors_df(cldta_feature_df)

cldta_feature_df <- convert_to_mcg_kg_min(cldta_feature_df)

cldta_feature_df <- cldta_feature_df %>%
  mutate(
    std_dose_mcg_kg_min = (phen_dose_mcg_kg_min / 10) +
      nepi_dose_mcg_kg_min +
      vaso_dose_mcg_kg_min * 2.5,
    
    n_pressors = (phen_dose_mcg_kg_min > 0) +
      (nepi_dose_mcg_kg_min > 0) +
      (vaso_dose_mcg_kg_min > 0)
  )



# Pre-op variables ----
cldta_feature_df <- get_preop_features(
  df_to_join_to       = cldta_feature_df,
  preop_vars_keep_na  = c(
    "bmi", "age", "bsa", "blrbn_pr", "bun_pr", "creat_pr",
    "hgb_pr", "lvef_pr", "tbill_pr", "mcv_pr", "pot_pr", "ptinr_pr",
    "wbc_pr", "gfr_pr", "meld", "wt", "iv_aocc"
  ), 
  preop_vars_na_to_0  = c("female",
                          "hx_chf", "hx_csurg", "afib_pr", "hx_copd", "hx_smoke", "hx_dm",
                          "hx_dyslp", "hx_endo", "hx_htn", "hx_mi", "hx_pad", "hx_cva",
                          "hx_aorta", "hx_arysu", "hx_cabg", "hx_CVInt", "hx_cngsu", "hx_icd",
                          "hx_ppm", "dial_pr", "varr_pr", "sp_valve", "sp_av", "sp_mv",
                          "sp_pv", "sp_tv", "sp_aorta", "sp_cabg", "sp_afib", "sp_asdpf",
                          "sp_chd", "sp_atrmy", "sp_peric", "sp_sepmy", "sp_lv", "sp_tx",
                          "sp_cend", "robot" 
  ),
  dummy_vars          = c("surg_yr", "nyha_pr", "race_comb", "diabetes_preop_status"),
  dummy_remove_first  = FALSE
)

cldta_feature_df <- cldta_feature_df %>%
  mutate(across(
    c("bmi", "age", "bsa", "female", "blrbn_pr", "bun_pr", "creat_pr",
      "hgb_pr", "lvef_pr", "tbill_pr", "mcv_pr", "pot_pr", "ptinr_pr",
      "wbc_pr", "gfr_pr", "meld", "precpb_pao2_mean", "precpb_gluc_io_mean", "precpb_hgb_io_mean"),
    ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x)
  ))


cldta_feature_df <- cldta_feature_df %>%
  mutate(across(
    .cols = contains("_dose_mcg_kg_min"),
    .fns = ~ .x * wt * 60 / 1000,
    .names = "{str_replace(.col, '_dose_mcg_kg_min', '_dose_mg_hr')}"
  ))

cldta_feature_df <- cldta_feature_df %>% 
  mutate(iv_aocc_hr = iv_aocc / 60,
         aocc_prop = 100 * pmin(iv_aocc_hr / cpb_duration, 1))




# Add baseline clusters in ----
# cluster_df <- get_baseline_clusters()
# 
# cldta_feature_df <- cldta_feature_df %>%
#   left_join(cluster_df, by = "surg_id")

summary(cldta_feature_df)
names(cldta_feature_df)

# Modeling ---- 

formula <- as.formula(
  lactate_lactate_end ~ lactate_lactate_start +
    cpb_duration +
    aocc_prop + 
    iv_evtof1_lactate_end + 
    do2_interp_m2_mean + 
    # ci_mean + 
    # hgb_io_mean + 
    # pao2_mean + 
    map_mean +
    svo2_mean +
    vtemp_mean + 
    ohgb_io_mean +
    paco2_mean + 
    std_dose_mg_hr + 
    epi_dose_mg_hr +
    milr_dose_mg_hr +
    nitp_dose_mg_hr + 
    precpb_phen_dose_mcg_kg +
    precpb_epi_dose_mcg_kg +
    precpb_nepi_dose_mcg_kg +
    precpb_vaso_dose_mcg_kg +
    precpb_nitp_dose_mcg_kg + 
    postcpb_phen_dose_mcg_kg_min +
    postcpb_epi_dose_mcg_kg_min +
    postcpb_nepi_dose_mcg_kg_min +
    postcpb_vaso_dose_mcg_kg_min +
    precpb_map_mean + 
    precpb_hgb_io_mean +
    precpb_pao2_mean + 
    # postcpb_map_mean +
    # postcpb_hgb_io_mean +
    # postcpb_pao2_mean + 
    gluc_io_mean +
    bmi + age + bsa + 
    wbc_pr + hgb_pr + creat_pr + lvef_pr + bun_pr + pot_pr + 
    tbill_pr + sp_cabg + hx_dm + robot + sp_aorta + hx_copd 
    # cluster
    # gfr_pr
    # mcv_pr + ptinr_pr +
    # hx_chf + hx_csurg + afib_pr + hx_copd + hx_smoke +
    # hx_dyslp + hx_endo + hx_htn + hx_mi + hx_pad + hx_cva +
    # hx_aorta + hx_arysu + hx_cabg + hx_CVInt + hx_cngsu + hx_icd +
    # hx_ppm + dial_pr + varr_pr +
    # sp_valve + sp_av + sp_mv + sp_pv + sp_tv + sp_aorta +
    # sp_afib + sp_asdpf + sp_chd + sp_atrmy + sp_peric + sp_sepmy +
    # sp_lv + sp_tx + sp_cend +
    # nyha_pr + race_comb
)


model_data <- na.omit(cldta_feature_df[, all.vars(formula)])
model_data <- as.data.frame(model_data)
dim(model_data)


# seed ----
set.seed(123)

# outcome and predictors
y <- model_data$lactate_lactate_end
X <- model_data[, setdiff(names(model_data), "lactate_lactate_end")]
k <- 10
params <- list(
  objective = "reg:squarederror",
  eval_metric = "rmse",
  eta = 0.02,
  max_depth = 5,
  subsample = 0.8,
  colsample_bytree = 0.8
)

xgboost_results <- k_fold_xgboost(X, y, k, params = params, 
                                  get_partial = TRUE, get_ale = TRUE,
                                  quantile_points = 100) 
print_rmse_mean_ci(xgboost_results$performance_list)
partial_list <- xgboost_results$partial_list
model_list <- xgboost_results$model_list
folds <- xgboost_results$folds


# Performance ----
plot_actual_vs_pred(X, y, model_list, folds, performance_list = xgboost_results$performance_list)

# VIMP ----
vimp_res <- plot_vimp(model_list, metric = "Gain", top_n = 40)
vimp_res$plot 

# Partial Plots ----
pps <- plot_pdp(partial_list, X, variables = NULL, numeric_var_order = vimp_res$var_order)
pps$numeric 

pps$factor

plot_pdp(xgboost_results, X, variables = vimp_res$var_order[1:25])

plot_pdp(partial_list, X, variables = c("ci_mean", "map_mean", "svo2_mean"), numeric_var_order = c("ci_mean", "map_mean", "svo2_mean"))
plot_pdp(partial_list, X,   variables = c(
  "phen_dose_mcg_kg_min",
  "epi_dose_mcg_kg_min",
  "nepi_dose_mcg_kg_min",
  "vaso_dose_mcg_kg_min",
  "milr_dose_mcg_kg_min",
  "precpb_phen_dose_mcg_kg",
  "precpb_epi_dose_mcg_kg",
  "precpb_nepi_dose_mcg_kg",
  "precpb_vaso_dose_mcg_kg"
), numeric_var_order = vimp_res$var_order)


plot_pdp(partial_list, X, variables = c("hx_aorta", "hx_dyslp", "sp_cabg", 
                                        "hx_pad"))


# Conditional PDP ----

## 2 independent variables ---- 

age_groups <- cut(
  X$age,
  breaks = c(0, 50, 70, 80, Inf),   # adjust bins as needed
  labels = c("<50", "50–69", "70–79", "≥80"),
  right = FALSE,
  ordered_result = TRUE
)

conditional_pdp_age_map <- compute_conditional_pdp(
  model_list = model_list,
  X = X,  
  group_var = age_groups,
  facet_var1 = NULL, facet_var2 = NULL,
  pred_var = "map_mean", 
  grid.resolution = 50)


plot_conditional_pdp(
  summary_df = conditional_pdp_age_map,
  X = X,
  group_var = age_groups,
  facet_var1 = NULL,
  facet_var2 = NULL,
  pred_var = "map_mean",
  pred_var_label = "Mean MAP (mmHg)",
  group_label = "Age groups"
)


map_groups <- cut(
  X$ci_mean,
  breaks = c(0, 1.6, 1.8, 1.9, 2.0, 2.2, Inf),
  labels = c("<1.6", "1.6–1.79", "1.8-1.89", "1.9–1.99", "2.0–2.19", "≥2.2"),
  right = FALSE,
  ordered_result = TRUE
)

conditional_pdp_age_map <- compute_conditional_pdp(
  model_list = model_list,
  X = X,  
  group_var = map_groups,
  facet_var1 = NULL, facet_var2 = NULL,
  pred_var = "wbc_pr", 
  grid.resolution = 50)


plot_conditional_pdp(
  summary_df = conditional_pdp_age_map,
  X = X,
  group_var = map_groups,
  facet_var1 = NULL,
  facet_var2 = NULL,
  pred_var = "wbc_pr",
  pred_var_label = "Pre-op WBCs",
  group_label = "Mean CI"
)

## 3 independent variables ----
cut_groups <- cut(
  X$phen_dose_mcg_kg,
  breaks = c(0, 20, 40, 60, 80, 100, 150, Inf),
  labels = c("<20", "20–39", "40–59", "60–79", "80–99", "100–149", "≥150"),
  right = FALSE,
  ordered_result = TRUE
)

facet_groups <- cut(
  X$ci_mean,
  breaks = c(0, 1.6, 1.8, 1.9, 2.0, 2.2, Inf),
  labels = c("<1.6", "1.6–1.79", "1.8-1.89", "1.9–1.99", "2.0–2.19", "≥2.2"),
  right = FALSE,
  ordered_result = TRUE
)


# Table of N in each combination 
data.frame(
  ci_group = facet_groups,
  phen_group = cut_groups
) %>%
  tbl_summary(
    by = ci_group,            
    statistic = list(all_categorical() ~ "{n}"), 
    missing = "no"
  ) %>%
  add_overall() 


conditional_pdp_ci_phen_map <- compute_conditional_pdp(
  model_list = model_list,
  X = X,  
  group_var = cut_groups,
  facet_var1 = facet_groups, facet_var2 = NULL,
  pred_var = "map_mean", 
  grid.resolution = 50)


plot_conditional_pdp(
  summary_df = conditional_pdp_ci_phen_map,
  X = X,
  group_var = cut_groups,
  facet_var1 = facet_groups,
  facet_var2 = NULL,
  pred_var = "map_mean",
  pred_var_label = "Mean MAP (mmHg)",
  group_label = "Phenylephrine groups",
  facet_label1 = "Mean CI"
)

## 4 independent variables ----
cut_groups <- cut(
  X$phen_dose_mcg_kg,
  breaks = c(0, 40, 80, 120, Inf),
  labels = c("<40", "40–79", "80–119", "≥120"),
  right = FALSE,
  ordered_result = TRUE
)

facet_groups1 <- cut(
  X$ci_mean,
  breaks = c(0, 1.8, 2.0, 2.2, Inf),
  labels = c("<1.8", "1.8–1.99", "2.0–2.19", "≥2.2"),
  right = FALSE,
  ordered_result = TRUE
)

facet_groups2 <- cut(
  X$age,
  breaks = c(0, 50, 70, 80, Inf),
  labels = c("<50", "50–69", "70–79", "≥80"),
  right = FALSE,
  ordered_result = TRUE
)

conditional_pdp_ci_phen_age_map <- compute_conditional_pdp(
  model_list = model_list,
  X = X,                     
  group_var = cut_groups,     
  facet_var1 = facet_groups1, 
  facet_var2 = facet_groups2, 
  pred_var = "map_mean",
  grid.resolution = 50
)

plot_conditional_pdp(
  summary_df = conditional_pdp_ci_phen_age_map,
  X = X,
  group_var = cut_groups,
  facet_var1 = facet_groups1,
  facet_var2 = facet_groups2,
  pred_var = "map_mean",
  pred_var_label = "Mean MAP (mmHg)",
  group_label = "Total phenylephrine dose (mcg/kg)",
  facet_label1 = "Mean CI (L/min/m²)",
  facet_label2 = "Age Group"
)




# Sepsis: WBC, & CI ----
cut_groups <- cut(
  X$map_mean,
  breaks = c(0, 65, 70, 75, Inf),
  labels = c("<65", "65–69", "70–74", "≥75"),
  right = FALSE,
  ordered_result = TRUE
)

facet_groups1 <- cut(
  X$wbc_pr,
  breaks = c(0, 4.5, 11, Inf),
  labels = c("<4.5", "4.5–10.9", "≥11"),
  right = FALSE,
  ordered_result = TRUE
)

facet_groups2 <- cut(
  X$ci_mean,
  breaks = c(0, 1.8, 2.0, 2.2, Inf),
  labels = c("<1.8", "1.8–1.99", "2.0–2.19", "≥2.2"),
  right = FALSE,
  ordered_result = TRUE
)

conditional_pdp_ci_phen_wbc_map <- compute_conditional_pdp(
  model_list = model_list,
  X = X,                     
  group_var = cut_groups,     
  facet_var1 = facet_groups1, 
  facet_var2 = facet_groups2, 
  pred_var = "phen_dose_mcg_kg_min",
  grid.resolution = 50
)

plot_conditional_pdp(
  summary_df = conditional_pdp_ci_phen_wbc_map,
  X = X,
  group_var = cut_groups,
  facet_var1 = facet_groups1,
  facet_var2 = facet_groups2,
  pred_var = "phen_dose_mcg_kg_min",
  group_label = "Mean MAP (mmHg)",
  pred_var_label = "Total phenylephrine dose (mcg/kg/min)",
  facet_label1 = "WBC Group (cells/µL)",
  facet_label2 = "Mean CI (L/min/m²)",
  x_lim_truncation = TRUE,
  lower_lim_truncation = 0,
  upper_lim_truncation = 0.95
)


# Optimal MAP: Age x BMI ----
cut_groups <- cut(
  X$map_mean,
  breaks = c(0, 65, 70, 75, Inf),
  labels = c("<65", "65–69", "70–74", "≥75"),
  right = FALSE,
  ordered_result = TRUE
)
summary(cut_groups)

facet_groups1 <- cut(
  X$bmi,
  breaks = c(0, 25, 30, Inf),
  labels = c("<25", "25–29.9", "≥30"),
  right = FALSE,
  ordered_result = TRUE
)
summary(facet_groups1)

facet_groups2 <- cut(
  X$age,
  breaks = c(0, 55, 65, 75, Inf),
  labels = c("<55", "55–64", "65–74", "≥75"),
  right = FALSE,
  ordered_result = TRUE
)
summary(facet_groups2)


conditional_pdp_ci_phen_wbc_map <- compute_conditional_pdp(
  model_list = model_list,
  X = X,                     
  group_var = cut_groups,     
  facet_var1 = facet_groups1, 
  facet_var2 = facet_groups2, 
  pred_var = "phen_dose_mcg_kg_min",
  grid.resolution = 50
)

plot_conditional_pdp(
  summary_df = conditional_pdp_ci_phen_wbc_map,
  X = X,
  group_var = cut_groups,
  facet_var1 = facet_groups1,
  facet_var2 = facet_groups2,
  pred_var = "phen_dose_mcg_kg_min",
  group_label = "Mean MAP (mmHg)",
  facet_label1 = "BMI",
  facet_label2 = "Age",
  x_lim_truncation = TRUE,
  lower_lim_truncation = 0,
  upper_lim_truncation = 0.95,
  switch_facets = FALSE
)


