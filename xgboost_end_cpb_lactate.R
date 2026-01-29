source("/studies/cardiac/support/cpb/adequate/analyses/fangm/utils.R")
source("/studies/cardiac/support/cpb/adequate/analyses/fangm/xgboost/utils_xgboost.R")
source("/studies/cardiac/support/cpb/adequate/analyses/fangm/xgboost/utils_3d.R")

# Define start and end lactate  ----
cldta_first_last <- get_first_last_lactate_df()

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
  get_time_summary(df_onpump_list[["map"]], "map", cutoff = 60, direction = "below") %>% 
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

## Pressors ---- 
cldta_feature_df <- get_pressor_df(cldta_feature_df)
cldta_feature_df <- get_pressor_df_mcg(cldta_feature_df)
cldta_feature_df <- get_precpb_pressors_df(cldta_feature_df)

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


summary(cldta_feature_df)
names(cldta_feature_df)

sapply(cldta_feature_df, function(x) {
  if(is.numeric(x)) {
    c(min = min(x, na.rm = TRUE),
      q1 = quantile(x, 0.25, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      mean = mean(x, na.rm = TRUE),
      q3 = quantile(x, 0.75, na.rm = TRUE),
      max = max(x, na.rm = TRUE),
      na_count = sum(is.na(x)))
  } else {
    table(x)
  }
})

cldta_feature_df <- cldta_feature_df %>% 
  mutate(gluc_max_above_mean = gluc_io_max - gluc_io_mean)

# Modeling ---- 

formula <- as.formula(lactate_lactate_end ~ 
  age +
  aocc_prop +
  bmi +
  bsa +
  bun_pr +
  cpb_duration +
  creat_pr +
  do2_interp_m2_mean +
  epi_dose_mg_hr +
  gfr_pr +
  gluc_io_mean +
  hgb_pr +
  hx_aorta +
  hx_csurg +
  hx_dm +
  lactate_lactate_start +
  lvef_pr +
  map_mean +
  mcv_pr +
  milr_dose_mg_hr +
  ohgb_io_mean +
  paco2_mean +
  pot_pr +
  precpb_epi_dose_mcg_kg +
  precpb_gluc_io_mean +
  precpb_map_mean +
  precpb_nepi_dose_mcg_kg +
  precpb_pao2_mean +
  precpb_phen_dose_mcg_kg +
  sp_atrmy +
  sp_cabg +
  std_dose_mg_hr +
  svo2_mean +
  tbill_pr +
  vtemp_mean +
  wbc_pr
)

# only keep people who were found in perfusionist database
# assume all others have missing phenylephrine dose 
surg_id_with_perfusionist_pressors <- df_onpump_list[["med_phen"]] %>% 
  filter(neobolus == 1) %>% 
  pull(surg_id) %>% 
  unique()

cldta_feature_df <- cldta_feature_df %>% 
  filter(surg_id %in% surg_id_with_perfusionist_pressors) 

model_data <- cldta_feature_df[, c(all.vars(formula), "surg_id")]
model_data_with_surg_id <- na.omit(model_data)
model_surg_id <- model_data$surg_id
model_data <- model_data_with_surg_id %>%  select(-surg_id)
dim(model_data)


## seed ----
set.seed(123)

# outcome and predictors
y <- model_data$lactate_lactate_end
X <- model_data[, setdiff(names(model_data), "lactate_lactate_end")]
k <- 10

params <- list(
  objective = "reg:squarederror",
  eval_metric = "rmse",
  eta = 0.05,
  max_depth = 5,
  subsample = 0.8,
  colsample_bytree = 0.8
)

xgboost_results <- k_fold_xgboost(X, y, k, params = params, get_partial = TRUE, get_ale = TRUE,
                                  quantile_points = 100) 
print_rmse_mean_ci(xgboost_results$performance_list)
partial_list <- xgboost_results$partial_list
model_list <- xgboost_results$model_list
folds <- xgboost_results$folds

dir.create("/studies/cardiac/support/cpb/adequate/analyses/fangm/xgboost/xgb_models", showWarnings = FALSE)

for (i in seq_along(model_list)) {
  xgb.save(
    model_list[[i]],
    fname = paste0("/studies/cardiac/support/cpb/adequate/analyses/fangm/xgboost/xgb_models/xgb_fold_", i, ".json")
  )
}


# Performance ----
plot_actual_vs_pred(X, y, model_list, folds, 
                    performance_list = xgboost_results$performance_list) + 
  coord_cartesian(xlim = c(0, 14), ylim = c(0, 14))


# VIMP ----
vimp_res <- plot_vimp(model_list, metric = "Gain")
plot_vimp(model_list, metric = "Cover", top_n = 10)
plot_vimp(model_list, metric = "Frequency", top_n = 10)


# Partial Plots ----
pps <- plot_pdp(xgboost_results, X, numeric_var_order = vimp_res$var_order)