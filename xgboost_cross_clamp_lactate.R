source("/studies/cardiac/support/cpb/adequate/analyses/fangm/utils.R")
source("/studies/cardiac/support/cpb/adequate/analyses/fangm/xgboost/utils_xgboost.R")
source("/studies/cardiac/support/cpb/adequate/analyses/fangm/xgboost/utils_3d.R")


# Find valid Cross-Clamp Times ----
has_aocc <- cldta_onpump_baseline_df %>% 
  filter(!is.na(dtmccon) & !is.na(dtmccof))
dim(has_aocc)


has_aocc <- has_aocc %>%
  mutate(iv_aocc_calc = as.numeric(difftime(dtmccof, dtmccon, units = "hours")),
         iv_aocc_hr   = iv_aocc / 60)

ggplot(has_aocc, aes(x = iv_aocc_calc, y = iv_aocc_hr)) +
  geom_point(alpha = 0.6, color = "blue") +
  labs(x = "STS Cross Clamp Time (hours)", y = "Calculated Cross Clamp Time (Time Off - On)") +
  theme_minimal()


has_aocc <- has_aocc %>%
  mutate(aocc_between_cpb = factor(dtmccon >= dtmbpon1 & dtmccof <= dtmbpof1), 
         aocc_end_less_than_start = factor(dtmccon > dtmbpof1))

summary(has_aocc$aocc_between_cpb)
summary(has_aocc$aocc_end_less_than_start)

has_aocc_filtered <- has_aocc %>%
  filter(
    aocc_end_less_than_start == FALSE,
    aocc_between_cpb == TRUE,
    abs(iv_aocc_calc - iv_aocc_hr) <= 5/60   
  ) 

dim(has_aocc_filtered)

cldta_with_aocc <- cldta %>% 
  filter(surg_id %in% has_aocc_filtered$surg_id)


cldta_with_aocc <- cldta_with_aocc %>%
  mutate(
    # convert iv_evth (hours) to seconds and add to dtmbpon1
    iv_abs_time = dtmbpon1 + iv_evth * 3600,  
    
    # hours since cross-clamp start
    iv_aocc_on = as.numeric(difftime(iv_abs_time, dtmccon, units = "hours")),
    
    # hours since cross-clamp end
    iv_aocc_of = as.numeric(difftime(iv_abs_time, dtmccof, units = "hours"))
  )


cldta_while_aocc <- cldta_with_aocc %>% 
  filter(iv_aocc_on >= 0 & iv_aocc_of <= 0)

rows_per_surg_id <- cldta_while_aocc %>%
  count(surg_id, name = "n_rows")

# Histogram
ggplot(rows_per_surg_id, aes(x = n_rows)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black", boundary = 0.5) +
  scale_x_continuous(breaks = 1:max(rows_per_surg_id$n_rows)) +
  labs(
    title = "Histogram of Rows per Patient",
    x = "Number of Lactate Measurements while Cross-Clamped",
    y = "Count of surg_id"
  ) +
  theme_minimal()


# Lactate modeling ----
get_first_last_lactate_aocc_df <- function(lactate_df, iv_evth_lower = -0.5, iv_evth_upper = 0.5, 
                                      iv_evtof1_lower = -0.25, iv_evtof1_upper = 0.25) {
  cldta_first_last <- lactate_df %>%
    group_by(surg_id) %>%
    filter(n() > 1) %>%
    ungroup() %>%
    mutate(lactate_type = case_when(
      iv_aocc_on   >= iv_evth_lower & iv_aocc_on   <= iv_evth_upper   ~ "lactate_start",
      iv_aocc_of >= iv_evtof1_lower & iv_aocc_of <= iv_evtof1_upper ~ "lactate_end",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(lactate_type)) %>%
    group_by(surg_id, lactate_type) %>%
    slice_min(
      order_by = case_when(
        lactate_type == "lactate_start" ~ abs(iv_aocc_on),
        lactate_type == "lactate_end"   ~ abs(iv_aocc_of)
      ),
      with_ties = FALSE
    ) %>%
    ungroup() %>%
    group_by(surg_id) %>%
    filter(n_distinct(lactate_type) == 2) %>%
    ungroup() %>%
    # keep one copy of cpb_duration separately
    group_by(surg_id) %>%
    mutate(cpb_duration = first(cpb_duration)) %>%
    ungroup() %>%
    pivot_wider(
      id_cols = c(surg_id, cpb_duration),   # keep cpb_duration as-is
      names_from = lactate_type,
      values_from = c(lactate, iv_aocc_on, iv_aocc_of, iv_evth),
      names_glue = "{.value}_{lactate_type}"
    )
  
  print(dim(cldta_first_last)) 
  print(get_n_surg(cldta_first_last)) 
  
  return(cldta_first_last)
}


cldta_with_aocc_first_last_lactate <- get_first_last_lactate_aocc_df(cldta_with_aocc)


# Feature engineering ---- 
df_onaocc_list <- list()  # initialize new list

for (nm in names(df_onpump_list)) {
  df_current <- df_onpump_list[[nm]]
  
  # Determine which columns to join
  join_cols <- c("surg_id", "dtmccon", "dtmccof")
  if (!"dtmbpon1" %in% colnames(df_current)) {
    join_cols <- c(join_cols, "dtmbpon1")
  }
  
  df_onaocc_list[[nm]] <- df_current %>%
    left_join(
      has_aocc_filtered %>% select(all_of(join_cols)),
      by = "surg_id"
    ) %>%
    mutate(
      # convert iv_evth (hours) to seconds and add to dtmbpon1
      iv_abs_time = dtmbpon1 + iv_evth * 3600,  
      
      # hours since cross-clamp start
      iv_aocc_on = as.numeric(difftime(iv_abs_time, dtmccon, units = "hours")),
      
      # hours since cross-clamp end
      iv_aocc_of = as.numeric(difftime(iv_abs_time, dtmccof, units = "hours"))
    ) %>%
    filter(iv_aocc_on >= 0 & iv_aocc_of <= 0)
}

ggplot(df_onaocc_list[["ci"]], aes(x = ci)) +
  geom_histogram(binwidth = 0.2, fill = "steelblue", color = "black") +
  labs(
    title = "Histogram of Pump Flow while Cross Clamped",
    x = "Pump Flow",
    y = "Frequency"
  ) +
  theme_minimal()

cldta_feature_df <- cldta_with_aocc_first_last_lactate

cldta_feature_df <- cldta_feature_df %>%
  get_mean_summary(df_onaocc_list[["ci"]], "ci") %>%
  get_mean_summary(df_onaocc_list[["do2_interp"]], "do2_interp_m2") %>%
  get_mean_summary(df_onaocc_list[["map"]], "map") %>%
  get_mean_summary(df_onaocc_list[["svo2"]], "svo2") %>%
  get_mean_summary(df_onaocc_list[["glucose"]], "gluc_io") %>%
  get_mean_summary(df_onaocc_list[["hgb"]], "hgb_io") %>%
  get_mean_summary(df_onaocc_list[["pao2"]], "pao2") %>% 
  get_mean_summary(df_onaocc_list[["temp"]], "temp") %>% 
  get_mean_summary(df_onaocc_list[["ohgb"]], "ohgb_io") %>% 
  get_mean_summary(df_onaocc_list[["paco2"]], "paco2") %>% 
  get_min_summary(df_onaocc_list[["map"]], "map") %>% 
  get_sd_summary(df_onaocc_list[["map"]], "map") %>% 
  get_sd_summary(df_onaocc_list[["ci"]], "ci")

atemp_omit_first_five <- df_onaocc_list[["art_temp"]] #%>% 
  # filter(iv_evth > 5/60)
vtemp_omit_first_five <- df_onaocc_list[["ven_temp"]] #%>% 
  # filter(iv_evth > 5/60)

cldta_feature_df <- cldta_feature_df %>%
  get_mean_summary(atemp_omit_first_five, "atemp") %>%
  get_mean_summary(vtemp_omit_first_five, "vtemp") 


cldta_feature_df <- cldta_feature_df %>%
  get_time_summary(df_onaocc_list[["do2_interp"]], "do2_interp", cutoff = 500, direction = "below") %>% 
  get_time_summary(df_onaocc_list[["do2_interp"]], "do2_interp", cutoff = 800, direction = "above") %>%
  get_time_summary(df_onaocc_list[["map"]], "map", cutoff = 70, direction = "below") %>%   
  get_time_summary(df_onaocc_list[["map"]], "map", cutoff = 60, direction = "below") %>% 
  get_time_summary(df_onaocc_list[["map"]], "map", cutoff = 85, direction = "above") %>% 
  get_time_summary(df_onaocc_list[["svo2"]], "svo2", cutoff = 75, direction = "below") %>% 
  get_time_summary(df_onaocc_list[["svo2"]], "svo2", cutoff = 85, direction = "above")

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

# Modeling ---- 

formula <- as.formula(
  lactate_lactate_end ~ lactate_lactate_start +
    iv_aocc_hr + 
    # cpb_duration +
    # aocc_prop + 
    # iv_aocc_hr + 
    do2_interp_m2_mean + 
    # ci_response_slope + corr + ci_response_prop + ci_sd + 
    # ci_mean + 
    # hgb_io_mean + 
    # pao2_mean + 
    map_mean +
    # map_time_below_70 +
    # map_time_above_85 +
    # map_time_below_70_frac +
    # map_time_below_60_frac +
    # map_time_above_85_frac +
    svo2_mean +
    # svo2_time_below_75 + 
    # svo2_time_above_85 +
    # svo2_time_below_75_frac +
    # svo2_time_above_85_frac +
    # phen_dose_mcg_kg +
    # epi_dose_mcg_kg +
    # nepi_dose_mcg_kg +
    # vaso_dose_mcg_kg +
    # milr_dose_mcg_kg +
    # temp_mean +
    # atemp_mean + 
    vtemp_mean + 
    ohgb_io_mean +
    paco2_mean + 
    std_dose_mg_hr + 
    # phen_dose_mcg_kg_min +
    epi_dose_mg_hr +
    # nepi_dose_mcg_kg_min +
    # vaso_dose_mcg_kg_min +
    milr_dose_mg_hr +
    nitp_dose_mg_hr + 
    precpb_phen_dose_mcg_kg +
    precpb_epi_dose_mcg_kg +
    precpb_nepi_dose_mcg_kg +
    precpb_vaso_dose_mcg_kg +
    precpb_nitp_dose_mcg_kg + 
    # precpb_milr_dose_mcg_kg +
    precpb_map_mean + 
    # precpb_hgb_io_mean +
    precpb_pao2_mean + 
    precpb_gluc_io_mean + 
    gluc_io_mean +
    bmi + age + bsa + 
    wbc_pr + hgb_pr + creat_pr + lvef_pr + bun_pr + pot_pr + 
    tbill_pr + sp_cabg + hx_dm + robot  
  # +   map_min + map_sd + ci_sd
  # meld +
  # # cluster     diabetes_preop_status
  # # gfr_pr  blrbn_pr
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
  # colsample_bynode = sqrt(dim(model_data)[2]) / dim(model_data)[2]
  colsample_bytree = 0.8
)

xgboost_results <- k_fold_xgboost(X, y, k, params = params, get_partial = TRUE, get_ale = FALSE,
                                  quantile_points = 50) 
print_rmse_mean_ci(xgboost_results$performance_list)
partial_list <- xgboost_results$partial_list
model_list <- xgboost_results$model_list
folds <- xgboost_results$folds

# Performance ----
plot_actual_vs_pred(X, y, model_list, folds, 
                    performance_list = xgboost_results$performance_list) + 
  coord_cartesian(xlim = c(0, 12), ylim = c(0, 12))


# VIMP ----
vimp_res <- plot_vimp(model_list, metric = "Gain")
plot_vimp(model_list, metric = "Gain", top_n = 50)
plot_vimp(model_list, metric = "Cover", top_n = 30)
plot_vimp(model_list, metric = "Frequency", top_n = 10)


# Partial Plots ----
pps <- plot_pdp(xgboost_results, X, numeric_var_order = vimp_res$var_order)
pps$numeric 
pps$factor

plot_pdp(xgboost_results, X, variables = vimp_res$var_order[1:25])


plot_pdp(xgboost_results, X, variables = vimp_res$var_order[3:11])

p1 <- plot_pdp(xgboost_results, X, variables = c("gluc_io_mean", "precpb_gluc_io_mean"))$numeric
p2 <- plot_pdp(xgboost_results, X, variables = c("hx_dm"))$factor
p1 | p2
y1 <- ggplot_build(p1)$layout$panel_params[[1]]$y.range
y2 <- ggplot_build(p2)$layout$panel_params[[1]]$y.range
shared_y <- range(c(y1, y2))
p1_fixed <- p1 + ylim(shared_y)
p2_fixed <- p2 + ylim(shared_y)
p1_fixed | p2_fixed



# End AOCC Lactate vs End CPB Lactate ----
cldta_first_last <- get_first_last_lactate_df()


lactate_compare <- cldta_with_aocc_first_last_lactate %>%
  select(surg_id, lactate_crossclamp_end = lactate_lactate_end, iv_evth_lactate_end, cpb_duration) %>%
  inner_join(
    cldta_first_last %>% select(surg_id, lactate_cpb_end = lactate_lactate_end, iv_evth_lactate_end_cpb = iv_evth_lactate_end),
    by = "surg_id"
  ) %>%
  mutate(is_same_lactate = factor(iv_evth_lactate_end == iv_evth_lactate_end_cpb)) %>%
  left_join(has_aocc_filtered %>% select(surg_id, iv_aocc_hr), by = "surg_id") %>% 
  mutate(cpb_duration_not_cross_clamped = cpb_duration - iv_aocc_hr)


dim(lactate_compare)

p1 <- ggplot(lactate_compare, aes(x = lactate_crossclamp_end, y = lactate_cpb_end, 
                                  color = iv_aocc_hr, shape = is_same_lactate)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue", size = 0.5) +
  labs(
    title = "End Cross-Clamp vs End CPB Lactate",
    x = "Lactate at end of cross-clamp",
    y = "Lactate at end of CPB",
    color = "Cross Clamp Duration (hours)",
    shape = "Is Same Lactate Value?"      # <<< ADD SHAPE LABEL
  ) +
  scale_shape_manual(values = c(`TRUE` = 25, `FALSE` = 16)) +
  theme_minimal() +
  coord_fixed(ratio = 1) + 
  scale_color_viridis_c(option = "plasma") +
  theme(legend.position = "bottom")


p2 <- ggplot(lactate_compare, aes(x = lactate_crossclamp_end, y = lactate_cpb_end, 
                                  color = cpb_duration_not_cross_clamped, shape = is_same_lactate)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue", size = 0.5) +
  labs(
    title = "End Cross-Clamp vs End CPB Lactate",
    x = "Lactate at end of cross-clamp",
    y = "Lactate at end of CPB",
    color = "CPB duration not cross clamped (hours)",
    shape = "Is Same Lactate Value?"      # <<< ADD SHAPE LABEL
  ) +
  scale_shape_manual(values = c(`TRUE` = 25, `FALSE` = 16)) +
  theme_minimal() +
  coord_fixed(ratio = 1) + 
  scale_color_viridis_c(option = "plasma") +
  theme(legend.position = "bottom")

p1 | p2


# Better understand AOCC vs CPB duration ----

short_duration <- cldta_onpump_baseline_df %>% 
  filter(cpb_duration < 1.75)
mid_duration <- cldta_onpump_baseline_df %>% 
  filter(cpb_duration >= 1.75 & cpb_duration < 2.1)
long_duration <- cldta_onpump_baseline_df %>% 
  filter(cpb_duration >= 2.1)

get_n_surg(short_duration)
get_n_surg(mid_duration)
get_n_surg(long_duration)

cldta_onpump_baseline_df <- cldta_onpump_baseline_df %>% 
  mutate(n_components = sp_afib + sp_aotic + spu_asdv + sp_asdpf + sp_atrmy + sp_avrpl +
           sp_avrpr + sp_cabg + sp_cend + sp_dscao + sp_lv +
           sp_mvrpl + sp_mvrpr + sp_peric + sp_pvrpl + sp_pvrpr +
           sp_sepmy + sp_tmr + sp_tvrpl + sp_tvrpr + sp_tx) 



combined_df <- cldta_onpump_baseline_df %>%
  mutate(duration_group = case_when(
    cpb_duration < 1.75             ~ "Short",
    cpb_duration >= 1.75 & cpb_duration < 2.1 ~ "Mid",
    cpb_duration >= 2.1             ~ "Long"
  )) %>%
  mutate(duration_group = factor(duration_group, levels = c("Short", "Mid", "Long"))) %>%
  select(duration_group, all_of(names(sp_labels)), n_components)

combined_df %>%
  tbl_summary(
    by = duration_group,                
    statistic = all_categorical() ~ "{n} ({p}%)",  
    label = sp_labels,
    missing = "no"                
  ) %>%
  add_p() %>%       # <-- adds p-value column
  as_gt()



