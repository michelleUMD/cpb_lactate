source("/studies/cardiac/support/cpb/adequate/analyses/fangm/utils.R")
source("/studies/cardiac/support/cpb/adequate/analyses/fangm/xgboost/utils_xgboost.R")
source("/studies/cardiac/support/cpb/adequate/analyses/fangm/xgboost/utils_3d.R")
pacman::p_load(ggeffects, gratia, plotly, lme4, mgcv)

# Data set up ----
feature_labels_to_english <- c(
  map_delta_adv = "ΔMAP (MAPₜ+₁ − MAPₜ)",
  pred = "ΔMAP (MAPₜ+₁ − MAPₜ)",
  map_lead = "MAPₜ+₁(mmHg)ₜ",
  map = "MAPₜ (mmHg)ₜ",
  map_lag   = "MAPₜ₋₁ (mmHg)",
  map_lag2   = "MAPₜ₋₂ (mmHg)",
  tpr = "TPRₜ", 
  tpr_lag   = "TPRₜ₋₁",
  tpr_lag3  = "TPRₜ₋₃",
  tpr_lag2  = "TPRₜ₋₂",
  ci_delta  = "ΔCI (CIₜ − CIₜ₋₁)",
  tpr_delta  = "TPR (TPRₜ − TPRₜ₋₁)",
  ci_lag    = "CIₜ₋₁"
)


map_feature_df <- df_tpr_interp_with_cross_clamp %>%
  arrange(surg_id, stretch_id, iv_evth) %>%
  group_by(surg_id, stretch_id) %>%     
  mutate(
    ci_delta  = ci  - lag(ci),
    map_delta = map - lag(map),
    
    map_delta_adv3 = lead(map, 3) - lead(map, 2),
    map_delta_adv2 = lead(map, 2) - lead(map),
    map_delta_adv = lead(map) - map,
    map_delta_lag = lag(map)  - lag(map, 2),
    map_lead = lead(map),
    map_lead2 = lead(map, 2),
    map_lead3 = lead(map, 3),
    
    ci_lead  = lead(ci),
    ci_lead2  = lead(ci, 2), 
    ci_lead3  = lead(ci, 3),
    ci_lag  = lag(ci),
    ci_lag2  = lag(ci, 2),
    ci_lag3  = lag(ci, 3),
    
    ci_delta_lag  = lag(ci)  - lag(ci, 2),
    ci_delta_lag2  = lag(ci, 2)  - lag(ci, 3),
    ci_delta_lag3  = lag(ci, 3)  - lag(ci, 4),
    tpr_delta = tpr - lag(tpr),
    tpr_delta_lag = lag(tpr)  - lag(tpr, 2),
    tpr_delta_lag2 = lag(tpr, 2)  - lag(tpr, 3),
    
    map_lag = lag(map),
    map_lag2 = lag(map, 2),
    
    
    tpr_lag  = lag(tpr),
    tpr_lag2 = lag(tpr, 2),
    tpr_lag3 = lag(tpr, 3)
  ) %>%
  ungroup()

dim(map_feature_df)
summary(map_feature_df)
get_n_surg(map_feature_df)


# Main figure  ----
map_feature_df_set_ci_delta <- map_feature_df %>% 
  filter(ci_delta >= 0.12) 
# plot_map_and_delta(map_feature_df_set_ci_delta, custom_title = "ΔCI ≥ 0.12")
map_ci_long <- get_df_long(map_feature_df_set_ci_delta) 
# compare_spline_models(map_ci_long)
p1 <- get_model_plot(map_ci_long, title = "ΔPump Flow Index ≥ +0.12 L·min⁻¹·m⁻²", 
                     ci_model_split_map_delta  = FALSE,
                     map_model_split_map_delta = TRUE,
                     map_ylim = c(62, 77))
p1
map_feature_df_set_ci_delta <- map_feature_df %>% 
  filter(ci_delta >= 0.09 & ci_delta <= 0.11) 
dim(map_feature_df_set_ci_delta)
# plot_map_and_delta(map_feature_df_set_ci_delta, custom_title = "ΔCI ~ +0.1")
map_ci_long <- get_df_long(map_feature_df_set_ci_delta) 
# compare_spline_models(map_ci_long)
p2 <- get_model_plot(map_ci_long, title = "ΔPump Flow Index +0.09 - 0.11 L·min⁻¹·m⁻²", 
                     ci_model_split_map_delta  = FALSE,
                     map_model_split_map_delta = TRUE,
                     map_ylim = c(62, 77)) 
p2

map_feature_df_set_ci_delta <- map_feature_df %>% 
  filter(ci_delta >= 0.04 & ci_delta <= 0.06) 
# plot_map_and_delta(map_feature_df_set_ci_delta, custom_title = "ΔCI ~ +0.05")
map_ci_long <- get_df_long(map_feature_df_set_ci_delta) 
# compare_spline_models(map_ci_long)
p3 <- get_model_plot(map_ci_long, title = "ΔPump Flow Index +0.04 - 0.06 L·min⁻¹·m⁻²", 
                     ci_model_split_map_delta  = FALSE,
                     map_model_split_map_delta = TRUE,
                     map_ylim = c(62, 77)) 
p3 | p2 | p1 
p4 <- p3 | p2 | p1 
save("map_ci_desc_knots.png", p4, 15, 10)

# Supplemental figure: stratify by CI values ----
## < 2 ----
map_feature_df_set_ci <- map_feature_df %>% 
  filter(ci_lag < 2.0)
dim(map_feature_df_set_ci)

map_feature_df_set_ci_delta <- map_feature_df_set_ci %>% 
  filter(ci_delta >= 0.12) 
# plot_map_and_delta(map_feature_df_set_ci_delta, custom_title = "ΔCI ≥ 0.12")
map_ci_long <- get_df_long(map_feature_df_set_ci_delta) 
# compare_spline_models(map_ci_long)
p1 <- get_model_plot(map_ci_long, title = "ΔCI ≥ 0.12", 
                     ci_model_split_map_delta  = FALSE,
                     map_model_split_map_delta = TRUE,
                     ci_ylim = c(1.85, 2.2),
                     map_ylim = c(62, 77))

map_feature_df_set_ci_delta <- map_feature_df_set_ci %>% 
  filter(ci_delta >= 0.09 & ci_delta <= 0.11) 
dim(map_feature_df_set_ci_delta)
# plot_map_and_delta(map_feature_df_set_ci_delta, custom_title = "ΔCI ~ +0.1")
map_ci_long <- get_df_long(map_feature_df_set_ci_delta) 
# compare_spline_models(map_ci_long)
p2 <- get_model_plot(map_ci_long, title = "ΔCI ~ +0.09 - 0.11", 
                     ci_model_split_map_delta  = FALSE,
                     map_model_split_map_delta = TRUE,
                     ci_ylim = c(1.85, 2.2),
                     map_ylim = c(62, 77)) 

map_feature_df_set_ci_delta <- map_feature_df_set_ci %>% 
  filter(ci_delta >= 0.04 & ci_delta <= 0.06) 
# plot_map_and_delta(map_feature_df_set_ci_delta, custom_title = "ΔCI ~ +0.05")
map_ci_long <- get_df_long(map_feature_df_set_ci_delta) 
# compare_spline_models(map_ci_long)
p3 <- get_model_plot(map_ci_long, title = "ΔCI ~ +0.04 - 0.06", 
                     ci_model_split_map_delta  = FALSE,
                     map_model_split_map_delta = TRUE,
                     ci_ylim = c(1.85, 2.2),
                     map_ylim = c(62, 77)) 
(p3 | p2 | p1) + 
  plot_annotation(
    title = "CI < 2 L·min⁻¹·m⁻²"
  )

## 2 - 2.2 ----
map_feature_df_set_ci <- map_feature_df %>% 
  filter(ci_lag >= 2.0 & ci_lag < 2.2)
dim(map_feature_df_set_ci)

map_feature_df_set_ci_delta <- map_feature_df_set_ci %>% 
  filter(ci_delta >= 0.12) 
# plot_map_and_delta(map_feature_df_set_ci_delta, custom_title = "ΔCI ≥ 0.12")
map_ci_long <- get_df_long(map_feature_df_set_ci_delta) 
# compare_spline_models(map_ci_long)
p1 <- get_model_plot(map_ci_long, title = "ΔCI ≥ 0.12", 
                     ci_model_split_map_delta  = FALSE,
                     map_model_split_map_delta = TRUE,
                     ci_ylim = c(1.95, 2.35),
                     map_ylim = c(62, 77))

map_feature_df_set_ci_delta <- map_feature_df_set_ci %>% 
  filter(ci_delta >= 0.09 & ci_delta <= 0.11) 
dim(map_feature_df_set_ci_delta)
# plot_map_and_delta(map_feature_df_set_ci_delta, custom_title = "ΔCI ~ +0.1")
map_ci_long <- get_df_long(map_feature_df_set_ci_delta) 
# compare_spline_models(map_ci_long)
p2 <- get_model_plot(map_ci_long, title = "ΔCI ~ +0.09 - 0.11",  
                     ci_model_split_map_delta  = FALSE,
                     map_model_split_map_delta = TRUE,
                     ci_ylim = c(1.95, 2.35),
                     map_ylim = c(62, 77)) 

map_feature_df_set_ci_delta <- map_feature_df_set_ci %>% 
  filter(ci_delta >= 0.04 & ci_delta <= 0.06) 
# plot_map_and_delta(map_feature_df_set_ci_delta, custom_title = "ΔCI ~ +0.05")
map_ci_long <- get_df_long(map_feature_df_set_ci_delta) 
# compare_spline_models(map_ci_long)
p3 <- get_model_plot(map_ci_long, title = "ΔCI ~ +0.04 - 0.06", 
                     ci_model_split_map_delta  = FALSE,
                     map_model_split_map_delta = TRUE,
                     ci_ylim = c(1.95, 2.35),
                     map_ylim = c(62, 77)) 
(p3 | p2 | p1) + 
  plot_annotation(
    title = "CI [2.0 - 2.2) L·min⁻¹·m⁻²"
  )

## 2.2 - 2.4 ----
map_feature_df_set_ci <- map_feature_df %>% 
  filter(ci_lag >= 2.2 & ci_lag < 2.4)
dim(map_feature_df_set_ci)

map_feature_df_set_ci_delta <- map_feature_df_set_ci %>% 
  filter(ci_delta >= 0.12) 
# plot_map_and_delta(map_feature_df_set_ci_delta, custom_title = "ΔCI ≥ 0.12")
map_ci_long <- get_df_long(map_feature_df_set_ci_delta) 
# compare_spline_models(map_ci_long)
p1 <- get_model_plot(map_ci_long, title = "ΔCI ≥ 0.12", 
                     ci_model_split_map_delta  = FALSE,
                     map_model_split_map_delta = TRUE,
                     ci_ylim = c(2.25, 2.6),
                     map_ylim = c(62, 77))

map_feature_df_set_ci_delta <- map_feature_df_set_ci %>% 
  filter(ci_delta >= 0.09 & ci_delta <= 0.11) 
dim(map_feature_df_set_ci_delta)
# plot_map_and_delta(map_feature_df_set_ci_delta, custom_title = "ΔCI ~ +0.1")
map_ci_long <- get_df_long(map_feature_df_set_ci_delta) 
# compare_spline_models(map_ci_long)
p2 <- get_model_plot(map_ci_long, title = "ΔCI ~ +0.09 - 0.11", 
                     ci_model_split_map_delta  = FALSE,
                     map_model_split_map_delta = TRUE,
                     ci_ylim = c(2.25, 2.6),
                     map_ylim = c(62, 77)) 

map_feature_df_set_ci_delta <- map_feature_df_set_ci %>% 
  filter(ci_delta >= 0.04 & ci_delta <= 0.06) 
# plot_map_and_delta(map_feature_df_set_ci_delta, custom_title = "ΔCI ~ +0.05")
map_ci_long <- get_df_long(map_feature_df_set_ci_delta) 
# compare_spline_models(map_ci_long)
p3 <- get_model_plot(map_ci_long, title = "ΔCI ~ +0.04 - 0.06", 
                     ci_model_split_map_delta  = FALSE,
                     map_model_split_map_delta = TRUE,
                     ci_ylim = c(2.25, 2.6),
                     map_ylim = c(62, 77)) 
(p3 | p2 | p1) + 
  plot_annotation(
    title = "CI [2.2, 2.4) L·min⁻¹·m⁻²"
  )


## >= 2.4 ----
map_feature_df_set_ci <- map_feature_df %>% 
  filter(ci_lag >= 2.4)
dim(map_feature_df_set_ci)

map_feature_df_set_ci_delta <- map_feature_df_set_ci %>% 
  filter(ci_delta >= 0.12) 
# plot_map_and_delta(map_feature_df_set_ci_delta, custom_title = "ΔCI ≥ 0.12")
map_ci_long <- get_df_long(map_feature_df_set_ci_delta) 
# compare_spline_models(map_ci_long)
p1 <- get_model_plot(map_ci_long, title = "ΔCI ≥ 0.12", 
                     ci_model_split_map_delta  = FALSE,
                     map_model_split_map_delta = TRUE,
                     ci_ylim = c(2.4, 2.8),
                     map_ylim = c(62, 77))

map_feature_df_set_ci_delta <- map_feature_df_set_ci %>% 
  filter(ci_delta >= 0.09 & ci_delta <= 0.11) 
dim(map_feature_df_set_ci_delta)
# plot_map_and_delta(map_feature_df_set_ci_delta, custom_title = "ΔCI ~ +0.1")
map_ci_long <- get_df_long(map_feature_df_set_ci_delta) 
# compare_spline_models(map_ci_long)
p2 <- get_model_plot(map_ci_long, title = "ΔCI ~ +0.09 - 0.11", 
                     ci_model_split_map_delta  = FALSE,
                     map_model_split_map_delta = TRUE,
                     ci_ylim = c(2.4, 2.8),
                     map_ylim = c(62, 77)) 

map_feature_df_set_ci_delta <- map_feature_df_set_ci %>% 
  filter(ci_delta >= 0.04 & ci_delta <= 0.06) 
# plot_map_and_delta(map_feature_df_set_ci_delta, custom_title = "ΔCI ~ +0.05")
map_ci_long <- get_df_long(map_feature_df_set_ci_delta) 
# compare_spline_models(map_ci_long)
p3 <- get_model_plot(map_ci_long, title = "ΔCI ~ +0.04 - 0.06", 
                     ci_model_split_map_delta  = FALSE,
                     map_model_split_map_delta = TRUE,
                     ci_ylim = c(2.4, 2.8),
                     map_ylim = c(62, 77)) 
(p3 | p2 | p1) + 
  plot_annotation(
    title = "CI ≥ 2.4 L·min⁻¹·m⁻²"
  )