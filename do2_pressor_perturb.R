source("/studies/cardiac/support/cpb/adequate/analyses/fangm/cluster/utils_cluster.R")


pred_new_lactate_compare_colored_y <- function(X, y, model_list, folds,
                                               perturb_std = TRUE, 
                                               perturb_do2 = TRUE) {
  k <- length(model_list)
  
  # Prepare datasets
  X_orig <- X
  X_perturbed <- X
  
  caption <- c()   # store caption components
  
  if (perturb_std) {
    delta_std <- sd(X$std_dose_mg_hr) / 5
    X_perturbed$std_dose_mg_hr <- pmax(0, X_perturbed$std_dose_mg_hr - delta_std)
    caption <- c(caption, paste0("Pressor dose decreased by ", round(delta_std, 2), " mg/hr"))
    cat("Subtracting", round(delta_std, 3), "from pressor dose\n")
  } else {
    X_perturbed$std_dose_mg_hr <- X$std_dose_mg_hr
  }
  
  if (perturb_do2) {
    delta_do2 <- sd(X$do2_interp_m2_mean) / 5
    X_perturbed$do2_interp_m2_mean <- X_perturbed$do2_interp_m2_mean + delta_do2
    caption <- c(caption, paste0("DO₂ increased by ", round(delta_do2, 2), " mL/min/m²"))
    cat("Adding", round(delta_do2, 3), "to DO₂\n")
  } else {
    X_perturbed$do2_interp_m2_mean <- X$do2_interp_m2_mean
  }
  
  caption <- paste(caption, collapse = ". ")
  
  X_mat_orig <- x_to_x_mat(X_orig)
  X_mat_pert <- x_to_x_mat(X_perturbed)
  
  # Collect predictions
  pred_df <- bind_rows(
    lapply(1:k, function(i) {
      valid_idx <- which(folds == i)
      fit <- model_list[[i]]
      
      pred_orig <- predict(fit, xgb.DMatrix(X_mat_orig[valid_idx, ]))
      pred_pert <- predict(fit, xgb.DMatrix(X_mat_pert[valid_idx, ]))
      
      data.frame(
        fold = factor(i),
        Actual = y[valid_idx],
        Predicted_Original = pred_orig,
        Predicted_Perturbed = pred_pert
      )
    })
  )
  
  # Plot: Original vs Perturbed, colored by true Y
  p <- ggplot(pred_df, aes(x = Predicted_Original, y = Predicted_Perturbed, color = Actual)) +
    geom_point(alpha = 0.7) +
    scale_color_viridis_c(option = "plasma") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    labs(
      x = "Original Predicted Lactate (mmol/L)",
      y = "Perturbed Predicted Lactate (mmol/L)",
      color = "True Lactate",
      title = "Predicted Lactate: Original vs Perturbed",
      subtitle = caption
    ) +
    theme_minimal()
  
  print(p)
  
  return(pred_df)
}


pred_df_do2_std_dose <- pred_new_lactate_compare_colored_y(X, y, model_list, folds)
head(pred_df_do2_std_dose)
pred_df_do2 <- pred_new_lactate_compare_colored_y(X, y, model_list, folds, perturb_std = FALSE)
pred_df_std <- pred_new_lactate_compare_colored_y(X, y, model_list, folds, perturb_do2 = FALSE)

head(pred_df_std)

pred_df_do2_std_dose$delta <- pred_df_do2_std_dose$Predicted_Perturbed -
  pred_df_do2_std_dose$Predicted_Original

pred_df_do2$delta <- pred_df_do2$Predicted_Perturbed -
  pred_df_do2$Predicted_Original

pred_df_std$delta <- pred_df_std$Predicted_Perturbed -
  pred_df_std$Predicted_Original

pred_df_do2_std_dose$scenario <- "Change Both"
pred_df_do2$scenario <- "Increase DO₂ only"
pred_df_std$scenario <- "Decrease Pressor only"

all_pred_df <- bind_rows(pred_df_do2_std_dose,
                         pred_df_do2,
                         pred_df_std)

summary_df <- all_pred_df %>%
  group_by(scenario) %>%
  summarize(
    mean_delta = mean(delta, na.rm = TRUE),
    median_delta = median(delta, na.rm = TRUE),
    sd_delta = sd(delta, na.rm = TRUE),
    # IQR_delta = IQR(delta, na.rm = TRUE),
    p15_delta = quantile(delta, 0.15, na.rm = TRUE),
    p85_delta = quantile(delta, 0.85, na.rm = TRUE)
    # min_delta = min(delta, na.rm = TRUE),
    # max_delta = max(delta, na.rm = TRUE)
  )
summary_t <- summary_df %>%
  pivot_longer(
    cols = -scenario,
    names_to = "statistic",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = scenario,
    values_from = value
  )


summary_t %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  gt() %>%
  tab_header(
    title = "Summary of Δ (Perturbed − Original) by Scenario"
  )

# Boxplot
ggplot(all_pred_df, aes(x = scenario, y = delta, fill = scenario)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Change in Predicted Lactate After Perturbation",
       x = "Perturbation Scenario",
       y = "Δ Predicted (Perturbed - Original)") +
  theme_minimal()

# Violin 
ggplot(all_pred_df, aes(x = scenario, y = delta, fill = scenario)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Change in Predicted Lactate After Perturbation",
    x = "Perturbation Scenario",
    y = "Δ Predicted (Perturbed - Original)"
  ) +
  theme_minimal()

ggplot(all_pred_df, aes(x = delta, fill = scenario)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ scenario, ncol = 1, scales = "free_y") +
  labs(title = "Density of Δ Predicted",
       x = "Δ Predicted (Perturbed - Original)",
       y = "Density", fill = "Scenario") +
  theme_minimal() +   theme(legend.position = "none")


# Dot connection plots ----
pred_new_lactate_connected <- function(X, y, model_list, folds) {
  k <- length(model_list)
  
  # Prepare datasets
  X_orig <- X
  X_perturbed <- X
  X_perturbed$std_dose_mg_hr <- pmax(0, X_perturbed$std_dose_mg_hr - (sd(X$std_dose_mg_hr) / 10))
  X_perturbed$do2_interp_m2_mean <- X_perturbed$do2_interp_m2_mean + (sd(X$do2_interp_m2_mean) / 10)
  
  X_mat_orig <- x_to_x_mat(X_orig)
  X_mat_pert <- x_to_x_mat(X_perturbed)
  
  # Collect predictions
  pred_df <- bind_rows(
    lapply(1:k, function(i) {
      valid_idx <- which(folds == i)
      fit <- model_list[[i]]
      
      pred_orig <- predict(fit, xgb.DMatrix(X_mat_orig[valid_idx, ]))
      pred_pert <- predict(fit, xgb.DMatrix(X_mat_pert[valid_idx, ]))
      
      data.frame(
        fold = factor(i),
        Sample = valid_idx,
        Actual = y[valid_idx],
        Predicted_Original = pred_orig,
        Predicted_Perturbed = pred_pert
      )
    })
  )
  
  # Reshape for plotting
  pred_long <- pred_df %>%
    pivot_longer(cols = c(Predicted_Original, Predicted_Perturbed),
                 names_to = "Scenario",
                 values_to = "Predicted") %>%
    mutate(Scenario = ifelse(Scenario == "Predicted_Original", "Original", "Perturbed"))
  
  # Connected dot plot
  p <- ggplot(pred_long, aes(x = Scenario, y = Predicted, group = Sample, color = Actual)) +
    geom_line(alpha = 0.5) +   # connect original -> perturbed
    geom_point(size = 2) +
    scale_color_viridis_c(option = "plasma") +
    labs(
      x = "",
      y = "Predicted Lactate",
      color = "True Lactate",
      title = "Predicted Lactate: Original vs Perturbed (Connected Dots)"
    ) +
    theme_minimal()
  
  print(p)
  
  return(pred_df)
}

pred_df <- pred_new_lactate_compare_colored_y(X, y, model_list, folds)
head(pred_df)

ggplot(pred_df, aes(color = Actual)) +
  geom_segment(
    aes(
      x = Predicted_Original,
      xend = Predicted_Perturbed,
      y = Actual,     # same y for both endpoints → horizontal line
      yend = Actual
    ),
    alpha = 0.4
  ) +
  geom_point(aes(x = Predicted_Original, y = Actual), alpha = 0.7) +
  geom_point(aes(x = Predicted_Perturbed, y = Actual), alpha = 0.7, shape = 17) +
  scale_color_viridis_c(option = "plasma") +
  labs(
    x = "Predicted Original",
    y = "True Lactate",
    color = "Actual Lactate",
    title = "Original vs Perturbed Predictions",
    subtitle = "Each segment connects a sample’s original and perturbed prediction"
  ) +
  theme_minimal()



a <- pred_df_do2_std_dose %>% 
  mutate(change_both_predicted = Predicted_Perturbed) %>%
  select(-delta, -scenario, -Predicted_Perturbed)
b <- pred_df_do2 %>% 
  mutate(inc_do2_predicted = Predicted_Perturbed) %>% 
  select(inc_do2_predicted)
c <- pred_df_std %>% 
  mutate(dec_pressor_predicted = Predicted_Perturbed) %>% 
  select(dec_pressor_predicted)

plot_df <- cbind(a, b, c)

plot_df_long <- plot_df %>%
  mutate(id = row_number()) %>%
  pivot_longer(
    cols = c(Predicted_Original, change_both_predicted, inc_do2_predicted, dec_pressor_predicted),
    names_to = "PredictionType",
    values_to = "PredictedValue"
  ) %>%
  mutate(
    Xcat = recode(PredictionType,
                  Predicted_Original = "Original",
                  change_both_predicted = "Change Both",
                  inc_do2_predicted = "Increase DO₂ only",
                  dec_pressor_predicted = "Decrease Pressor only"),
    Xcat = factor(Xcat, levels = c("Increase DO₂ only", 
                                   "Original", 
                                   "Decrease Pressor only",
                                   "Change Both"))
  )

ggplot(plot_df_long, aes(x = Xcat, y = PredictedValue, group = id, color = Actual)) +
  geom_line(alpha = 0.4) +
  geom_point(size = 2) +
  scale_color_viridis_c(option = "plasma") +
  labs(
    title = "Predicted Lactate: Original vs Perturbed Scenarios",
    x = "",
    y = "Predicted Lactate (mmol/L)",
    color = "Actual Lactate (mmol/L)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 



plot_df_features <- plot_df %>%
  mutate(
    do2_change     = inc_do2_predicted     - Predicted_Original,
    pressor_change = dec_pressor_predicted - Predicted_Original,
    both_change    = change_both_predicted - Predicted_Original,
    
    do2_resp     = do2_change < 0,
    pressor_resp = pressor_change < 0,
    both_resp    = both_change < 0,
    
    synergistic = do2_resp & pressor_resp &
      both_change < do2_change &
      both_change < pressor_change,
    
    response_category = factor(case_when(
      do2_resp & !pressor_resp ~ "DO2 responsive only",
      pressor_resp & !do2_resp ~ "Pressor responsive only",
      do2_resp & pressor_resp & !synergistic ~ "Dual responsive (less than one change alone)",
      synergistic ~ "Dual responsive (more than one change alone)",
      TRUE ~ "No response"
    ),       
    levels = c(
      "No response",
      "DO2 responsive only",
      "Pressor responsive only",
      "Dual responsive (less than one change alone)",
      "Dual responsive (more than one change alone)"
    ))
  ) %>%
  mutate(
    cluster = case_when(
      response_category == "No response" ~ 1,
      response_category == "DO2 responsive only" ~ 2,
      response_category == "Pressor responsive only" ~ 3,
      response_category == "Dual responsive (less than one change alone)" ~ 4,
      response_category == "Dual responsive (more than one change alone)" ~ 5,
    )
  )

summary(plot_df_features$response_category)
summary(factor(plot_df_features$cluster))


plot_df_features$surg_id <- model_surg_id 

cluster_df <- plot_df_features %>% 
  select(surg_id, cluster)

get_tables_comparing_clusters(cluster_df)


plot_df_features_with_hx <- plot_df_features %>%  left_join(
  cldta_onpump_baseline_df %>% 
    select(surg_id, hx_copd, hx_dm, hx_dmtrt), by = "surg_id") 

plot_df_features_with_hx %>%
  select(response_category, hx_copd, hx_dm, hx_dmtrt) %>%
  mutate(across(c(hx_copd, hx_dm, hx_dmtrt), as.factor)) %>%  # ensure TRUE/FALSE are factors
  tbl_summary(
    by = response_category,
    type = all_categorical() ~ "categorical",   
    missing = "no",
    statistic = all_categorical() ~ "{n} ({p}%)",
    percent = "row",                        # row-wise percentages
    label = baseline_var_to_english
  ) %>%
  add_p(
    test = everything() ~ "chisq.test"     
  )
