pacman::p_load(ggh4x, splines)

# Plot violin boxplots of MAP vs time 
plot_map_and_delta <- function(df, custom_title = "") {
  
  # --- MAP long format ---
  map_plot_df <- df %>%
    select(
      map_lag2,
      map_lag,
      map,
      map_lead,
      map_lead2,
      map_lead3
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = "map_time",
      values_to = "map_value"
    ) %>% 
    mutate(
      map_time = factor(
        map_time,
        levels = c(
          "map_lag2",
          "map_lag",
          "map",
          "map_lead",
          "map_lead2",
          "map_lead3"
        )
      )
    )
  
  map_n_df <- map_plot_df %>%
    group_by(map_time) %>%
    summarise(
      n = sum(!is.na(map_value)),
      y = median(map_value, na.rm = TRUE)
    )
  
  map_plot <- ggplot(map_plot_df, aes(x = map_time, y = map_value)) +
    geom_violin(trim = TRUE, alpha = 0.6) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
    geom_text(
      data = map_n_df,
      aes(x = map_time, y = y + 5, label = paste0("N=", n)),
      inherit.aes = FALSE,
      size = 3.5
    ) +
    theme_bw() +
    labs(
      x = "Time of MAP",
      y = "MAP",
      title = "MAP over time"
    ) +
    coord_cartesian(ylim = c(60, 80))
  
  # --- ΔMAP long format ---
  delta_plot_df <- df %>%
    select(
      map_delta_lag,
      map_delta,
      map_delta_adv,
      map_delta_adv2,
      map_delta_adv3
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = "map_time",
      values_to = "delta_value"
    ) %>% 
    mutate(
      map_time = factor(
        map_time,
        levels = c(
          "map_delta_lag",
          "map_delta",
          "map_delta_adv",
          "map_delta_adv2",
          "map_delta_adv3"
        )
      )
    )
  
  delta_n_df <- delta_plot_df %>%
    group_by(map_time) %>%
    summarise(
      n = sum(!is.na(delta_value)),
      y = median(delta_value, na.rm = TRUE)
    )
  
  delta_plot <- ggplot(delta_plot_df, aes(x = map_time, y = delta_value)) +
    geom_violin(trim = TRUE, alpha = 0.6) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
    geom_text(
      data = delta_n_df,
      aes(x = map_time, y = y + 1, label = paste0("N=", n)),
      inherit.aes = FALSE,
      size = 3.5
    ) +
    theme_bw() +
    labs(
      x = "Time of MAP delta",
      y = "ΔMAP",
      title = "MAP Δ over time"
    ) +
    coord_cartesian(ylim = c(-10, 10))
  
  
  combined_plot <- map_plot / delta_plot +
    plot_annotation(
      title = custom_title
    )
  
  return(combined_plot)
}


# Compare spline numbers between mixed effect models 
compare_spline_models <- function(data, outcomes = c("map", "ci"), spline_dfs = 3:6) {
  
  # helper function to fit a single model
  fit_spline_model <- function(outcome, df, data) {
    formula <- as.formula(paste0(outcome,
                                 " ~ ns(time_offset, df = ", df, 
                                 ") + (1 | surg_id) + (1 | surg_id:iv_evth)"))
    lmer(formula, data = data, REML = FALSE)
  }
  
  # create a plain list of models
  model_list <- pmap(
    list(outcome = rep(outcomes, each = length(spline_dfs)),
         df = rep(spline_dfs, times = length(outcomes))),
    fit_spline_model,
    data = data
  )
  
  # store in tibble
  models <- tibble(
    outcome = rep(outcomes, each = length(spline_dfs)),
    df = rep(spline_dfs, times = length(outcomes)),
    model = model_list
  )
  
  # extract numeric stats
  models <- models %>%
    mutate(
      logLik = map_dbl(model, ~ as.numeric(logLik(.x))),
      AIC    = map_dbl(model, ~ as.numeric(AIC(.x))),
      BIC    = map_dbl(model, ~ as.numeric(BIC(.x)))
    )
  
  # run LRTs for each outcome
  run_lrt <- function(models_df) {
    n <- nrow(models_df)
    pvals <- rep(NA_real_, n)
    chisq <- rep(NA_real_, n)
    df_diff <- rep(NA_real_, n)
    
    for(i in 2:n) {
      an <- anova(models_df$model[[i-1]], models_df$model[[i]])
      chisq[i] <- an$Chisq[2]
      df_diff[i] <- an$`Df`[2]
      pvals[i] <- an$`Pr(>Chisq)`[2]
    }
    
    models_df %>%
      mutate(chisq = chisq,
             df_diff = df_diff,
             p_value = pvals)
  }
  
  model_stats_lrt <- models %>%
    group_by(outcome) %>%
    arrange(df) %>%
    group_modify(~ run_lrt(.x)) %>%
    ungroup()
  
  return(model_stats_lrt)
}

get_df_long <- function(map_feature_df_set_ci_delta) {
  map_ci_long <- map_feature_df_set_ci_delta %>%
    select(
      surg_id,
      iv_evth, map_delta,
      map_lag2, map_lag, map, map_lead, map_lead2, map_lead3,
      ci_lag2,  ci_lag,  ci,  ci_lead,  ci_lead2,  ci_lead3
    ) %>%
    mutate(mdelta = map_delta) %>% 
    select(-map_delta) %>% 
    pivot_longer(
      cols = matches("^(map|ci)"),
      names_to = c(".value", "map_version"),
      names_pattern = "^(map|ci)(?:_(.*))?$"
    ) %>%
    mutate(
      map_version = replace_na(map_version, "0"),
      time_offset = case_when(
        map_version == "lag2"  ~ -2,
        map_version == "lag"   ~ -1,
        map_version == ""     ~  0,
        map_version == "lead"  ~  1,
        map_version == "lead2" ~  2,
        map_version == "lead3" ~  3
      ),
      surg_id_iv_evth = paste(surg_id, iv_evth)
    ) %>% 
    mutate(map_delta = mdelta) 
  
  return(map_ci_long)
}

# Fit final mixed effects models and plot 
get_model_plot <- function(
    data,
    spline_df = 5,
    time_range = c(-2, 3),
    title = "",
    ci_ylim  = c(1.95, 2.3),
    map_ylim = c(65, 75),
    ci_model_split_map_delta  = FALSE,
    map_model_split_map_delta = FALSE,
    show_knots = FALSE,
    add_delta_annotation = TRUE   # new parameter
) {
  
  data_pos <- filter(data, map_delta >= 0)
  data_neg <- filter(data, map_delta < 0)
  
  time_grid <- tibble(time_offset = seq(time_range[1], time_range[2], by = 0.05))
  eval_times <- c(-1, 0, 1, 2)
  
  # MAP ----
  if (!map_model_split_map_delta) {
    map_model <- lmer(
      map ~ ns(time_offset, df = spline_df) + (1 | surg_id) + (1 | surg_id:iv_evth),
      data = data, REML = FALSE,
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    )
    map_pred <- time_grid %>%
      mutate(fitted = predict(map_model, newdata = ., re.form = NA),
             group = "All",
             outcome = "MAP (mmHg)",
             alpha = nrow(data),
             linetype = "combined")
    
  } else {
    n_pos <- nrow(data_pos)
    n_neg <- nrow(data_neg)
    
    map_model_all <- lmer(
      map ~ ns(time_offset, df = spline_df) + (1 | surg_id) + (1 | surg_id:iv_evth),
      data = data, REML = FALSE,
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    )
    map_model_pos <- lmer(
      map ~ ns(time_offset, df = spline_df) + (1 | surg_id) + (1 | surg_id:iv_evth),
      data = data_pos, REML = FALSE,
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    )
    map_model_neg <- lmer(
      map ~ ns(time_offset, df = spline_df) + (1 | surg_id) + (1 | surg_id:iv_evth),
      data = data_neg, REML = FALSE,
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    )
    
    map_pred <- bind_rows(
      time_grid %>% mutate(fitted = predict(map_model_all, newdata = ., re.form = NA),
                           group = "All", alpha = nrow(data), linetype = "combined"),
      time_grid %>% mutate(fitted = predict(map_model_pos, newdata = ., re.form = NA),
                           group = paste0("ΔMAP ≥ 0 (N=", n_pos, ")"), alpha = n_pos, linetype = "split"),
      time_grid %>% mutate(fitted = predict(map_model_neg, newdata = ., re.form = NA),
                           group = paste0("ΔMAP < 0 (N=", n_neg, ")"), alpha = n_neg, linetype = "split")
    ) %>% mutate(outcome = "MAP (mmHg)")
  }
  
  # CI ----
  if (!ci_model_split_map_delta) {
    ci_model <- lmer(
      ci ~ ns(time_offset, df = spline_df) + (1 | surg_id) + (1 | surg_id:iv_evth),
      data = data, REML = FALSE,
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    )
    ci_pred <- time_grid %>%
      mutate(fitted = predict(ci_model, newdata = ., re.form = NA),
             group = "All",
             outcome = "CI (L·min⁻¹·m⁻²)",
             alpha = nrow(data),
             linetype = "combined")
    
  } else {
    n_pos <- nrow(data_pos)
    n_neg <- nrow(data_neg)
    
    ci_model_pos <- lmer(
      ci ~ ns(time_offset, df = spline_df) + (1 | surg_id) + (1 | surg_id:iv_evth),
      data = data_pos, REML = FALSE,
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    )
    ci_model_neg <- lmer(
      ci ~ ns(time_offset, df = spline_df) + (1 | surg_id) + (1 | surg_id:iv_evth),
      data = data_neg, REML = FALSE,
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    )
    
    ci_pred <- bind_rows(
      time_grid %>% mutate(fitted = predict(ci_model_pos, newdata = ., re.form = NA),
                           group = paste0("ΔMAP ≥ 0 (N=", n_pos, ")"), alpha = n_pos, linetype = "combined"),
      time_grid %>% mutate(fitted = predict(ci_model_neg, newdata = ., re.form = NA),
                           group = paste0("ΔMAP < 0 (N=", n_neg, ")"), alpha = n_neg, linetype = "combined")
    ) %>% mutate(outcome = "CI (L·min⁻¹·m⁻²)")
  }
  
  plot_df <- bind_rows(map_pred, ci_pred) %>%
    mutate(outcome = factor(outcome, levels = c("CI (L·min⁻¹·m⁻²)", "MAP (mmHg)")))
  
  if (map_model_split_map_delta && add_delta_annotation) {
    map_eval <- plot_df %>%
      filter(outcome == "MAP (mmHg)", group != "All", time_offset %in% eval_times) %>%
      group_by(group, time_offset) %>%
      summarise(fitted = mean(fitted), .groups = "drop")
    
    delta_df <- map_eval %>%
      group_by(group) %>%
      group_modify(~{
        df <- arrange(.x, time_offset)
        n <- nrow(df)
        if (n < 2) return(tibble())
        
        comb <- expand.grid(i = seq_len(n), j = seq_len(n)) %>%
          filter(i < j) %>%
          mutate(t1 = df$time_offset[i], t2 = df$time_offset[j],
                 f1 = df$fitted[i], f2 = df$fitted[j],
                 delta = f2 - f1) %>%
          filter(delta > 0)
        if (nrow(comb) == 0) return(tibble())
        comb[which.max(comb$delta), ]
      }) %>% ungroup()
    
    delta_points <- bind_rows(
      delta_df %>% transmute(group, time_offset = t1, fitted = f1, outcome = "MAP (mmHg)"),
      delta_df %>% transmute(group, time_offset = t2, fitted = f2, outcome = "MAP (mmHg)")
    )
    
    delta_labels <- delta_df %>%
      transmute(group, time_offset = t2, fitted = f2,
                label = paste0("Δ = +", round(delta, 1), " mmHg"),
                outcome = "MAP (mmHg)") %>%
      mutate(nudge_y = case_when(
        grepl("ΔMAP ≥ 0", group) ~ 2,
        grepl("ΔMAP < 0", group)  ~ -2,
        TRUE ~ 0
      ))
  }
  
  group_colors <- c("All" = "black")
  split_groups <- unique(plot_df$group[plot_df$group != "All"])
  split_colors <- scales::hue_pal()(length(split_groups))
  names(split_colors) <- split_groups
  all_colors <- c(group_colors, split_colors)
  
  p <- ggplot(plot_df, aes(time_offset, fitted, color = group, linetype = linetype, alpha = alpha)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = all_colors) +
    scale_linetype_manual(values = c(combined = "solid", split = "dashed"), guide = "none") +
    scale_alpha_continuous(range = c(0.4, 1), guide = "none") +
    facet_wrap(~ outcome, scales = "free_y", ncol = 1) +
    scale_x_continuous(breaks = seq(time_range[1], time_range[2], by = 1)) +
    ggh4x::facetted_pos_scales(
      y = list(
        outcome == "MAP (mmHg)" ~ scale_y_continuous(limits = map_ylim),
        outcome == "CI (L·min⁻¹·m⁻²)" ~ scale_y_continuous(limits = ci_ylim)
      )
    ) +
    labs(x = "Time offset (6-min window)", y = "Fitted value", title = title, color = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  if (map_model_split_map_delta && add_delta_annotation && nrow(delta_df) > 0) {
    p <- p +
      geom_point(data = delta_points, aes(time_offset, fitted, color = group),
                 inherit.aes = FALSE, size = 3) +
      geom_text(data = delta_labels, aes(time_offset, fitted, label = label, color = group),
                inherit.aes = FALSE, nudge_y = delta_labels$nudge_y, size = 4)
  }
  
  if (show_knots) {
    knots <- c(attr(ns(data$time_offset, df = spline_df), "Boundary.knots"),
               attr(ns(data$time_offset, df = spline_df), "knots"))
    p <- p + geom_vline(xintercept = knots, linetype = "dotted", color = "grey")
  }
  
  p
}

