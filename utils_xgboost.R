pacman::p_load(fastDummies, xgboost, Matrix, dplyr, iml, plotly, htmlwidgets, ks)

# Define start and end lactate  ----
get_first_last_lactate_df <- function(iv_evth_lower = -0.5, iv_evth_upper = 0.5, 
                                      iv_evtof1_lower = -0.25, iv_evtof1_upper = 0.25) {
  cldta_first_last <- cldta %>%
    group_by(surg_id) %>%
    filter(n() > 1) %>%
    ungroup() %>%
    mutate(lactate_type = case_when(
      iv_evth   >= iv_evth_lower & iv_evth   <= iv_evth_upper   ~ "lactate_start",
      iv_evtof1 >= iv_evtof1_lower & iv_evtof1 <= iv_evtof1_upper ~ "lactate_end",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(lactate_type)) %>%
    group_by(surg_id, lactate_type) %>%
    slice_min(
      order_by = case_when(
        lactate_type == "lactate_start" ~ abs(iv_evth),
        lactate_type == "lactate_end"   ~ abs(iv_evtof1)
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
      values_from = c(lactate, iv_evth, iv_evtof1),
      names_glue = "{.value}_{lactate_type}"
    )
  
  print(dim(cldta_first_last)) 
  print(get_n_surg(cldta_first_last)) 
  
  return(cldta_first_last)
}

get_first_last_lactate_plots <- function(cldta_first_last, 
                                         facet_start = "Starting Lactate",
                                         facet_end = "Ending lactate",
                                         binwidth = 0.2) {
  
  df_start <- cldta_first_last %>%
    mutate(
      group = ifelse(iv_evth_lactate_start <= 0, "Pre-CPB", "On CPB"),
      facet = facet_start
    ) %>%
    select(lactate = lactate_lactate_end, group, facet)

  df_end <- cldta_first_last %>%
    mutate(
      group = ifelse(iv_evtof1_lactate_end <= 0, "On CPB", "Post-CPB"),
      facet = facet_end
    ) %>%
    select(lactate = lactate_lactate_end, group, facet)
  
  df_all <- bind_rows(df_start, df_end) %>%
    mutate(
      group = factor(group, levels = c("Pre-CPB", "On CPB", "Post-CPB")),
      facet = factor(facet, levels = c(facet_start, facet_end))
    )
  
  p <- ggplot(df_all, aes(x = lactate, fill = group)) +
    geom_histogram(alpha = 0.4, position = "identity", binwidth = binwidth) +
    facet_wrap(~facet, ncol = 2, scales = "fixed") +
    scale_fill_manual(values = c(
      "Pre-CPB" = "#1B9E77",   # dark teal
      "On CPB"  = "#D95F02",   # rust orange
      "Post-CPB"= "#7570B3"    # indigo
    )) + 
    labs(
      title = paste0("Starting and Ending lactate distribution (N = ", get_n_surg(cldta_first_last), ")"),
      x = "Lactate (mmol/L)",
      y = "Count",
      fill = "Group"
    ) +
    theme_minimal()
  
  return(p)
}

plot_lactate_over_time_quartiles <- function(df,
                                             start_time = "iv_evth_lactate_start",
                                             end_time   = "iv_evtof1_lactate_end",
                                             start_lact = "lactate_lactate_start",
                                             end_lact   = "lactate_lactate_end",
                                             n_bins = 20) {
  
  # Create time bins for averaging
  df <- df %>%
    mutate(
      start_bin = ntile(.data[[start_time]], n_bins),
      end_bin   = ntile(.data[[end_time]], n_bins)
    )
  
  summarize_ci <- function(x) {
    n <- sum(!is.na(x))
    mean_val <- mean(x, na.rm = TRUE)
    sd_val <- sd(x, na.rm = TRUE)
    error <- 1.96 * sd_val / sqrt(n)
    tibble(mean = mean_val, lower = mean_val - error, upper = mean_val + error)
  }
  
  # Summarize start lactate
  avg_start <- df %>%
    group_by(start_bin) %>%
    summarise(time = mean(.data[[start_time]], na.rm = TRUE),
              summarize_ci(.data[[start_lact]]), .groups = "drop")
  
  # Summarize end lactate
  avg_end <- df %>%
    group_by(end_bin) %>%
    summarise(time = mean(.data[[end_time]], na.rm = TRUE),
              summarize_ci(.data[[end_lact]]), .groups = "drop")
  
  # Start Lactate plot
  p_start <- ggplot(avg_start, aes(x = time, y = mean)) +
    geom_line(color = "steelblue", size = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.2) +
    geom_point(color = "steelblue", size = 2) +
    labs(x = "Time from CPB start (hours)", y = "Baseline Lactate (mmol/L)") +
    ylim(0, 6) +
    theme_minimal(base_size = 14)
  
  # Peak Lactate plot
  p_end <- ggplot(avg_end, aes(x = time, y = mean)) +
    geom_line(color = "salmon", size = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "salmon", alpha = 0.2) +
    geom_point(color = "salmon", size = 2) +
    labs(x = "Time from CPB end (hours)", y = "Peak Lactate (mmol/L)") +
    ylim(0, 6) +
    theme_minimal(base_size = 14)
  
  # Combine with patchwork and add global title
  combined_plot <- p_start + p_end + 
    plot_layout(ncol = 2, guides = "collect") +
    plot_annotation(title = paste0(n_bins, "-Quantile Average Baseline and Peak Lactate Over Time"),
                    theme = theme(plot.title = element_text(size = 18)))
  
  return(combined_plot)
}


plot_time_hist <- function(df,
                           start_time = "iv_evth_lactate_start",
                           end_time   = "iv_evtof1_lactate_end",
                           bins = 30) {
  
  # Start time histogram
  p_start <- ggplot(df, aes(x = .data[[start_time]])) +
    geom_histogram(bins = bins, fill = "steelblue", color = "black", alpha = 0.7) +
    labs(x = "Time from CPB start (hours)", y = "Count") +
    theme_minimal()
  
  # End time histogram
  p_end <- ggplot(df, aes(x = .data[[end_time]])) +
    geom_histogram(bins = bins, fill = "salmon", color = "black", alpha = 0.7) +
    labs(x = "Time from CPB end (hours)", y = "Count") +
    theme_minimal()
  
  # Combine with patchwork
  combined_plot <- p_start + p_end + plot_layout(ncol = 2, guides = "collect")
  
  return(combined_plot)
}

# Peak lactate ----
get_first_peak_lactate_df <- function(iv_evth_lower = -1, iv_evth_upper = 0.5, 
                                      iv_evtof1_lower = -0.25, iv_evtof1_upper = 12) {
  
  cldta_first_last <- cldta %>%
    group_by(surg_id) %>%
    filter(n() > 1) %>%
    ungroup() %>%
    mutate(lactate_type = case_when(
      iv_evth >= iv_evth_lower & iv_evth <= iv_evth_upper ~ "lactate_start",
      iv_evtof1 >= iv_evtof1_lower & iv_evtof1 <= iv_evtof1_upper ~ "lactate_end",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(lactate_type)) %>%
    group_by(surg_id, lactate_type) %>%
    slice_min(
      order_by = case_when(
        lactate_type == "lactate_start" ~ abs(iv_evth),
        lactate_type == "lactate_end"   ~ -lactate   # pick the largest lactate value
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
      values_from = c(lactate, iv_evth, iv_evtof1),
      names_glue = "{.value}_{lactate_type}"
    )
  
  print(dim(cldta_first_last)) 
  print(get_n_surg(cldta_first_last)) 
  
  return(cldta_first_last)
}


# Feature engineering ---- 

## Longitudinal vars ----
get_mean_summary <- function(df_to_join_to, df, var) {
  var_sym <- sym(var)
  
  summary_df <- df %>%
    arrange(surg_id, iv_evth) %>%
    group_by(surg_id) %>%
    summarise(
      mean_val = mean(!!var_sym, na.rm = TRUE),
      .groups = "drop"
    )
  
  df_to_join_to <- df_to_join_to %>%
    left_join(
      summary_df %>%
        rename(!!paste0(var, "_mean") := mean_val),
      by = "surg_id"
    )
  
  return(df_to_join_to)
}

get_mean_summary_segments <- function(df_to_join_to, df, var) {
  var_sym <- sym(var)
  
  summary_df <- df %>%
    arrange(surg_id, iv_evth) %>%
    group_by(surg_id) %>%
    mutate(
      segment = case_when(
        iv_evth <= 5/60 ~ "first_five",
        iv_evtof1 >= -5/60   ~ "last_five",
        TRUE            ~ "middle_cpb"
      )
    ) %>%
    group_by(surg_id, segment) %>%
    summarise(
      mean_val = mean(!!var_sym, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = segment,
      values_from = mean_val
    )
  
  # Rename columns AFTER pivot_wider
  summary_df <- summary_df %>% 
    rename_with(~ paste0(var, "_", .x, "_mean"), -surg_id)
  
  df_to_join_to <- df_to_join_to %>%
    left_join(summary_df, by = "surg_id")
  
  return(df_to_join_to)
}

get_min_summary <- function(df_to_join_to, df, var) {
  var_sym <- sym(var)
  
  summary_df <- df %>%
    group_by(surg_id) %>%
    summarise(
      min_val = min(!!var_sym, na.rm = TRUE),
      .groups = "drop"
    )
  
  df_to_join_to <- df_to_join_to %>%
    left_join(
      summary_df %>%
        rename(!!paste0(var, "_min") := min_val),
      by = "surg_id"
    )
  
  return(df_to_join_to)
}

get_max_summary <- function(df_to_join_to, df, var) {
  var_sym <- sym(var)
  
  summary_df <- df %>%
    group_by(surg_id) %>%
    summarise(
      max_val = max(!!var_sym, na.rm = TRUE),
      .groups = "drop"
    )
  
  df_to_join_to <- df_to_join_to %>%
    left_join(
      summary_df %>%
        rename(!!paste0(var, "_max") := max_val),
      by = "surg_id"
    )
  
  return(df_to_join_to)
}

get_sd_summary <- function(df_to_join_to, df, var) {
  var_sym <- sym(var)
  
  summary_df <- df %>%
    group_by(surg_id) %>%
    summarise(
      sd_val = sd(!!var_sym, na.rm = TRUE),
      .groups = "drop"
    )
  
  df_to_join_to <- df_to_join_to %>%
    left_join(
      summary_df %>%
        rename(!!paste0(var, "_sd") := sd_val),
      by = "surg_id"
    )
  
  return(df_to_join_to)
}

get_precpb_mean_summary <- function(df_to_join_to, df, var, within_x_hours = 24) {
  var_sym <- sym(var)

  summary_df <- df %>%
    filter(iv_evth < 0 & iv_evth > -within_x_hours) %>% 
    group_by(surg_id) %>%
    summarise(
      mean_val = mean(!!var_sym, na.rm = TRUE),
      .groups = "drop"
    )
  
  df_to_join_to <- df_to_join_to %>%
    left_join(
      summary_df %>%
        rename(!!paste0("precpb_", var, "_mean") := mean_val),
      by = "surg_id"
    )
  
  return(df_to_join_to)
}

get_postcpb_mean_summary <- function(df_to_join_to, df, var) {
  var_sym <- sym(var)
  
  # bring iv_evth_lactate_end into df
  df <- df %>%
    left_join(
      df_to_join_to %>% select(surg_id, iv_evth_lactate_end),
      by = "surg_id"
    )
  
  # filter by each patient’s lactate end time
  summary_df <- df %>%
    filter(iv_evtof1 > 0 & iv_evth < iv_evth_lactate_end) %>%
    group_by(surg_id) %>%
    summarise(
      mean_val = mean(!!var_sym, na.rm = TRUE),
      .groups = "drop"
    )
  
  # join back
  df_to_join_to <- df_to_join_to %>%
    left_join(
      summary_df %>%
        rename(!!paste0("postcpb_", var, "_mean") := mean_val),
      by = "surg_id"
    )
  
  return(df_to_join_to)
}

get_time_summary <- function(df_to_join_to, df, var, cutoff, direction = c("below", "above"),
                             window_cutoff = 10/60) {
  direction <- match.arg(direction)
  var_sym <- sym(var)

  summary_df <- df %>%
    arrange(surg_id, iv_evth) %>%
    group_by(surg_id) %>%
    mutate(dt = lead(iv_evth) - iv_evth) %>%
    summarise(
      time_hr = sum(
        if_else(
          (direction == "below" & (!!var_sym < cutoff)) |
            (direction == "above" & (!!var_sym > cutoff)),
          dt, 0
        ),
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    mutate(
      flag = ifelse(is.na(time_hr), NA_integer_,
                    ifelse(time_hr >= window_cutoff, 1L, 0L))
    )

  df_to_join_to <- df_to_join_to %>%
    left_join(
      summary_df %>%
        rename(
          !!paste0(var, "_time_", direction, "_", cutoff) := time_hr,
          !!paste0(var, "_at_least_", round(window_cutoff*60), "_min_", direction, "_", cutoff) := flag
        ),
      by = "surg_id"
    )

  return(df_to_join_to)
}


convert_to_cpb_frac <- function(df) {
  # match columns that contain _time_above_ or _time_below_ anywhere
  time_cols <- grep("_time_(above|below)_", names(df), value = TRUE)

  df %>%
    mutate(across(
      all_of(time_cols),
      ~ ifelse(!is.na(cpb_duration) & cpb_duration > 0, .x * 100 / cpb_duration, NA_real_),
      .names = "{.col}_frac"
    ))
}

convert_to_mcg_kg_min <- function(df) {
  # match columns that contain _time_above_ or _time_below_ anywhere
  time_cols <- grep("_dose_mcg_kg$", names(df), value = TRUE)
  
  df %>%
    mutate(across(
      all_of(time_cols),
      ~ ifelse(!is.na(cpb_duration) & cpb_duration > 0, .x / (60 * cpb_duration), NA_real_),
      .names = "{.col}_min"
    ))
}

## Pressors ----
get_pressor_df <- function(df_to_join_to) {
  
  summed_list <- map2(names(med_cols), med_cols, function(df_name, col_name) {
    df_onpump_list[[df_name]] %>%
      group_by(surg_id) %>%
      summarise(!!col_name := sum(.data[[col_name]], na.rm = TRUE), .groups = "drop")
  }) %>%
    compact()
  
  summed_df <- reduce(summed_list, full_join, by = "surg_id") 
  
  df_to_join_to <- df_to_join_to %>%
    left_join(summed_df, by = "surg_id") %>%
    mutate(across(all_of(setdiff(names(summed_df)[-1], "phen_dose_mcg_kg")), ~replace_na(.x, 0)))
  
  return(df_to_join_to)
}

get_pressor_df_mcg <- function(df_to_join_to) {
  
  summed_list <- map2(names(med_cols_mcg), med_cols_mcg, function(df_name, col_name) {
    df_onpump_list[[df_name]] %>%
      group_by(surg_id) %>%
      summarise(!!col_name := sum(.data[[col_name]], na.rm = TRUE), .groups = "drop")
  }) %>%
    compact()
  
  summed_df <- reduce(summed_list, full_join, by = "surg_id") 
  
  df_to_join_to <- df_to_join_to %>%
    left_join(summed_df, by = "surg_id") %>%
    mutate(across(all_of(setdiff(names(summed_df)[-1], "phen_dose_mcg")), ~replace_na(.x, 0)))
  
  return(df_to_join_to)
}

get_precpb_pressors_df <- function(df_to_join_to, within_x_hours = 24) {
  summed_list <- map2(names(med_cols), med_cols, function(df_name, col_name) {
    if (df_name == "med_phen") {
      df_list[[df_name]] %>%
        filter(neobolus == 0 & precpb == 1  & iv_evth > -within_x_hours) %>% 
        group_by(surg_id) %>%
        summarise(!!col_name := sum(.data[[col_name]], na.rm = TRUE), .groups = "drop")
    } else {
      df_list[[df_name]] %>%
        filter(precpb == 1 & iv_evth > -within_x_hours) %>% 
        group_by(surg_id) %>%
        summarise(!!col_name := sum(.data[[col_name]], na.rm = TRUE), .groups = "drop")
    }
  }) %>%
    compact() 
  
  summed_df <- reduce(summed_list, full_join, by = "surg_id")
  summed_df <- summed_df %>%
    rename_with(~paste0("precpb_", .), .cols = -surg_id)
  df_to_join_to <- df_to_join_to %>%
    left_join(summed_df, by = "surg_id") %>%
    mutate(across(matches("^precpb_(phen|epi|nepi|vaso|milr|nitp)"), ~replace_na(.x, 0)))
  
  return(df_to_join_to)
}

get_postcpb_pressors_df <- function(df_to_join_to) {
  summed_list <- map2(names(med_cols), med_cols, function(df_name, col_name) {
    
    med_df <- df_list[[df_name]] %>%
      # join the lactate end time and replace negative values with 0
      left_join(df_to_join_to %>% 
                  select(surg_id, iv_evth_lactate_end) %>%
                  mutate(iv_evth_lactate_end = pmax(iv_evth_lactate_end, 0.001)), # If peak before CPB end then should be 0 
                by = "surg_id")
    
    if (df_name == "med_phen") {
      med_df <- med_df %>%
        filter(neobolus == 0 & iv_evtof1 > 0 & iv_evth < iv_evth_lactate_end)
    } else {
      med_df <- med_df %>%
        filter(iv_evtof1 > 0 & iv_evth < iv_evth_lactate_end)
    }
    
    med_df %>%
      group_by(surg_id, iv_evth_lactate_end) %>%
      summarise(sum_dose = sum(.data[[col_name]], na.rm = TRUE), .groups = "drop") %>%
      mutate(rate_per_min = ifelse(iv_evth_lactate_end > 0, sum_dose / (60 * iv_evth_lactate_end), 0)) %>%
      select(surg_id, rate_per_min) %>%
      rename(!!col_name := rate_per_min)
    
  }) %>%
    compact()
  
  summed_df <- reduce(summed_list, full_join, by = "surg_id") %>%
    rename_with(~paste0("postcpb_", ., "_min"), .cols = -surg_id)
  
  df_to_join_to <- df_to_join_to %>%
    left_join(summed_df, by = "surg_id") %>%
    mutate(across(matches("^postcpb_(phen|epi|nepi|vaso|milr)"), ~replace_na(.x, 0)))
  
  return(df_to_join_to)
}


## Preop features ----
get_preop_features <- function(df_to_join_to,
                               preop_vars_keep_na, 
                               preop_vars_na_to_0,        # variables to NA→0 numeric
                               dummy_vars,        # variables to expand into dummies
                               dummy_remove_first = FALSE) {
  # These vars keep NAs 
  preop_hx_vars <- cldta_onpump_baseline_df %>%
    select(surg_id, all_of(preop_vars_keep_na))
  df_to_join_to <- df_to_join_to %>%
    left_join(preop_hx_vars, by = "surg_id") 
  
  # These vars NAs --> 0 
  preop_hx_vars <- cldta_onpump_baseline_df %>%
    select(surg_id, all_of(preop_vars_na_to_0))
  df_to_join_to <- df_to_join_to %>%
    left_join(preop_hx_vars, by = "surg_id") %>%
    mutate(across(all_of(preop_vars_na_to_0),
                  ~ as.numeric(replace_na(.x, 0))))
  
  # These vars can keep original + compute dummies 
  df_to_join_to <- df_to_join_to %>%
    left_join(cldta_onpump_baseline_df %>% select(surg_id, all_of(dummy_vars)), by = "surg_id")
  
  for (v in dummy_vars) {
    df_to_join_to <- df_to_join_to %>%
      mutate("{v}" := factor(.data[[v]], levels = sort(unique(.data[[v]])), labels = sort(unique(.data[[v]])))) %>%
      mutate("{v}" := fct_explicit_na(.data[[v]], na_level = "0")) %>%
      fastDummies::dummy_cols(select_columns = v,
                              remove_first_dummy = dummy_remove_first,
                              remove_selected_columns = FALSE)
  }
  
  return(df_to_join_to)
}

# Clustering ----
pacman::p_load(pheatmap)

get_baseline_clusters <- function(k = 8) {
  # --- Select baseline vars ---
  baseline_var_to_english_combined <- c(
    # age       = "Age",
    # bmi       = "Body Mass Index",
    # bsa       = "Body Surface Area",
    sex       = "Sex",
    race_comb = "Race",
    hx_csurg  = "History of Cardiac Surgery",
    hx_cvd    = "History of Cerebrovascular Disease",
    hx_copd   = "History of COPD",
    # hx_dm     = "History of Diabetes Mellitus",
    hx_dyslp  = "History of Dyslipidemia",
    hx_endo   = "History of Endocarditis",
    hx_chf    = "History of Congestive Heart Failure",
    hx_htn    = "History of Hypertension",
    hx_mi     = "History of Myocardial Infarction",
    hx_pad    = "History of Peripheral Artery Disease",
    hx_smoke  = "History of Smoking",
    hx_cva    = "History of Stroke",
    hx_aorta  = "History of Aortic Disease",
    hx_arysu  = "History of Aortic Aneurysm",
    hx_cabg   = "History of CABG",
    hx_CVInt  = "History of Cardiovascular Intervention",
    hx_cngsu  = "History of Congenital Heart Surgery",
    hx_icd    = "Implanted Cardioverter Defibrillator",
    hx_ppm    = "Permanent Pacemaker",
    afib_pr   = "Preop Atrial Fibrillation",
    dial_pr   = "Preop Dialysis",
    chb_pr    = "Preop Complete Heart Block",
    varr_pr   = "Preop Aortic Regurgitation",
    blrbn_pr  = "Preop Bilirubin", 
    # bun_pr    = "Preop BUN",
    creat_pr  = "Preop Creatinine",
    # hgb_pr    = "Preop Hemoglobin",
    # lvef_pr   = "Preop LVEF",
    bnpl_pr   = "Preop BNP",
    tbill_pr  = "Preop Total Bilirubin", # STOP HERE 
    chol_pr   = "Preop Cholesterol",
    esr_pr    = "Preop ESR",
    hdl_pr    = "Preop HDL",
    hgbac_pr  = "Preop Hemoglobin A1c",
    ldl_pr    = "Preop LDL",
    lymp_pr   = "Preop Lymphocyte Count",
    # mcv_pr    = "Preop MCV",
    pot_pr    = "Preop Potassium",
    ptinr_pr  = "Preop INR",
    # wbc_pr    = "Preop WBC",
    # gfr_pr    = "Preop GFR",
    sp_valve  = "Any Valve",
    sp_av     = "Any Aortic Valve",
    sp_mv     = "Any Mitral Valve",
    sp_pv     = "Any Pulmonic Valve",
    sp_tv     = "Any Tricuspid Valve Procedure",
    sp_aorta  = "Any Aortic Procedure",
    # sp_cabg   = "Coronary Artery Bypass Grafting",
    sp_afib   = "Atrial Fibrillation",
    sp_asdpf  = "ASD / PFO Closure",
    sp_chd    = "Congenital Heart Defect Repair",
    sp_atrmy  = "Atrial Myxoma Resection",
    sp_peric  = "Pericardiectomy",
    sp_sepmy  = "Septal Myectomy",
    sp_lv     = "Any Left Ventricular",
    sp_tx     = "Heart Transplant",
    sp_cend   = "Carotid Endarterectomy"
  )
  
  vars_to_keep <- intersect(names(baseline_var_to_english_combined), names(cldta))
  
  cldta_baseline <- cldta_onpump_baseline_df %>%
    select(all_of(c("surg_id", vars_to_keep))) %>%
    group_by(surg_id) %>%
    dplyr::slice(1) %>%
    ungroup()
  
  # --- Filter out vars with too many NAs ---
  na_threshold <- 8000
  na_counts <- colSums(is.na(cldta_baseline))
  cols_to_keep <- names(na_counts)[na_counts <= na_threshold]
  
  cldta_baseline_filtered <- cldta_baseline %>%
    select(all_of(cols_to_keep))
  
  # --- Impute missing values ---
  cldta_baseline_imputed <- cldta_baseline_filtered %>%
    mutate(across(where(is.logical), ~ ifelse(is.na(.), FALSE, .))) %>%
    mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
  
  cldta_baseline_imputed$female <- ifelse(tolower(cldta_baseline_imputed$sex) == "m", 0,
                                          ifelse(tolower(cldta_baseline_imputed$sex) == "f", 1, NA))
  cldta_baseline_imputed$race_comb <- as.numeric(factor(cldta_baseline_imputed$race_comb,
                                                        levels = c("White", "Black", "Other"))) - 1
  
  cldta_baseline_imputed <- cldta_baseline_imputed %>%
    mutate(across(where(is.logical), ~ as.integer(.)))
  
  surg_id_vector <- cldta_baseline_imputed$surg_id
  cldta_baseline_imputed <- cldta_baseline_imputed %>%
    select(-surg_id, -sex)
  
  # --- Winsorize function ---
  winsorize <- function(x, lower = 0.01, upper = 0.99) {
    q <- quantile(x, probs = c(lower, upper), na.rm = TRUE)
    pmin(pmax(x, q[1]), q[2])
  }
  
  # --- Winsorize numeric cols (exclude binary) ---
  numeric_cols <- cldta_baseline_imputed %>%
    select(where(is.numeric)) %>%
    select(where(~ length(unique(.)) > 3)) %>%
    names()
  
  cldta_baseline_winsor <- cldta_baseline_imputed %>%
    mutate(across(all_of(numeric_cols), ~ winsorize(.x, 0.01, 0.99)))
  
  # --- Scale numeric cols ---
  cldta_baseline_scaled <- cldta_baseline_winsor %>%
    mutate(across(all_of(numeric_cols), ~ scale(.x)))
  
  # --- Run pheatmap (silent) to get tree ---
  cldta_baseline_scaled_df <- as.data.frame(cldta_baseline_scaled)
  rownames(cldta_baseline_scaled_df) <- cldta_baseline$surg_id
  
  p <- pheatmap(cldta_baseline_scaled_df,
                clustering_method = "ward.D2",
                show_rownames = FALSE,
                show_colnames = FALSE,
                scale = "none",
                silent = TRUE)
  
  # --- Cut dendrogram into clusters ---
  row_clusters <- cutree(p$tree_row, k = k)
  
  cluster_df <- data.frame(
    surg_id = names(row_clusters),
    cluster = as.factor(row_clusters),
    stringsAsFactors = FALSE
  )
  
  return(cluster_df)
}


# Modeling ---- 
# x_to_x_mat <- function(X) {
#   binary_vars <- names(X)[sapply(X, function(col) {
#     is.factor(col) && length(levels(col)) == 2 ||
#       all(unique(col) %in% c(0, 1))
#   })]
#   
#   multi_level_factors <- names(X)[sapply(X, is.factor) & !(names(X) %in% binary_vars)]
#   
#   X_dummy <- fastDummies::dummy_cols(
#     X,
#     select_columns = multi_level_factors,
#     remove_first_dummy = FALSE,
#     remove_selected_columns = TRUE
#   )
#   
#   X_mat <- as.matrix(X_dummy)
#   
#   return(X_mat)
# }

x_to_x_mat <- function(X) {
  # identify binary vars (factor with 2 levels OR 0/1 numeric/integer)
  binary_vars <- names(X)[sapply(X, function(col) {
    (is.factor(col) && length(levels(col)) == 2) ||
      (all(unique(col) %in% c(0, 1)))
  })]
  
  # identify multi-level factors (factor vars not in binary_vars)
  factor_vars <- names(X)[sapply(X, is.factor)]
  multi_level_factors <- setdiff(factor_vars, binary_vars)
  
  # dummy encode only multi-level factors
  if (length(multi_level_factors) > 0) {
    X_dummy <- fastDummies::dummy_cols(
      X,
      select_columns = multi_level_factors,
      remove_first_dummy = FALSE,
      remove_selected_columns = TRUE
    )
  } else {
    X_dummy <- X
  }
  
  X_dummy[] <- lapply(X_dummy, function(col) as.numeric(as.character(col)))
  
  X_mat <- as.matrix(X_dummy)
  return(X_mat)
}


get_partial_grid <- function(X, var_name, quantile_points = 50) {
  # Use only unique non-NA values
  vals <- unique(X[[var_name]][!is.na(X[[var_name]])])
  
  # Compute percentiles across truncated unique values
  q_probs <- seq(0, 1, length.out = quantile_points)
  grid_q  <- quantile(vals, probs = q_probs, type = 7)
  
  # Return as data.frame
  df <- data.frame(var = sort(unique(grid_q)))
  colnames(df) <- var_name
  return(df)
}

k_fold_xgboost <- function(X, y, k = 10, params, 
                           nrounds = 1200, 
                           quantile_points = 50,
                           early_stopping_rounds = 20, 
                           seed = 123, 
                           verbose = 0, 
                           get_partial = TRUE,
                           get_ale = FALSE,
                           loss_function = NULL,
                           eval_function = NULL, 
                           is_classification = FALSE,
                           is_multiclass = FALSE) {   
  cat("Performing", k, "fold cross-validation of XGBoost model\n")
  
  set.seed(seed)
  
  X_mat <- x_to_x_mat(X)
  
  # mapping original vars -> actual columns in X_mat
  factor_map <- list()
  for (var in names(X)) {
    if (is.factor(X[[var]])) {
      factor_map[[var]] <- grep(paste0("^", var), colnames(X_mat), value = TRUE)
    } else {
      factor_map[[var]] <- var
    }
  }
  
  # build grids for PDP
  grid_list <- lapply(names(X), function(var_name) {
    if (is.factor(X[[var_name]])) {
      dummies <- factor_map[[var_name]]
      df <- diag(1, nrow = length(dummies), ncol = length(dummies))
      df <- as.data.frame(df)
      colnames(df) <- dummies
      df
    } else if (all(unique(X[[var_name]]) %in% c(0, 1))) {
      df <- data.frame(var_name = c(0, 1))
      colnames(df) <- var_name
      df
    } else {
      get_partial_grid(X, var_name, quantile_points)
    }
  })
  names(grid_list) <- names(X)
  
  # K-fold splits
  folds <- sample(rep(1:k, length.out = nrow(X)))
  
  # metric labels
  metric_name <- params$eval_metric
  train_col <- paste0(metric_name, "_train")
  valid_col <- paste0(metric_name, "_valid")
  
  # storage
  partial_list <- list()
  ale_list <- list()
  model_list <- list()
  performance_list <- data.frame(fold = integer())
  
  # Determine model type
  # objective <- params$objective
  # is_classification <- grepl("multi:|binary:", objective)
  num_classes <- ifelse(is_multiclass, params$num_class, 1)
  
  # CV loop
  for (i in 1:k) {
    train_idx <- which(folds != i)
    valid_idx <- which(folds == i)
    
    dtrain <- xgb.DMatrix(data = X_mat[train_idx, ], label = y[train_idx])
    dvalid <- xgb.DMatrix(data = X_mat[valid_idx, ], label = y[valid_idx])
    
    if(!is.null(loss_function)) { 
      
      fit <- xgb.train(
        params = params,
        data = dtrain,
        nrounds = nrounds,
        watchlist = list(train = dtrain, eval = dvalid),
        early_stopping_rounds = early_stopping_rounds,
        verbose = verbose, 
        obj = loss_function,
        feval = eval_function, 
        maximize = FALSE
      )
    }
    else {
      fit <- xgb.train(
        params = params,
        data = dtrain,
        nrounds = nrounds,
        watchlist = list(train = dtrain, eval = dvalid),
        early_stopping_rounds = early_stopping_rounds,
        verbose = verbose
      )
    }

    model_list[[i]] <- fit
    
    best_iter <- fit$best_iteration
    evaluation_log <- fit$evaluation_log
    metric_cols <- colnames(evaluation_log)
    train_col_name <- metric_cols[grepl("^train", metric_cols)][1]
    valid_col_name <- metric_cols[grepl("^eval", metric_cols)][1]
    train_val <- evaluation_log[[train_col_name]][best_iter]
    valid_val <- evaluation_log[[valid_col_name]][best_iter]
    curr_row <- data.frame(fold = i, train_metric = train_val, valid_metric = valid_val)
    performance_list <- rbind(performance_list, curr_row)
    print(curr_row)
    
    # --- Partial Dependence ---
    if (get_partial) {
      for (var in names(X)) {
        pred_cols <- factor_map[[var]]
        
        if (!is_classification) {
          pd <- pdp::partial(
            object = fit,
            pred.var = pred_cols,
            train = as.data.frame(X_mat[train_idx, , drop = FALSE]),
            pred.grid = grid_list[[var]],
            progress = "none"
          )
          pd$fold <- i
          pd$variable <- var
          partial_list[[length(partial_list) + 1]] <- pd
        } else {
          for (cls in 1:num_classes) {
            pd <- pdp::partial(
              object = fit,
              pred.var = pred_cols,
              train = as.data.frame(X_mat[train_idx, , drop = FALSE]),
              pred.grid = grid_list[[var]],
              type = "classification",
              which.class = cls,
              progress = "none",
              prob = TRUE
            )
            pd$fold <- i
            pd$variable <- var
            pd$class <- cls
            partial_list[[length(partial_list) + 1]] <- pd
          }
        }
      }
    }
    
    # --- ALE via iml ---
    if (get_ale) {

      predictor <- Predictor$new(
        model = fit,
        data = as.data.frame(X_mat[train_idx, , drop = FALSE]),
        y = y[train_idx],
        predict.function = function(m, newdata) {
          pred <- predict(m, as.matrix(newdata))
          if (is_classification) {
            if (num_classes > 1) {
              mat <- matrix(pred, ncol = num_classes, byrow = TRUE)
              return(mat)
            } else {
              return(pred)
            }
          } else {
            return(pred)
          }
        }
      )
      
      for (var in names(X)) {
        if (!is_classification) {
          
          grid_points_numeric <- grid_list[[var]][,1]
          ale_obj <- FeatureEffect$new(predictor, 
                                       feature = var, 
                                       method = "ale",
                                       grid.points = grid_points_numeric)

          df_ale <- ale_obj$results
          
          # standardize to match PDP output
          df_ale <- data.frame(
            x = df_ale[[var]],
            yhat = df_ale$.value,
            variable = var,
            fold = i
          )
          colnames(df_ale)[1] <- var   # rename "x" column to actual var name
          ale_list[[length(ale_list) + 1]] <- df_ale
          
        } else {
          for (cls in 1:num_classes) {
            ale_obj <- FeatureEffect$new(
              predictor,
              feature = var,
              method = "ale",
              class = cls
            )
            df_ale <- ale_obj$results
            
            df_ale <- data.frame(
              x = df_ale[[var]],
              yhat = df_ale$.value,
              variable = var,
              class = cls,
              fold = i
            )
            colnames(df_ale)[1] <- var
            ale_list[[length(ale_list) + 1]] <- df_ale
          }
        }
      }
    }
  } 
  
  return(list(
    model_list = model_list,
    performance_list = performance_list,
    partial_list = partial_list,
    ale_list = ale_list,
    folds = folds
  ))
}



print_rmse_mean_ci <- function(performance_list) {
  rmse_vals <- performance_list$valid_metric
  
  mean_rmse <- mean(rmse_vals)
  ci_rmse   <- quantile(rmse_vals, probs = c(0.025, 0.975))
  
  cat("Mean RMSE for validation set:", round(mean_rmse, 3),
      " (95% CI:", round(ci_rmse[1], 3), "–", round(ci_rmse[2], 3), ")\n")
}


print_xgb_performance <- function(xgb_results, y_numeric,
                                  metric = "mlogloss") {
  perf <- xgb_results$performance_list
  n_folds <- nrow(perf)
  
  # Compute baseline mlogloss from class frequencies
  y_factor <- factor(y_numeric)
  class_probs <- table(y_factor) / length(y_factor)
  mlogloss_baseline <- -mean(log(class_probs[as.character(y_factor)]))
  
  cat("XGBoost 10-fold CV performance:\n")
  
  # Function to compute mean and 95% CI
  mean_ci <- function(x) {
    x <- x[!is.na(x)]
    n <- length(x)
    m <- mean(x)
    se <- sd(x)/sqrt(n)
    ci_lower <- m - 1.96*se
    ci_upper <- m + 1.96*se
    list(mean = m, lower = ci_lower, upper = ci_upper)
  }
  
  # Print mlogloss if it exists
  if (metric == "mlogloss") {
    ci <- mean_ci(perf$valid_metric)
    cat(sprintf("- Mean mlogloss (validation): %.3f (95%% CI: %.3f - %.3f)\n",
                ci$mean, ci$lower, ci$upper))
    cat(sprintf("- Baseline mlogloss (naive class freq): %.3f\n", mlogloss_baseline))

  }
  
  # Print merror if it exists
  if (metric == "merror") {
    ci <- mean_ci(perf$valid_metric)
    cat(sprintf("- Mean misclassification rate (merror, validation): %.3f (95%% CI: %.3f - %.3f)\n",
                ci$mean, ci$lower, ci$upper))
    cat(sprintf("- Mean accuracy (1 - merror): %.3f (95%% CI: %.3f - %.3f)\n",
                1-ci$mean, 1-ci$upper, 1-ci$lower)) # flip CI for accuracy
  }
  
  # Return values invisibly for further use
  invisible(list(
    mean_mlogloss = if ("mlogloss_valid" %in% names(perf)) ci$mean else NULL,
    baseline_mlogloss = if ("mlogloss_valid" %in% names(perf)) mlogloss_baseline else NULL,
    mean_merror = if ("merror_valid" %in% names(perf)) ci$mean else NULL
  ))
}

# Model performance ---- 
plot_actual_vs_pred <- function(X, y, model_list, folds, performance_list = NULL, outcome_label = "Lactate (mmol/L)") {
  k <- length(model_list)
  
  X_mat <- x_to_x_mat(X)
  
  # Collect predictions per fold
  pred_df <- bind_rows(
    lapply(1:k, function(i) {
      valid_idx <- which(folds == i)
      fit <- model_list[[i]]
      
      dvalid <- xgb.DMatrix(data = X_mat[valid_idx, ], label = y[valid_idx])
      pred_valid <- predict(fit, dvalid)
      
      data.frame(
        fold = factor(i),
        Actual = y[valid_idx],
        Predicted = pred_valid
      )
    })
  )
  
  # Plot
  p <- ggplot(pred_df, aes(x = Actual, y = Predicted, color = fold)) +
    geom_point(alpha = 0.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(
      title = paste("XGBoost performance over", k, "fold Cross Validation"),
      x = paste("True", outcome_label),
      y = paste("Predicted", outcome_label),
      color = "Fold"
    ) +
    theme_minimal()
  
  if (!is.null(performance_list)) {
    rmse_vals <- xgboost_results$performance_list$valid_metric
    mean_rmse <- mean(rmse_vals)
    ci_rmse   <- quantile(rmse_vals, probs = c(0.025, 0.975))
    
    rmse_text <- paste0(
      "RMSE = ", sprintf("%.2f", mean_rmse),
      " (95% CI: ", sprintf("%.2f", ci_rmse[1]), "–", sprintf("%.2f", ci_rmse[2]), ")"
    )
    
    p <- p + annotate(
      "text",
      x = Inf, y = -Inf,
      hjust = 1, vjust = -1.2,   # bottom-right
      label = rmse_text,
      size = 4, color = "black"
    )
    
  }
  return(p)
}

# VIMP ---- 
plot_vimp <- function(model_list, metric = "Gain", top_n = NULL,
                      feature_labels = formula_var_name_to_english) {
  # feature_labels <- formula_var_name_to_english

  # Extract importance from each model
  all_importance <- lapply(model_list, function(mod) {
    imp <- xgb.importance(model = mod)
    imp$Feature <- as.character(imp$Feature)
    imp
  })

  # Combine and summarize
  importance_stats <- bind_rows(all_importance) %>%
    group_by(Feature) %>%
    summarise(
      mean_val = mean(.data[[metric]], na.rm = TRUE),
      sd_val   = sd(.data[[metric]], na.rm = TRUE),
      n        = n(),
      .groups  = "drop"
    ) %>%
    mutate(
      se = sd_val / sqrt(n),
      ci_lower = mean_val - 1.96 * se,
      ci_upper = mean_val + 1.96 * se
    )

  # English mapping of feature names
  if (!is.null(feature_labels)) {
    importance_stats$Feature_label <- feature_labels[importance_stats$Feature]
    importance_stats$Feature_label <- ifelse(
      is.na(importance_stats$Feature_label),
      importance_stats$Feature,
      importance_stats$Feature_label
    )
  } else {
    importance_stats$Feature_label <- importance_stats$Feature
  }

  # Order variables by mean importance
  importance_stats <- importance_stats %>% arrange(desc(mean_val))

  # Keep only top_n if specified
  if (!is.null(top_n)) {
    importance_stats <- importance_stats[1:min(top_n, nrow(importance_stats)), ]
  }

  var_order <- importance_stats$Feature_label

  # Plot
  p <- ggplot(importance_stats, aes(x = reorder(Feature_label, mean_val), y = mean_val)) +
    geom_col(fill = "royalblue1") +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.3, color = "black") +
    coord_flip() +
    labs(
      x = "Feature",
      y = paste("Average", metric, "(95% CI)"),
      title = "XGBoost Variable Importance" # paste("XGBoost Variable Importance (", metric, ")", sep = "")
    ) +
    theme_minimal()

  print(p)
  return(list(plot = p, var_order = names(var_order)))
}


# Partial Plots ----
plot_pdp_aats <- function(xgboost_results, X, variables = NULL, 
                     numeric_var_order = NULL,
                     x_lim_truncation = TRUE, 
                     is_ALE_plot = FALSE,
                     ncol = NULL) {  
  
  if(is_ALE_plot) {
    partial_list <- xgboost_results$ale_list
  } else {
    partial_list <- xgboost_results$partial_list
  }
  
  k <- length(unique(unlist(lapply(partial_list, function(df) unique(df$fold)))))
  feature_labels <- formula_var_name_to_english
  
  # Identify categorical vars (factors or binary 0/1)
  factor_vars <- names(X)[sapply(names(X), function(var) {
    is.factor(X[[var]]) || all(unique(X[[var]]) %in% c(0, 1))
  })]
  
  partial_list_clean <- lapply(partial_list, function(df) {
    varname <- unique(df$variable)
    matching_cols <- grep(paste0("^", varname, "_"), names(df), value = TRUE)
    
    if (length(matching_cols) > 0) {
      idx <- apply(df[matching_cols], 1, function(row) which(row == 1)[1])
      labels <- gsub(paste0(varname, "_"), "", matching_cols)
      df$x <- factor(labels[idx], levels = labels)
      
    } else if (!is.null(df[[varname]])) {
      x_col <- df[[varname]]
      if (varname %in% factor_vars) {
        if (is.factor(X[[varname]])) {
          df$x <- factor(x_col, levels = levels(X[[varname]]))
        } else {
          df$x <- factor(x_col, levels = c(0, 1), labels = c("No", "Yes"))
        }
      } else {
        df$x <- as.numeric(x_col)
      }
      
    } else {
      df$x <- as.numeric(df[[1]])
    }
    
    df[, c("variable", "x", "yhat", "fold")]
  })
  
  partial_df <- do.call(rbind, partial_list_clean)
  
  if (!is.null(variables)) {
    partial_df <- partial_df %>% filter(variable %in% variables)
  }
  
  factor_df <- partial_df %>% filter(variable %in% factor_vars)
  numeric_df <- partial_df %>% filter(!variable %in% factor_vars)
  
  plots <- list()
  
  # Prepare numeric rug
  numeric_rug_df <- lapply(setdiff(names(X), factor_vars), function(var) {
    data.frame(variable = var, x = X[[var]])
  }) %>% bind_rows() %>%
    filter(variable %in% unique(numeric_df$variable))
  
  # 99% range 
  x_limits <- numeric_rug_df %>%
    group_by(variable) %>%
    summarise(
      x_lim_lower = quantile(x, 0.005, na.rm = TRUE),
      x_lim_upper = quantile(x, 0.995, na.rm = TRUE),
      .groups = "drop"
    )
  
  if(x_lim_truncation) {
    numeric_rug_df <- numeric_rug_df %>%
      left_join(x_limits, by = "variable") %>%
      filter(x >= x_lim_lower & x <= x_lim_upper) %>%
      select(-x_lim_lower, -x_lim_upper)
  }
  
  # Numeric PDPs
  if (nrow(numeric_df) > 0) {
    
    pdp_summary <- numeric_df %>%
      group_by(variable, x) %>%
      summarise(
        mean_y = mean(yhat, na.rm = TRUE),
        se_y   = sd(yhat, na.rm = TRUE)/sqrt(n()),
        .groups = "drop"
      )
    
    if(x_lim_truncation) {
      pdp_summary <- pdp_summary %>%
        left_join(x_limits, by = "variable") %>%
        filter(as.numeric(x) >= x_lim_lower & as.numeric(x) <= x_lim_upper) %>%
        select(-x_lim_lower, -x_lim_upper)
    }
    
    if(!is.null(variables) & is.null(numeric_var_order)) {
      ordered_vars <- variables
    } else if (!is.null(numeric_var_order)) {
      ordered_vars <- numeric_var_order[numeric_var_order %in% pdp_summary$variable]
      extra_vars <- setdiff(unique(pdp_summary$variable), ordered_vars)
      ordered_vars <- c(ordered_vars, extra_vars)
    } else {
      ordered_vars <- unique(pdp_summary$variable)
    }
    
    label_map <- sapply(ordered_vars, function(var) {
      if (!is.null(feature_labels[var]) && !is.na(feature_labels[var])) feature_labels[var] else var
    })
    names(label_map) <- ordered_vars
    
    pdp_summary$variable_label_display <- sapply(pdp_summary$variable, function(var) {
      if (!is.null(label_map[var]) && !is.na(label_map[var])) label_map[var] else var
    })
    
    pdp_summary$variable_label_display <- factor(
      pdp_summary$variable_label_display,
      levels = label_map
    )
    
    numeric_rug_df$variable_label_display <- sapply(numeric_rug_df$variable, function(var) {
      if (!is.null(label_map[var]) && !is.na(label_map[var])) label_map[var] else var
    })
    numeric_rug_df$variable_label_display <- factor(
      numeric_rug_df$variable_label_display,
      levels = label_map
    )
    
    pdp_summary <- pdp_summary %>%
      mutate(
        lower = mean_y - 2*se_y,
        upper = mean_y + 2*se_y
      )
    
    plots$numeric <- ggplot(pdp_summary, aes(x = as.numeric(x), y = mean_y,
                                             color = variable_label_display,
                                             fill = variable_label_display)) +
      geom_line() +
      geom_point(size = 0.5) +
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = variable_label_display),
                  alpha = 0.2, color = NA) +
      # geom_rug(data = numeric_rug_df, aes(x = x), inherit.aes = FALSE,
      #          sides = "b", alpha = 0.1) +
      facet_wrap(~ variable_label_display, scales = "free_x", ncol = ncol) +
      labs(
        # x = "Predictor value",
        x = NULL, 
        y = ifelse(is_ALE_plot, "ALE Lactate", "Predicted Lactate (mmol/L)")
      ) +
      theme_minimal() +
      guides(color = "none", fill = "none") +
      theme(
        strip.text = element_text(size = 12),
        panel.spacing = unit(1, "lines")
      ) +
      # 🔹 Skip the first color in the default ggplot hue palette
      scale_color_hue(h.start = 30) +  # start hue wheel at 30° (avoids red)
      scale_fill_hue(h.start = 30)
    
      # guides(color = "none", fill = "none")  +
      # theme(strip.text = element_text(size = 12))
  }
  
  # Factor PDPs
  n_levels <- length(unique(factor_df$x))
  if (nrow(factor_df) > 0) {
    factor_df$variable_label_display <- sapply(factor_df$variable, function(var) {
      if (!is.null(feature_labels[var]) && !is.na(feature_labels[var])) feature_labels[var] else var
    })
    
    if (is.null(numeric_var_order) & !is.null(variables)) {
      factor_levels <- feature_labels[variables]
    } else if (!is.null(numeric_var_order)) {
      ordered_factor_vars <- numeric_var_order[numeric_var_order %in% factor_df$variable]
      extra_vars <- setdiff(unique(factor_df$variable), ordered_factor_vars)
      ordered_factor_vars <- c(ordered_factor_vars, extra_vars)
      factor_levels <- sapply(ordered_factor_vars, function(var) {
        if (!is.null(feature_labels[var]) && !is.na(feature_labels[var])) feature_labels[var] else var
      })
    } else {
      factor_levels <- unique(factor_df$variable_label_display)
    }
    
    factor_df$variable_label_display <- factor(factor_df$variable_label_display, levels = factor_levels)
    
    n_levels <- length(unique(factor_df$x))
    
    p <- ggplot(factor_df, aes(x = factor(x), y = yhat, fill = x)) +
      geom_boxplot(outlier.size = 0.8) +
      facet_wrap(~ variable_label_display, scales = "free_x", ncol = ncol) +  # 🔹 Added ncol
      labs(
        x = NULL,
        y = ifelse(is_ALE_plot, "ALE Lactate", "Predicted Lactate (mmol/L)")
      ) +
      theme_minimal() +
      theme(legend.position = "none")
    
    if (n_levels <= 8) {
      p <- p + scale_fill_brewer(palette = "Set2")
    } else {
      p <- p + scale_fill_viridis_d(option = "turbo")
    }
    
    plots$factor <- p
  }
  
  return(plots)
}

plot_pdp <- function(xgboost_results, X, variables = NULL,
                     numeric_var_order = NULL,
                     x_lim_truncation = TRUE,
                     is_ALE_plot = FALSE, outcome_label = "Lactate (mmol/L)",
                     feature_labels = formula_var_name_to_english) {

  if(is_ALE_plot) {
    partial_list <- xgboost_results$ale_list
  }
  else {
    partial_list <- xgboost_results$partial_list
  }

  k <- length(unique(unlist(lapply(partial_list, function(df) unique(df$fold)))))
  # feature_labels <- formula_var_name_to_english

  # Identify categorical vars (factors or binary 0/1)
  factor_vars <- names(X)[sapply(names(X), function(var) {
    is.factor(X[[var]]) || all(unique(X[[var]]) %in% c(0, 1))
  })]

  partial_list_clean <- lapply(partial_list, function(df) {
    varname <- unique(df$variable)

    # Find all columns that start with this varname
    matching_cols <- grep(paste0("^", varname, "_"), names(df), value = TRUE)

    if (length(matching_cols) > 0) {
      # Collapse one-hot encoding into a single categorical label
      idx <- apply(df[matching_cols], 1, function(row) which(row == 1)[1])
      labels <- gsub(paste0(varname, "_"), "", matching_cols)
      df$x <- factor(labels[idx], levels = labels)

    } else if (!is.null(df[[varname]])) {
      # Standard case: variable exists as a column
      x_col <- df[[varname]]
      if (varname %in% factor_vars) {
        if (is.factor(X[[varname]])) {
          df$x <- factor(x_col, levels = levels(X[[varname]]))
        } else {
          df$x <- factor(x_col, levels = c(0, 1), labels = c("No", "Yes"))
        }
      } else {
        df$x <- as.numeric(x_col)
      }

    } else {
      # Fallback: first column
      df$x <- as.numeric(df[[1]])
    }

    df[, c("variable", "x", "yhat", "fold")]
  })

  partial_df <- do.call(rbind, partial_list_clean)

  # Subset variables if requested
  if (!is.null(variables)) {
    partial_df <- partial_df %>% filter(variable %in% variables)
  }

  # Split numeric vs factor
  factor_df <- partial_df %>% filter(variable %in% factor_vars)
  numeric_df <- partial_df %>% filter(!variable %in% factor_vars)

  plots <- list()

  # Prepare numeric rug
  numeric_rug_df <- lapply(setdiff(names(X), factor_vars), function(var) {
    data.frame(variable = var, x = X[[var]])
  }) %>% bind_rows() %>%
    filter(variable %in% unique(numeric_df$variable))


  # 99% range
  x_limits <- numeric_rug_df %>%
    group_by(variable) %>%
    summarise(
      x_lim_lower = quantile(x, 0.005, na.rm = TRUE),
      x_lim_upper = quantile(x, 0.995, na.rm = TRUE),
      .groups = "drop"
    )

  if(x_lim_truncation) {
    numeric_rug_df <- numeric_rug_df %>%
      left_join(x_limits, by = "variable") %>%
      filter(x >= x_lim_lower & x <= x_lim_upper) %>%
      select(-x_lim_lower, -x_lim_upper)
  }

  # Numeric PDPs
  if (nrow(numeric_df) > 0) {

    pdp_summary <- numeric_df %>%
      group_by(variable, x) %>%
      summarise(
        mean_y = mean(yhat, na.rm = TRUE),
        se_y   = sd(yhat, na.rm = TRUE)/sqrt(n()),
        .groups = "drop"
      )

    if(x_lim_truncation) {
      pdp_summary <- pdp_summary %>%
        left_join(x_limits, by = "variable") %>%
        filter(as.numeric(x) >= x_lim_lower & as.numeric(x) <= x_lim_upper) %>%
        select(-x_lim_lower, -x_lim_upper)
    }

    if(!is.null(variables) & is.null(numeric_var_order)) {
      ordered_vars <- variables
    }
    else if (!is.null(numeric_var_order)) {
      # Keep only variables actually present
      ordered_vars <- numeric_var_order[numeric_var_order %in% pdp_summary$variable]

      # Add any extra variables from pdp_summary that weren't in numeric_var_order
      # For example if var was so unimportant it never made it into VIMP
      extra_vars <- setdiff(unique(pdp_summary$variable), ordered_vars)
      ordered_vars <- c(ordered_vars, extra_vars)
    } else {
      ordered_vars <- unique(pdp_summary$variable)
    }

    # Map human-readable labels
    label_map <- sapply(ordered_vars, function(var) {
      if (!is.null(feature_labels[var]) && !is.na(feature_labels[var])) feature_labels[var] else var
    })
    names(label_map) <- ordered_vars


    pdp_summary$variable_label_display <- sapply(pdp_summary$variable, function(var) {
      if (!is.null(label_map[var]) && !is.na(label_map[var])) label_map[var] else var
    })

    print(    label_map[duplicated(label_map)])
    pdp_summary$variable_label_display <- factor(
      pdp_summary$variable_label_display,
      levels = label_map
    )

    numeric_rug_df$variable_label_display <- sapply(numeric_rug_df$variable, function(var) {
      if (!is.null(label_map[var]) && !is.na(label_map[var])) label_map[var] else var
    })
    numeric_rug_df$variable_label_display <- factor(
      numeric_rug_df$variable_label_display,
      levels = label_map
    )

    # Compute ribbons
    pdp_summary <- pdp_summary %>%
      mutate(
        lower = mean_y - 2*se_y,
        upper = mean_y + 2*se_y
      )


    plots$numeric <- ggplot(pdp_summary, aes(x = as.numeric(x), y = mean_y, color = variable_label_display, fill = variable_label_display)) +
      geom_line() +
      geom_point(size = 0.5) +
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = variable_label_display), alpha = 0.2, color = NA) +
      geom_rug(data = numeric_rug_df, aes(x = x), inherit.aes = FALSE, sides = "b", alpha = 0.1) +
      facet_wrap(~ variable_label_display, scales = "free_x") +
      labs(
        # title = paste("Partial Dependence with", k, "fold CV"),
           x = "Predictor value",
           y = ifelse(is_ALE_plot, paste("ALE", outcome_label), paste("Predicted", outcome_label))) +
      theme_minimal() +
      guides(color = "none", fill = "none")  +
      theme(
        strip.text = element_text(size = 12)   # 🔹 facet label font size
      )
  }

  # Factor PDPs
  # Determine number of levels in the factor for color scale
  n_levels <- length(unique(factor_df$x))
  if (nrow(factor_df) > 0) {
    # Map human-readable labels
    factor_df$variable_label_display <- sapply(factor_df$variable, function(var) {
      if (!is.null(feature_labels[var]) && !is.na(feature_labels[var])) feature_labels[var] else var
    })


    # Order factor variables
    if (is.null(numeric_var_order) & !is.null(variables)) {
      factor_levels <- feature_labels[variables]
    }
    else if (!is.null(numeric_var_order)) {
      ordered_factor_vars <- numeric_var_order[numeric_var_order %in% factor_df$variable]
      extra_vars <- setdiff(unique(factor_df$variable), ordered_factor_vars)
      ordered_factor_vars <- c(ordered_factor_vars, extra_vars)

      # Apply human-readable labels
      factor_levels <- sapply(ordered_factor_vars, function(var) {
        if (!is.null(feature_labels[var]) && !is.na(feature_labels[var])) feature_labels[var] else var
      })
    } else {
      factor_levels <- unique(factor_df$variable_label_display)
    }

    factor_df$variable_label_display <- factor(factor_df$variable_label_display, levels = factor_levels)

    # Determine number of levels in the factor for color scale
    n_levels <- length(unique(factor_df$x))

    p <- ggplot(factor_df, aes(x = factor(x), y = yhat, fill = x)) +
      geom_boxplot(outlier.size = 0.8) +
      facet_wrap(~ variable_label_display, scales = "free_x") +
      labs(
        # title = paste("Partial Dependence with", k, "fold CV"),
           x = NULL,
          y = ifelse(is_ALE_plot, paste("ALE", outcome_label), paste("Predicted", outcome_label))) +
      theme_minimal() +
      theme(legend.position = "none")

    if (n_levels <= 8) {
      p <- p + scale_fill_brewer(palette = "Set2")
    } else {
      p <- p + scale_fill_viridis_d(option = "turbo")
    }

    plots$factor <- p
  }
  return(plots)
}



# Conditional PDP  ----
# Compute PDP results (slow) 

compute_conditional_pdp <- function(model_list, X, group_var,
                                    facet_var1 = NULL, facet_var2 = NULL,
                                    pred_var, grid.resolution = 50) {
  # Convert to factors and preserve levels
  group_var <- factor(group_var, levels = levels(group_var) %||% unique(group_var))
  facet_var1 <- if (!is.null(facet_var1)) factor(facet_var1, levels = levels(facet_var1) %||% unique(facet_var1)) else NULL
  facet_var2 <- if (!is.null(facet_var2)) factor(facet_var2, levels = levels(facet_var2) %||% unique(facet_var2)) else NULL
  
  group_levels <- levels(group_var)
  facet1_levels <- if (!is.null(facet_var1)) levels(facet_var1) else NULL
  facet2_levels <- if (!is.null(facet_var2)) levels(facet_var2) else NULL
  
  k <- length(model_list)
  X_mat <- model.matrix(~ . - 1, data = X)
  
  conditional_partial_list <- list()
  
  for (i in seq_len(k)) {
    fit <- model_list[[i]]
    
    for (gr in group_levels) {
      for (fc1 in facet1_levels %||% NA) {
        for (fc2 in facet2_levels %||% NA) {
          
          cond1 <- if (is.null(facet_var1)) {
            rep(TRUE, length(group_var))  # all rows match
          } else {
            facet_var1 == fc1
          }
          
          cond2 <- if (is.null(facet_var2)) {
            rep(TRUE, length(group_var))
          } else {
            facet_var2 == fc2
          }
          
          idx <- which(group_var == gr & cond1 & cond2)
          
          if (length(idx) == 0) next
          message(sprintf("Fold %d / %d | Group: %s | Facet1: %s | Facet2: %s | Rows: %d",
                          i, k, gr,
                          ifelse(is.null(fc1), "NA", fc1),
                          ifelse(is.null(fc2), "NA", fc2),
                          length(idx)))
          
          pd <- pdp::partial(
            object = fit,
            pred.var = pred_var,
            train = as.data.frame(X_mat[idx, , drop = FALSE]),
            pred.grid = get_partial_grid(X[idx, , drop = FALSE], pred_var, grid.resolution)
          )
          
          pd$group_level <- factor(gr, levels = group_levels)
          if (!is.null(facet_var1)) pd$facet_var1 <- facet_var1[idx][1]
          if (!is.null(facet_var2)) pd$facet_var2 <- facet_var2[idx][1]
          pd$fold <- i
          
          conditional_partial_list[[length(conditional_partial_list) + 1]] <- pd
        }
      }
    }
  }
  
  conditional_all <- bind_rows(conditional_partial_list)
  
  # Summarize across folds
  summary_df <- conditional_all %>%
    group_by_at(vars(all_of(pred_var), group_level, facet_var1, facet_var2)) %>%
    summarise(
      yhat_mean = mean(yhat, na.rm = TRUE),
      yhat_sd   = sd(yhat, na.rm = TRUE),
      n_folds   = n(),
      yhat_se   = yhat_sd / sqrt(n_folds),
      yhat_lower = yhat_mean - 1.96 * yhat_se,
      yhat_upper = yhat_mean + 1.96 * yhat_se,
      .groups = "drop"
    )
  
  return(summary_df)
}


# Plot PDP results (fast) 
plot_conditional_pdp <- function(summary_df, X, group_var,
                                 facet_var1 = NULL, facet_var2 = NULL,
                                 pred_var, 
                                 group_label = "Group",
                                 facet_label1 = "Facet 1",
                                 facet_label2 = "Facet 2",
                                 x_lim_truncation = TRUE, 
                                 lower_lim_truncation = 0.005,
                                 upper_lim_truncation = 0.995,
                                 switch_facets = FALSE) {
  
  # Prepare rug data
  rug_df <- X %>%
    mutate(
      group_level = group_var,
      facet_var1  = facet_var1,
      facet_var2  = facet_var2
    ) %>%
    select(all_of(pred_var), group_level, facet_var1, facet_var2)
  
  # 99% range truncation
  x_lim_lower <- quantile(X[[pred_var]], lower_lim_truncation, na.rm = TRUE)
  x_lim_upper <- quantile(X[[pred_var]], upper_lim_truncation, na.rm = TRUE)
  
  if(x_lim_truncation) {
    rug_df <- rug_df %>%
      filter(.data[[pred_var]] >= x_lim_lower & .data[[pred_var]] <= x_lim_upper)
    summary_df <- summary_df %>% 
      filter(.data[[pred_var]] >= x_lim_lower & .data[[pred_var]] <= x_lim_upper)
  }
  
  # Custom facet label function
  facet_labeller <- function(labels) {
    if ("facet_var1" %in% names(labels)) {
      labels$facet_var1 <- paste0(facet_label1, ": ", labels$facet_var1)
    }
    if ("facet_var2" %in% names(labels)) {
      labels$facet_var2 <- paste0(facet_label2, ": ", labels$facet_var2)
    }
    labels
  }
  
  # Base plot
  p <- ggplot(summary_df, aes(x = !!sym(pred_var), y = yhat_mean, color = group_level)) +
    geom_line() +
    geom_ribbon(aes(ymin = yhat_lower, ymax = yhat_upper, fill = group_level),
                alpha = 0.2, color = NA) +
    geom_rug(data = rug_df, aes(x = !!sym(pred_var), color = group_level),
             inherit.aes = FALSE, sides = "b", alpha = 0.2) +
    labs(
      x = formula_var_name_to_english[[pred_var]],
      y = "Predicted Lactate (mmol/L)",
      color = group_label,
      fill  = group_label
    ) +
    guides(color = guide_legend(override.aes = list(fill = NA))) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Facet handling with switch_facets
  if (!is.null(facet_var1) & !is.null(facet_var2)) {
    if (switch_facets) {
      p <- p + facet_grid(facet_var2 ~ facet_var1, labeller = facet_labeller, drop = FALSE)
    } else {
      p <- p + facet_grid(facet_var1 ~ facet_var2, labeller = facet_labeller, drop = FALSE)
    }
  } else if (!is.null(facet_var1)) {
    p <- p + facet_wrap(~facet_var1, labeller = facet_labeller, drop = FALSE)
  } else if (!is.null(facet_var2)) {
    p <- p + facet_wrap(~facet_var2, labeller = facet_labeller, drop = FALSE)
  }
  
  return(p)
}


# Group labels ----
get_labels_from_breaks <- function(breaks, include_lower = FALSE) {
  n <- length(breaks)
  labels <- character(n - 1)
  
  # Determine max decimals
  get_decimals <- function(x) {
    if (is.infinite(x)) return(0)
    x_str <- sub(".*\\.", "", format(x, scientific = FALSE, trim = TRUE))
    ifelse(grepl("\\.", format(x, scientific = FALSE)), nchar(x_str), 0)
  }
  
  max_dec <- max(sapply(breaks[is.finite(breaks)], get_decimals))
  step <- 10^(-(max_dec + 1))  # one more decimal place
  
  format_num <- function(x) {
    # Format with max_dec + 1 decimal places, then remove trailing zeros if needed
    format(round(x, max_dec + 1), nsmall = max_dec + 1, trim = TRUE, scientific = FALSE)
  }
  
  for (i in seq_len(n - 1)) {
    lower <- breaks[i]
    upper <- breaks[i + 1]
    
    lower_fmt <- format_num(lower)
    
    if (is.infinite(upper)) {
      labels[i] <- paste0("≥", lower_fmt)
    } else if (i == 1 && !include_lower) {
      labels[i] <- paste0("<", format_num(upper))
    } else {
      upper_adj <- upper - step
      labels[i] <- paste0(lower_fmt, "–", format_num(upper_adj))
    }
  }
  
  labels
}

