pacman::p_load(lme4, umap, dbscan, randomForestSRC, Rtsne)

get_kmeans_cluster_df <- function(lactate_matrix, k, remap_clusters = TRUE) {
  set.seed(123)
  kmeans_result <- kmeans(lactate_matrix, centers = k, nstart = 20)
  
  cluster_df <- tibble(
    surg_id = surg_ids,
    cluster = kmeans_result$cluster
  )
  
  # Change cluster numbers to be in ascending order of average lactate values 
  if(remap_clusters) {
    # Compute mean lactate per cluster
    cluster_means <- cluster_df %>%
      left_join(
        tibble(surg_id = surg_ids, lactate = rowMeans(lactate_matrix, na.rm = TRUE)),
        by = "surg_id"
      ) %>%
      group_by(cluster) %>%
      summarize(mean_lactate = mean(lactate, na.rm = TRUE), .groups = "drop") %>%
      arrange(mean_lactate) %>%
      mutate(new_cluster = row_number())  # 1 = lowest mean, k = highest mean
    
    # Remap cluster IDs
    cluster_df <- cluster_df %>%
      left_join(cluster_means %>% select(cluster, new_cluster), by = "cluster") %>%
      mutate(cluster = new_cluster) %>%
      select(-new_cluster)
  }
  
  return(cluster_df)
}

get_pressor_table <- function(cluster_df) {
  
  # ---- Individual pressor tables ----
  med_dose_sum_summary_list <- purrr::map(names(med_cols), function(med_name) {
    df <- df_onpump_list[[med_name]]
    dose_var <- med_cols[[med_name]]
    if (!(dose_var %in% names(df))) return(NULL)
    
    df_sum <- df %>%
      group_by(surg_id) %>%
      summarize(
        !!dose_var := sum(.data[[dose_var]], na.rm = TRUE),
        cpb_duration = first(cpb_duration),
        .groups = "drop"
      )
    
    df_merged <- df_sum %>%
      right_join(cluster_df, by = "surg_id") %>%
      left_join(cldta_onpump_baseline_df %>% select(surg_id), by = "surg_id") %>%
      mutate(
        !!paste0(dose_var, "_norm")     := .data[[dose_var]] / (60 * .data[["cpb_duration"]]),
        !!paste0(dose_var, "_std")      := standardize_dose(med_name, .data[[dose_var]]),
        !!paste0(dose_var, "_norm_std") := standardize_dose(med_name, .data[[dose_var]]) / (60 * .data[["cpb_duration"]])
      )
    
    # labels
    dose_label          <- var_name_to_english_with_units[[dose_var]]
    dose_norm_var       <- paste0(dose_var, "_norm")
    dose_norm_label     <- paste0(dose_label, " per min")
    
    df_merged %>%
      select(cluster, all_of(c(dose_var, dose_norm_var))) %>%
      add_n(location = "column") %>%
      tbl_summary(
        by = cluster,
        statistic = all_continuous() ~ "{median} ({p15}, {p85})",
        label = setNames(
          list(dose_label, dose_norm_label),
          c(dose_var, dose_norm_var)
        ),
        missing = "no"
      ) %>%
      add_p() %>%
      bold_labels()
  })
  
  med_dose_sum_summary_list_clean <- med_dose_sum_summary_list[!sapply(med_dose_sum_summary_list, is.null)]
  
  # ---- Composite summary across all pressors in NE equivalents ----
  standardized_dose_list <- purrr::map(names(med_cols), function(med_name) {
    df <- df_onpump_list[[med_name]]
    dose_var <- med_cols[[med_name]]
    if (!(dose_var %in% names(df))) return(NULL)
    
    df %>%
      group_by(surg_id) %>%
      summarize(
        std_dose = sum(standardize_dose(med_name, .data[[dose_var]]), na.rm = TRUE),
        cpb_duration = first(cpb_duration),
        .groups = "drop"
      )
  })
  
  standardized_dose_df <- bind_rows(standardized_dose_list) %>%
    group_by(surg_id) %>%
    summarize(
      total_std_dose = sum(std_dose, na.rm = TRUE),
      cpb_duration = first(cpb_duration),
      .groups = "drop"
    )
  
  pressor_summary_df <- cluster_df %>%
    left_join(standardized_dose_df, by = "surg_id") %>%
    mutate(std_dose_per_min = total_std_dose / (60 * cpb_duration))
  
  pressor_composite_table <- pressor_summary_df %>%
    select(cluster, total_std_dose, std_dose_per_min) %>%
    tbl_summary(
      by = cluster,
      statistic = all_continuous() ~ "{median} ({p15}, {p85})",
      label = list(
        total_std_dose ~ "Total NE-equivalent dose",
        std_dose_per_min ~ "NE-equivalent dose per min"
      ),
      missing = "no"
    ) %>%
    add_n(location = "column") %>%
    add_p() %>%
    bold_labels()  %>%
    adjust_tbl_p(method = "fdr")
  
  # ---- Stack all tables together ----
  pressor_table <- tbl_stack(
    c(med_dose_sum_summary_list_clean, list(Composite = pressor_composite_table))
  )
  
  return(pressor_table)
}

get_pressor_table_mcg_hr <- function(cluster_df) {
  
  all_phen_surg_id <- df_onpump_list[["med_phen"]] %>% 
    filter(neobolus == 1) %>% 
    pull(surg_id) %>% 
    unique()
  
  
  # ---- Individual pressor tables ----
  med_dose_sum_summary_list <- purrr::map(names(med_cols_mcg), function(med_name) {
    df <- df_onpump_list[[med_name]]
    dose_var <- med_cols_mcg[[med_name]]
    if (!(dose_var %in% names(df))) return(NULL)
    
    df_sum <- df %>%
      group_by(surg_id) %>%
      summarize(
        !!dose_var := sum(.data[[dose_var]], na.rm = TRUE),
        cpb_duration = first(cpb_duration),
        .groups = "drop"
      )
    
    df_merged <- df_sum %>%
      right_join(cluster_df, by = "surg_id") %>%
      left_join(cldta_onpump_baseline_df %>% select(surg_id), by = "surg_id") %>%
      mutate(
        !!paste0(dose_var, "_norm")     := .data[[dose_var]] / (.data[["cpb_duration"]]),
        !!paste0(dose_var, "_std")      := standardize_dose(med_name, .data[[dose_var]]),
        !!paste0(dose_var, "_norm_std") := standardize_dose(med_name, .data[[dose_var]]) / (.data[["cpb_duration"]])
      )
    
    dose_label          <- var_name_to_english_with_units[[dose_var]]
    dose_norm_var       <- paste0(dose_var, "_norm")
    dose_norm_label     <- paste0(dose_label, " per hour")
    
    df_merged %>%
      select(cluster, all_of(c(dose_var, dose_norm_var))) %>%
      tbl_summary(
        by = cluster,
        statistic = all_continuous() ~ "{median} ({p15}, {p85})",
        label = setNames(
          list(dose_label, dose_norm_label),
          c(dose_var, dose_norm_var)
        ),
        missing = "no"
      ) %>%
      add_n(location = "column") %>%
      add_p() %>%
      bold_labels()  %>%
      adjust_tbl_p(method = "fdr")
  })
  
  med_dose_sum_summary_list_clean <- med_dose_sum_summary_list[!sapply(med_dose_sum_summary_list, is.null)]
  
  # ---- Composite summary across all pressors in NE equivalents ----
  standardized_dose_list <- purrr::map(names(med_cols_mcg), function(med_name) {
    df <- df_onpump_list[[med_name]]
    dose_var <- med_cols_mcg[[med_name]]
    if (!(dose_var %in% names(df))) return(NULL)
    
    df %>%
      group_by(surg_id) %>%
      summarize(
        std_dose = sum(standardize_dose(med_name, .data[[dose_var]]), na.rm = TRUE),
        cpb_duration = first(cpb_duration),
        .groups = "drop"
      )
  })
  
  standardized_dose_df <- bind_rows(standardized_dose_list) %>%
    filter(surg_id %in% all_phen_surg_id) %>% 
    group_by(surg_id) %>%
    summarize(
      total_std_dose = sum(std_dose, na.rm = TRUE),
      cpb_duration = first(cpb_duration),
      .groups = "drop"
    )
  
  pressor_summary_df <- cluster_df %>%
    left_join(standardized_dose_df, by = "surg_id") %>%
    mutate(std_dose_per_min = total_std_dose / (cpb_duration))
  
  pressor_composite_table <- pressor_summary_df %>%
    select(cluster, total_std_dose, std_dose_per_min) %>%
    tbl_summary(
      by = cluster,
      statistic = all_continuous() ~ "{median} ({p15}, {p85})",
      label = list(
        total_std_dose ~ "Total NE-equivalent dose",
        std_dose_per_min ~ "NE-equivalent dose per hour"
      ),
      missing = "no"
    ) %>%
    add_n(location = "column") %>%
    add_p() %>%
    bold_labels()
  
  # ---- Stack all tables together ----
  pressor_table <- tbl_stack(
    c(med_dose_sum_summary_list_clean, list(Composite = pressor_composite_table))
  )
  
  return(pressor_table)
}

adjust_tbl_p <- function(tbl, method = "fdr") {
  tbl %>%
    modify_table_body(
      ~ .x %>%
        dplyr::mutate(
          p.value = ifelse(!is.na(p.value),
                           p.adjust(p.value, method = method),
                           NA)
        )
    ) %>%
    modify_header(p.value ~ "**p-value (adjusted)**")
}


# cluster_df must have column cluster 
get_tables_comparing_clusters <- function(
    cluster_df,
    vals_to_get_summary_for = c(
      "lactate", "do2_interp", "hgb_io", "gluc_io", "ci", "map",
      "ohgb_io", "pao2", "svo2", "pvo2", "cao2_interp",
      "cvo2_interp", "o2_consume_interp"
    ),
    df_long = NULL, 
    cluster_order = NULL,
    tables_to_return = c("baseline", "surgery", "longitudinal", "pressor", "complications") 
) {
  # Reorder clusters if order is provided
  if (!is.null(cluster_order)) {
    cluster_df <- cluster_df %>%
      mutate(cluster = factor(cluster, levels = cluster_order))
  }
  
  out <- list()  # Store results
  
  # Baseline characteristics ----
  if ("baseline" %in% tables_to_return) {
    cldta_onpump_with_cluster <- cldta_onpump_baseline_df %>%
      left_join(cluster_df, by = "surg_id")
    
    summary_df <- cldta_onpump_with_cluster %>%
      select(cluster, age, bmi, sex, race_comb, cpb_duration, wbc_pr, starts_with("hx_"))
    
    baseline_table <- summary_df %>%
      tbl_summary(
        by = cluster,
        statistic = list(
          age ~ "{mean} ({sd})",
          bmi ~ "{mean} ({sd})",
          cpb_duration ~ "{median} ({p15}, {p85})"
        ),
        label = baseline_var_to_english,
        missing = "no"
      ) %>%
      add_p() %>%
      bold_labels()  %>%
      adjust_tbl_p(method = "fdr")
    
    out$baseline_table <- baseline_table
  }
  
  # Surgery type ----
  if ("surgery" %in% tables_to_return) {
    surgery_df <- cldta_onpump_with_cluster %>%
      select(cluster, all_of(names(sp_labels)))
    
    surgery_table <- surgery_df %>%
      tbl_summary(
        by = cluster,
        label = sp_labels,
        missing = "no"
      ) %>%
      add_p() %>%  
      adjust_tbl_p(method = "fdr") %>% 
      bold_labels()
    
    out$surgery_table <- surgery_table
  }
  
  # Longitudinal variables summaries ----
  if ("longitudinal" %in% tables_to_return) {
    
    # Map of dataset name -> actual column names in that dataset
    df_colname_map <- list(
      hct = "hct_io",
      hgb = "hgb_io",
      ohgb = "ohgb_io",
      glucose = "glucose_io" # add others as needed
    )
    
    cpb_var_list <- map2(df_onpump_list, names(df_onpump_list), function(df, df_name) {
      
      # Use the actual column names for this dataset
      vars_to_use <- intersect(names(df), df_colname_map[[df_name]])
      if (length(vars_to_use) == 0) return(NULL)
      
      # Keep your original averaging logic
      df_avg <- df %>%
        group_by(surg_id) %>%
        summarize(
          across(all_of(vars_to_use), ~ mean(.x, na.rm = TRUE)),
          .groups = "drop"
        ) %>%
        left_join(cluster_df, by = "surg_id")
      
      # Compute median, p15, p85 per cluster
      df_summary <- df_avg %>%
        group_by(cluster) %>%
        summarize(across(all_of(vars_to_use), list(
          median = ~ median(.x, na.rm = TRUE),
          p15    = ~ quantile(.x, 0.15, na.rm = TRUE),
          p85    = ~ quantile(.x, 0.85, na.rm = TRUE)
        ), .names = "{.col}_{.fn}"), .groups = "drop")
      
      # Preformat median/p15/p85 to 3 significant figures
      df_formatted <- df_summary
      for (var in vars_to_use) {
        df_formatted[[var]] <- paste0(
          style_sigfig(df_summary[[paste0(var, "_median")]], 3), " (",
          style_sigfig(df_summary[[paste0(var, "_p15")]], 3), ", ",
          style_sigfig(df_summary[[paste0(var, "_p85")]], 3), ")"
        )
      }
      
      df_formatted <- df_formatted %>%
        select(cluster, all_of(vars_to_use))
      
      # Return gtsummary table
      tbl_summary(
        df_formatted,
        by = cluster,
        statistic = all_continuous() ~ "{value}",
        label = as.list(var_name_to_english_with_units),
        missing = "no"
      ) %>%
        add_n(location = "column") %>%
        add_p() %>%
        bold_labels() %>%
        adjust_tbl_p(method = "fdr")
    })
    
    cpb_var_list_clean <- cpb_var_list[!sapply(cpb_var_list, is.null)]
    
    # Lactate table (same 3 sig figs formatting)
    lactate_df <- cldta_onpump_lactate_mean %>%
      group_by(surg_id) %>%
      summarize(mean_lactate = mean(mean_lactate, na.rm = TRUE), .groups = "drop") %>%
      left_join(cluster_df, by = "surg_id")
    
    lactate_summary <- lactate_df %>%
      group_by(cluster) %>%
      summarize(
        median = median(mean_lactate, na.rm = TRUE),
        p15    = quantile(mean_lactate, 0.15, na.rm = TRUE),
        p85    = quantile(mean_lactate, 0.85, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(`Mean Lactate (mmol/L)` = paste0(
        style_sigfig(median, 3), " (",
        style_sigfig(p15, 3), ", ",
        style_sigfig(p85, 3), ")"
      )) %>%
      select(cluster, `Mean Lactate (mmol/L)`)
    
    lactate_table <- tbl_summary(
      lactate_summary,
      by = cluster,
      statistic = all_continuous() ~ "{value}",
      label = list(`Mean Lactate (mmol/L)` ~ "Mean Lactate (mmol/L)"),
      missing = "no"
    ) %>%
      add_n(location = "column") %>%
      add_p() %>%
      bold_labels() %>%
      adjust_tbl_p(method = "fdr")
    
    # Stack tables
    longitudinal_table <- tbl_stack(c(list(`Mean Lactate (mmol/L)` = lactate_table), cpb_var_list_clean))
    out$longitudinal_table <- longitudinal_table
  }

  # Pressors ----
  if ("pressor" %in% tables_to_return) {
    out$pressor_table <- get_pressor_table_mcg_hr(cluster_df) 
  }
  
  if ("complications" %in% tables_to_return) {
    
    po_labs <- c("creat_po", "ast_po", "alt_po")
    
    immediate_postop_labs <- purrr::map(
      po_labs,
      ~ df_list[[.x]] %>%
        filter(iv_evtof1 < 1) %>%
        arrange(iv_evtof1) %>%
        group_by(surg_id) %>%
        slice(1) %>%
        ungroup() %>%
        select(surg_id, !!paste0(.x, "_immediate") := all_of(.x))
    ) %>%
      purrr::reduce(full_join, by = "surg_id")
    
    peak_48h_labs <- purrr::map(
      po_labs,
      ~ df_list[[.x]] %>%
        filter(iv_evtof1 > 0 & iv_evtof1 <= 48) %>%
        group_by(surg_id) %>%
        slice_max(.data[[.x]], n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        select(surg_id, !!paste0(.x, "_peak48h") := all_of(.x))
    ) %>%
      purrr::reduce(full_join, by = "surg_id")
    
    complications_with_labs_df <- cldta_onpump_baseline_df %>%
      left_join(cluster_df, by = "surg_id") %>%
      select(surg_id, cluster, all_of(names(complication_labels))) %>%
      left_join(immediate_postop_labs, by = "surg_id") %>%
      left_join(peak_48h_labs, by = "surg_id")
    
    lab_labels <- c(
      creat_po_immediate = "Creatinine (immediate post-op)",
      ast_po_immediate   = "AST (immediate post-op)",
      alt_po_immediate   = "ALT (immediate post-op)",
      creat_po_peak48h   = "Peak creatinine (0–48h)",
      ast_po_peak48h     = "Peak AST (0–48h)",
      alt_po_peak48h     = "Peak ALT (0–48h)"
    )
    
    all_labels <- c(complication_labels, lab_labels)
    
    complications_table <- complications_with_labs_df %>%
      select(-surg_id) %>%
      tbl_summary(
        by = cluster,
        statistic = all_continuous() ~ "{median} ({p15}, {p85})",
        label = all_labels,
        missing = "no"
      ) %>%
      add_p(
        test = list(
          c(
            "COtAFib", "CPVntLng", "bldiprod", "bldirbc",
            "bldicry", "bldiffp", "bldiplt", "dead30d"
          ) ~ "chisq.test",
          c(
            "CRenFail", "CRenDial", "CNStrokP",
            "hdeath", "CISeptic", "CIStDeep"
          ) ~ "fisher.test",
          all_continuous() ~ "kruskal.test"
        )
      ) %>%
      bold_labels() %>%
      italicize_levels() %>%
      adjust_tbl_p(method = "fdr")
    
    out$complications_table <- complications_table
  }
  
  # if ("complications" %in% tables_to_return) {
  #   complications_df <- cldta_onpump_baseline_df %>%
  #     left_join(cluster_df, by = "surg_id") %>%
  #     select(cluster, all_of(names(complication_labels)))
  #   
  #   complications_table <- complications_df %>%
  #     tbl_summary(
  #       by = cluster,
  #       statistic = all_continuous() ~ "{median} ({p15}, {p85})",
  #       label = complication_labels,
  #       missing = "no"
  #     ) %>%
  #     add_p(
  #       test = list(
  #         c("COtAFib", "CPVntLng", "bldiprod", "bldirbc", "bldicry", "bldiffp", "bldiplt", "dead30d") ~ "chisq.test",
  #         c("CRenFail", "CRenDial", "CNStrokP",
  #           "hdeath", "CISeptic", "CIStDeep") ~ "fisher.test",
  #         all_continuous() ~ "kruskal.test"
  #       )
  #     ) %>%
  #     bold_labels() %>%
  #     italicize_levels()  %>%
  #     adjust_tbl_p(method = "fdr")
  #   
  #   out$complications_table <- complications_table
  #   
  #   po_labs <- c("creat_po", "ast_po", "alt_po")
  #   
  #   immediate_postop_df <- df_list[[po_col]] %>%
  #     filter(iv_evtof1 < 1) %>%
  #     select(
  #       surg_id,
  #       lab_po = !!sym(po_col)
  #     )
  #   peak_creat_df <- df_list[["creat_po"]] %>% 
  #     filter(iv_evtof1 > 0 & iv_evtof1 <= 48) %>% 
  #     group_by(surg_id) %>% 
  #     slice_max(creat_po, n = 1, with_ties = FALSE) %>% 
  #     summarize(
  #       peak_creat = creat_po,
  #       iv_evtof1  = iv_evtof1,
  #       .groups = "drop"
  #     )
  #   
  # }
  
  return(out)
}


get_spaghetti_lactate <- function(cluster_df, cluster_order = NULL, standardize_y = TRUE, x_lim = NULL, y_lim = NULL,
                                  xlab = "Time since CPB start (hours)", ylab = "Lactate (mmol/L)") {
  cldta_onpump_clustered <- cldta_onpump %>%
    inner_join(cluster_df %>% select(surg_id, cluster), by = "surg_id")
  
  # Count unique patients per cluster
  cluster_counts <- cldta_onpump_clustered %>%
    distinct(surg_id, cluster) %>%
    count(cluster, name = "N")
  
  # Add labeled cluster info
  cldta_onpump_clustered <- cldta_onpump_clustered %>%
    left_join(cluster_counts, by = "cluster") %>%
    mutate(cluster_label = paste0("Cluster ", cluster, " (N=", N, ")"))
  
  # Apply cluster order if provided
  if (!is.null(cluster_order)) {
    cldta_onpump_clustered <- cldta_onpump_clustered %>%
      mutate(cluster_label = factor(cluster_label,
                                    levels = paste0("Cluster ", cluster_order, " (N=",
                                                    cluster_counts$N[match(cluster_order, cluster_counts$cluster)], ")")))
  }
  
  # Base plot
  p <- ggplot(cldta_onpump_clustered, aes(x = iv_evth, y = lactate, group = surg_id)) +
    geom_line(alpha = 0.3, linewidth = 0.2) +
    facet_wrap(~ cluster_label, scales = ifelse(standardize_y, "fixed", "free_y"), ncol = 4) +
    labs(
      title = "Lactate vs Time by CPB Cluster",
      x = xlab,
      y = ylab
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 12) 
    )
  
  if (!is.null(y_lim) || !is.null(x_lim)) {
    p <- p + coord_cartesian(
      ylim = if (!is.null(y_lim)) c(0, y_lim) else NULL,
      xlim = if (!is.null(x_lim)) c(0, x_lim) else NULL
    )
  }
  
  return(p)
}


get_spaghetti_lactate_pre30 <- function(cldta_filtered, 
                                        cluster_df, 
                                        title = "Pre-CPB Lactate trends by Cluster", 
                                        cluster_order = NULL, 
                                        standardize_y = TRUE, 
                                        x_lim = -3, y_lim = NULL) {
  cldta_precpb_clustered <- cldta_filtered %>% 
    inner_join(cluster_df %>% select(surg_id, cluster), by = "surg_id")

  # Count unique patients per cluster
  cluster_counts <- cldta_precpb_clustered %>%
    distinct(surg_id, cluster) %>%
    count(cluster, name = "N")
  
  # Add labeled cluster info
  cldta_precpb_clustered <- cldta_precpb_clustered %>%
    left_join(cluster_counts, by = "cluster") %>%
    mutate(cluster_label = paste0("Cluster ", cluster, " (N=", N, ")"))
  
  # Apply cluster order if provided
  if (!is.null(cluster_order)) {
    cldta_precpb_clustered <- cldta_precpb_clustered %>%
      mutate(cluster_label = factor(cluster_label,
                                    levels = paste0("Cluster ", cluster_order, " (N=",
                                                    cluster_counts$N[match(cluster_order, cluster_counts$cluster)], ")")))
  }
  

  p <- ggplot(cldta_precpb_clustered, aes(x = iv_evth, y = lactate, group = surg_id)) +
    geom_line(alpha = 0.3, linewidth = 0.2) +
    facet_wrap(~ cluster_label, scales = ifelse(standardize_y, "fixed", "free_y"), ncol = 4) +
    labs(
      title = title,
      x = "Time since CPB start (hours)",
      y = "Lactate (mmol/L)"
    ) +
    theme_minimal() +
    theme(strip.text = element_text(size = 12)) +
    coord_cartesian(xlim = c(x_lim, 0), ylim = if (!is.null(y_lim)) c(0, y_lim) else NULL)
  
  return(p)
}

# Dimentionality reduction plots ---- 
get_pca_plot <- function(mat, cluster_vector = NULL) {
  pca_result <- prcomp(mat, scale. = TRUE)
  df <- as_tibble(pca_result$x[, 1:2]) %>%
    rename(PC1 = 1, PC2 = 2)
  if (!is.null(cluster_vector)) {
    df <- df %>% mutate(cluster = factor(cluster_vector))
  }
  
  ggplot(df, aes(x = PC1, y = PC2, color = cluster)) +
    geom_point(size = 2, alpha = 0.7) +
    labs(title = "PCA Plot", x = "PC1", y = "PC2") +
    theme_minimal()
}

get_umap_plot <- function(mat, cluster_vector = NULL, umap_config = umap.defaults) {
  umap_result <- umap(mat, config = umap_config)
  df <- as_tibble(umap_result$layout) %>%
    rename(UMAP1 = V1, UMAP2 = V2)
  if (!is.null(cluster_vector)) {
    df <- df %>% mutate(cluster = factor(cluster_vector))
  }
  
  ggplot(df, aes(x = UMAP1, y = UMAP2, color = cluster)) +
    geom_point(size = 2, alpha = 0.6) +
    labs(title = "UMAP Plot", x = "UMAP1", y = "UMAP2") +
    theme_minimal()
}

get_tsne_plot <- function(mat, cluster_vector = NULL, perplexity = 30, max_iter = 1000, seed = 42) {
  set.seed(seed)
  tsne_result <- Rtsne(mat, dims = 2, perplexity = perplexity, max_iter = max_iter, verbose = FALSE)
  df <- as_tibble(tsne_result$Y) %>%
    rename(tSNE1 = V1, tSNE2 = V2)
  if (!is.null(cluster_vector)) {
    df <- df %>% mutate(cluster = factor(cluster_vector))
  }
  
  ggplot(df, aes(x = tSNE1, y = tSNE2, color = cluster)) +
    geom_point(size = 2, alpha = 0.6) +
    labs(title = "t-SNE Plot", x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
    theme_minimal()
}


get_tsne_plot <- function(mat, cluster_vector = NULL,
                          perplexity = 30, max_iter = 1000, seed = 42) {
  # Check input
  if (!is.matrix(mat) && !is.data.frame(mat)) {
    stop("Input must be a matrix or data frame.")
  }
  mat <- as.matrix(mat)
  
  # Identify and remove duplicate rows
  dedup_idx <- !duplicated(mat)
  mat_unique <- mat[dedup_idx, ]
  
  # Subset cluster vector if provided
  if (!is.null(cluster_vector)) {
    if (length(cluster_vector) != nrow(mat)) {
      stop("cluster_vector must be the same length as the number of rows in mat.")
    }
    cluster_vector <- cluster_vector[dedup_idx]
  }
  
  # Run t-SNE on unique data
  set.seed(seed)
  tsne_result <- Rtsne(mat_unique, dims = 2, perplexity = perplexity,
                       max_iter = max_iter, verbose = TRUE)
  
  # Create t-SNE plot dataframe
  tsne_df <- as_tibble(tsne_result$Y) %>%
    rename(tSNE1 = V1, tSNE2 = V2)
  
  if (!is.null(cluster_vector)) {
    tsne_df <- tsne_df %>% mutate(cluster = factor(cluster_vector))
  }
  
  # Plot
  p <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2)) +
    geom_point(size = 2, alpha = 0.7,
               aes(color = if (!is.null(cluster_vector)) cluster else NULL)) +
    labs(
      title = "t-SNE on Random Effects (duplicates removed)",
      x = "t-SNE 1", y = "t-SNE 2",
      color = if (!is.null(cluster_vector)) "Cluster" else NULL
    ) +
    theme_minimal()
  
  return(p)
}
