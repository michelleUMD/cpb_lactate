pacman::p_load(plotly, htmlwidgets, ks, akima, roxygen2)

# Save html figure ----
# figure_dir is already defined 
save <- function(filename, plot_obj = NULL, width = 25, height = 15, selfcontained = TRUE) {
  
  if(is.null(plot_obj)) {
    if (!grepl("\\.png$", filename, ignore.case = TRUE)) {
      filename <- paste0(filename, ".png")
    }
    ggsave(
      filename = file.path(figure_dir, filename),
      plot = last_plot(),   # grabs the last ggplot
      width = width,
      height = height,
      dpi = 300
    )
  } else if (inherits(plot_obj, "ggplot")) {
    if (!grepl("\\.png$", filename, ignore.case = TRUE)) {
      filename <- paste0(filename, ".jpg")
    }
    ggsave(
      filename = file.path(figure_dir, filename),
      plot = plot_obj,
      width = width,
      height = height,
      dpi = 300
    )
    message("ggplot object saved to: ", normalizePath(file.path(figure_dir, filename)))
    
    # If an htmlwidget (e.g., plotly, leaflet)
  } else if (inherits(plot_obj, "htmlwidget")) {
    if (!grepl("\\.html?$", filename, ignore.case = TRUE)) {
      filename <- paste0(filename, ".html")
    }
    outfile <- file.path(figure_dir, filename)
    htmlwidgets::saveWidget(plot_obj, file = outfile, selfcontained = selfcontained)
    message("HTML plot saved to: ", normalizePath(outfile))
    
  } else {
    print("ERROR: Not a supported type")
  }
}


# 2D + Color ----

#' Plot 2D Partial Dependence
#'
#' This function computes and plots the 2-dimensional partial dependence of a fitted model
#' with respect to two predictor variables. The resulting plot shows the predicted response
#' as a color gradient across the grid of the two variables.
#'
#' @param fit A fitted model object (e.g., from `xgboost`, `randomForest`, or any model
#'   compatible with `pdp::partial`).
#' @param X A data frame containing the predictor variables used for the model.
#' @param var1 Character. The first predictor variable for the 2D partial dependence plot.
#' @param var2 Character. The second predictor variable for the 2D partial dependence plot.
#' @param quantile_points Integer. Number of grid points per variable (used to create the evaluation grid). Default is 50.
#' @param xlim Numeric vector of length 2. Limits for the x-axis. Default is `NULL` (automatic).
#' @param ylim Numeric vector of length 2. Limits for the y-axis. Default is `NULL` (automatic).
#' @param title Character. Optional plot title. Default is automatically generated from variable names.
#'
#' @return A `ggplot` object showing the 2D partial dependence of the model with respect
#'   to `var1` and `var2`.
#'
#' @details
#' The function first creates a grid of values for `var1` and `var2` using `get_partial_grid()`.
#' It then computes the partial dependence over this grid using `pdp::partial()`. The resulting
#' predicted values are plotted using `ggplot2` with a color gradient representing the predicted response.
#'
#' @export
plot_2d_pdp <- function(fit, X, var1, var2, 
                        quantile_points = 50,
                        xlim = NULL, ylim = NULL,
                        title = NULL) {
  
  # grids for each variable
  grid1 <- get_partial_grid(X, var1, quantile_points = quantile_points)
  grid2 <- get_partial_grid(X, var2, quantile_points = quantile_points)
  
  # combined grid
  grid_2d <- expand.grid(
    grid1[[var1]],
    grid2[[var2]]
  )
  colnames(grid_2d) <- c(var1, var2)
  
  # partial dependence
  pdp_2d <- pdp::partial(
    object = fit,
    pred.var = c(var1, var2),
    train = X,
    pred.grid = grid_2d,
    progress = "none"
  )
  
  # plot
  p <- ggplot(pdp_2d, aes_string(x = var1, y = var2, color = "yhat")) +
    geom_point(size = 2, alpha = 0.7) +    # smaller points, semi-transparent
    scale_color_viridis_c(option = "plasma") +
    theme_minimal() +
    labs(
      title = ifelse(is.null(title), 
                     paste("2D Partial Dependence:", var1, "x", var2), 
                     title),
      x = var1,
      y = var2,
      color = "Predicted Response"
    )
  
  if (!is.null(xlim) | !is.null(ylim)) {
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  }
  
  return(p)
}

plot_2d_ale <- function(fit, X, var1, var2, 
                        xlim = NULL, ylim = NULL,
                        show_counts = TRUE,
                        title = NULL) {
  
  # Custom predict function for xgboost
  pred_fun <- function(model, newdata) {
    newdata_mat <- as.matrix(newdata)
    predict(model, newdata_mat)
  }
  
  # Create Predictor object
  predictor <- Predictor$new(
    model = fit,
    data = X,
    y = NULL,
    predict.fun = pred_fun
  )
  
  # Compute 2D ALE
  ale_eff <- FeatureEffect$new(
    predictor,
    feature = c(var1, var2),
    method = "ale"
  )
  
  ale_df <- ale_eff$results
  
  # Optional: overlay patient counts per bin
  if (show_counts) {
    df_counts <- X %>%
      mutate(
        bin1 = cut(.data[[var1]], breaks = unique(c(ale_df$.left, ale_df$.right)), include.lowest = TRUE),
        bin2 = cut(.data[[var2]], breaks = unique(c(ale_df$.bottom, ale_df$.top)), include.lowest = TRUE)
      ) %>%
      group_by(bin1, bin2) %>%
      summarise(n_patients = n(), .groups = "drop")
    
    ale_df$text <- paste0("Patients: ", rep(df_counts$n_patients, length.out = nrow(ale_df)))
  } else {
    ale_df$text <- NULL
  }
  
  # Plot using the rectangle boundaries
  p <- ggplot(
    ale_df,
    aes(xmin = .left, xmax = .right, ymin = .bottom, ymax = .top, fill = .ale, text = text)
  ) +
    geom_rect(color = NA) +  # ALE heatmap grid
    scale_fill_viridis_c(option = "plasma") +
    theme_minimal() +
    labs(
      title = ifelse(is.null(title), paste("2D ALE:", var1, "x", var2), title),
      x = var1,
      y = var2,
      fill = "ALE"
    )
  
  if (!is.null(xlim) | !is.null(ylim)) {
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  }
  
  return(p)
}

# plot_2d_ale <- function(fit, X, var1, var2, 
#                         # grid_points = 50,
#                         xlim = NULL, ylim = NULL,
#                         show_counts = TRUE,
#                         title = NULL) {
#   
#   # Custom predict function for xgboost
#   pred_fun <- function(model, newdata) {
#     newdata_mat <- as.matrix(newdata)
#     predict(model, newdata_mat)
#   }
#   
#   # Create Predictor object
#   predictor <- Predictor$new(
#     model = fit,
#     data = X,
#     y = NULL,
#     predict.fun = pred_fun
#   )
#   
#   # Compute 2D ALE
#   ale_eff <- FeatureEffect$new(
#     predictor,
#     feature = c(var1, var2),
#     method = "ale"
#     # grid.size = grid_points
#   )
#   
#   ale_df <- ale_eff$results
#   
#   View(ale_df)
#   
#   # Optional: overlay patient counts per bin
#   if (show_counts) {
#     df_counts <- X %>%
#       mutate(
#         bin1 = cut(.data[[var1]], breaks = unique(ale_df[[var1]]), include.lowest = TRUE),
#         bin2 = cut(.data[[var2]], breaks = unique(ale_df[[var2]]), include.lowest = TRUE)
#       ) %>%
#       group_by(bin1, bin2) %>%
#       summarise(n_patients = n(), .groups = "drop")
#     
#     ale_df$text <- paste0("Patients: ", rep(df_counts$n_patients, length.out = nrow(ale_df)))
#   } else {
#     ale_df$text <- NULL
#   }
#   
#   p <- ggplot(ale_df, aes_string(x = var1, y = var2, color = ".ale", text = "text")) +
#     geom_point(size = 2, alpha = 0.8) +   # smaller, semi-transparent points
#     scale_color_viridis_c(option = "plasma") +
#     theme_minimal() +
#     labs(
#       title = ifelse(is.null(title), paste("2D ALE:", var1, "x", var2), title),
#       x = var1,
#       y = var2,
#       color = "ALE"
#     )
#   
#   if (!is.null(xlim) | !is.null(ylim)) {
#     p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
#   }
#   
#   return(p)
# }

plot_2d_ale_k_fold <- function(model_list, var1, var2, X,
                               n_tiles1 = 30, n_tiles2 = 30, 
                               xlim = NULL, ylim = NULL) {
  ale_results <- list()
  
  pred_fun <- function(model, newdata) {
    newdata_mat <- as.matrix(newdata)
    predict(model, newdata_mat)
  }

  for (fold in seq_along(model_list)) {
    predictor <- Predictor$new(model_list[[fold]], data = X, predict.fun = pred_fun)
    ale_eff <- FeatureEffect$new(predictor, feature = c(var1, var2), method = "ale")
    ale_results[[fold]] <- ale_eff$results
  }
  # grid1 <- get_partial_grid(X, var1)[[var1]]
  # grid2 <- get_partial_grid(X, var2)[[var2]]
  
  grid1 <- seq(min(X[[var1]]), max(X[[var1]]), length.out = n_tiles1)
  grid2 <- seq(min(X[[var2]]), max(X[[var2]]), length.out = n_tiles2)
  common_grid <- expand.grid(var1 = grid1, var2 = grid2)
  
  
  interpolated_list <- lapply(ale_results, function(df) {
    # Use bin centers as x,y
    x_center <- (df$.left + df$.right) / 2
    y_center <- (df$.bottom + df$.top) / 2
    z <- df$.ale
    
    interp_out <- with(df, akima::interp(
      x = x_center, y = y_center, z = z,
      xo = grid1, yo = grid2, linear = TRUE, extrap = FALSE
    ))
    
    expand.grid(var1 = interp_out$x, var2 = interp_out$y) %>%
      mutate(ale = as.vector(interp_out$z))
  })
  avg_df <- Reduce(function(x, y) full_join(x, y, by = c("var1","var2")),
                   interpolated_list) %>%
    rowwise() %>%
    mutate(ale_mean = mean(c_across(starts_with("ale")), na.rm = TRUE)) %>% 
    filter(!is.na(ale_mean))
  
  ggplot(avg_df, aes(x = var1, y = var2, fill = ale_mean)) +
    geom_tile() +
    scale_fill_viridis_c(option = "plasma") +
    theme_minimal() +
    labs(title = paste("Mean 2D ALE across folds"),
         x = var1, y = var2, fill = "ALE")
  
  p <- ggplot(avg_df, aes(x = var1, y = var2, fill = ale_mean)) +
    geom_tile() +
    # overlay raw patient data
    geom_density_2d(data = X, aes(x = .data[[var1]], y = .data[[var2]]),
                    inherit.aes = FALSE, color = "white", alpha = 0.25) + 
    scale_fill_viridis_c(option = "plasma") +
    theme_minimal() +
    labs(
      title = paste("Mean 2D ALE across folds"),
      x = formula_var_name_to_english[[var1]],
      y = formula_var_name_to_english[[var2]],
      fill = "ALE"
    ) + theme (legend.position = "bottom")
  
    if (!is.null(xlim) | !is.null(ylim)) {
      p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
    }
  
  print(p)
  return(p)
}

# 3D mesh ----
# Filters using kernel based density 
plot_filtered_pdp <- function(model_data, pdp_4d_summary, density_quantile = 0.1) {
  # predictor columns = everything except those starting with "yhat"
  predictor_cols <- setdiff(names(pdp_4d_summary), grep("^yhat", names(pdp_4d_summary), value = TRUE))
  
  # extract observed variables
  obs <- model_data[, predictor_cols, drop = FALSE]
  
  # fit multivariate KDE on observed data
  kde_fit <- ks::kde(x = obs)
  
  # evaluate density at PDP grid points
  grid_points <- as.matrix(pdp_4d_summary[, predictor_cols, drop = FALSE])
  dens_vals <- predict(kde_fit, x = grid_points)
  
  # filter: keep points above chosen quantile
  cutoff <- quantile(dens_vals, density_quantile)
  valid_mask <- dens_vals > cutoff
  pdp_filtered <- pdp_4d_summary[valid_mask, ]
  
  if (length(predictor_cols) != 3) {
    stop("Currently this plotting supports exactly 3 predictor variables (for 3D).")
  }
  
  p <- plot_ly(
    data = pdp_filtered,
    x = as.formula(paste0("~", predictor_cols[1])),
    y = as.formula(paste0("~", predictor_cols[2])),
    z = as.formula(paste0("~", predictor_cols[3])),
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = 2,
      opacity = 0.6,
      color = ~yhat_mean,
      colorscale = "Viridis",
      colorbar = list(title = "Predicted lactate (mmol/L)")
    )
  ) %>%
    layout(scene = list(
      xaxis = list(title = formula_var_name_to_english[[ predictor_cols[1] ]]),
      yaxis = list(title = formula_var_name_to_english[[ predictor_cols[2] ]]),
      zaxis = list(title = formula_var_name_to_english[[ predictor_cols[3] ]])
    ))
  
  print(p)
  return(p)
}

# Surface ----

get_pdp_truncated <- function(pdp_df, X, x_var, y_var, slider_var, percent_truncation = 0.95) {
  # Central proportion to keep (default 95%)
  lower <- (1 - percent_truncation) / 2
  upper <- 1 - lower
  
  # 1. Compute quantiles for truncation
  x_low <- quantile(X[[x_var]], lower, na.rm = TRUE)
  x_high <- quantile(X[[x_var]], upper, na.rm = TRUE)
  
  y_low <- quantile(X[[y_var]], lower, na.rm = TRUE)
  y_high <- quantile(X[[y_var]], upper, na.rm = TRUE)
  
  slider_low <- quantile(X[[slider_var]], lower, na.rm = TRUE)
  slider_high <- quantile(X[[slider_var]], upper, na.rm = TRUE)
  
  # 2. Filter pdp_df to keep only values within range
  pdp_trunc <- pdp_df %>%
    filter(
      .data[[x_var]] >= x_low & .data[[x_var]] <= x_high,
      .data[[y_var]] >= y_low & .data[[y_var]] <= y_high,
      .data[[slider_var]] >= slider_low & .data[[slider_var]] <= slider_high
    )
  
  return(pdp_trunc)
}



plot_conditional_surface <- function(pdp_df, X, x_var, y_var, z_var, slider_var,
                                     percent_truncation = 0.95) {
  
  # Map axis labels
  x_label <- formula_var_name_to_english[[x_var]] %||% x_var
  y_label <- formula_var_name_to_english[[y_var]] %||% y_var
  z_label <- formula_var_name_to_english[[z_var]] %||% z_var
  slider_label <- formula_var_name_to_english[[slider_var]] %||% slider_var
  
  pdp_df <- get_pdp_truncated(
    pdp_df = pdp_df,
    X = X,
    x_var = x_var,
    y_var = y_var,
    slider_var = slider_var,
    percent_truncation = percent_truncation
  )
  
  x_vals <- sort(unique(pdp_df[[x_var]]))
  y_vals <- sort(unique(pdp_df[[y_var]]))
  slider_vals <- sort(unique(pdp_df[[slider_var]]))
  
  
  
  # Build frames directly from pdp_df
  frames <- lapply(slider_vals, function(sv) {
    df_slice <- pdp_df[pdp_df[[slider_var]] == sv, ]
    
    # Order by x then y exactly as in your data
    df_slice <- df_slice[order(df_slice[[x_var]], df_slice[[y_var]]), ]
    
    # Reshape z-values into matrix
    # x_vals number of rows, so each represents a different y_val
    z_matrix <- t(matrix(df_slice[[z_var]],
                         nrow = length(x_vals),
                         ncol = length(y_vals),
                         byrow = TRUE))
    
    list(name = as.character(round(sv, 2)),
         data = list(list(z = z_matrix)))
  })
  
  # Initial surface (first slider value)
  first_slice <- pdp_df[pdp_df[[slider_var]] == slider_vals[1], ]
  first_slice <- first_slice[order(first_slice[[x_var]], first_slice[[y_var]]), ]
  z_matrix_first <- t(matrix(first_slice[[z_var]],
                             nrow = length(x_vals),
                             ncol = length(y_vals),
                             byrow = TRUE))
  
  p <- plot_ly(
    x = x_vals,
    y = y_vals,
    z = z_matrix_first,
    type = "surface",
    colorscale = "Viridis"
  ) %>%
    layout(
      scene = list(
        xaxis = list(title = x_label, range = range(x_vals)),
        yaxis = list(title = y_label, range = range(y_vals)),
        zaxis = list(title = z_label, range = range(pdp_df[[z_var]]))
      )
    )
  
  # Slider
  p$x$frames <- frames
  steps <- lapply(seq_along(frames), function(i) {
    list(label = as.character(round(slider_vals[i], 2)),
         method = "animate",
         args = list(list(frames[[i]]$name),
                     list(mode = "immediate",
                          frame = list(duration = 0, redraw = TRUE),
                          transition = list(duration = 0))))
  })
  
  p <- p %>% layout(
    sliders = list(list(
      active = 0,
      currentvalue = list(prefix = paste0(slider_label, ": ")),
      steps = steps
    ))
  )
  
  print(p)
  return(p)
}


# Stratified surfaces ----
compute_stratified_surface_pdp <- function(
    X,
    model_list,
    grid_var1,       # e.g. "do2_interp_m2_mean"
    grid_var2,       # e.g. "map_mean" or "std_dose_mg_hr"
    strata_var,       # variable used for cut (e.g. "std_dose_mg_hr")
    strata_breaks,    # numeric breaks for cut()
    pred_name = "yhat_mean"
) {

  grid1 <- get_partial_grid(X, grid_var1)
  grid2 <- get_partial_grid(X, grid_var2)
  
  custom_grid <- tidyr::crossing(
    !!grid_var1 := grid1[[grid_var1]],
    !!grid_var2 := grid2[[grid_var2]]
  )
  

  strata_vals <- cut(
    X[[strata_var]],
    breaks = strata_breaks,
    labels = get_labels_from_breaks(strata_breaks), 
    right = FALSE,
    ordered_result = TRUE
  )
  
  strata_levels <- levels(strata_vals)
  
  # Build matrix version of X for xgb.DMatrix
  X_mat <- x_to_x_mat(X)
  
  # Loop over each group
  pdp_list_strata <- purrr::map(strata_levels, function(strat) {
    idx <- which(strata_vals == strat)
    train_sub <- as.data.frame(X_mat[idx, , drop = FALSE])
    
    # Compute PDP for each fold
    pdp_fold_list <- purrr::map(model_list, function(mod) {
      pdp::partial(
        object    = mod,
        pred.var  = c(grid_var1, grid_var2),
        train     = train_sub,
        pred.grid = custom_grid,
        type      = "regression",
        .f        = function(object, newdata) {
          predict(object, xgb.DMatrix(data.matrix(newdata)))
        }
      )
    })
    
    # Average PDP folds
    dplyr::bind_rows(pdp_fold_list) %>%
      dplyr::group_by(
        .data[[grid_var1]],
        .data[[grid_var2]]
      ) %>%
      dplyr::summarise(
        !!pred_name := mean(yhat, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(!!strata_var := strat)
  })
  
  dplyr::bind_rows(pdp_list_strata)
}


plot_stratified_surface_pdp <- function(
    df,
    x_var,
    y_var,
    z_var,
    group_col,
    x_title = NULL,
    y_title = NULL,
    z_title = NULL,
    strata_title = NULL, 
    strata_title_offset = 0.15,
    opacity = 0.8,
    label_size = 12,
    label_color = "black",
    x_offset = 0,
    y_offset = 0,
    z_offset = 0
) {
  # Default axis titles
  if (is.null(x_title)) x_title <- x_var
  if (is.null(y_title)) y_title <- y_var
  if (is.null(z_title)) z_title <- z_var
  
  p <- plotly::plot_ly()
  groups <- unique(df[[group_col]])
  
  max_z <- 0
  
  for (g in groups) {
    df_g <- df[df[[group_col]] == g, ]
    
    # sorted unique values
    x_vals <- sort(unique(df_g[[x_var]]))
    y_vals <- sort(unique(df_g[[y_var]]))
    
    # order by x then y
    df_g <- df_g[order(df_g[[x_var]], df_g[[y_var]]), ]
    
    # reshape z into matrix: rows = length(y_vals), byrow = FALSE
    z_mat <- matrix(df_g[[z_var]], nrow = length(y_vals), byrow = FALSE)
    
    # add surface
    p <- p %>%
      add_surface(
        x = x_vals,
        y = y_vals,
        z = z_mat,
        name = as.character(g),
        opacity = opacity,
        showscale = FALSE
      )
    
    # corner label (min x, max y)
    corner_x <- min(x_vals) + x_offset
    corner_y <- max(y_vals) + y_offset
    ix <- which.min(x_vals)
    iy <- which.max(y_vals)
    corner_z <- z_mat[iy, ix] - z_offset
    
    max_z <- max(max_z, corner_z)
    
    p <- p %>%
      add_text(
        x = corner_x,
        y = corner_y,
        z = corner_z,
        text = g,
        showlegend = FALSE,
        textposition = "top center",
        textfont = list(size = label_size, color = label_color)
      )
  }
  
  # axis labels
  p %>%
    layout(
      scene = list(
        xaxis = list(title = x_title),
        yaxis = list(title = y_title),
        zaxis = list(title = z_title)
      )
    )
  if(!is.null(strata_title)) {
    p <- p %>%
      add_text(
        x = corner_x,
        y = corner_y,
        z = max_z + strata_title_offset,
        text = strata_title, 
        showlegend = FALSE,
        textposition = "top center",
        textfont = list(size = 14, color = "black")
      )
  }
  return(p)
}


