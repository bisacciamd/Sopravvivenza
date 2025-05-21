# ===========================================================================
# R/visualization.R
# ===========================================================================

#' Create Kaplan-Meier Plot Comparing Two Groups
#'
#' This function creates and optionally saves a Kaplan-Meier plot comparing
#' survival curves for two or more groups.
#'
#' @param combined_data Combined dataset with both groups
#' @param group_var Name of the group variable in combined_data
#' @param time_var Name of the time variable in combined_data (default: "time")
#' @param status_var Name of the status variable in combined_data (default: "status")
#' @param title Plot title (default: "Kaplan-Meier Plot")
#' @param group_labels Vector of labels for the groups (default: NULL, will use group names)
#' @param colors Vector of colors for the groups (default: c("#2E9FDF", "#E7B800"))
#' @param filename Output file name without extension (default: "km_plot")
#' @param output_dir Output directory (default: NULL, no file saved)
#' @param file_format File format for saving: "pdf", "png", or "both" (default: "both")
#' @param show_risk_table Whether to show risk table (default: TRUE)
#' @param show_pvalue Whether to show p-value (default: TRUE)
#' @param show_conf_int Whether to show confidence intervals (default: TRUE)
#' @param xlab X-axis label (default: "Time")
#' @param ylab Y-axis label (default: "Survival Probability")
#' @param risk_table_height Height of risk table relative to plot (default: 0.25)
#'
#' @return A list with plot object and, if requested, paths to saved files
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' create_km_plot(
#'   combined_data,
#'   group_var = "group",
#'   title = "Treatment vs. Control"
#' )
#' 
#' # Customized plot
#' create_km_plot(
#'   combined_data,
#'   group_var = "group",
#'   group_labels = c("Control Group", "Treatment Group"),
#'   colors = c("blue", "red"),
#'   title = "Survival Comparison",
#'   output_dir = "figures",
#'   show_pvalue = FALSE
#' )
#' }
create_km_plot <- function(combined_data, group_var = "group", 
                          time_var = "time", status_var = "status",
                          title = "Kaplan-Meier Plot", 
                          group_labels = NULL, colors = c("#2E9FDF", "#E7B800"),
                          filename = "km_plot", output_dir = NULL,
                          file_format = "both", show_risk_table = TRUE,
                          show_pvalue = TRUE, show_conf_int = TRUE,
                          xlab = "Time", ylab = "Survival Probability",
                          risk_table_height = 0.25) {
  
  # Check if required packages are available
  if (!requireNamespace("survminer", quietly = TRUE)) {
    stop("Package 'survminer' is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  # Check if output directory exists, create if requested and not existing
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # Create Surv formula
  surv_formula <- stats::as.formula(paste("survival::Surv(", time_var, ",", status_var, ") ~", group_var))
  
  # Create the survival fit
  fit <- survival::survfit(surv_formula, data = combined_data)
  
  # Set default group labels if not provided
  if (is.null(group_labels)) {
    group_labels <- levels(combined_data[[group_var]])
  }
  
  # Create plot using survminer
  km_plot <- survminer::ggsurvplot(
    fit,
    data = combined_data,
    title = title,
    pval = show_pvalue,
    conf.int = show_conf_int,
    risk.table = show_risk_table,
    legend.labs = group_labels,
    palette = colors,
    risk.table.height = risk_table_height,
    ggtheme = ggplot2::theme_minimal(),
    xlab = xlab,
    ylab = ylab
  )
  
  # Save as requested format(s)
  output_files <- list()
  
  if (!is.null(output_dir)) {
    if (file_format %in% c("pdf", "both")) {
      pdf_path <- file.path(output_dir, paste0(filename, ".pdf"))
      grDevices::pdf(pdf_path, width = 8, height = 6)
      print(km_plot)
      grDevices::dev.off()
      output_files$pdf <- pdf_path
    }
    
    if (file_format %in% c("png", "both")) {
      png_path <- file.path(output_dir, paste0(filename, ".png"))
      grDevices::png(png_path, width = 800, height = 600, res = 100)
      print(km_plot)
      grDevices::dev.off()
      output_files$png <- png_path
    }
  }
  
  # Return results
  results <- list(
    plot = km_plot,
    fit = fit,
    output_files = output_files
  )
  
  return(results)
}

#' Create Forest Plot for Meta-Analysis
#'
#' This function creates and optionally saves a forest plot visualizing
#' meta-analysis results.
#'
#' @param meta_results Results from perform_fixed_effects_meta or perform_random_effects_meta
#' @param title Plot title (default: "Forest Plot: Meta-Analysis Results")
#' @param xlab X-axis label (default: "Hazard Ratio")
#' @param filename Output file name without extension (default: "forest_plot")
#' @param output_dir Output directory (default: NULL, no file saved)
#' @param file_format File format for saving: "pdf", "png", or "both" (default: "both")
#' @param log_scale Whether to use log scale for x-axis (default: TRUE)
#' @param show_weights Whether to show study weights (default: TRUE)
#' @param show_stats Whether to show heterogeneity statistics (default: TRUE)
#' @param null_line_value Value for reference line (default: 1 for hazard ratios)
#' @param text_size Base text size for the plot (default: 12)
#'
#' @return A list with plot object and, if requested, paths to saved files
#' @export
#'
#' @examples
#' \dontrun{
#' # Perform meta-analysis
#' meta_data <- data.frame(
#'   study = c("Study 1", "Study 2", "Study 3"),
#'   hr = c(1.5, 2.1, 1.8),
#'   lower = c(1.1, 1.6, 1.2),
#'   upper = c(2.0, 2.8, 2.7)
#' )
#' 
#' meta_results <- perform_fixed_effects_meta(meta_data)
#' 
#' # Create forest plot
#' create_forest_plot(
#'   meta_results,
#'   title = "Forest Plot: Hazard Ratios",
#'   output_dir = "figures"
#' )
#' }
create_forest_plot <- function(meta_results, 
                              title = "Forest Plot: Meta-Analysis Results",
                              xlab = "Hazard Ratio",
                              filename = "forest_plot", 
                              output_dir = NULL,
                              file_format = "both",
                              log_scale = TRUE,
                              show_weights = TRUE,
                              show_stats = TRUE,
                              null_line_value = 1,
                              text_size = 12) {
  
  # Extract meta data and results
  meta_data <- meta_results$meta_data
  
  # Check if weights are available
  has_weights <- "weight_percent" %in% names(meta_data)
  if (show_weights && !has_weights) {
    warning("Weight percentages not found in meta_data. Setting show_weights to FALSE.")
    show_weights <- FALSE
  }
  
  # Ensure meta_data has study column
  if (!"study" %in% names(meta_data)) {
    stop("meta_data must contain a 'study' column")
  }
  
  # Check if output directory exists, create if requested and not existing
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # Set margins
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  graphics::par(mar = c(5, 10, 4, 4))
  
  # Determine x-axis range based on HR values
  max_hr <- max(c(meta_data$upper, meta_results$ci_upper)) * 1.1
  min_hr <- min(c(meta_data$lower, meta_results$ci_lower)) * 0.9
  min_hr <- max(0.01, min_hr)  # Ensure not too close to 0 for log scale
  
  # Create forest plot function
  create_plot <- function() {
    # Create empty plot
    graphics::plot(NA, NA, xlim = c(min_hr, max_hr), 
                 ylim = c(0, nrow(meta_data) + 3), 
                 xlab = xlab, ylab = "", yaxt = "n", 
                 log = ifelse(log_scale, "x", ""),
                 cex.lab = text_size/12, cex.axis = text_size/12)
    
    # Add reference line
    graphics::abline(v = null_line_value, lty = 2)
    
    # Add study names and results
    for (i in 1:nrow(meta_data)) {
      y_pos <- nrow(meta_data) - i + 1
      
      # Study name
      graphics::text(min_hr * 0.9, y_pos, meta_data$study[i], pos = 4, cex = text_size/12)
      
      # Draw point for HR
      graphics::points(meta_data$hr[i], y_pos, pch = 15, cex = text_size/10)
      
      # Draw confidence interval
      graphics::segments(meta_data$lower[i], y_pos, meta_data$upper[i], y_pos, lwd = 2)
      
      # Add HR and CI as text
      hr_text <- sprintf("%.2f (%.2f-%.2f)", 
                        meta_data$hr[i], meta_data$lower[i], meta_data$upper[i])
      graphics::text(max_hr * 1.1, y_pos, hr_text, pos = 2, cex = text_size/12)
      
      # Add weight if available
      if (show_weights && has_weights) {
        weight_text <- sprintf("%.1f%%", meta_data$weight_percent[i])
        graphics::text(max_hr * 1.3, y_pos, weight_text, pos = 2, cex = text_size/12)
      }
    }
    
    # Add pooled result
    rect_y_bottom <- 0.2
    rect_y_top <- 0.8
    graphics::rect(meta_results$ci_lower, rect_y_bottom, 
                 meta_results$ci_upper, rect_y_top, 
                 col = "gray80", border = NA)
    graphics::points(meta_results$pooled_hr, (rect_y_bottom + rect_y_top)/2, pch = 18, cex = text_size/8)
    graphics::text(min_hr * 0.9, (rect_y_bottom + rect_y_top)/2, "Pooled", font = 2, pos = 4, cex = text_size/12)
    
    pooled_text <- sprintf("%.2f (%.2f-%.2f)", 
                         meta_results$pooled_hr, meta_results$ci_lower, meta_results$ci_upper)
    graphics::text(max_hr * 1.1, (rect_y_bottom + rect_y_top)/2, 
                 pooled_text, pos = 2, cex = text_size/12, font = 2)
    
    if (show_weights) {
      p_text <- sprintf("p = %.3f", meta_results$p_value)
      graphics::text(max_hr * 1.3, (rect_y_bottom + rect_y_top)/2, 
                   p_text, pos = 2, cex = text_size/12, font = 2)
    }
    
    # Add column headers
    graphics::text(min_hr * 0.9, nrow(meta_data) + 1, "Study", 
                 pos = 4, font = 2, cex = text_size/12)
    graphics::text(max_hr * 1.1, nrow(meta_data) + 1, "HR (95% CI)", 
                 pos = 2, font = 2, cex = text_size/12)
    
    if (show_weights) {
      graphics::text(max_hr * 1.3, nrow(meta_data) + 1, "Weight", 
                   pos = 2, font = 2, cex = text_size/12)
    }
    
    # Add titles
    graphics::title(main = title, cex.main = text_size/10)
    
    # Add heterogeneity statistics if requested
    if (show_stats && !is.null(meta_results$heterogeneity)) {
      het <- meta_results$heterogeneity
      het_text <- sprintf("Heterogeneity: Q = %.2f, df = %d, p = %.3f, I² = %.1f%%",
                        het$q, het$df, het$p, het$i_squared)
      graphics::mtext(het_text, side = 3, line = -1.5, cex = 0.8 * text_size/12)
    }
    
    # Add interpretation
    favors_left <- "Favors Reference"
    favors_right <- "Favors Comparison"
    graphics::mtext(paste(favors_left, "    ", favors_right), 
                  side = 3, line = 0, cex = 0.9 * text_size/12)
  }
  
  # Save plots if requested
  output_files <- list()
  
  if (!is.null(output_dir)) {
    if (file_format %in% c("pdf", "both")) {
      pdf_path <- file.path(output_dir, paste0(filename, ".pdf"))
      grDevices::pdf(pdf_path, width = 10, height = 8)
      create_plot()
      grDevices::dev.off()
      output_files$pdf <- pdf_path
    }
    
    if (file_format %in% c("png", "both")) {
      png_path <- file.path(output_dir, paste0(filename, ".png"))
      grDevices::png(png_path, width = 1000, height = 800, res = 100)
      create_plot()
      grDevices::dev.off()
      output_files$png <- png_path
    }
  }
  
  # Create plot in current device
  create_plot()
  
  # Return results
  results <- list(
    meta_results = meta_results,
    output_files = output_files
  )
  
  return(results)
}

#' Create Dose-Response Plot for Meta-Regression
#'
#' This function creates and optionally saves a dose-response plot visualizing
#' the relationship between a moderator variable and effect size.
#'
#' @param spline_results Results from perform_spline_meta_regression
#' @param title Plot title (default: NULL, auto-generated)
#' @param xlab X-axis label (default: NULL, uses moderator name)
#' @param ylab Y-axis label (default: "Hazard Ratio")
#' @param filename Output file name without extension (default: "dose_response_plot")
#' @param output_dir Output directory (default: NULL, no file saved)
#' @param file_format File format for saving: "pdf", "png", or "both" (default: "both")
#' @param show_points Whether to show original study points (default: TRUE)
#' @param show_error_bars Whether to show error bars for study points (default: TRUE)
#' @param show_labels Whether to show study labels (default: FALSE)
#' @param null_line_value Value for reference line (default: 1 for hazard ratios)
#' @param ribbon_alpha Transparency of confidence interval ribbon (default: 0.2)
#' @param point_size Size of study points (default: 3)
#' @param line_size Size of spline curve line (default: 1)
#'
#' @return A list with plot object and, if requested, paths to saved files
#' @export
#'
#' @examples
#' \dontrun{
#' # Perform spline meta-regression
#' meta_data <- data.frame(
#'   study = c("Study 1", "Study 2", "Study 3"),
#'   hr = c(1.5, 2.1, 1.8),
#'   lower = c(1.1, 1.6, 1.2),
#'   upper = c(2.0, 2.8, 2.7),
#'   gls = c(-16, -14, -12)  # Moderator variable
#' )
#' 
#' spline_results <- perform_spline_meta_regression(
#'   meta_data, 
#'   moderator_name = "gls"
#' )
#' 
#' # Create dose-response plot
#' create_dose_response_plot(
#'   spline_results,
#'   output_dir = "figures"
#' )
#' }
create_dose_response_plot <- function(spline_results, 
                                     title = NULL,
                                     xlab = NULL,
                                     ylab = "Hazard Ratio",
                                     filename = "dose_response_plot", 
                                     output_dir = NULL,
                                     file_format = "both",
                                     show_points = TRUE,
                                     show_error_bars = TRUE,
                                     show_labels = FALSE,
                                     null_line_value = 1,
                                     ribbon_alpha = 0.2,
                                     point_size = 3,
                                     line_size = 1) {
  
  # Check if required packages are available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  # Extract data from spline_results
  meta_data <- spline_results$meta_data
  pred_grid <- spline_results$predictions
  moderator_name <- spline_results$moderator_name
  
  # Set default axis labels if not provided
  if (is.null(xlab)) {
    xlab <- moderator_name
  }
  
  # Set default title if not provided
  if (is.null(title)) {
    title <- paste("Dose-Response Relationship:", moderator_name, "vs. Hazard Ratio")
  }
  
  # Check if output directory exists, create if requested and not existing
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # Create ggplot
  p <- ggplot2::ggplot()
  
  # Add confidence interval ribbon
  p <- p + ggplot2::geom_ribbon(
    data = pred_grid,
    ggplot2::aes(
      x = .data[[moderator_name]],
      ymin = hr_lower,
      ymax = hr_upper
    ),
    alpha = ribbon_alpha,
    fill = "grey70"
  )
  
  # Add spline curve
  p <- p + ggplot2::geom_line(
    data = pred_grid,
    ggplot2::aes(
      x = .data[[moderator_name]],
      y = hr
    ),
    size = line_size
  )
  
  # Add study points if requested
  if (show_points) {
    p <- p + ggplot2::geom_point(
      data = meta_data,
      ggplot2::aes(
        x = .data[[moderator_name]],
        y = hr
      ),
      size = point_size
    )
  }
  
  # Add error bars if requested
  if (show_error_bars && show_points) {
    p <- p + ggplot2::geom_errorbar(
      data = meta_data,
      ggplot2::aes(
        x = .data[[moderator_name]],
        ymin = lower,
        ymax = upper
      ),
      width = 0.2
    )
  }
  
  # Add study labels if requested
  if (show_labels && "study" %in% names(meta_data)) {
    p <- p + ggplot2::geom_text(
      data = meta_data,
      ggplot2::aes(
        x = .data[[moderator_name]],
        y = hr,
        label = study
      ),
      nudge_y = 0.1,
      size = 3
    )
  }
  
  # Add reference line
  p <- p + ggplot2::geom_hline(
    yintercept = null_line_value,
    linetype = "dashed"
  )
  
  # Add titles and labels
  p <- p + ggplot2::labs(
    title = title,
    x = xlab,
    y = ylab
  )
  
  # Use minimal theme
  p <- p + ggplot2::theme_minimal() + 
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "none"
    )
  
  # Save plots if requested
  output_files <- list()
  
  if (!is.null(output_dir)) {
    if (file_format %in% c("pdf", "both")) {
      pdf_path <- file.path(output_dir, paste0(filename, ".pdf"))
      ggplot2::ggsave(pdf_path, p, width = 10, height = 8)
      output_files$pdf <- pdf_path
    }
    
    if (file_format %in% c("png", "both")) {
      png_path <- file.path(output_dir, paste0(filename, ".png"))
      ggplot2::ggsave(png_path, p, width = 10, height = 8, dpi = 100)
      output_files$png <- png_path
    }
  }
  
  # Return results
  results <- list(
    plot = p,
    spline_results = spline_results,
    output_files = output_files
  )
  
  return(results)
}

#' Create Funnel Plot for Publication Bias Assessment
#'
#' This function creates and optionally saves a funnel plot to assess
#' publication bias in meta-analysis.
#'
#' @param meta_results Results from perform_fixed_effects_meta or perform_random_effects_meta
#' @param precision_measure Measure for y-axis: "se" (standard error, default), "precision" (1/se), or "sample" (sample size)
#' @param title Plot title (default: "Funnel Plot")
#' @param xlab X-axis label (default: "Effect Size")
#' @param ylab Y-axis label (default: NULL, auto-generated)
#' @param filename Output file name without extension (default: "funnel_plot")
#' @param output_dir Output directory (default: NULL, no file saved)
#' @param file_format File format for saving: "pdf", "png", or "both" (default: "both")
#' @param show_contours Whether to show contour lines for significance (default: TRUE)
#' @param show_se_lines Whether to show standard error guidance lines (default: TRUE)
#' @param inverted Whether to invert the y-axis (default: TRUE for "se", FALSE otherwise)
#' @param eggers_line Whether to show Egger's regression line (default: TRUE)
#'
#' @return A list with plot object and, if requested, paths to saved files
#' @export
#'
#' @examples
#' \dontrun{
#' # Perform meta-analysis
#' meta_data <- data.frame(
#'   study = c("Study 1", "Study 2", "Study 3"),
#'   hr = c(1.5, 2.1, 1.8),
#'   lower = c(1.1, 1.6, 1.2),
#'   upper = c(2.0, 2.8, 2.7)
#' )
#' 
#' meta_results <- perform_fixed_effects_meta(meta_data)
#' 
#' # Create funnel plot
#' create_funnel_plot(
#'   meta_results,
#'   output_dir = "figures"
#' )
#' }
create_funnel_plot <- function(meta_results, 
                              precision_measure = "se",
                              title = "Funnel Plot",
                              xlab = "Effect Size",
                              ylab = NULL,
                              filename = "funnel_plot", 
                              output_dir = NULL,
                              file_format = "both",
                              show_contours = TRUE,
                              show_se_lines = TRUE,
                              inverted = ifelse(precision_measure == "se", TRUE, FALSE),
                              eggers_line = TRUE) {
  
  # Extract meta data
  meta_data <- meta_results$meta_data
  
  # Check if required columns exist
  if (!all(c("loghr", "se") %in% names(meta_data))) {
    stop("meta_data must contain 'loghr' and 'se' columns")
  }
  
  # Calculate precision measure for y-axis
  if (precision_measure == "se") {
    y_values <- meta_data$se
    if (is.null(ylab)) ylab <- "Standard Error"
  } else if (precision_measure == "precision") {
    y_values <- 1 / meta_data$se
    if (is.null(ylab)) ylab <- "Precision (1/SE)"
  } else if (precision_measure == "sample") {
    if (!"sample_size" %in% names(meta_data)) {
      stop("meta_data must contain 'sample_size' column when precision_measure = 'sample'")
    }
    y_values <- sqrt(meta_data$sample_size)
    if (is.null(ylab)) ylab <- "Study Size (sqrt)"
  } else {
    stop("precision_measure must be 'se', 'precision', or 'sample'")
  }
  
  # Add precision measure to meta_data
  meta_data$precision_measure <- y_values
  
  # Perform Egger's test if requested
  if (eggers_line) {
    if (precision_measure == "se") {
      # For standard error, regress effect against se
      eggers_model <- stats::lm(meta_data$loghr ~ meta_data$se)
    } else {
      # For other measures, use the inverse relationship
      eggers_model <- stats::lm(meta_data$loghr ~ I(1/meta_data$precision_measure))
    }
    eggers_summary <- summary(eggers_model)
    eggers_p <- eggers_summary$coefficients[2, "Pr(>|t|)"]
  }
  
  # Check if output directory exists, create if requested and not existing
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # Set up plot parameters
  pooled_effect <- meta_results$loghr  # Log scale for symmetry
  if (is.null(pooled_effect) && !is.null(meta_results$pooled_hr)) {
    pooled_effect <- log(meta_results$pooled_hr)
  }
  
  # Define x-axis range
  x_range <- range(c(meta_data$loghr, pooled_effect)) * 1.2
  
  # Define y-axis range and limits
  if (precision_measure == "se") {
    y_max <- max(meta_data$se) * 1.2
    y_range <- c(0, y_max)
  } else {
    y_min <- min(meta_data$precision_measure) * 0.8
    y_max <- max(meta_data$precision_measure) * 1.2
    y_range <- c(y_min, y_max)
  }
  
  # Create plot function
  create_plot <- function() {
    # Create empty plot
    if (inverted) {
      graphics::plot(meta_data$loghr, meta_data$precision_measure, 
                   xlim = x_range, ylim = rev(y_range),
                   xlab = xlab, ylab = ylab, type = "n",
                   main = title)
    } else {
      graphics::plot(meta_data$loghr, meta_data$precision_measure, 
                   xlim = x_range, ylim = y_range,
                   xlab = xlab, ylab = ylab, type = "n",
                   main = title)
    }
    
    # Add vertical line at pooled effect
    graphics::abline(v = pooled_effect, lty = 2)
    
    # Add contour lines if requested
    if (show_contours && precision_measure == "se") {
      # Draw contours for p = 0.01, 0.05, 0.1
      z_values <- c(2.58, 1.96, 1.64)  # Z-scores for p = 0.01, 0.05, 0.1
      max_se <- max(meta_data$se)
      se_range <- seq(0, max_se, length.out = 100)
      
      for (z in z_values) {
        # Left side
        graphics::lines(pooled_effect - z * se_range, se_range, lty = 3, col = "gray60")
        # Right side
        graphics::lines(pooled_effect + z * se_range, se_range, lty = 3, col = "gray60")
      }
      
      # Add legend for contour lines
      graphics::legend("bottomright", 
                     legend = c("p = 0.01", "p = 0.05", "p = 0.10"),
                     lty = 3, col = "gray60", bty = "n")
    }
    
    # Add Egger's regression line if requested
    if (eggers_line) {
      if (precision_measure == "se") {
        # Draw regression line for effect ~ se
        se_range <- seq(0, max(meta_data$se), length.out = 100)
        pred_effects <- stats::predict(eggers_model, 
                                     newdata = data.frame(se = se_range))
        graphics::lines(pred_effects, se_range, col = "red", lty = 2)
        
        # Add Egger's test p-value
        p_text <- sprintf("Egger's test: p = %.3f", eggers_p)
        graphics::mtext(p_text, side = 3, line = 0.5, cex = 0.8)
      } else {
        # Draw regression line for other precision measures
        prec_range <- seq(min(meta_data$precision_measure), 
                         max(meta_data$precision_measure), 
                         length.out = 100)
        pred_effects <- stats::predict(eggers_model, 
                                     newdata = data.frame(
                                       precision_measure = prec_range))
        graphics::lines(pred_effects, prec_range, col = "red", lty = 2)
        
        # Add Egger's test p-value
        p_text <- sprintf("Egger's test: p = %.3f", eggers_p)
        graphics::mtext(p_text, side = 3, line = 0.5, cex = 0.8)
      }
    }
    
    # Add study points
    graphics::points(meta_data$loghr, meta_data$precision_measure, pch = 19)
    
    # Add grid for readability
    graphics::grid()
  }
  
  # Save plots if requested
  output_files <- list()
  
  if (!is.null(output_dir)) {
    if (file_format %in% c("pdf", "both")) {
      pdf_path <- file.path(output_dir, paste0(filename, ".pdf"))
      grDevices::pdf(pdf_path, width = 8, height = 8)
      create_plot()
      grDevices::dev.off()
      output_files$pdf <- pdf_path
    }
    
    if (file_format %in% c("png", "both")) {
      png_path <- file.path(output_dir, paste0(filename, ".png"))
      grDevices::png(png_path, width = 800, height = 800, res = 100)
      create_plot()
      grDevices::dev.off()
      output_files$png <- png_path
    }
  }
  
  # Create plot in current device
  create_plot()
  
  # Return results
  results <- list(
    meta_results = meta_results,
    output_files = output_files,
    eggers_p = if (eggers_line) eggers_p else NULL
  )
  
  return(results)
}