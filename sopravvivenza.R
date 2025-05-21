# KM-Digitizer: A Template for Reconstructing IPD from Digitized Survival Curves
# Author: Giandomenico Bisaccia
# License: MIT
# GitHub: https://github.com/bisacciamd/km-digitizer

#=====================================================================
# DOCUMENTATION
#=====================================================================
# 
# This script provides a comprehensive template for reconstructing 
# Individual Patient Data (IPD) from digitized Kaplan-Meier survival curves
# using the method described by Guyot et al.
#
# Workflow:
# 1. Digitize curves using WebPlotDigitizer or similar software
# 2. Save digitized points as CSV files (time, survival)
# 3. Load these files into R using this script
# 4. Reconstruct IPD and validate the reconstruction
# 5. Calculate hazard ratios and perform meta-analysis
#
# Usage:
# - Modify the 'study_files' list to point to your digitized curve files
# - Adjust parameters like initial number at risk and total events
# - Run the script to generate reconstructed IPD, validation plots, and analyses
#

#=====================================================================
# DEPENDENCIES
#=====================================================================

# Basic packages
if (!require("survival")) install.packages("survival")
if (!require("survminer")) install.packages("survminer")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("readr")) install.packages("readr")
if (!require("writexl")) install.packages("writexl")
if (!require("splines")) install.packages("splines")

# Load packages
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(readr)
library(writexl)
library(splines)

# Optional packages for meta-analysis
# Uncomment if needed:
# if (!require("meta")) install.packages("meta")
# if (!require("metafor")) install.packages("metafor")
# library(meta)
# library(metafor)

#=====================================================================
# DIRECTORY SETUP
#=====================================================================

# Create output directories
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

#=====================================================================
# CORE FUNCTIONS
#=====================================================================

#' Read and Preprocess Digitized Curve Data
#'
#' @param file_path Path to CSV file with digitized curve data
#' @return Preprocessed dataframe with time and survival probabilities
#'
read_curve_data <- function(file_path) {
  # Read the CSV file
  data <- tryCatch({
    read_csv(file_path, col_names = FALSE, show_col_types = FALSE)
  }, error = function(e) {
    cat("Error reading file:", file_path, "\n")
    cat("Error message:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(data)) return(NULL)
  
  # Check if data has at least 2 columns
  if (ncol(data) < 2) {
    cat("Warning: File", file_path, "does not have at least 2 columns\n")
    return(NULL)
  }
  
  # Rename columns for clarity
  names(data)[1:2] <- c("time", "surv")
  
  # Keep only the first two columns
  data <- data[, c("time", "surv")]
  
  # Ensure survival probability is between 0 and 1
  # For some studies values might be in percentages
  if (max(data$surv, na.rm = TRUE) > 1.1) {
    data$surv <- data$surv / 100
  }
  
  # Ensure time is in months (assume it is if max time is > 100)
  # Otherwise, convert from years to months if likely in years
  if (max(data$time, na.rm = TRUE) < 100) {
    # Check if likely in years (max < 10)
    if (max(data$time, na.rm = TRUE) < 10) {
      data$time <- data$time * 12
    }
  }
  
  # Sort by time
  data <- data[order(data$time), ]
  
  # Ensure the curve starts at time 0 and surv = 1
  if (min(data$time, na.rm = TRUE) > 0) {
    data <- rbind(data.frame(time = 0, surv = 1), data)
  } else {
    # If there's already a time 0, ensure surv = 1
    data$surv[data$time == 0] <- 1
  }
  
  return(data)
}

#' Reconstruct Individual Patient Data from Digitized Survival Curve
#'
#' @param time_points Vector of time points from digitized curve
#' @param surv_probs Vector of survival probabilities from digitized curve
#' @param n_risk_initial Initial number at risk (if known)
#' @param total_events Total number of events (if known)
#' @param study_name Name of the study (for output files)
#' @param arm_name Name of the study arm/group (for output files)
#' @param risk_times Optional vector of time points for number at risk
#' @param risk_numbers Optional vector of numbers at risk at risk_times
#' @param output_dir Directory for output files
#' @param make_plot Whether to create validation plot
#' @return List containing reconstructed IPD and fit object
#'
reconstruct_ipd <- function(time_points, surv_probs, n_risk_initial = 100, 
                           total_events = NULL, study_name = "study", 
                           arm_name = "arm", risk_times = NULL, 
                           risk_numbers = NULL, output_dir = "results",
                           make_plot = TRUE) {
  
  # Create directory for study outputs
  study_dir <- file.path(output_dir, study_name)
  dir.create(study_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Ensure time points and survival probabilities are sorted
  idx <- order(time_points)
  time_points <- time_points[idx]
  surv_probs <- surv_probs[idx]
  
  # Ensure start at time 0 with survival 1
  if (time_points[1] > 0) {
    time_points <- c(0, time_points)
    surv_probs <- c(1, surv_probs)
  } else if (surv_probs[1] != 1) {
    surv_probs[1] <- 1
  }
  
  # Initialize variables
  n <- n_risk_initial
  times <- numeric(0)
  events <- numeric(0)
  
  # Process risk data if available
  has_risk_data <- !is.null(risk_times) && !is.null(risk_numbers) && 
                  length(risk_times) == length(risk_numbers)
  
  # Generate IPD for each interval
  for (i in 1:(length(time_points) - 1)) {
    t1 <- time_points[i]
    t2 <- time_points[i + 1]
    s1 <- surv_probs[i]
    s2 <- surv_probs[i + 1]
    
    # If risk data is available, update n
    if (has_risk_data) {
      # Find closest risk time point
      if (i > 1) {  # Skip first interval
        closest_idx <- which.min(abs(risk_times - t1))
        if (length(closest_idx) > 0) {
          n <- risk_numbers[closest_idx]
        }
      }
    }
    
    # Calculate hazard for this interval
    if (s1 > 0) {
      h <- (s1 - s2) / s1
    } else {
      h <- 0
    }
    
    # Calculate expected number of events
    d <- round(n * h)
    
    # Generate event times (spread throughout interval)
    if (d > 0) {
      event_time_points <- seq(t1, t2, length.out = d + 2)
      event_time_points <- event_time_points[2:(length(event_time_points) - 1)]
      
      times <- c(times, event_time_points)
      events <- c(events, rep(1, length(event_time_points)))
    }
    
    # Update number at risk
    n <- n - d
  }
  
  # Create IPD data frame
  ipd <- data.frame(
    time = times,
    status = events
  )
  
  # Adjust for total events if specified
  if (!is.null(total_events) && length(events) > 0 && sum(events) != total_events) {
    current_events <- sum(events)
    
    if (current_events < total_events) {
      # Need more events - add random events
      more_needed <- total_events - current_events
      # Generate random times within range
      random_times <- runif(more_needed, min(time_points), max(time_points))
      new_events <- data.frame(
        time = random_times,
        status = rep(1, more_needed)
      )
      ipd <- rbind(ipd, new_events)
    } else if (current_events > total_events) {
      # Too many events - randomly remove some
      to_remove <- current_events - total_events
      event_indices <- which(ipd$status == 1)
      if (length(event_indices) >= to_remove) {
        remove_indices <- sample(event_indices, to_remove)
        ipd <- ipd[-remove_indices, ]
      }
    }
  }
  
  # Sort by time
  ipd <- ipd[order(ipd$time), ]
  
  # Create a validation plot
  if (make_plot) {
    fit <- survfit(Surv(time, status) ~ 1, data = ipd)
    
    # Plot to PDF file
    pdf(file.path(study_dir, paste0(study_name, "_", arm_name, "_validation.pdf")), 
        width = 8, height = 6)
    
    plot(fit, main = paste("Validation -", study_name, "-", arm_name),
         xlab = "Time", ylab = "Survival Probability", 
         conf.int = FALSE, mark.time = FALSE)
    points(time_points, surv_probs, col = "red", pch = 19)
    legend("topright", c("Reconstructed", "Original Digitized"), 
           col = c("black", "red"), lty = c(1, NA), pch = c(NA, 19))
    
    dev.off()
    
    # Also save a PNG version
    png(file.path(study_dir, paste0(study_name, "_", arm_name, "_validation.png")), 
        width = 800, height = 600, res = 100)
    
    plot(fit, main = paste("Validation -", study_name, "-", arm_name),
         xlab = "Time", ylab = "Survival Probability", 
         conf.int = FALSE, mark.time = FALSE)
    points(time_points, surv_probs, col = "red", pch = 19)
    legend("topright", c("Reconstructed", "Original Digitized"), 
           col = c("black", "red"), lty = c(1, NA), pch = c(NA, 19))
    
    dev.off()
  } else {
    fit <- survfit(Surv(time, status) ~ 1, data = ipd)
  }
  
  # Save IPD to CSV file
  write_csv(ipd, file.path(study_dir, paste0(study_name, "_", arm_name, "_ipd.csv")))
  
  return(list(
    ipd = ipd,
    fit = fit,
    original_data = data.frame(time = time_points, surv = surv_probs)
  ))
}

#' Calculate Hazard Ratio Between Two Groups
#'
#' @param data1 IPD dataframe for group 1
#' @param data2 IPD dataframe for group 2
#' @param group1_name Name for group 1
#' @param group2_name Name for group 2
#' @param covariates Optional vector of covariate names to adjust for
#' @return List with hazard ratio and related statistics
#'
calculate_hazard_ratio <- function(data1, data2, group1_name = "group1", 
                                  group2_name = "group2", covariates = NULL) {
  
  # Add group indicators
  data1$group <- group1_name
  data2$group <- group2_name
  
  # Combine datasets
  combined_data <- rbind(data1, data2)
  
  # Set group as a factor with group1 as reference
  combined_data$group <- factor(combined_data$group, levels = c(group1_name, group2_name))
  
  # Create formula for Cox model
  if (is.null(covariates)) {
    formula <- as.formula("Surv(time, status) ~ group")
  } else {
    formula <- as.formula(paste("Surv(time, status) ~ group +", 
                               paste(covariates, collapse = " + ")))
  }
  
  # Fit Cox model
  cox_model <- coxph(formula, data = combined_data)
  summary_cox <- summary(cox_model)
  
  # Extract hazard ratio and CI
  hr <- exp(summary_cox$coefficients[1, "coef"])
  hr_lower <- exp(summary_cox$coefficients[1, "coef"] - 1.96 * summary_cox$coefficients[1, "se(coef)"])
  hr_upper <- exp(summary_cox$coefficients[1, "coef"] + 1.96 * summary_cox$coefficients[1, "se(coef)"])
  p_value <- summary_cox$coefficients[1, "Pr(>|z|)"]
  
  return(list(
    hr = hr,
    lower = hr_lower,
    upper = hr_upper,
    p_value = p_value,
    cox_model = cox_model,
    combined_data = combined_data
  ))
}

#' Create Kaplan-Meier Plot Comparing Two Groups
#'
#' @param combined_data Combined dataset with both groups
#' @param group_var Name of the group variable
#' @param title Plot title
#' @param group_labels Vector of labels for the groups
#' @param colors Vector of colors for the groups
#' @param filename Output file name
#' @param output_dir Output directory
#' @return NULL (creates plot files)
#'
create_km_plot <- function(combined_data, group_var = "group", 
                          title = "Kaplan-Meier Plot", 
                          group_labels = NULL, colors = c("#2E9FDF", "#E7B800"),
                          filename = "km_plot", output_dir = "figures") {
  
  # Check if output directory exists, create if not
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Create the survival fit
  fit <- survfit(as.formula(paste("Surv(time, status) ~", group_var)), 
                data = combined_data)
  
  # Set default group labels if not provided
  if (is.null(group_labels)) {
    group_labels <- levels(combined_data[[group_var]])
  }
  
  # Create plot using survminer
  km_plot <- ggsurvplot(
    fit,
    data = combined_data,
    title = title,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    legend.labs = group_labels,
    palette = colors,
    risk.table.height = 0.25,
    ggtheme = theme_minimal(),
    xlab = "Time",
    ylab = "Survival Probability"
  )
  
  # Save as PDF
  pdf(file.path(output_dir, paste0(filename, ".pdf")), width = 8, height = 6)
  print(km_plot)
  dev.off()
  
  # Save as PNG
  png(file.path(output_dir, paste0(filename, ".png")), 
      width = 800, height = 600, res = 100)
  print(km_plot)
  dev.off()
}

#=====================================================================
# META-ANALYSIS FUNCTIONS
#=====================================================================

#' Perform Fixed-Effects Meta-Analysis
#'
#' @param meta_data Data frame with columns: study, hr, lower, upper
#' @param output_dir Directory for output files
#' @return List with meta-analysis results
#'
perform_fixed_effects_meta <- function(meta_data, output_dir = "results") {
  
  # Calculate log hazard ratios and standard errors
  meta_data$loghr <- log(meta_data$hr)
  meta_data$se <- (log(meta_data$upper) - log(meta_data$lower)) / (2 * 1.96)
  
  # Calculate weights
  meta_data$weight <- 1 / meta_data$se^2
  meta_data$weight_percent <- 100 * meta_data$weight / sum(meta_data$weight)
  
  # Calculate pooled estimate
  weighted_mean <- sum(meta_data$loghr * meta_data$weight) / sum(meta_data$weight)
  weighted_se <- sqrt(1 / sum(meta_data$weight))
  
  # Calculate 95% CI
  ci_lower <- exp(weighted_mean - 1.96 * weighted_se)
  ci_upper <- exp(weighted_mean + 1.96 * weighted_se)
  pooled_hr <- exp(weighted_mean)
  
  # Calculate z-score and p-value
  z_score <- weighted_mean / weighted_se
  p_value <- 2 * (1 - pnorm(abs(z_score)))
  
  # Create forest plot
  create_forest_plot(meta_data, pooled_hr, ci_lower, ci_upper, 
                    p_value, output_dir)
  
  # Save results
  results <- list(
    meta_data = meta_data,
    pooled_hr = pooled_hr,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    p_value = p_value,
    z_score = z_score
  )
  
  # Save summary to text file
  sink(file.path(output_dir, "meta_analysis_summary.txt"))
  cat("Fixed-Effects Meta-Analysis Results\n")
  cat("==================================\n\n")
  cat("Pooled Hazard Ratio: ", round(pooled_hr, 3), 
      " (95% CI: ", round(ci_lower, 3), " - ", round(ci_upper, 3), ")\n")
  cat("Z-score: ", round(z_score, 3), "\n")
  cat("P-value: ", format.pval(p_value, digits = 3), "\n\n")
  
  cat("Individual Study Results:\n")
  cat("------------------------\n\n")
  for (i in 1:nrow(meta_data)) {
    cat(meta_data$study[i], ": HR = ", round(meta_data$hr[i], 2), 
        " (95% CI: ", round(meta_data$lower[i], 2), " - ", 
        round(meta_data$upper[i], 2), "), Weight = ", 
        round(meta_data$weight_percent[i], 1), "%\n")
  }
  sink()
  
  return(results)
}

#' Create Forest Plot for Meta-Analysis
#'
#' @param meta_data Data frame with study results
#' @param pooled_hr Pooled hazard ratio
#' @param ci_lower Lower confidence interval
#' @param ci_upper Upper confidence interval
#' @param p_value P-value for pooled result
#' @param output_dir Output directory
#' @return NULL (creates plot files)
#'
create_forest_plot <- function(meta_data, pooled_hr, ci_lower, ci_upper, 
                              p_value, output_dir = "results") {
  
  # Set up the PDF file
  pdf(file.path(output_dir, "forest_plot.pdf"), width = 10, height = 8)
  
  # Set margins
  par(mar = c(5, 10, 4, 4))
  
  # Determine x-axis range based on HR values
  max_hr <- max(c(meta_data$upper, ci_upper)) * 1.1
  min_hr <- min(c(meta_data$lower, ci_lower)) * 0.9
  min_hr <- max(0, min_hr)  # Ensure not negative
  
  # Create empty plot
  plot(NA, NA, xlim = c(min_hr, max_hr), 
       ylim = c(0, nrow(meta_data) + 3), 
       xlab = "Hazard Ratio", ylab = "", yaxt = "n", 
       log = "x")  # Log scale for HR
  
  # Add reference line at HR = 1
  abline(v = 1, lty = 2)
  
  # Add study names and results
  for (i in 1:nrow(meta_data)) {
    y_pos <- nrow(meta_data) - i + 1
    
    # Study name
    text(min_hr * 0.9, y_pos, meta_data$study[i], pos = 4, cex = 0.9)
    
    # Draw point for HR
    points(meta_data$hr[i], y_pos, pch = 15, cex = 1.2)
    
    # Draw confidence interval
    segments(meta_data$lower[i], y_pos, meta_data$upper[i], y_pos, lwd = 2)
    
    # Add HR and CI as text
    text(max_hr * 1.1, y_pos, 
         sprintf("%.2f (%.2f-%.2f)", meta_data$hr[i], meta_data$lower[i], meta_data$upper[i]), 
         pos = 2, cex = 0.8)
    
    # Add weight
    text(max_hr * 1.3, y_pos, 
         sprintf("%.1f%%", meta_data$weight_percent[i]), 
         pos = 2, cex = 0.8)
  }
  
  # Add pooled result
  rect(ci_lower, 0.2, ci_upper, 0.8, col = "gray80", border = NA)
  points(pooled_hr, 0.5, pch = 18, cex = 1.5)
  text(min_hr * 0.9, 0.5, "Pooled", font = 2, pos = 4)
  text(max_hr * 1.1, 0.5, 
       sprintf("%.2f (%.2f-%.2f)", pooled_hr, ci_lower, ci_upper), 
       pos = 2, cex = 0.8, font = 2)
  text(max_hr * 1.3, 0.5, 
       sprintf("p = %.3f", p_value), 
       pos = 2, cex = 0.8, font = 2)
  
  # Add column headers
  text(min_hr * 0.9, nrow(meta_data) + 1, "Study", pos = 4, font = 2)
  text(max_hr * 1.1, nrow(meta_data) + 1, "HR (95% CI)", pos = 2, font = 2)
  text(max_hr * 1.3, nrow(meta_data) + 1, "Weight", pos = 2, font = 2)
  
  # Add titles
  title(main = "Forest Plot: Meta-Analysis Results", 
        line = 1)
  
  # Add interpretation
  mtext("Favors Group 1     Favors Group 2", side = 3, line = 0)
  
  dev.off()
  
  # Also create a PNG version
  png(file.path(output_dir, "forest_plot.png"), 
      width = 1000, height = 800, res = 100)
  
  # Repeat the plotting code for PNG
  par(mar = c(5, 10, 4, 4))
  plot(NA, NA, xlim = c(min_hr, max_hr), 
       ylim = c(0, nrow(meta_data) + 3), 
       xlab = "Hazard Ratio", ylab = "", yaxt = "n", 
       log = "x")
  abline(v = 1, lty = 2)
  for (i in 1:nrow(meta_data)) {
    y_pos <- nrow(meta_data) - i + 1
    text(min_hr * 0.9, y_pos, meta_data$study[i], pos = 4, cex = 0.9)
    points(meta_data$hr[i], y_pos, pch = 15, cex = 1.2)
    segments(meta_data$lower[i], y_pos, meta_data$upper[i], y_pos, lwd = 2)
    text(max_hr * 1.1, y_pos, 
         sprintf("%.2f (%.2f-%.2f)", meta_data$hr[i], meta_data$lower[i], meta_data$upper[i]), 
         pos = 2, cex = 0.8)
    text(max_hr * 1.3, y_pos, 
         sprintf("%.1f%%", meta_data$weight_percent[i]), 
         pos = 2, cex = 0.8)
  }
  rect(ci_lower, 0.2, ci_upper, 0.8, col = "gray80", border = NA)
  points(pooled_hr, 0.5, pch = 18, cex = 1.5)
  text(min_hr * 0.9, 0.5, "Pooled", font = 2, pos = 4)
  text(max_hr * 1.1, 0.5, 
       sprintf("%.2f (%.2f-%.2f)", pooled_hr, ci_lower, ci_upper), 
       pos = 2, cex = 0.8, font = 2)
  text(max_hr * 1.3, 0.5, 
       sprintf("p = %.3f", p_value), 
       pos = 2, cex = 0.8, font = 2)
  text(min_hr * 0.9, nrow(meta_data) + 1, "Study", pos = 4, font = 2)
  text(max_hr * 1.1, nrow(meta_data) + 1, "HR (95% CI)", pos = 2, font = 2)
  text(max_hr * 1.3, nrow(meta_data) + 1, "Weight", pos = 2, font = 2)
  title(main = "Forest Plot: Meta-Analysis Results", 
        line = 1)
  mtext("Favors Group 1     Favors Group 2", side = 3, line = 0)
  
  dev.off()
}

#' Perform Spline Meta-Regression
#'
#' @param meta_data Data frame with columns: study, hr, lower, upper, moderator
#' @param moderator_name Name of the moderator variable
#' @param df Degrees of freedom for spline (default: 2)
#' @param output_dir Directory for output files
#' @return List with meta-regression results
#'
perform_spline_meta_regression <- function(meta_data, moderator_name, 
                                          df = 2, output_dir = "results") {
  
  # Calculate log hazard ratios and standard errors
  meta_data$loghr <- log(meta_data$hr)
  meta_data$se <- (log(meta_data$upper) - log(meta_data$lower)) / (2 * 1.96)
  
  # Check if moderator exists in data
  if (!moderator_name %in% names(meta_data)) {
    stop(paste("Moderator variable", moderator_name, "not found in meta_data"))
  }
  
  # Extract moderator values
  moderator_values <- meta_data[[moderator_name]]
  
  # Create spline basis
  spline_basis <- ns(moderator_values, df = df)
  
  # Fit spline model
  spline_model <- lm(loghr ~ spline_basis, data = meta_data, weights = 1/se^2)
  
  # Create prediction grid
  pred_grid <- data.frame(
    moderator = seq(
      min(moderator_values),
      max(moderator_values),
      length.out = 100
    )
  )
  names(pred_grid)[1] <- moderator_name
  
  # Create design matrix for predictions
  pred_basis <- predict(ns(moderator_values, df = df), 
                       newx = pred_grid[[moderator_name]])
  pred_matrix <- matrix(pred_basis, ncol = df)
  
  # Predict using spline model
  pred_fit <- predict(spline_model, 
                     newdata = data.frame(spline_basis = I(pred_matrix)), 
                     se.fit = TRUE)
  
  pred_grid$loghr <- pred_fit$fit
  pred_grid$hr <- exp(pred_grid$loghr)
  
  # Calculate confidence intervals
  pred_grid$lower <- exp(pred_grid$loghr - 1.96 * pred_fit$se.fit)
  pred_grid$upper <- exp(pred_grid$loghr + 1.96 * pred_fit$se.fit)
  
  # Create dose-response plot
  create_dose_response_plot(meta_data, pred_grid, moderator_name, output_dir)
  
  # Save results
  results <- list(
    meta_data = meta_data,
    pred_grid = pred_grid,
    spline_model = spline_model,
    moderator_name = moderator_name
  )
  
  # Calculate specific values at particular moderator points
  specific_points <- seq(
    min(moderator_values),
    max(moderator_values),
    length.out = 5
  )
  
  specific_results <- data.frame(
    moderator = specific_points,
    hr = numeric(length(specific_points)),
    lower = numeric(length(specific_points)),
    upper = numeric(length(specific_points))
  )
  names(specific_results)[1] <- moderator_name
  
  for (i in 1:length(specific_points)) {
    point <- specific_points[i]
    idx <- which.min(abs(pred_grid[[moderator_name]] - point))
    specific_results$hr[i] <- pred_grid$hr[idx]
    specific_results$lower[i] <- pred_grid$lower[idx]
    specific_results$upper[i] <- pred_grid$upper[idx]
  }
  
  results$specific_results <- specific_results
  
  # Save to CSV
  write_csv(pred_grid, file.path(output_dir, "spline_predictions.csv"))
  write_csv(specific_results, file.path(output_dir, "specific_point_results.csv"))
  
  # Save summary to text file
  sink(file.path(output_dir, "spline_regression_summary.txt"))
  cat("Spline Meta-Regression Results\n")
  cat("=============================\n\n")
  cat("Moderator variable:", moderator_name, "\n")
  cat("Degrees of freedom:", df, "\n\n")
  
  cat("Model summary:\n")
  print(summary(spline_model))
  
  cat("\nHazard ratios at specific points:\n")
  for (i in 1:nrow(specific_results)) {
    cat(moderator_name, "=", round(specific_results[[moderator_name]][i], 2),
        ": HR =", round(specific_results$hr[i], 2),
        "(95% CI:", round(specific_results$lower[i], 2), "-",
        round(specific_results$upper[i], 2), ")\n")
  }
  sink()
  
  return(results)
}

#' Create Dose-Response Plot
#'
#' @param meta_data Original meta-analysis data
#' @param pred_grid Prediction grid with fitted values
#' @param moderator_name Name of the moderator variable
#' @param output_dir Output directory
#' @return NULL (creates plot files)
#'
create_dose_response_plot <- function(meta_data, pred_grid, moderator_name, 
                                     output_dir = "results") {
  
  # Set up the PDF file
  pdf(file.path(output_dir, "dose_response_plot.pdf"), width = 10, height = 8)
  
  # Set margins
  par(mar = c(5, 5, 4, 2))
  
  # Determine y-axis range
  max_y <- max(c(meta_data$upper, pred_grid$upper), na.rm = TRUE) * 1.1
  min_y <- min(c(meta_data$lower, pred_grid$lower), na.rm = TRUE) * 0.9
  min_y <- max(0, min_y)  # Ensure not negative
  
  # Plot the curve
  plot(pred_grid[[moderator_name]], pred_grid$hr, type = "l", lwd = 2,
       xlim = range(pred_grid[[moderator_name]]),
       ylim = c(min_y, max_y),
       xlab = moderator_name,
       ylab = "Hazard Ratio",
       main = paste("Dose-Response Relationship:", moderator_name, "vs. Hazard Ratio"))
  
  # Add confidence interval
  polygon(c(pred_grid[[moderator_name]], rev(pred_grid[[moderator_name]])),
          c(pred_grid$lower, rev(pred_grid$upper)),
          col = rgb(0.5, 0.5, 0.5, 0.2), border = NA)
  
  # Add data points with error bars
  points(meta_data[[moderator_name]], meta_data$hr, pch = 19, cex = 1.2)
  segments(meta_data[[moderator_name]], meta_data$lower,
           meta_data[[moderator_name]], meta_data$upper,
           lwd = 2)
  
  # Add reference line
  abline(h = 1, lty = 2)
  
  # Add grid
  grid()
  
  # Add study labels
  if ("study" %in% names(meta_data)) {
    text(meta_data[[moderator_name]], meta_data$hr, 
         labels = meta_data$study, 
         pos = 3, offset = 0.5, cex = 0.7)
  }
  
  dev.off()
  
  # Also create a PNG version
  png(file.path(output_dir, "dose_response_plot.png"), 
      width = 1000, height = 800, res = 100)
  
  # Repeat plotting code for PNG
  par(mar = c(5, 5, 4, 2))
  plot(pred_grid[[moderator_name]], pred_grid$hr, type = "l", lwd = 2,
       xlim = range(pred_grid[[moderator_name]]),
       ylim = c(min_y, max_y),
       xlab = moderator_name,
       ylab = "Hazard Ratio",
       main = paste("Dose-Response Relationship:", moderator_name, "vs. Hazard Ratio"))
  polygon(c(pred_grid[[moderator_name]], rev(pred_grid[[moderator_name]])),
          c(pred_grid$lower, rev(pred_grid$upper)),
          col = rgb(0.5, 0.5, 0.5, 0.2), border = NA)
  points(meta_data[[moderator_name]], meta_data$hr, pch = 19, cex = 1.2)
  segments(meta_data[[moderator_name]], meta_data$lower,
           meta_data[[moderator_name]], meta_data$upper,
           lwd = 2)
  abline(h = 1, lty = 2)
  grid()
  if ("study" %in% names(meta_data)) {
    text(meta_data[[moderator_name]], meta_data$hr, 
         labels = meta_data$study, 
         pos = 3, offset = 0.5, cex = 0.7)
  }
  
  dev.off()
  
  # Create a ggplot2 version (often better for publication)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot() +
      geom_ribbon(data = pred_grid, 
                 aes(x = .data[[moderator_name]], 
                     ymin = lower, ymax = upper), 
                 alpha = 0.2) +
      geom_line(data = pred_grid, 
               aes(x = .data[[moderator_name]], y = hr), 
               size = 1) +
      geom_point(data = meta_data, 
                aes(x = .data[[moderator_name]], y = hr), 
                size = 3) +
      geom_errorbar(data = meta_data, 
                   aes(x = .data[[moderator_name]], 
                       ymin = lower, ymax = upper), 
                   width = 0.2) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      labs(title = paste("Dose-Response Relationship"),
           subtitle = paste(moderator_name, "vs. Hazard Ratio"),
           x = moderator_name,
           y = "Hazard Ratio") +
      theme_minimal() +
      theme(
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )
    
    # Add study labels if available
    if ("study" %in% names(meta_data)) {
      p <- p + geom_text(data = meta_data, 
                        aes(x = .data[[moderator_name]], 
                            y = hr, 
                            label = study), 
                        nudge_y = 0.1, size = 3)
    }
    
    # Save plot
    ggsave(file.path(output_dir, "dose_response_ggplot.pdf"), 
           p, width = 10, height = 8)
    ggsave(file.path(output_dir, "dose_response_ggplot.png"), 
           p, width = 10, height = 8, dpi = 100)
  }
}

#=====================================================================
# EXAMPLE USAGE
#=====================================================================

# This section demonstrates how to use the functions above
# with example data. Uncomment and modify for your specific needs.

#' # Example: Load digitized curve data
#' example_data <- read_curve_data("path/to/digitized_curve.csv")
#' 
#' # Example: Reconstruct IPD
#' example_ipd <- reconstruct_ipd(
#'   example_data$time,
#'   example_data$surv,
#'   n_risk_initial = 100,
#'   total_events = 30,
#'   study_name = "Example Study",
#'   arm_name = "Group A"
#' )
#' 
#' # Example: Calculate hazard ratio between two groups
#' group1_ipd <- reconstruct_ipd(
#'   time_points = c(0, 6, 12, 18, 24),
#'   surv_probs = c(1, 0.9, 0.8, 0.7, 0.6),
#'   n_risk_initial = 100,
#'   total_events = 40,
#'   study_name = "Example",
#'   arm_name = "group1"
#' )$ipd
#' 
#' group2_ipd <- reconstruct_ipd(
#'   time_points = c(0, 6, 12, 18, 24),
#'   surv_probs = c(1, 0.8, 0.6, 0.5, 0.4),
#'   n_risk_initial = 100,
#'   total_events = 60,
#'   study_name = "Example",
#'   arm_name = "group2"
#' )$ipd
#' 
#' hr_result <- calculate_hazard_ratio(
#'   group1_ipd, 
#'   group2_ipd, 
#'   "Control", 
#'   "Treatment"
#' )
#' 
#' # Example: Create meta-analysis data frame
#' meta_data <- data.frame(
#'   study = c("Study 1", "Study 2", "Study 3"),
#'   hr = c(1.5, 2.1, 1.8),
#'   lower = c(1.1, 1.6, 1.2),
#'   upper = c(2.0, 2.8, 2.7),
#'   moderator = c(-16, -14, -12)  # e.g., GLS values
#' )
#' 
#' # Example: Perform meta-analysis
#' meta_results <- perform_fixed_effects_meta(meta_data)
#' 
#' # Example: Perform spline meta-regression
#' spline_results <- perform_spline_meta_regression(
#'   meta_data, 
#'   moderator_name = "moderator"
#' )

#=====================================================================
# STEP-BY-STEP WORKFLOW FOR YOUR PROJECT
#=====================================================================

# 1. Define your studies and digitized curve files
# Replace this with your actual study data

#' study_files <- list(
#'   study1 = list(
#'     reduced = "data/study1_reduced.csv",
#'     preserved = "data/study1_preserved.csv",
#'     n_risk_reduced = 50,
#'     n_risk_preserved = 52,
#'     events_reduced = 22,
#'     events_preserved = 17
#'   ),
#'   study2 = list(
#'     reduced = "data/study2_reduced.csv",
#'     preserved = "data/study2_preserved.csv",
#'     n_risk_reduced = 48,
#'     n_risk_preserved = 49,
#'     events_reduced = 20,
#'     events_preserved = 15
#'   )
#' )
#' 
#' # 2. Process each study
#' results <- list()
#' meta_data <- data.frame(
#'   study = character(),
#'   hr = numeric(),
#'   lower = numeric(),
#'   upper = numeric(),
#'   p_value = numeric(),
#'   moderator = numeric(),  # e.g., GLS threshold
#'   stringsAsFactors = FALSE
#' )
#' 
#' # Process binary comparison studies
#' for (study_name in names(study_files)) {
#'   study_data <- study_files[[study_name]]
#'   
#'   # Read curve data
#'   reduced_data <- read_curve_data(study_data$reduced)
#'   preserved_data <- read_curve_data(study_data$preserved)
#'   
#'   # Reconstruct IPD
#'   reduced_result <- reconstruct_ipd(
#'     reduced_data$time,
#'     reduced_data$surv,
#'     study_data$n_risk_reduced,
#'     study_data$events_reduced,
#'     study_name,
#'     "reduced"
#'   )
#'   
#'   preserved_result <- reconstruct_ipd(
#'     preserved_data$time,
#'     preserved_data$surv,
#'     study_data$n_risk_preserved,
#'     study_data$events_preserved,
#'     study_name,
#'     "preserved"
#'   )
#'   
#'   # Calculate hazard ratio
#'   hr_result <- calculate_hazard_ratio(
#'     preserved_result$ipd,  # Use preserved as reference
#'     reduced_result$ipd,
#'     "preserved",
#'     "reduced"
#'   )
#'   
#'   # Create KM plot
#'   create_km_plot(
#'     hr_result$combined_data,
#'     group_var = "group",
#'     title = paste("Survival Comparison -", study_name),
#'     group_labels = c("Preserved GLS", "Reduced GLS"),
#'     filename = paste0(study_name, "_km_comparison")
#'   )
#'   
#'   # Store results
#'   results[[study_name]] <- list(
#'     reduced = reduced_result,
#'     preserved = preserved_result,
#'     hr_result = hr_result
#'   )
#'   
#'   # Add to meta-analysis data
#'   meta_data <- rbind(meta_data, data.frame(
#'     study = study_name,
#'     hr = hr_result$hr,
#'     lower = hr_result$lower,
#'     upper = hr_result$upper,
#'     p_value = hr_result$p_value,
#'     moderator = -16,  # Replace with actual GLS threshold for this study
#'     stringsAsFactors = FALSE
#'   ))
#' }
#' 
#' # 3. Perform meta-analysis
#' meta_results <- perform_fixed_effects_meta(meta_data)
#' 
#' # 4. Perform spline meta-regression
#' spline_results <- perform_spline_meta_regression(
#'   meta_data,
#'   moderator_name = "moderator"  # e.g., GLS values
#' )
#' 
#' # 5. Summary of results
#' cat("\nMeta-Analysis Results:\n")
#' cat("Pooled HR: ", round(meta_results$pooled_hr, 2),
#'     " (95% CI: ", round(meta_results$ci_lower, 2),
#'     " - ", round(meta_results$ci_upper, 2), ")\n")
#' cat("P-value: ", format.pval(meta_results$p_value, digits = 3), "\n\n")
#' 
#' cat("Dose-Response Analysis - HR at specific points:\n")
#' specific_results <- spline_results$specific_results
#' for (i in 1:nrow(specific_results)) {
#'   cat("GLS = ", round(specific_results$moderator[i], 1),
#'       "%: HR = ", round(specific_results$hr[i], 2),
#'       " (95% CI: ", round(specific_results$lower[i], 2),
#'       " - ", round(specific_results$upper[i], 2), ")\n")
#' }
