# ===========================================================================
# R/utils.R - Utility functions for sopravvivenza package
# ===========================================================================

#' Check if Required Packages are Installed
#'
#' This function checks if the required packages are installed and offers
#' to install them if they are missing.
#'
#' @param packages Character vector of package names to check
#' @param prompt Whether to prompt user before installing (default: TRUE)
#'
#' @return Logical vector indicating which packages are installed
#' @keywords internal
check_packages <- function(packages, prompt = TRUE) {
  installed <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  
  if (all(installed)) {
    return(invisible(installed))
  }
  
  missing_pkgs <- packages[!installed]
  msg <- sprintf(
    "The following packages are required but not installed: %s\n",
    paste(missing_pkgs, collapse = ", ")
  )
  
  if (prompt) {
    msg <- paste0(msg, "Would you like to install them now? (Y/n)")
    answer <- readline(msg)
    
    if (tolower(substr(answer, 1, 1)) != "n") {
      install.packages(missing_pkgs)
      installed <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
    }
  } else {
    warning(msg)
  }
  
  return(invisible(installed))
}

#' Extract Numbers at Risk from Kaplan-Meier Plot Table
#'
#' This function extracts number-at-risk data from a digitized risk table
#' commonly found below Kaplan-Meier plots.
#'
#' @param risk_data Data frame with risk data (time and numbers)
#' @param time_points Vector of time points to match with risk data
#' @param n_initial Initial number at risk (for time 0)
#' @param interpolate Whether to interpolate missing values (default: TRUE)
#'
#' @return Vector of numbers at risk at each time point
#' @export
extract_risk_numbers <- function(risk_data, time_points, n_initial, interpolate = TRUE) {
  # Ensure risk_data has time and n_risk columns
  if (!all(c("time", "n_risk") %in% names(risk_data))) {
    stop("risk_data must contain 'time' and 'n_risk' columns")
  }
  
  # Add time 0 if not present
  if (min(risk_data$time) > 0) {
    risk_data <- rbind(data.frame(time = 0, n_risk = n_initial), risk_data)
  }
  
  # Sort by time
  risk_data <- risk_data[order(risk_data$time), ]
  
  # Initialize output vector
  risk_numbers <- numeric(length(time_points))
  
  # Fill in known values
  for (i in seq_along(time_points)) {
    t <- time_points[i]
    
    # Find exact match
    exact_match <- which(risk_data$time == t)
    if (length(exact_match) > 0) {
      risk_numbers[i] <- risk_data$n_risk[exact_match[1]]
      next
    }
    
    # If interpolation is requested, find closest points and interpolate
    if (interpolate) {
      lower_idx <- max(which(risk_data$time <= t), na.rm = TRUE)
      upper_idx <- min(which(risk_data$time >= t), na.rm = TRUE)
      
      if (length(lower_idx) > 0 && length(upper_idx) > 0 && lower_idx != upper_idx) {
        # Linear interpolation
        t_lower <- risk_data$time[lower_idx]
        t_upper <- risk_data$time[upper_idx]
        n_lower <- risk_data$n_risk[lower_idx]
        n_upper <- risk_data$n_risk[upper_idx]
        
        # Interpolate
        risk_numbers[i] <- n_lower + (n_upper - n_lower) * (t - t_lower) / (t_upper - t_lower)
        risk_numbers[i] <- round(risk_numbers[i])
      } else if (length(lower_idx) > 0) {
        # Use lower bound
        risk_numbers[i] <- risk_data$n_risk[lower_idx]
      } else if (length(upper_idx) > 0) {
        # Use upper bound
        risk_numbers[i] <- risk_data$n_risk[upper_idx]
      }
    }
  }
  
  return(risk_numbers)
}

#' Calculate Survival Probabilities from IPD
#'
#' This function calculates survival probabilities at specific time points
#' from reconstructed individual patient data.
#'
#' @param ipd Individual patient data with time and status columns
#' @param time_points Vector of time points at which to calculate survival probabilities
#' @param time_var Name of time variable in ipd (default: "time")
#' @param status_var Name of status variable in ipd (default: "status")
#'
#' @return Data frame with time points and corresponding survival probabilities
#' @export
calculate_survival_from_ipd <- function(ipd, time_points, 
                                       time_var = "time", status_var = "status") {
  
  # Fit survival object
  surv_formula <- stats::as.formula(paste("survival::Surv(", time_var, ",", status_var, ") ~ 1"))
  fit <- survival::survfit(surv_formula, data = ipd)
  
  # Extract survival probabilities at specified time points
  surv_summary <- summary(fit, times = time_points)
  
  # Create data frame with results
  result <- data.frame(
    time = time_points,
    surv = surv_summary$surv
  )
  
  # Handle time points beyond the last event
  na_idx <- is.na(result$surv)
  if (any(na_idx)) {
    # For times beyond the last event, use the last available survival probability
    last_surv <- fit$surv[length(fit$surv)]
    result$surv[na_idx] <- last_surv
  }
  
  return(result)
}

#' Calculate Confidence Intervals for Reconstructed Curves
#'
#' This function calculates confidence intervals for survival probabilities
#' from reconstructed individual patient data.
#'
#' @param ipd Individual patient data with time and status columns
#' @param time_points Vector of time points at which to calculate confidence intervals
#' @param conf_level Confidence level (default: 0.95)
#' @param time_var Name of time variable in ipd (default: "time")
#' @param status_var Name of status variable in ipd (default: "status")
#' @param interval_type Type of confidence interval: "log" (default), "log-log", "plain", or "arcsin"
#'
#' @return Data frame with time points, survival probabilities, and confidence limits
#' @export
calculate_survival_ci <- function(ipd, time_points, conf_level = 0.95,
                                 time_var = "time", status_var = "status",
                                 interval_type = "log") {
  
  # Fit survival object
  surv_formula <- stats::as.formula(paste("survival::Surv(", time_var, ",", status_var, ") ~ 1"))
  fit <- survival::survfit(surv_formula, data = ipd, conf.type = interval_type,
                         conf.int = conf_level)
  
  # Extract survival probabilities and confidence intervals at specified time points
  surv_summary <- summary(fit, times = time_points)
  
  # Create data frame with results
  result <- data.frame(
    time = time_points,
    surv = surv_summary$surv,
    lower = surv_summary$lower,
    upper = surv_summary$upper
  )
  
  # Handle time points beyond the last event
  na_idx <- is.na(result$surv)
  if (any(na_idx)) {
    # For times beyond the last event, use the last available values
    last_idx <- max(which(!is.na(result$surv)))
    if (length(last_idx) > 0) {
      result$surv[na_idx] <- result$surv[last_idx]
      result$lower[na_idx] <- result$lower[last_idx]
      result$upper[na_idx] <- result$upper[last_idx]
    }
  }
  
  return(result)
}

#' Estimate Total Events from Survival Curve
#'
#' This function estimates the total number of events from a digitized survival curve
#' and initial number at risk.
#'
#' @param time_points Vector of time points from digitized curve
#' @param surv_probs Vector of survival probabilities from digitized curve
#' @param n_risk_initial Initial number at risk
#'
#' @return Estimated total number of events
#' @export
estimate_total_events <- function(time_points, surv_probs, n_risk_initial) {
  
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
  
  # Calculate estimated events in each interval
  total_events <- 0
  n <- n_risk_initial
  
  for (i in 1:(length(time_points) - 1)) {
    s1 <- surv_probs[i]
    s2 <- surv_probs[i + 1]
    
    if (s1 > 0) {
      h <- (s1 - s2) / s1
    } else {
      h <- 0
    }
    
    events <- round(n * h)
    total_events <- total_events + events
    
    # Update number at risk
    n <- n - events
  }
  
  return(total_events)
}

#' Compare Original and Reconstructed Survival Curves
#'
#' This function compares the original digitized survival curve with
#' the reconstructed curve from IPD to assess reconstruction quality.
#'
#' @param original_data Data frame with original time and surv columns
#' @param reconstructed_data Result from reconstruct_ipd function
#' @param time_var Name of time variable in IPD (default: "time")
#' @param status_var Name of status variable in IPD (default: "status")
#'
#' @return Data frame with comparison metrics
#' @export
#'
#' @examples
#' \dontrun{
#' # Read digitized curve data
#' curve_data <- read_curve_data("curve.csv")
#' 
#' # Reconstruct IPD
#' result <- reconstruct_ipd(
#'   curve_data$time,
#'   curve_data$surv,
#'   n_risk_initial = 100
#' )
#' 
#' # Compare original and reconstructed curves
#' comparison <- compare_curves(
#'   data.frame(time = curve_data$time, surv = curve_data$surv),
#'   result
#' )
#' }
compare_curves <- function(original_data, reconstructed_data,
                          time_var = "time", status_var = "status") {
  
  # Extract original time points and survival probabilities
  orig_time <- original_data$time
  orig_surv <- original_data$surv
  
  # Extract IPD
  ipd <- reconstructed_data$ipd
  
  # Calculate survival probabilities from IPD at original time points
  surv_formula <- stats::as.formula(paste("survival::Surv(", time_var, ",", status_var, ") ~ 1"))
  fit <- survival::survfit(surv_formula, data = ipd)
  
  recon_surv <- summary(fit, times = orig_time)$surv
  
  # Handle NAs in reconstructed survival (for times beyond last event)
  na_idx <- is.na(recon_surv)
  if (any(na_idx)) {
    last_surv <- fit$surv[length(fit$surv)]
    recon_surv[na_idx] <- last_surv
  }
  
  # Calculate absolute differences
  abs_diff <- abs(orig_surv - recon_surv)
  
  # Calculate root mean squared error
  rmse <- sqrt(mean(abs_diff^2))
  
  # Calculate maximum absolute difference
  max_diff <- max(abs_diff)
  
  # Calculate mean absolute difference
  mean_diff <- mean(abs_diff)
  
  # Create data frame with comparison results
  comparison <- data.frame(
    time = orig_time,
    original_surv = orig_surv,
    reconstructed_surv = recon_surv,
    absolute_diff = abs_diff
  )
  
  # Add summary metrics
  attr(comparison, "rmse") <- rmse
  attr(comparison, "max_diff") <- max_diff
  attr(comparison, "mean_diff") <- mean_diff
  
  return(comparison)
}

#' Validate Reconstruction Quality
#'
#' This function validates the quality of IPD reconstruction by comparing
#' the original and reconstructed survival curves and returns quality metrics.
#'
#' @param original_data Data frame with original time and surv columns
#' @param reconstructed_data Result from reconstruct_ipd function
#' @param time_var Name of time variable in IPD (default: "time")
#' @param status_var Name of status variable in IPD (default: "status")
#' @param threshold RMSE threshold for acceptable reconstruction (default: 0.02)
#'
#' @return List with validation results and quality metrics
#' @export
#'
#' @examples
#' \dontrun{
#' # Read digitized curve data
#' curve_data <- read_curve_data("curve.csv")
#' 
#' # Reconstruct IPD
#' result <- reconstruct_ipd(
#'   curve_data$time,
#'   curve_data$surv,
#'   n_risk_initial = 100
#' )
#' 
#' # Validate reconstruction quality
#' validation <- validate_reconstruction(
#'   data.frame(time = curve_data$time, surv = curve_data$surv),
#'   result
#' )
#' }
validate_reconstruction <- function(original_data, reconstructed_data,
                                   time_var = "time", status_var = "status",
                                   threshold = 0.02) {
  
  # Compare curves
  comparison <- compare_curves(original_data, reconstructed_data,
                             time_var = time_var, status_var = status_var)
  
  # Extract quality metrics
  rmse <- attr(comparison, "rmse")
  max_diff <- attr(comparison, "max_diff")
  mean_diff <- attr(comparison, "mean_diff")
  
  # Determine if reconstruction quality is acceptable
  acceptable <- rmse <= threshold
  
  # Create validation results
  result <- list(
    acceptable = acceptable,
    rmse = rmse,
    max_diff = max_diff,
    mean_diff = mean_diff,
    comparison = comparison,
    threshold = threshold
  )
  
  # Add quality grade
  if (rmse <= 0.01) {
    result$quality <- "Excellent"
  } else if (rmse <= 0.02) {
    result$quality <- "Good"
  } else if (rmse <= 0.05) {
    result$quality <- "Fair"
  } else {
    result$quality <- "Poor"
  }
  
  return(result)
}

#' Extract Event Times from Reconstructed IPD
#'
#' This function extracts the times of events from reconstructed IPD,
#' which can be useful for further analysis or visualization.
#'
#' @param ipd Individual patient data with time and status columns
#' @param time_var Name of time variable in IPD (default: "time")
#' @param status_var Name of status variable in IPD (default: "status")
#'
#' @return Vector of event times
#' @export
#'
#' @examples
#' \dontrun{
#' # Read digitized curve data
#' curve_data <- read_curve_data("curve.csv")
#' 
#' # Reconstruct IPD
#' result <- reconstruct_ipd(
#'   curve_data$time,
#'   curve_data$surv,
#'   n_risk_initial = 100
#' )
#' 
#' # Extract event times
#' event_times <- extract_event_times(result$ipd)
#' }
extract_event_times <- function(ipd, time_var = "time", status_var = "status") {
  # Ensure time_var and status_var exist in ipd
  if (!all(c(time_var, status_var) %in% names(ipd))) {
    stop("IPD must contain '", time_var, "' and '", status_var, "' columns")
  }
  
  # Extract event times
  event_times <- ipd[[time_var]][ipd[[status_var]] == 1]
  
  # Sort event times
  event_times <- sort(event_times)
  
  return(event_times)
}

#' Convert Between Time Units
#'
#' This function converts time values between different units,
#' which can be useful when dealing with curves using different time scales.
#'
#' @param time Numeric vector of time values to convert
#' @param from_unit Original time unit ("days", "weeks", "months", "years")
#' @param to_unit Target time unit ("days", "weeks", "months", "years")
#'
#' @return Numeric vector of converted time values
#' @export
#'
#' @examples
#' # Convert 24 months to years
#' convert_time_units(24, from_unit = "months", to_unit = "years")
#' 
#' # Convert 5 years to months
#' convert_time_units(5, from_unit = "years", to_unit = "months")
convert_time_units <- function(time, from_unit, to_unit) {
  # Define conversion factors (to days)
  to_days <- list(
    days = 1,
    weeks = 7,
    months = 30.44,  # Average month length
    years = 365.25   # Average year length including leap years
  )
  
  # Check if units are valid
  if (!from_unit %in% names(to_days) || !to_unit %in% names(to_days)) {
    stop("Units must be one of: days, weeks, months, years")
  }
  
  # Convert to days first
  days <- time * to_days[[from_unit]]
  
  # Convert from days to target unit
  result <- days / to_days[[to_unit]]
  
  return(result)
}

#' Check Data Consistency
#'
#' This function checks if digitized survival data is consistent and valid.
#'
#' @param time_points Vector of time points from digitized curve
#' @param surv_probs Vector of survival probabilities from digitized curve
#'
#' @return List with validation results and messages
#' @export
#'
#' @examples
#' # Check data consistency
#' time_points <- c(0, 6, 12, 18, 24)
#' surv_probs <- c(1, 0.9, 0.8, 0.75, 0.7)
#' check_data_consistency(time_points, surv_probs)
check_data_consistency <- function(time_points, surv_probs) {
  results <- list(
    valid = TRUE,
    messages = character(0)
  )
  
  # Check for equal length
  if (length(time_points) != length(surv_probs)) {
    results$valid <- FALSE
    results$messages <- c(results$messages, 
                        "Time points and survival probabilities must have the same length")
  }
  
  # Check for missing values
  if (any(is.na(time_points)) || any(is.na(surv_probs))) {
    results$valid <- FALSE
    results$messages <- c(results$messages, 
                        "Data contains missing values")
  }
  
  # Check for valid survival probabilities (between 0 and 1)
  if (any(surv_probs < 0 | surv_probs > 1, na.rm = TRUE)) {
    results$valid <- FALSE
    results$messages <- c(results$messages, 
                        "Survival probabilities must be between 0 and 1")
  }
  
  # Check for time point ordering
  if (!all(diff(time_points) >= 0, na.rm = TRUE)) {
    results$valid <- FALSE
    results$messages <- c(results$messages, 
                        "Time points must be in non-decreasing order")
  }
  
  # Check for survival probability non-increasing
  if (!all(diff(surv_probs) <= 0, na.rm = TRUE)) {
    results$valid <- FALSE
    results$messages <- c(results$messages, 
                        "Survival probabilities must be non-increasing")
  }
  
  # Check if curve starts with survival = 1
  if (surv_probs[1] != 1) {
    results$messages <- c(results$messages, 
                        "Warning: First survival probability is not 1; will be adjusted")
  }
  
  # Check if first time point is 0
  if (time_points[1] != 0) {
    results$messages <- c(results$messages, 
                        "Warning: First time point is not 0; will be adjusted")
  }
  
  return(results)
}

#' Format Hazard Ratio with Confidence Interval
#'
#' This function formats a hazard ratio and its confidence interval
#' as a string, suitable for reporting in publications.
#'
#' @param hr Hazard ratio value
#' @param lower Lower confidence limit
#' @param upper Upper confidence limit
#' @param digits Number of decimal places (default: 2)
#' @param separator String separating HR and CI (default: " ")
#' @param ci_format Format for confidence interval (default: "(95% CI: %s-%s)")
#'
#' @return Formatted string with HR and CI
#' @export
#'
#' @examples
#' # Format hazard ratio
#' format_hr(1.75, 1.25, 2.45)
#' # "1.75 (95% CI: 1.25-2.45)"
#' 
#' # Custom format
#' format_hr(1.75, 1.25, 2.45, digits = 3, ci_format = "[%s, %s]")
#' # "1.750 [1.250, 2.450]"
format_hr <- function(hr, lower, upper, digits = 2, 
                     separator = " ", ci_format = "(95% CI: %s-%s)") {
  
  # Format numbers with specified decimal places
  hr_str <- sprintf(paste0("%.", digits, "f"), hr)
  lower_str <- sprintf(paste0("%.", digits, "f"), lower)
  upper_str <- sprintf(paste0("%.", digits, "f"), upper)
  
  # Format confidence interval
  ci_str <- sprintf(ci_format, lower_str, upper_str)
  
  # Combine HR and CI
  result <- paste0(hr_str, separator, ci_str)
  
  return(result)
}

#' Extract Median Survival Time
#'
#' This function extracts the median survival time from digitized curve data
#' or from reconstructed IPD.
#'
#' @param time_points Vector of time points from digitized curve
#' @param surv_probs Vector of survival probabilities from digitized curve
#' @param ipd Optional individual patient data (if provided, time_points and surv_probs are ignored)
#' @param time_var Name of time variable in IPD (default: "time")
#' @param status_var Name of status variable in IPD (default: "status")
#' @param interpolate Whether to interpolate median survival time (default: TRUE)
#'
#' @return Median survival time
#' @export
#'
#' @examples
#' \dontrun{
#' # From digitized curve data
#' time_points <- c(0, 6, 12, 18, 24, 30, 36)
#' surv_probs <- c(1, 0.9, 0.8, 0.7, 0.5, 0.4, 0.3)
#' median_survival <- extract_median_survival(time_points, surv_probs)
#' 
#' # From reconstructed IPD
#' median_survival <- extract_median_survival(ipd = result$ipd)
#' }
extract_median_survival <- function(time_points = NULL, surv_probs = NULL, ipd = NULL,
                                   time_var = "time", status_var = "status",
                                   interpolate = TRUE) {
  
  # Check if IPD is provided
  if (!is.null(ipd)) {
    # Calculate median survival from IPD
    surv_formula <- stats::as.formula(paste("survival::Surv(", time_var, ",", status_var, ") ~ 1"))
    fit <- survival::survfit(surv_formula, data = ipd)
    
    # Extract median survival time
    median_survival <- stats::median(fit)
    
    return(median_survival)
  }
  
  # Check if digitized curve data is provided
  if (is.null(time_points) || is.null(surv_probs)) {
    stop("Either IPD or both time_points and surv_probs must be provided")
  }
  
  # Ensure time points and survival probabilities are sorted
  idx <- order(time_points)
  time_points <- time_points[idx]
  surv_probs <- surv_probs[idx]
  
  # Find median survival
  median_idx <- which(surv_probs <= 0.5)[1]
  
  # Check if median can be found
  if (is.na(median_idx)) {
    # Median survival is beyond the curve's range
    return(NA)
  }
  
  if (median_idx == 1) {
    # Median survival is at or before the first time point
    return(time_points[1])
  }
  
  # Get time points and survival probabilities around median
  t_before <- time_points[median_idx - 1]
  t_after <- time_points[median_idx]
  s_before <- surv_probs[median_idx - 1]
  s_after <- surv_probs[median_idx]
  
  # Check if survival is exactly 0.5 at the time point
  if (s_after == 0.5) {
    return(t_after)
  }
  
  # Interpolate median survival time if requested
  if (interpolate && s_before > 0.5 && s_after < 0.5) {
    # Linear interpolation
    median_time <- t_before + (t_after - t_before) * (0.5 - s_before) / (s_after - s_before)
    return(median_time)
  }
  
  # Otherwise, return the time point where survival first drops below 0.5
  return(t_after)
}