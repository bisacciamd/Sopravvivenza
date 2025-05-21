# ===========================================================================
# R/read_data.R
# ===========================================================================

#' Read and Preprocess Digitized Curve Data
#'
#' This function reads CSV files containing digitized Kaplan-Meier curve data
#' and preprocesses them for subsequent analysis.
#'
#' @param file_path Path to CSV file with digitized curve data
#' @param time_col Column number or name for time values (default: 1)
#' @param surv_col Column number or name for survival probabilities (default: 2)
#' @param header Whether file has a header row (default: FALSE)
#' @param convert_percentage Whether to convert survival percentages to proportions (default: TRUE)
#' @param time_scale Character string indicating time scale. If "auto" (default), 
#'        will attempt to detect scale. Other options: "months", "years", "days", "weeks"
#'
#' @return A dataframe with preprocessed time and survival probability data
#' @export
#'
#' @examples
#' \dontrun{
#' curve_data <- read_curve_data("path/to/digitized_curve.csv")
#' curve_data <- read_curve_data("path/to/digitized_curve.csv", 
#'                              time_col = "time", 
#'                              surv_col = "survival",
#'                              header = TRUE)
#' }
read_curve_data <- function(file_path, time_col = 1, surv_col = 2, header = FALSE,
                           convert_percentage = TRUE, time_scale = "auto") {
  # Check if file exists
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  # Read the CSV file
  data <- tryCatch({
    if (header) {
      readr::read_csv(file_path, show_col_types = FALSE)
    } else {
      readr::read_csv(file_path, col_names = FALSE, show_col_types = FALSE)
    }
  }, error = function(e) {
    stop("Error reading file: ", file_path, "\nError message: ", e$message)
  })
  
  # Check if data has sufficient columns
  if (ncol(data) < max(time_col, surv_col)) {
    stop("File does not have enough columns. Required: ", max(time_col, surv_col),
         ", Found: ", ncol(data))
  }
  
  # Extract time and survival data
  if (is.numeric(time_col) && is.numeric(surv_col)) {
    # Using column indices
    time_data <- data[[time_col]]
    surv_data <- data[[surv_col]]
  } else {
    # Using column names
    time_data <- data[[time_col]]
    surv_data <- data[[surv_col]]
  }
  
  # Create a new dataframe with only time and survival data
  processed_data <- data.frame(time = time_data, surv = surv_data)
  
  # Ensure survival probability is between 0 and 1
  if (convert_percentage && max(processed_data$surv, na.rm = TRUE) > 1.1) {
    processed_data$surv <- processed_data$surv / 100
  }
  
  # Handle time scale
  if (time_scale == "auto") {
    # Attempt to detect time scale
    if (max(processed_data$time, na.rm = TRUE) < 10) {
      # Likely years, convert to months
      processed_data$time <- processed_data$time * 12
    }
  } else if (time_scale == "years") {
    # Convert years to months
    processed_data$time <- processed_data$time * 12
  } else if (time_scale == "days") {
    # Convert days to months (approximate)
    processed_data$time <- processed_data$time / 30.44
  } else if (time_scale == "weeks") {
    # Convert weeks to months (approximate)
    processed_data$time <- processed_data$time / 4.35
  }
  
  # Sort by time
  processed_data <- processed_data[order(processed_data$time), ]
  
  # Ensure the curve starts at time 0 and surv = 1
  if (min(processed_data$time, na.rm = TRUE) > 0) {
    processed_data <- rbind(data.frame(time = 0, surv = 1), processed_data)
  } else {
    # If there's already a time 0, ensure surv = 1
    processed_data$surv[processed_data$time == 0] <- 1
  }
  
  return(processed_data)
}

#' Read Number-at-Risk Data
#'
#' This function reads data on the number of patients at risk at specific time points,
#' which is often presented in tables below Kaplan-Meier curves.
#'
#' @param file_path Path to CSV file with risk data
#' @param time_col Column number or name for time values (default: 1)
#' @param risk_col Column number or name for number at risk values (default: 2)
#' @param header Whether file has a header row (default: FALSE)
#' @param time_scale Character string indicating time scale. If "auto" (default), 
#'        will attempt to detect scale. Other options: "months", "years", "days", "weeks"
#'
#' @return A dataframe with time points and corresponding numbers at risk
#' @export
#'
#' @examples
#' \dontrun{
#' risk_data <- read_risk_data("path/to/risk_data.csv")
#' }
read_risk_data <- function(file_path, time_col = 1, risk_col = 2, header = FALSE,
                          time_scale = "auto") {
  # Check if file exists
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  # Read the CSV file
  data <- tryCatch({
    if (header) {
      readr::read_csv(file_path, show_col_types = FALSE)
    } else {
      readr::read_csv(file_path, col_names = FALSE, show_col_types = FALSE)
    }
  }, error = function(e) {
    stop("Error reading file: ", file_path, "\nError message: ", e$message)
  })
  
  # Check if data has sufficient columns
  if (ncol(data) < max(time_col, risk_col)) {
    stop("File does not have enough columns. Required: ", max(time_col, risk_col),
         ", Found: ", ncol(data))
  }
  
  # Extract time and risk data
  if (is.numeric(time_col) && is.numeric(risk_col)) {
    # Using column indices
    time_data <- data[[time_col]]
    risk_data <- data[[risk_col]]
  } else {
    # Using column names
    time_data <- data[[time_col]]
    risk_data <- data[[risk_col]]
  }
  
  # Create a new dataframe with only time and risk data
  processed_data <- data.frame(time = time_data, n_risk = risk_data)
  
  # Handle time scale
  if (time_scale == "auto") {
    # Attempt to detect time scale
    if (max(processed_data$time, na.rm = TRUE) < 10) {
      # Likely years, convert to months
      processed_data$time <- processed_data$time * 12
    }
  } else if (time_scale == "years") {
    # Convert years to months
    processed_data$time <- processed_data$time * 12
  } else if (time_scale == "days") {
    # Convert days to months (approximate)
    processed_data$time <- processed_data$time / 30.44
  } else if (time_scale == "weeks") {
    # Convert weeks to months (approximate)
    processed_data$time <- processed_data$time / 4.35
  }
  
  # Sort by time
  processed_data <- processed_data[order(processed_data$time), ]
  
  return(processed_data)
}
