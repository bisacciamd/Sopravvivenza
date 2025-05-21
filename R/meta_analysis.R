# ===========================================================================
# R/meta_analysis.R
# ===========================================================================

#' Perform Fixed-Effects Meta-Analysis
#'
#' This function performs a fixed-effects meta-analysis using the inverse variance method
#' on hazard ratios from multiple studies.
#'
#' @param meta_data Data frame with columns: study, hr, lower, upper. Can also include
#'        additional columns that will be preserved in the output.
#' @param log_transform Whether to log-transform hazard ratios for analysis (default: TRUE)
#' @param conf_level Confidence level for intervals (default: 0.95)
#' @param add_weights Whether to add weight columns to output data (default: TRUE)
#'
#' @return List with meta-analysis results, including pooled estimate and data with weights
#' @export
#'
#' @examples
#' \dontrun{
#' meta_data <- data.frame(
#'   study = c("Study 1", "Study 2", "Study 3"),
#'   hr = c(1.5, 2.1, 1.8),
#'   lower = c(1.1, 1.6, 1.2),
#'   upper = c(2.0, 2.8, 2.7)
#' )
#' 
#' meta_results <- perform_fixed_effects_meta(meta_data)
#' print(meta_results$pooled_hr)  # Pooled hazard ratio
#' }
perform_fixed_effects_meta <- function(meta_data, log_transform = TRUE, 
                                      conf_level = 0.95, add_weights = TRUE) {
  
  # Validate input data
  required_cols <- c("study", "hr", "lower", "upper")
  missing_cols <- setdiff(required_cols, names(meta_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in meta_data: ", paste(missing_cols, collapse = ", "))
  }
  
  # Ensure numeric types for calculations
  meta_data$hr <- as.numeric(meta_data$hr)
  meta_data$lower <- as.numeric(meta_data$lower)
  meta_data$upper <- as.numeric(meta_data$upper)
  
  # Calculate log hazard ratios and standard errors
  if (log_transform) {
    meta_data$loghr <- log(meta_data$hr)
    meta_data$se <- (log(meta_data$upper) - log(meta_data$lower)) / (2 * stats::qnorm((1 + conf_level) / 2))
  } else {
    meta_data$loghr <- meta_data$hr
    meta_data$se <- (meta_data$upper - meta_data$lower) / (2 * stats::qnorm((1 + conf_level) / 2))
  }
  
  # Calculate weights
  if (add_weights) {
    meta_data$weight <- 1 / meta_data$se^2
    meta_data$weight_percent <- 100 * meta_data$weight / sum(meta_data$weight)
  }
  
  # Calculate pooled estimate
  weighted_mean <- sum(meta_data$loghr * meta_data$weight) / sum(meta_data$weight)
  weighted_se <- sqrt(1 / sum(meta_data$weight))
  
  # Calculate confidence interval for pooled estimate
  alpha <- 1 - conf_level
  z_value <- stats::qnorm(1 - alpha/2)
  ci_lower <- weighted_mean - z_value * weighted_se
  ci_upper <- weighted_mean + z_value * weighted_se
  
  # Calculate z-score and p-value
  z_score <- weighted_mean / weighted_se
  p_value <- 2 * (1 - stats::pnorm(abs(z_score)))
  
  # Transform back if necessary
  if (log_transform) {
    pooled_hr <- exp(weighted_mean)
    pooled_ci_lower <- exp(ci_lower)
    pooled_ci_upper <- exp(ci_upper)
  } else {
    pooled_hr <- weighted_mean
    pooled_ci_lower <- ci_lower
    pooled_ci_upper <- ci_upper
  }
  
  # Calculate heterogeneity statistics
  q_stat <- sum(meta_data$weight * (meta_data$loghr - weighted_mean)^2)
  df <- nrow(meta_data) - 1
  p_heterogeneity <- 1 - stats::pchisq(q_stat, df)
  i_squared <- max(0, 100 * (q_stat - df) / q_stat)
  if (q_stat <= df) i_squared <- 0
  
  # Prepare results
  results <- list(
    meta_data = meta_data,
    pooled_hr = pooled_hr,
    ci_lower = pooled_ci_lower,
    ci_upper = pooled_ci_upper,
    p_value = p_value,
    z_score = z_score,
    heterogeneity = list(
      q = q_stat,
      df = df,
      p = p_heterogeneity,
      i_squared = i_squared
    ),
    log_transform = log_transform,
    conf_level = conf_level
  )
  
  return(results)
}

#' Perform Random-Effects Meta-Analysis
#'
#' This function performs a random-effects meta-analysis using the DerSimonian and Laird
#' method on hazard ratios from multiple studies.
#'
#' @param meta_data Data frame with columns: study, hr, lower, upper. Can also include
#'        additional columns that will be preserved in the output.
#' @param log_transform Whether to log-transform hazard ratios for analysis (default: TRUE)
#' @param conf_level Confidence level for intervals (default: 0.95)
#' @param add_weights Whether to add weight columns to output data (default: TRUE)
#'
#' @return List with meta-analysis results, including pooled estimate and data with weights
#' @export
#'
#' @examples
#' \dontrun{
#' meta_data <- data.frame(
#'   study = c("Study 1", "Study 2", "Study 3"),
#'   hr = c(1.5, 2.1, 1.8),
#'   lower = c(1.1, 1.6, 1.2),
#'   upper = c(2.0, 2.8, 2.7)
#' )
#' 
#' meta_results <- perform_random_effects_meta(meta_data)
#' print(meta_results$pooled_hr)  # Pooled hazard ratio
#' }
perform_random_effects_meta <- function(meta_data, log_transform = TRUE, 
                                       conf_level = 0.95, add_weights = TRUE) {
  
  # Validate input data
  required_cols <- c("study", "hr", "lower", "upper")
  missing_cols <- setdiff(required_cols, names(meta_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in meta_data: ", paste(missing_cols, collapse = ", "))
  }
  
  # Ensure numeric types for calculations
  meta_data$hr <- as.numeric(meta_data$hr)
  meta_data$lower <- as.numeric(meta_data$lower)
  meta_data$upper <- as.numeric(meta_data$upper)
  
  # Calculate log hazard ratios and standard errors
  if (log_transform) {
    meta_data$loghr <- log(meta_data$hr)
    meta_data$se <- (log(meta_data$upper) - log(meta_data$lower)) / (2 * stats::qnorm((1 + conf_level) / 2))
  } else {
    meta_data$loghr <- meta_data$hr
    meta_data$se <- (meta_data$upper - meta_data$lower) / (2 * stats::qnorm((1 + conf_level) / 2))
  }
  
  # Calculate initial weights (fixed-effects)
  meta_data$weight_fixed <- 1 / meta_data$se^2
  
  # Calculate fixed-effects pooled estimate (for heterogeneity calculation)
  weighted_mean_fixed <- sum(meta_data$loghr * meta_data$weight_fixed) / sum(meta_data$weight_fixed)
  
  # Calculate heterogeneity statistics
  q_stat <- sum(meta_data$weight_fixed * (meta_data$loghr - weighted_mean_fixed)^2)
  df <- nrow(meta_data) - 1
  p_heterogeneity <- 1 - stats::pchisq(q_stat, df)
  
  # Calculate between-study variance (tau-squared) using DerSimonian-Laird method
  tau_squared <- max(0, (q_stat - df) / (sum(meta_data$weight_fixed) - sum(meta_data$weight_fixed^2) / sum(meta_data$weight_fixed)))
  
  # Calculate I-squared
  i_squared <- max(0, 100 * (q_stat - df) / q_stat)
  if (q_stat <= df) i_squared <- 0
  
  # Calculate random-effects weights
  meta_data$weight <- 1 / (meta_data$se^2 + tau_squared)
  
  if (add_weights) {
    meta_data$weight_percent <- 100 * meta_data$weight / sum(meta_data$weight)
  }
  
  # Calculate random-effects pooled estimate
  weighted_mean <- sum(meta_data$loghr * meta_data$weight) / sum(meta_data$weight)
  weighted_se <- sqrt(1 / sum(meta_data$weight))
  
  # Calculate confidence interval for pooled estimate
  alpha <- 1 - conf_level
  z_value <- stats::qnorm(1 - alpha/2)
  ci_lower <- weighted_mean - z_value * weighted_se
  ci_upper <- weighted_mean + z_value * weighted_se
  
  # Calculate z-score and p-value
  z_score <- weighted_mean / weighted_se
  p_value <- 2 * (1 - stats::pnorm(abs(z_score)))
  
  # Transform back if necessary
  if (log_transform) {
    pooled_hr <- exp(weighted_mean)
    pooled_ci_lower <- exp(ci_lower)
    pooled_ci_upper <- exp(ci_upper)
  } else {
    pooled_hr <- weighted_mean
    pooled_ci_lower <- ci_lower
    pooled_ci_upper <- ci_upper
  }
  
  # Prepare results
  results <- list(
    meta_data = meta_data,
    pooled_hr = pooled_hr,
    ci_lower = pooled_ci_lower,
    ci_upper = pooled_ci_upper,
    p_value = p_value,
    z_score = z_score,
    heterogeneity = list(
      q = q_stat,
      df = df,
      p = p_heterogeneity,
      i_squared = i_squared,
      tau_squared = tau_squared
    ),
    log_transform = log_transform,
    conf_level = conf_level
  )
  
  return(results)
}

#' Perform Spline Meta-Regression
#'
#' This function performs a spline meta-regression to analyze non-linear relationships
#' between a moderator variable and effect sizes (e.g., hazard ratios).
#'
#' @param meta_data Data frame with columns: hr, lower, upper, and the moderator variable
#' @param moderator_name Name of the moderator variable in meta_data
#' @param df Degrees of freedom for spline (default: 2)
#' @param log_transform Whether to log-transform hazard ratios for analysis (default: TRUE)
#' @param conf_level Confidence level for intervals (default: 0.95)
#' @param method Estimation method: "REML" (default), "ML", "DL", or "FE"
#' @param n_pred Number of prediction points for smooth curve (default: 100)
#'
#' @return List with meta-regression results, including model, predictions, and specific points
#' @export
#'
#' @examples
#' \dontrun{
#' # Meta-regression with GLS values as moderator
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
#' }
perform_spline_meta_regression <- function(meta_data, moderator_name, df = 2, 
                                          log_transform = TRUE, conf_level = 0.95,
                                          method = "REML", n_pred = 100) {
  
  # Check if required packages are available
  if (!requireNamespace("splines", quietly = TRUE)) {
    stop("Package 'splines' is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  # Validate input data
  required_cols <- c("hr", "lower", "upper", moderator_name)
  missing_cols <- setdiff(required_cols, names(meta_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in meta_data: ", paste(missing_cols, collapse = ", "))
  }
  
  # Ensure numeric types for calculations
  meta_data$hr <- as.numeric(meta_data$hr)
  meta_data$lower <- as.numeric(meta_data$lower)
  meta_data$upper <- as.numeric(meta_data$upper)
  meta_data[[moderator_name]] <- as.numeric(meta_data[[moderator_name]])
  
  # Calculate log hazard ratios and standard errors
  if (log_transform) {
    meta_data$loghr <- log(meta_data$hr)
    meta_data$se <- (log(meta_data$upper) - log(meta_data$lower)) / (2 * stats::qnorm((1 + conf_level) / 2))
  } else {
    meta_data$loghr <- meta_data$hr
    meta_data$se <- (meta_data$upper - meta_data$lower) / (2 * stats::qnorm((1 + conf_level) / 2))
  }
  
  # Extract moderator values
  moderator_values <- meta_data[[moderator_name]]
  
  # Create spline basis
  spline_basis <- splines::ns(moderator_values, df = df)
  
  # Fit weighted regression model
  lm_model <- stats::lm(loghr ~ spline_basis, data = meta_data, weights = 1/se^2)
  
  # Create prediction grid
  pred_grid <- data.frame(
    moderator = seq(
      min(moderator_values, na.rm = TRUE),
      max(moderator_values, na.rm = TRUE),
      length.out = n_pred
    )
  )
  names(pred_grid)[1] <- moderator_name
  
  # Create design matrix for predictions
  pred_basis <- predict(splines::ns(moderator_values, df = df), 
                       newx = pred_grid[[moderator_name]])
  
  # Prepare matrix for prediction
  pred_matrix <- matrix(pred_basis, ncol = df)
  colnames(pred_matrix) <- colnames(spline_basis)
  
  # Predict using model
  pred_fit <- stats::predict(lm_model, 
                            newdata = data.frame(spline_basis = pred_matrix), 
                            se.fit = TRUE)
  
  # Extract predictions and standard errors
  pred_grid$loghr <- pred_fit$fit
  pred_grid$se_loghr <- pred_fit$se.fit
  
  # Calculate confidence intervals
  alpha <- 1 - conf_level
  z_value <- stats::qnorm(1 - alpha/2)
  pred_grid$loghr_lower <- pred_grid$loghr - z_value * pred_grid$se_loghr
  pred_grid$loghr_upper <- pred_grid$loghr + z_value * pred_grid$se_loghr
  
  # Transform back if requested
  if (log_transform) {
    pred_grid$hr <- exp(pred_grid$loghr)
    pred_grid$hr_lower <- exp(pred_grid$loghr_lower)
    pred_grid$hr_upper <- exp(pred_grid$loghr_upper)
  } else {
    pred_grid$hr <- pred_grid$loghr
    pred_grid$hr_lower <- pred_grid$loghr_lower
    pred_grid$hr_upper <- pred_grid$loghr_upper
  }
  
  # Calculate specific values at particular moderator points
  specific_points <- seq(
    min(moderator_values, na.rm = TRUE),
    max(moderator_values, na.rm = TRUE),
    length.out = 5
  )
  
  specific_results <- data.frame(
    moderator = specific_points,
    hr = numeric(length(specific_points)),
    hr_lower = numeric(length(specific_points)),
    hr_upper = numeric(length(specific_points))
  )
  names(specific_results)[1] <- moderator_name
  
  for (i in 1:length(specific_points)) {
    point <- specific_points[i]
    idx <- which.min(abs(pred_grid[[moderator_name]] - point))
    specific_results$hr[i] <- pred_grid$hr[idx]
    specific_results$hr_lower[i] <- pred_grid$hr_lower[idx]
    specific_results$hr_upper[i] <- pred_grid$hr_upper[idx]
  }
  
  # Get model summary for output
  model_summary <- summary(lm_model)
  
  # Prepare results
  results <- list(
    meta_data = meta_data,
    model = lm_model,
    model_summary = model_summary,
    predictions = pred_grid,
    specific_results = specific_results,
    moderator_name = moderator_name,
    df = df,
    log_transform = log_transform,
    conf_level = conf_level
  )
  
  return(results)
}

#' Test for Publication Bias (Egger's Test)
#'
#' This function performs Egger's test for publication bias in meta-analysis.
#'
#' @param meta_data Data frame with columns: hr, lower, upper, and optionally sample_size
#' @param log_transform Whether to log-transform hazard ratios for analysis (default: TRUE)
#' @param precision_measure Method to calculate precision: "se" (standard error, default) or "sample" (sample size)
#'
#' @return List with Egger's test results
#' @export
#'
#' @examples
#' \dontrun{
#' meta_data <- data.frame(
#'   study = c("Study 1", "Study 2", "Study 3"),
#'   hr = c(1.5, 2.1, 1.8),
#'   lower = c(1.1, 1.6, 1.2),
#'   upper = c(2.0, 2.8, 2.7)
#' )
#' 
#' eggers_result <- test_publication_bias(meta_data)
#' print(eggers_result$p_value)
#' }
test_publication_bias <- function(meta_data, log_transform = TRUE, 
                                 precision_measure = "se") {
  
  # Validate input data
  required_cols <- c("hr", "lower", "upper")
  missing_cols <- setdiff(required_cols, names(meta_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in meta_data: ", paste(missing_cols, collapse = ", "))
  }
  
  # Calculate log hazard ratios and standard errors
  if (log_transform) {
    meta_data$loghr <- log(meta_data$hr)
    meta_data$se <- (log(meta_data$upper) - log(meta_data$lower)) / (2 * 1.96)
  } else {
    meta_data$loghr <- meta_data$hr
    meta_data$se <- (meta_data$upper - meta_data$lower) / (2 * 1.96)
  }
  
  # Calculate precision (1/se or sample size)
  if (precision_measure == "se") {
    meta_data$precision <- 1 / meta_data$se
  } else if (precision_measure == "sample") {
    if (!"sample_size" %in% names(meta_data)) {
      stop("'sample_size' column is required when precision_measure = 'sample'")
    }
    meta_data$precision <- sqrt(meta_data$sample_size)
  } else {
    stop("precision_measure must be 'se' or 'sample'")
  }
  
  # Perform Egger's test (regress effect estimate against precision)
  eggers_model <- stats::lm(meta_data$loghr ~ meta_data$precision)
  eggers_summary <- summary(eggers_model)
  
  # Extract test statistics
  intercept <- eggers_summary$coefficients[1, "Estimate"]
  intercept_se <- eggers_summary$coefficients[1, "Std. Error"]
  t_value <- eggers_summary$coefficients[1, "t value"]
  p_value <- eggers_summary$coefficients[1, "Pr(>|t|)"]
  
  # Prepare results
  results <- list(
    intercept = intercept,
    intercept_se = intercept_se,
    t_value = t_value,
    p_value = p_value,
    model = eggers_model,
    log_transform = log_transform,
    precision_measure = precision_measure
  )
  
  return(results)
}