
# ===========================================================================
# R/hazard_ratio.R
# ===========================================================================

#' Calculate Hazard Ratio Between Two Groups
#'
#' This function calculates the hazard ratio and confidence intervals between
#' two groups using Cox proportional hazards model.
#'
#' @param data1 IPD dataframe for group 1 (reference group)
#' @param data2 IPD dataframe for group 2
#' @param group1_name Name for group 1 (reference)
#' @param group2_name Name for group 2
#' @param covariates Optional vector of covariate names to adjust for
#' @param conf_level Confidence level for intervals (default: 0.95)
#'
#' @return List with hazard ratio statistics and model results
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate hazard ratio between two groups
#' hr_result <- calculate_hazard_ratio(
#'   group1_ipd, group2_ipd, 
#'   "Control", "Treatment"
#' )
#' 
#' # With covariates
#' hr_result <- calculate_hazard_ratio(
#'   group1_ipd, group2_ipd, 
#'   "Control", "Treatment",
#'   covariates = c("age", "sex")
#' )
#' }
calculate_hazard_ratio <- function(data1, data2, group1_name = "group1", 
                                  group2_name = "group2", covariates = NULL,
                                  conf_level = 0.95) {
  
  # Add group indicators
  data1$group <- group1_name
  data2$group <- group2_name
  
  # Combine datasets
  combined_data <- rbind(data1, data2)
  
  # Set group as a factor with group1 as reference
  combined_data$group <- factor(combined_data$group, levels = c(group1_name, group2_name))
  
  # Create formula for Cox model
  if (is.null(covariates)) {
    formula <- stats::as.formula("survival::Surv(time, status) ~ group")
  } else {
    formula <- stats::as.formula(paste("survival::Surv(time, status) ~ group +", 
                               paste(covariates, collapse = " + ")))
  }
  
  # Fit Cox model
  cox_model <- survival::coxph(formula, data = combined_data)
  summary_cox <- summary(cox_model)
  
  # Extract hazard ratio and CI
  alpha <- 1 - conf_level
  z_value <- stats::qnorm(1 - alpha/2)
  
  hr <- exp(summary_cox$coefficients[1, "coef"])
  hr_se <- summary_cox$coefficients[1, "se(coef)"]
  hr_lower <- exp(summary_cox$coefficients[1, "coef"] - z_value * hr_se)
  hr_upper <- exp(summary_cox$coefficients[1, "coef"] + z_value * hr_se)
  p_value <- summary_cox$coefficients[1, "Pr(>|z|)"]
  
  return(list(
    hr = hr,
    lower = hr_lower,
    upper = hr_upper,
    se = hr_se,
    p_value = p_value,
    cox_model = cox_model,
    combined_data = combined_data,
    reference_group = group1_name,
    comparison_group = group2_name
  ))
}