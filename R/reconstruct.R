
# ===========================================================================
# R/reconstruct.R
# ===========================================================================

#' Reconstruct Individual Patient Data from Digitized Survival Curve
#'
#' This function implements the algorithm described by Guyot et al. to reconstruct
#' individual patient data (IPD) from digitized Kaplan-Meier curves.
#'
#' @param time_points Vector of time points from digitized curve
#' @param surv_probs Vector of survival probabilities from digitized curve
#' @param n_risk_initial Initial number of patients at risk
#' @param total_events Total number of events reported (if known)
#' @param study_name Name of the study (for labeling)
#' @param arm_name Name of the study arm/group (for labeling)
#' @param risk_times Optional vector of time points for number at risk data
#' @param risk_numbers Optional vector of numbers at risk at risk_times
#' @param output_dir Directory for saving output files (NULL for no saving)
#' @param make_plot Whether to create validation plot
#' @param adjust_for_total_events Whether to adjust IPD to match total events
#'
#' @return List containing reconstructed IPD and other results
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' result <- reconstruct_ipd(
#'   time_points = c(0, 6, 12, 18, 24),
#'   surv_probs = c(1, 0.9, 0.8, 0.7, 0.6),
#'   n_risk_initial = 100
#' )
#' 
#' # With number at risk data
#' result <- reconstruct_ipd(
#'   time_points = c(0, 6, 12, 18, 24),
#'   surv_probs = c(1, 0.9, 0.8, 0.7, 0.6),
#'   n_risk_initial = 100,
#'   risk_times = c(0, 12, 24),
#'   risk_numbers = c(100, 75, 50)
#' )
#' }
reconstruct_ipd <- function(time_points, surv_probs, n_risk_initial = 100, 
                           total_events = NULL, study_name = "study", 
                           arm_name = "arm", risk_times = NULL, 
                           risk_numbers = NULL, output_dir = NULL,
                           make_plot = TRUE, adjust_for_total_events = TRUE) {
  
  # Create directory for study outputs
  if (!is.null(output_dir)) {
    study_dir <- file.path(output_dir, study_name)
    dir.create(study_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
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
  if (adjust_for_total_events && !is.null(total_events) && length(events) > 0 && sum(events) != total_events) {
    current_events <- sum(events)
    
    if (current_events < total_events) {
      # Need more events - add random events
      more_needed <- total_events - current_events
      # Generate random times within range
      random_times <- stats::runif(more_needed, min(time_points), max(time_points))
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
    fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = ipd)
    
    if (!is.null(output_dir)) {
      # Save plot as PDF
      grDevices::pdf(file.path(study_dir, paste0(study_name, "_", arm_name, "_validation.pdf")), 
          width = 8, height = 6)
      
      graphics::plot(fit, main = paste("Validation -", study_name, "-", arm_name),
           xlab = "Time", ylab = "Survival Probability", 
           conf.int = FALSE, mark.time = FALSE)
      graphics::points(time_points, surv_probs, col = "red", pch = 19)
      graphics::legend("topright", c("Reconstructed", "Original Digitized"), 
             col = c("black", "red"), lty = c(1, NA), pch = c(NA, 19))
      
      grDevices::dev.off()
      
      # Also save as PNG
      grDevices::png(file.path(study_dir, paste0(study_name, "_", arm_name, "_validation.png")), 
          width = 800, height = 600, res = 100)
      
      graphics::plot(fit, main = paste("Validation -", study_name, "-", arm_name),
           xlab = "Time", ylab = "Survival Probability", 
           conf.int = FALSE, mark.time = FALSE)
      graphics::points(time_points, surv_probs, col = "red", pch = 19)
      graphics::legend("topright", c("Reconstructed", "Original Digitized"), 
             col = c("black", "red"), lty = c(1, NA), pch = c(NA, 19))
      
      grDevices::dev.off()
      
      # Save IPD to CSV file
      readr::write_csv(ipd, file.path(study_dir, paste0(study_name, "_", arm_name, "_ipd.csv")))
    }
  } else {
    fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = ipd)
  }
  
  return(list(
    ipd = ipd,
    fit = fit,
    original_data = data.frame(time = time_points, surv = surv_probs),
    n_risk_initial = n_risk_initial,
    total_events = sum(ipd$status),
    study_name = study_name,
    arm_name = arm_name
  ))
}
