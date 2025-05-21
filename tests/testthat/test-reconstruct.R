test_that("reconstruct_ipd works correctly", {
  Basic reconstruction test
  time_points <- c(0, 6, 12, 18, 24)
  surv_probs <- c(1, 0.9, 0.8, 0.7, 0.6)
  
  result <- reconstruct_ipd(
    time_points = time_points,
    surv_probs = surv_probs,
    n_risk_initial = 100,
    make_plot = FALSE
  )
  
  # Check basic properties of result
  expect_true(is.data.frame(result$ipd))
  expect_equal(ncol(result$ipd), 2)
  expect_true(all(c("time", "status") %in% names(result$ipd)))
  
  # Verify some expected properties
  expect_true(all(result$ipd$status %in% c(0, 1)))  # Status should be 0 or 1
  expect_true(all(result$ipd$time >= 0))            # Times should be non-negative
  expect_true(min(result$ipd$time) >= 0)            # Min time should be 0 or greater
  expect_true(max(result$ipd$time) <= max(time_points))  # Max time shouldn't exceed data
  
  # Check total events is reasonable (should be around 40 for this example)
  total_events <- sum(result$ipd$status)
  expect_true(total_events > 30 && total_events < 50)
  
  # Verify reconstruction matches original curve
  fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = result$ipd)
  surv_at_times <- summary(fit, times = time_points)$surv
  
  # Allow for small numerical differences due to discretization
  for (i in 1:length(surv_at_times)) {
    expect_true(abs(surv_at_times[i] - surv_probs[i]) < 0.05)
  }
})