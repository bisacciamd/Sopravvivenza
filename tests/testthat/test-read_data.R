test_that("read_curve_data works correctly", {
  # Test with proper input file
  test_file <- tempfile(fileext = ".csv")
  write.csv(data.frame(time = c(0, 6, 12), 
                       surv = c(1, 0.8, 0.6)), 
            test_file, row.names = FALSE)
  
  result <- read_curve_data(test_file, header = TRUE)
  expect_equal(nrow(result), 3)
  expect_equal(result$time[1], 0)
  expect_equal(result$surv[1], 1)
  
 # Test with percentage survival
  test_file_pct <- tempfile(fileext = ".csv")
  write.csv(data.frame(time = c(0, 6, 12), 
                       surv = c(100, 80, 60)), 
            test_file_pct, row.names = FALSE)
  
  result_pct <- read_curve_data(test_file_pct, header = TRUE)
  expect_equal(result_pct$surv[2], 0.8)
  
 # Test with year-based time
  test_file_year <- tempfile(fileext = ".csv")
  write.csv(data.frame(time = c(0, 0.5, 1), 
                       surv = c(1, 0.8, 0.6)), 
            test_file_year, row.names = FALSE)
  
  result_year <- read_curve_data(test_file_year, header = TRUE)
  expect_equal(result_year$time[2], 6)  0.5 years = 6 months
})