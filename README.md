# sopravvivenza: Reconstruct Individual Patient Data from Digitized Kaplan-Meier Curves

Please note this is under construction.

# Overview

The `sopravvivenza` package provides a comprehensive toolkit for reconstructing individual patient data (IPD) from digitized Kaplan-Meier survival curves using the method described by Guyot et al. This approach is particularly valuable for meta-analyses where only published survival curves are available rather than the original IPD.

# Installation

```r
# Install the development version from GitHub
install.packages("devtools")
devtools::install_github("bisacciamd/sopravvivenza")
```

# Workflow

The typical workflow consists of the following steps:

1. Digitize survival curves using software like WebPlotDigitizer
2. Save digitized points as CSV files (time, survival)
3. Load these files into R using functions in this package
4. Reconstruct IPD and validate the reconstruction
5. Calculate hazard ratios and perform meta-analysis

# Key Functions

- `read_curve_data()`: Read and preprocess digitized curve data
- `read_risk_data()`: Read number-at-risk data
- `reconstruct_ipd()`: Reconstruct individual patient data from digitized curves
- `calculate_hazard_ratio()`: Calculate hazard ratio between two groups
- `create_km_plot()`: Create Kaplan-Meier plots comparing survival curves
- `perform_fixed_effects_meta()`: Perform fixed-effects meta-analysis
- `perform_spline_meta_regression()`: Conduct spline meta-regression analysis

# Example

```r
library(sopravvivenza)

# Read digitized curve data
treatment_data <- read_curve_data("treatment_curve.csv")
control_data <- read_curve_data("control_curve.csv")

# Reconstruct IPD
treatment_ipd <- reconstruct_ipd(
  treatment_data$time, 
  treatment_data$surv,
  n_risk_initial = 100,
  study_name = "Example Study",
  arm_name = "Treatment"
)

control_ipd <- reconstruct_ipd(
  control_data$time, 
  control_data$surv,
  n_risk_initial = 100,
  study_name = "Example Study",
  arm_name = "Control"
)

# Calculate hazard ratio
hr_result <- calculate_hazard_ratio(
  control_ipd$ipd,  Reference group
  treatment_ipd$ipd,
  "Control",
  "Treatment"
)

Print results
print(paste0("HR: ", round(hr_result$hr, 2), 
            " (95% CI: ", round(hr_result$lower, 2), 
            "-", round(hr_result$upper, 2), 
            "), p = ", format.pval(hr_result$p_value, digits = 3)))

# Create Kaplan-Meier plot
create_km_plot(
  hr_result$combined_data,
  group_var = "group",
  title = "Treatment vs. Control",
  group_labels = c("Control", "Treatment")
)
```

# Citation

If you use this package, please cite both the package and the original Guyot et al. paper that describes the method:

```
Bisaccia G (2025). sopravvivenza: Reconstruct Individual Patient Data from Digitized Kaplan-Meier Curves. R package version 0.1.0.

Guyot P, Ades AE, Ouwens MJ, Welton NJ (2012). Enhanced secondary analysis of survival data: reconstructing the data from published Kaplan-Meier survival curves. BMC Medical Research Methodology, 12:9.
```

# License

This package is released under the MIT License.
