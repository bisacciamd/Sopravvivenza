% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{calculate_survival_from_ipd}
\alias{calculate_survival_from_ipd}
\title{Calculate Survival Probabilities from IPD}
\usage{
calculate_survival_from_ipd(
  ipd,
  time_points,
  time_var = "time",
  status_var = "status"
)
}
\arguments{
\item{ipd}{Individual patient data with time and status columns}

\item{time_points}{Vector of time points at which to calculate survival probabilities}

\item{time_var}{Name of time variable in ipd (default: "time")}

\item{status_var}{Name of status variable in ipd (default: "status")}
}
\value{
Data frame with time points and corresponding survival probabilities
}
\description{
This function calculates survival probabilities at specific time points
from reconstructed individual patient data.
}
