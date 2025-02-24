% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{which.lag}
\alias{which.lag}
\title{Find optimal lag of VAR according to certain criterion}
\usage{
which.lag(data, lag.max = 12, criterion = "AIC")
}
\arguments{
\item{lag.max}{a integer indicating max lag to test}

\item{criterion}{"AIC" or "BIC"}

\item{Y}{input time series data matrix}
}
\value{
optimal lag based on AIC or BIC
}
\description{
Find optimal lag of VAR according to certain criterion
}
