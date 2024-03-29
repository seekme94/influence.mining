% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph_util.R
\name{normalize_data}
\alias{normalize_data}
\title{Normalize a set of numeric variables in a dataset between 0 and 1. Non-numeric data will be skipped}
\usage{
normalize_data(data, columns)
}
\arguments{
\item{data}{is data frame to be normalized}

\item{columns}{is list of columns to be normalized}
}
\value{
normalized data frame
}
\description{
Normalize a set of numeric variables in a dataset between 0 and 1. Non-numeric data will be skipped
}
