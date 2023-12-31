% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/individual_covcor.R
\name{individual_covcor}
\alias{individual_covcor}
\title{Calculate group-specific covariance and correlation matrices}
\usage{
individual_covcor(X, y)
}
\arguments{
\item{X}{the data matrix}

\item{y}{the label}
}
\value{
COR.list is a list, each element contains a group-specific correlation matrix

COV.list is a list, each element contains a group-specific covariance matrix

COV.scale.list is a list, each element contains a standardized group-specific covariance matrix
}
\description{
Calculate group-specific covariance and correlation matrices
}
