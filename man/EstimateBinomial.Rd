% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Binomial.R
\name{EstimateBinomial}
\alias{EstimateBinomial}
\title{Binomial Component}
\usage{
EstimateBinomial(data, prior = NULL)
}
\arguments{
\item{data}{vector of flat data}

\item{prior}{an optional list of hyperparameters for prior distributions (See details in \code{\link{TheBayesteApproach}} how to define this list). Default is NULL.}
}
\value{
MAPE for the parameter of the binomial component
}
\description{
This function calculates the MAPE for the parameter of the binomial component
}
\keyword{Beta}
\keyword{Binomial}
\keyword{Likelihood}
\keyword{MAPE,}
\keyword{prior,}
