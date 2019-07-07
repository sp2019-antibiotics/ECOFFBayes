#' Binomial Component
#'
#' This function calculates the MAPE for the parameter of the binomial component
#' @param data vector of flat data
#' @param prior an optional list of hyperparameters for prior distributions (See details in \code{\link{TheBayesteApproach}} how to define this list). Default is NULL.
#' @keywords MAPE, Beta prior, Binomial Likelihood
#' @return MAPE for the parameter of the binomial component
EstimateBinomial <- function(data, prior = NULL){
  bin.y <- ifelse(data == 6, 1, 0)

  #Use hyperparameter if given
  ifelse(!is.null(prior$a), a.prior <- prior$a, a.prior <- 1)
  ifelse(!is.null(prior$b), b.prior <- prior$b, b.prior <- 1)

  #Calculating the posterior distribution
  a.posterior <- sum(bin.y) + a.prior
  b.posterior <- length(bin.y) - sum(bin.y) + b.prior

  #Beta estimator (mode) for density estimation
  mode.posterior <- (a.posterior - 1) / (a.posterior + b.posterior - 2)

  #Return parameter estimator for density estimation
  return(mode.posterior)
}
