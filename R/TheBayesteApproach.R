#' Bayesian estimation of distribution components from antibiotic resistance data
#'
#' This function performs BNN... mixture density estimation of binned antibiotic resistance data with Bayesian modeling
#'
#' This function can be used to estimate a BNN... mixture density from binned antibiotic resistance. It therefore uses a Bayesian set up
#' where it models the resistant antibiotics (ZD = 6) with a binomial component and all other non-resistant antibiotics with a mixture
#' of normal components. The data is assumed to be binned to integer values which minimum 6 and maximum 50. Additionally to the data, the arguments
#' prior and k may be specified. Thereby, the argument prior must be a list containing hyperparameters for the prior distribution and k must be
#' a value defining the number of modeled normal components. The elements of the prior list need to be named as the following:
#' \itemize{
#'  \item{"a"}{    shape parameter \eqn{\alpha} of the beta prior for \eqn{\theta} of the binomial component}
#'  \item{"b"}{    shape parameter \eqn{\beta} of the beta prior for \eqn{\theta} of the binomial component}
#'  \item{"b0"}{   hyperparameter \eqn{\mu} of the normal prior for \eqn{\mu} of the normal components}
#'  \item{"B0"}{   hyperparameters \eqn{\tau} (precision) of the normal prior for \eqn{\mu} of the normal components}
#'  \item{"c0"}{   hyper shape parameter \eqn{\gamma} of the gamma prior for \eqn{\tau} of the normal components}
#'  \item{"C0"}{   hyper lambda parameter \eqn{\lambda} of the gamma prior for \eqn{\tau} of the normal components}
#'  \item{"e0"}{   hyper lambda parameter \eqn{e_0} of the dirichlet prior for \eqn{\pi}}
#' }
#' The length of each hyperparameter can only be one, thus if one wants to use hyperparameters for the priors of \eqn{\mu} and \eqn{\tau},
#' they are used for all normal components. If the hyperparameters are not named in the that way the function still runs but ignores the
#' wrongly labeled elements of the list. In the first step the functions estimates the MAPE for the binomial component. If any
#' hyperparameter for the beta prior is defined as noted above, it will be used. Otherwise a and b are set to 1. Jags is used for Gibbs sampling
#' from the posterior distribution. Finally, even though hyperparameters and the prior itself are defined for \eqn{\tau}, the MAPE for \eqn{\sigma^2}
#'  will be estimated and returned the MAPE, besides the MAPE for the \eqn{\mu}.
#'
#' If the argument k is not defined by the user, the function assumes 20 components, applies Gibbs sampling and then post-processes
#' the outcome to determine the number of resulting components k that are filled. If the argument k is given, Gibbs sampling is
#' applied for the given number of components k.
#'
#' Finally, the function returns an object of the class 'TheBayesteApproach' which includes the estimated number of components, the MAPE
#' parameters for the \eqn{\mu}'s and \eqn{\sigma^2}'s of all normal components, the MAPE of the binomial component and the input data.
#' @param data vector of flat data or table of counts with groups as names
#' @param prior an optional list of hyperparameters for prior distributions (See details how to define this list). Default is NULL.
#' @param k optional value of how many normal components should be modeled. Default is NULL.
#' @keywords Bayesian modeling, Antibiotics, density estimation
#' @return An object of the S3 class 'TheBayesteApproach which includes a list including these elements:
#' \itemize{
#'  \item{k: }{    number of components}
#'  \item{normal.params: }{vector with the MAPE normal component parameters}
#'  \item{binom.res: }{ estimated MAPE binomial parameter}
#'  \item{data: }{ the flat data as vector}
#' }
#' @export
#' @seealso \code{\link{plot.TheBayesteApproach}}
#' @examples
#' data("Antibiotics")
#' bayes.density <- TheBayesteApproach(Antibiotics)
#' print(bayes.density)
#' plot(bayes.density)
TheBayesteApproach <- function(data, prior = NULL, k = NULL) {

  #In case of count data, converting it to flat data
  if (!is.null(names(data))) {
    data <- rep(6:(length(data)+5), data)
  }

  if (! (is.numeric(data) & is.vector(data))) {
    stop("Argument data must be numeric vector")
  }

  if(!all(data == floor(data))) {
    stop("data must only contain integer values")
  }

  if (!is.null(prior)){
    if (! is.list(prior)){
      stop("Argument prior must be a list")
    }
  }

  if (any(sapply(prior, length) != 1) | any(!sapply(prior, is.numeric)) ){
    stop("All arguments in prior list must be of length 1 with numeric values")
  }

  if (!is.null(k)){
    if (! (is.numeric(k) & length(k) == 1) ){
      stop("k must be a single number")
    }

    if(k < 0 | k > 10){
      stop("k must be an integer between 1 and 10")
    }
  }


  #Binomial Component density estimation
  binom.res <- EstimateBinomial(data, prior)


  #if we only have resistent data we return immediately & throw a warning
  if(all(data == 6)){
    warning("No normal component was fitted because only resistant data is given")

    ret.result <- list(k = NA,normal.params = NA, binom.res = binom.res, data = data)
    class(ret.result) <- "TheBayesteApproach"
    return(ret.result)
  }

  #too little data
  if(length(which(data>6)) < 16){
    stop("Too little non-resistant data!")
  }


  #Extract the non-resistent bacterias that are relevant for normal components
  data.resistent <- data[data != 6]

  #Gibbs Sampling
  output.gibbs.sampler <- CallJags(data.resistent, prior = prior, k = k)

  #Postprocessing: Get number of components
  if (is.null(k)){
    res.comp <- GetNumberComponents(output.gibbs.sampler)
    k <- res.comp$k
    output.gibbs.sampler <- res.comp$updated.gibbs
  }

  #Postprocessing: Estimate parameters using mode estimator
  result.estimation <- EstimateNormal(output.gibbs.sampler, k, y = data.resistent, prior = prior)

  ret.result <- list(k = k,normal.params = result.estimation$params, binom.res = binom.res, data = data)
  class(ret.result) <- "TheBayesteApproach"
  return(ret.result)
}



