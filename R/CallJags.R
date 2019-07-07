#' Gibbs Sampler
#'
#' This function can be used to perform gibbs sampling using JAGS and the R-package rjags.
#'
#' This function is part of the TheBayesteApproach function of the ECOFFBayes package. It can be used to performed Gibbs sampling to sample
#' the \eqn{\mu}'s and \eqn{\tau} from the posterior distribution. Thereby, the prior of the \eqn{\mu}'s are assumed to be normal and the
#' prior of \eqn{\tau} is assumed to be a gamma distribution. To account for the binned data, y is drawn from a truncated normal distribution.
#' The data argument must be a flat vector of the binned non-resistant bacteria observations. The prior argument must be a list and specified like
#' described in the function \code{\link{TheBayesteApproach}}. k is an optional parameter indicating the number of normal components.If the argument
#' k is not defined by the user, the function assumes 20 components. In total 10000 draws are done where only each 10th draw is really taken (thinning = 10).
#' A burning period of 1000 iterations is automatically set. Finally, \eqn{\tau} is converted to \eqn{\sigma}^2 and
#' returned is a matrix of the Gibbs sampling draws with the \eqn{\mu}, \eqn{\pi}, \eqn{\sigma}^2 and z.
#' @param y data
#' @param prior an optional list of hyperparameters for prior distributions (See details how to define this list). Default is NULL.
#' @param k optional value of how many normal components should be modeled. Default is NULL.
#' @export
#' @keywords Gibbs sampling, JAGS, posterior distribution
#' @return Returns a matrix which contains the gibbs sampler iterations after thining and deleting the burning period.
#' Columns include the estimated normal components parameters, the component probabilities and the classifications z for all observations.
#'
CallJags <- function(y, prior = NULL, k = NULL){

  #Input checks
  if (! (is.numeric(y) & is.vector(y))) {
    stop("Argument y must be numeric vector")
  }

  if(!all(y == floor(y))) {
    stop("'y' must only contain integer values")
  }

  if(any(y < 7) | any(y > 50)){
    stop("y must only be between 7 and 50")
  }

  if(length(y) < 16){
    stop("Too little non-resistant data!")
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

  # Initialize k if not defined
  if(is.null(k)){
        k <- 20
        e0.prior <- rep(0.0001,20)
  }else{
      ifelse(!is.null(prior$e0), e0.prior <- rep(prior$e0,k), e0.prior <- rep(0.0001,k))
  }

  #prior parameters given the user has not used prior parameters
  ifelse(!is.null(prior$b0), b0.prior <- prior$b0, b0.prior <- 0)
  ifelse(!is.null(prior$B0), B0.prior <- prior$B0, B0.prior <- 0.0001)
  ifelse(!is.null(prior$c0), c0.prior <- prior$c0, c0.prior <- 1)
  ifelse(!is.null(prior$C0), C0.prior <- prior$C0, C0.prior <- 10)


  #the model in the BUGS language
  model.string <- 'model {
    #data model
    for (i in 1:n) {
      y0[i] ~ dnorm(mu[z[i]], tau[z[i]])
      y[i] ~ dinterval(y0[i], lower)
      z[i] ~ dcat(pi)
    }

    #priors
    pi ~ ddirich(e0)
    for (k in 1:K) {
      mu[k] ~ dnorm(b0, B0)
      tau[k] ~ dgamma(c0, C0)
    }
  } '


  model <- textConnection(model.string)

  #data preparation
  y <- y - 6
  c <- seq(0.5,44.5,by=1)
  data <- list(y = y, n = length(y), K = k, lower = c,
               e0 = e0.prior, b0 = b0.prior,
               B0 = B0.prior, c0 = c0.prior, C0 = C0.prior)

  #intialize the model
  inits <- list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 1)
  model <- rjags::jags.model(model, data = data, inits = inits)
  model

  #this is the burn-in
  update(model, n.iter = 1000)
  #we want to get mu, tau, pi and the classification z out of the model
  vars <- c("mu", "tau", "pi","z")
  #the actual sampling, thining is set to 10
  samples <- rjags::coda.samples(model, vars, n.iter = 10000, thin = 10)[[1]]

  #transform tau into sigma^2 and renaming it and retransform the location of mu
  samples[,grepl("tau",colnames(samples))] <- 1/samples[,grepl("tau",colnames(samples))]
  colnames(samples)[grepl("tau",colnames(samples))] <- paste("sigma2[", 1:k, "]", sep = "")

  samples[,grepl("mu",colnames(samples))] <- samples[,grepl("mu",colnames(samples))]+rep(6,k)

  return(samples)
}
