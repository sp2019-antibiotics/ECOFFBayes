#' Normal Components
#'
#' This internal function calculates the log likelihood for each step in the gibbs sampling
#' @param gibbs matrix with parameters sampled with gibbs algorithm
#' @param k Number of components
#' @param y vector of flat data
#' @param prior an optional list of hyperparameters for prior distributions
#' @keywords logLikelihood
#' @return logLikelihood
LogLikelihood <- function(gibbs, k, y, prior) {

  #relevant parameter columns
  theta <- gibbs[, 1:(3*k)]

  #prior settings
  if(is.null(prior)){
    b0.prior <- 0
    B0.prior <- 0.0001
    c0.prior <- 1
    C0.prior <- 10
  }else{
    b0.prior <- prior$b0
    B0.prior <- prior$B0
    c0.prior <- prior$c0
    C0.prior <- prior$C0
  }


  #density estimator
  EstDensity <- function(j, theta, k){
    z <- vector("numeric",k)

    for(i in 1:k){
      z[i] <- pnorm(j, theta[i], sd = sqrt(theta[(k*2+i)]))
      z[i] <- z[i]*theta[k+i]
    }
    return(sum(z))
  }
  #boundaries of the bins
  bj <- seq(7.5, 50.5, by = 1)
  aj <- seq(6.5, 49.5, by = 1)

  res <- vector("numeric",nrow(gibbs))

  #helper function, computes the log-likelihood
  logL <- function(theta) {

    y0 <- c(7:50)
    #computes the bin frequencies
    nj <- merge(data.frame(y = y0),as.data.frame(table(y)),by="y",all.x = TRUE)
    nj[is.na(nj[,2]),2] <- 0

    #computes the likelihood part
    bin <- vector("numeric",length = length(nj))
    for(i in 1:nrow(nj)){
      l <- nj[i, 2] * log(EstDensity(bj[i],theta = theta, k = k) -
                                 EstDensity(aj[i],theta = theta, k = k))
      bin[i] <- ifelse(is.nan(l),0,l)
    }

    #computes the priors
    prior <- vector("numeric",length = length(k))
    for(i in 1:k){
      prior[i] <- log(dnorm(theta[i], b0.prior, sqrt(1/B0.prior))) +
        log(dgamma(1/theta[k * 2 + i], shape = c0.prior, rate = C0.prior)) +
        log(dbeta(theta[k + i], 1, 1))
    }
    #res[s] <- sum(bin+sum(prior))
    return(sum(bin + sum(prior)))
  }


  res <- apply(theta, 1, logL)
  return(res)
}

#' Normal Components
#'
#' Calculates the log likelihood for each step in the gibbs sampling.
#' @param gibbs Matrix with parameters sampled with gibbs algorithm
#' @param k Number components
#' @param prior Prior list
#' @param y Data
#' @import stats
EstimateNormal <- function(gibbs, k, y, prior) {

  # Estimates the parameter of the normal distributions by maximizing the log likelihood.
  logL <- LogLikelihood(gibbs, k, y, prior)
  estimated.params <- gibbs[which.max(logL),1:(3*k)]
  return(list(params = estimated.params))
}

