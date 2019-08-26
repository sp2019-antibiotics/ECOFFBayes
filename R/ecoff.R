#' ECOFF Finder
#'
#' This function calculates the ECOFF for the wild type of an 'TheBayesteApproach' object
#'
#' This function takes an object of the class TheBayesteApproach and calculates the ECOFF. The ECOFF is the
#' epidemiological cut-off value that separates bacteria without (wild type) and with acquired resistance
#' mechanisms (non-wild type) to the antibiotics in question. In order to calculate the ECOFF, first the wild type
#' component is selected by finding the component where the sum of its \eqn{\pi_k} and  the \eqn{\pi}'s of all components
#' with a higher \eqn{\mu_k} exceed a certain threshold (by default: 0.3). Once this component is found a certain
#' quantile is taken as the ECOFF value (by default: 0.01)
#'
#'
#' @param obj object of class TheBayesteApproach
#' @param pi.level The level of summed up component probabilities from the distributions coming from the "right" side that needs to be exceeded (default is 0.3)
#' @param quantile The quantile from the component where the pi.level is exceeded (default is 0.01)
#' @keywords ECOFF, wild type, antibiotics
#' @export
#' @importFrom methods is
#' @examples
#'
#' \donttest{data("Antibiotics")
#' bayes.density <- TheBayesteApproach(Antibiotics)
#' FindECOFF(bayes.density, pi.level = 0.3, quantile = 0.02)
#' plot(bayes.density, ECOFF = TRUE)}
#'
FindECOFF <- function(obj, pi.level = 0.3, quantile = 0.01){
  #input check
  if(!is(obj, "TheBayesteApproach")){
    stop("Needs object of type TheBayesteApproach")
  }

  if(any(is.na(obj$normal.params))){
    return(NA)
  }

  if (length(pi.level) != 1 || pi.level > 1 || pi.level < 0){
    stop("Argument pi.level must be one number between 0 and 1")
  }

  if (length(quantile) != 1 || quantile > 1 || quantile < 0){
    stop("Argument quantile must be one number between 0 and 1")
  }

  np <- obj$normal.params
  k <- obj$k
  order.mu <- order(np[1:obj$k], decreasing = TRUE)
  pi <- np[grepl("pi",names(np))]
  pi.ordered <- pi[order.mu]

  #Find normal component which exceeds the pi.level threshold
  exceed <- TRUE
  i <- sum <-  0
  while(exceed){
    i <- i + 1
    sum <- sum + pi.ordered[i]
    if(sum > pi.level){
      exceed = FALSE
    }
  }

  #Calculate the quantile of that distribution
  mu.ordered <- np[1:k][order.mu]
  sigma2.ordered <- np[(2*k+1):(3*k)][order.mu]
  ECOFF <- round(qnorm(quantile, mean = mu.ordered[i], sd = sqrt(sigma2.ordered[i])))

  if(ECOFF < 7){
    ECOFF <- 7
  }

  return(ECOFF)
}
