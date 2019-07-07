#' Number of Components
#'
#' This internal function postprocesses the gibbs sampling output and returns the resulting number of components and MAPE estimates
#' @param gibbs Matrix with parameters sampled (output of function CallJags)
#' @keywords Gibbs sampler
#' @return Number of components and MAPE estimates of the parameters of the normal components

GetNumberComponents <- function(gibbs) {

  x <- as.data.frame(gibbs[,grepl("z",colnames(gibbs))])
  dist.k <- apply(x, 1,function(x) length(unique(x)))

  # number of groups
  updated.k <- as.numeric(names(which.max(table(dist.k))))
  print(paste('Number of components:', updated.k))

  # all rows which have new.k groups
  updated.gibbs <- gibbs[dist.k == updated.k, ,drop=FALSE]
  x.updated <- x[dist.k == updated.k, ,drop=FALSE]

  #extract the correct columns
  unique.comb <- matrix(t(apply(x.updated, 1, function(x) sort(unique(x), decreasing = TRUE))),ncol = 1*updated.k)
  labels <- t(apply(t(unique.comb), 2, function(x){c(paste0("mu[",x,"]"),
                                                     paste0("pi[",x,"]"),
                                                     paste0("sigma2[",x,"]"))}))

  #column indices for the used mu, pi, sigma columns
  ind <- t(apply(labels, 1, function(x){match(x, colnames(updated.gibbs))}))

  #create a new matrix with samples from the posterior
  mat <- matrix(0,ncol = 3*updated.k, nrow = nrow(updated.gibbs))
  for(i in 1: nrow(ind)){
      mat[i, 1:(updated.k*3)] <- updated.gibbs[i, ind[i,]]
  }

  #Setting
  label.comp <- c(paste0("mu[",1:updated.k,"]"), paste0("pi[",1:updated.k,"]"),paste0("sigma2[",1:updated.k,"]"))
  colnames(mat) <- c(label.comp)

  return(list(k = updated.k, updated.gibbs = mat))

}



