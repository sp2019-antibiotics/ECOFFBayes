#' Plot TheBayesteApproach
#'
#' This function plots a histogram of the data and its estimated mixture density.
#'
#' This function plots a histogram of the data and its estimated mixture density. Additionally to the plotted density,
#' an ECOFF value can be added in the plot by setting the ECOFF argument to true. Further, the arguments pi.level and
#' qunatile can be defined to calculate the ECOFF in a certain way (see \code{\link{FindECOFF}} for Details). Further
#' plot arguments may be provided.
#'
#'
#' @param x An object of class TheBayesteApproach
#' @param xlim limits of x axis
#' @param xlab label of x axis
#' @param main the plot title
#' @param ... additional plot parameters
#' @param ECOFF logical value indicating if the ECOFF should be computed and added to the plot
#' @param quantile quantile of the component with the highest mean that is taken for the ECOFF. Default is 0.01
#' @param pi.level The level of summed up pi from the distributions coming from the "right" side that needs to be exceeded, ECOFF argument (Default: 0.3)
#' @keywords plot.TheBayesteApproach, plotting density, antibiotics
#' @export
#' @import graphics
#' @examples
#' \donttest{data("Antibiotics")
#' bayes.density <- TheBayesteApproach(Antibiotics)
#' print(bayes.density)
#' plot(bayes.density)}
plot.TheBayesteApproach <- function(x, xlim = c(6,50), xlab = "ZD",main = "Density estimator", ECOFF = FALSE, pi.level = 0.3, quantile = 0.01, ...){

  x0 <- seq(7,50,0.001)

  if(any(is.na(x$normal.params))){

    plot.hist <- hist(x$data,breaks  = seq(min(x$data),max(x$data)+1,by=1)-0.01, plot = FALSE)
    plot(plot.hist, freq = FALSE, xlim = xlim, xlab = xlab, main = main,...)
    lines(c(6.5,6.5), c(0,x$binom.res), col = "red",lwd = 2)

  }else{

    z <- matrix(rep(rep(0,length(x0)), x$`k`), nrow = x$`k`)

    for(i in 1:x$`k`){
      z[i,] <- dnorm(x0, x$normal.params[i], sd = sqrt(x$normal.params[(x$`k`*2+i)]))
      z[i,] <- z[i,]*x$normal.params[x$`k`+i]
    }


    z.total <- (1-x$binom.res)*colSums(z)

    plot.hist <- graphics::hist(x$data,breaks  = seq(min(x$data),max(x$data)+1,by=1)-0.01, plot = FALSE)
    y_max <- max(as.vector(z),plot.hist$density)
    plot(plot.hist, freq = FALSE, xlim = xlim, xlab = xlab, main = main,
         ylim = c(0,y_max),...)

    for(i in 1:x$`k`){
      lines(x0,z[i,], col = "green",lty=2)
    }

    lines(x0,z.total, col = "red", lwd = 2)


    lines(c(6.5,6.5), c(0,x$binom.res), col = "red",lwd = 2)

    if(ECOFF){
      ecoff <- FindECOFF(x, pi.level = pi.level, quantile = quantile)
      lines(x = c(ecoff,ecoff), y = c(0,y_max), col = "black", lwd = 3)
      text(x = ecoff, y = y_max+(y_max/20), labels = "ECOFF", xpd=NA)
    }
  }


  legend("topright", legend = c("Weighted estimated density", "Unweighted normal components"),
         col = c("red", "green"), lty = c(1, 2))

}
