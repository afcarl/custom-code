compareDistributions <- function (Ds) {
  #chisq require(distr)
  require(fitdistrplus)
  dens1 <- density(Ds,from=0)
  range_x1 <- range(dens1$x)
  range_y1 <- range(dens1$y)
  fit <- fitdist(Ds,"lnorm")
  dens2 <- density(rlnorm(1000, meanlog=fit$estimate[1], sdlog=fit$estimate[2]),from=0)
  range_x2 <- range(dens2$x)
  range_y2 <- range(dens2$y)
  plot(dens1, xlim=c(0,max(range_x1,range_x2)*1.25), ylim=c(0,max(range_y1,range_y2)*1.1))
  # chisq lines(seq(0,max(range_x1,range_x2)*1.25,0.01),dchisq(seq(0,max(range_x1,range_x2)*1.25,0.01),df=mean(Ds)), col="green")
  # chisq lines(seq(0,max(range_x1,range_x2)*1.25,0.01),dchisq(seq(0,max(range_x1,range_x2)*1.25,0.01),df=igamma(mean(log(Ds/2)))*2), col="darkgreen",lwd=1.5)
  lines(seq(0,max(range_x1,range_x2)*1.25,0.01), dlnorm(seq(0,max(range_x1,range_x2)*1.25,0.01), meanlog=fit$estimate[1], sdlog=fit$estimate[2]), col="darkgreen",lwd=1.5)
}

compareDistributions <- function (Ds) {
  #chisq require(distr)
  require(fitdistrplus)
  dens1 <- density(Ds,from=0)
  range_x1 <- range(dens1$x)
  range_y1 <- range(dens1$y)
  fit <- fitdist(Ds,"lnorm")
  plot(dens1, xlim=c(0,max(range_x1)*1.25), ylim=c(0,max(range_y1)*1.1))
  # chisq lines(seq(0,max(range_x1,range_x2)*1.25,0.01),dchisq(seq(0,max(range_x1,range_x2)*1.25,0.01),df=mean(Ds)), col="green")
  # chisq lines(seq(0,max(range_x1,range_x2)*1.25,0.01),dchisq(seq(0,max(range_x1,range_x2)*1.25,0.01),df=igamma(mean(log(Ds/2)))*2), col="darkgreen",lwd=1.5)
  lines(seq(0,max(range_x1)*1.25,0.01), dlnorm(seq(0,max(range_x1)*1.25,0.01), meanlog=fit$estimate[1], sdlog=fit$estimate[2]), col="darkgreen",lwd=1.5)
}
