pltBvarPosterior <- function(bvarObj,lag=1){

  # preliminary calculations
  nLength <- length(bvarObj$varnames)
  pltIntercepts <- list()
  pltPosterior <- list()

  # If Model has intercept, remove intercept and store posterior for intercept in different variable
  if(bvarObj$intercept == TRUE){
    #
    betaFinalDf <- bvarObj$betadraws[-c(1),,]
    bvarIntercept <- bvarObj$betadraws[1,,]


    for(ii in 1:nLength){
      tempDf <- data.frame(Intercept = bvarIntercept[ii,])
      p1 <- ggplot(tempDf,aes(Intercept)) + geom_density()
      p1 <- p1 + xlab(paste("Intercept on",bvarObj$varnames[ii]))
      pltIntercepts[[ii]] <- p1
    }

    nc <- sqrt(nLength)
    do.call("grid.arrange",c(pltIntercepts,ncol=nc))

  }
  else{
    betaFinalDF <- bvarObj$betadraws
  }

  # Posterior densities for VAR-Coefficients
  for(ii in 1:nLength){

    for(jj in 1:nLength){

      tempDf <- data.frame(xx = betaFinalDf[ii + (lag - 1) * nLength,jj,])
      p1 <- ggplot(tempDf,aes(xx)) + geom_density()
      p1 <- p1 + xlab(paste("Posterior density for",bvarObj$varnames[ii],"on",bvarObj$varnames[2],sep=" "))
      pltPosterior[[(jj-1)*nLength+ii]] <- p1

    }

  }

  do.call("grid.arrange",c(pltPosterior,ncol=nLength))

}
