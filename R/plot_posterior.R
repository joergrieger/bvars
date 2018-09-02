pltBvarPosterior <- function(bvarObj,lag = 1){

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

pltTvarPosterior <- function(tvarObj, lag = 1){

  # preliminary calculations
  nLength <- length(tvarObj$varnames)
  pltIntercepts1 <- list()
  pltPosterior1  <- list()
  pltIntercepts2 <- list()
  pltPosterior2  <- list()

  # Plot density of threshold variable
  thDf <- data.frame(threshold = tvarObj$tardraws)
  p1 <- ggplot(thDf,aes(threshold)) + geom_density()
  print(p1)
  readline("Press [Enter] to continue")

  # Plot Histogram of delay

  delDf <- data.frame(delay = as.integer(tv$deldraws))
  p1 <- ggplot(delDf,aes(delay)) + geom_histogram()
  print(p1)
  readline("Press [Enter] to continue")

  # Plot time series of regime
  regLength <- dim(tvarObj$regimes)
  regMean   <- array(0,dim=c(regLength))
  regUpper  <- array(0,dim=c(regLength))
  regLower  <- array(0,dim=c(regLength))
  for(ii in 1:regLength){
    regMean[ii]   <- mean(tvarObj$regime[ii,])
    regUpper[ii]  <- quantile(tvarObj$regime[ii,],probs = 0.95)
    regLower[ii]  <- quantile(tvarObj$regime[ii,],probs = 0.05)
  }
  regDf <- data.frame(x=seq(1:regLength),regime = as.vector(regMean), Upper = as.vector(regUpper), Lower = as.vector(regLower))
  print(regDf)
  readline("Press Enter")
  p1 <- ggplot(regDf) + geom_line(mapping=aes(x=x,y=regime)) + geom_line(mapping=aes(x=x,y=Upper)) + geom_line(mapping=aes(x=x,y=Lower))
  print(p1)
  readline("Press [Enter] to continue")

  # Plot the posterior draws of the intercept (if available)
  if(tvarObj$Intercept == TRUE){

    betaIntercept1 <- tvarObj$Alphadraws[1,,1,]
    betaIntercept2 <- tvarObj$Alphadraws[1,,2,]

    betaFinal1    <- tvarObj$Alphadraws[-c(1),,1,]
    betaFinal2    <- tvarObj$Alphadraws[-c(1),,2,]


    pltPostIntercept1 <- list()
    pltPostIntercept2 <- list()

    for(ii in 1:nLength){

      tempDf1 <- data.frame(beta = betaIntercept1[ii,])
      tempDf2 <- data.frame(beta = betaIntercept1[ii,])

      p1 <- ggplot(tempDf1,aes(beta)) + geom_density()
      p2 <- ggplot(tempDf2,aes(beta)) + geom_density()

      pltPostIntercept1[[ii]] <- p1
      pltPostIntercept2[[ii]] <- p2
    }

    # Plot intercepts
    do.call("grid.arrange",c(pltPostIntercept1,ncol=sqrt(nLength)))
    readline("Press [Enter] to continue")

    do.call("grid.arrange",c(pltPostIntercept2,ncol=sqrt(nLength)))
    readline("Press [Enter] to continue")

  }
  else{

    betaFinal1    <- tvarObj$Alphadraws[,,1,]
    betaFinal2    <- tvarObj$Alphadraws[,,2,]

  }

  pltPosterior1 <- list()
  pltPosterior2 <- list()

  for(lag in 1:tvarObj$NoLags){

    for(ii in 1:nLength){

      for(jj in 1:nLength){

        tempDf <- data.frame(xx = betaFinal1[ii + (lag - 1) * nLength,jj,])
        p1 <- ggplot(tempDf,aes(xx)) + geom_density()
        p1 <- p1 + xlab(paste("Posterior density for",tvarObj$varnames[ii],"on",tvarObj$varnames[2],sep=" "))
        pltPosterior1[[(jj-1)*nLength+ii]] <- p1


        tempDf <- data.frame(xx = betaFinal2[ii + (lag - 1) * nLength,jj,])
        p1 <- ggplot(tempDf,aes(xx)) + geom_density()
        p1 <- p1 + xlab(paste("Posterior density for",tvarObj$varnames[ii],"on",tvarObj$varnames[2],sep=" "))
        pltPosterior2[[(jj-1)*nLength+ii]] <- p1

      }

    }

    do.call("grid.arrange",c(pltPosterior1,ncol=nLength))
    readline("Press [Enter] to continue.")

    do.call("grid.arrange",c(pltPosterior2,ncol=nLength))
    readline("Press [Enter] to continue.")
  }

}
