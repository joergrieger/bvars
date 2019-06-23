irf.bvar  <- function(bvObj, nhor = 12, ncores = 1, irfquantiles = c(0.05,0.95), ident = 1, restrictions=NULL){

  intercept <- bvObj$intercept
  betadraws <- bvObj$betadraws
  sigmadraws <- bvObj$sigmadraws
  NoLags <- bvObj$NoLags
  nreps <- dim(betadraws)[3]
  k <- dim(sigmadraws)[1]
  irfdraws <- array(0,dim=c(k,k,nhor,nreps))
  irffinal <- array(0,dim=c(k,k,nhor,3))

  if(ncores>1 && !require(parallel)){

    stop("The parallel package has to be installed")

  }

  if(ncores == 1){ # No parallelization

    for(ii in 1:nreps){

      Alpha <- betadraws[,,ii]
      Sigma <- sigmadraws[,,ii]

      irf <- compirf(A=Alpha,Sigma=Sigma,NoLags=NoLags,intercept=intercept,nhor = nhor)
      irfdraws[,,,ii] <- irf

    }

  }
  else{ # parallel version

    # Register workers
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    xtmp <- foreach(ii = 1:nreps) %dopar% {
      print(ii)

      Alpha <- betadraws[,,ii]
      Sigma <- sigmadraws[,,ii]

      irf <- compirf(A=Alpha,Sigma=Sigma,NoLags=NoLags,intercept=intercept,nhor = nhor)

    }

    # stop workers
    parallel::stopCluster(cl)

    for(ii in 1:nreps){

      irfdraws[,,,ii] <- xtmp[[ii]]

    }

  }

  irflower <- min(irfquantiles)
  irfupper <- max(irfquantiles)
  for(jj in 1:k){
    for(kk in 1:k){
      for(ll in 1:nhor){

        irffinal[jj,kk,ll,1] <- median(irfdraws[jj,kk,ll,])
        irffinal[jj,kk,ll,2] <- quantile(irfdraws[jj,kk,ll,],probs=irflower)
        irffinal[jj,kk,ll,3] <- quantile(irfdraws[jj,kk,ll,],probs=irfupper)

      }
    }
  }

  relist <- structure(list(irfdraws=irffinal,irfhorizon=nhor,varnames=bvObj$varnames),class = "bvirf")

  return(relist)

}

irf.tvar  <- function(tvObj, nhor = 12, ncores = 1, irfquantiles = c(0.05,0.95), ident = 1, restrictions =NULL, bootrep=50){

  nLength    <- dim(tvObj$Alphadraws)[4]
  Alphadraws <- tvObj$Alphadraws
  Sigmadraws <- tvObj$Sigmadraws
  thVar      <- tvObj$thVar
  tardraws   <- tvObj$tardraws
  deldraws   <- tvObj$deldraws
  NoLags     <- tvObj$NoLags
  Intercept  <- tvObj$Intercept


  K <- dim(Sigmadraws[,,,1])[1]

  Irfdraws <- array(0,dim=c(K,K,nhor,2,nLength))

  if(ncores == 1){

    for(jj in 1:nLength){


      Alpha <- Alphadraws[,,,jj]
      Sigma <- Sigmadraws[,,,jj]
      tart <- tardraws[jj]
      thDelay <- deldraws[jj]

      xsplit <- splitVariables(y = tvObj$mydata,lags = NoLags, thDelay = thDelay, thresh = thVar, tart = tart, intercept = Intercept)

      for(ii in 1:K){
        if(ident==1){

          xx <- tirf(xsplit$ystar,xsplit$ytest,Alpha[,,1],Alpha[,,2],Sigma[,,1],Sigma[,,2],tart,thVar,thDelay,NoLags,nhor,Intercept,shockvar=ii,bootrep)

        }
        else if(ident==2){

          xx <- tirfsign(xsplit$ystar,xsplit$ytest,Alpha[,,1],Alpha[,,2],Sigma[,,1],Sigma[,,2],tart,thVar,thDelay,NoLags,nhor,Intercept,shockvar=ii,bootrep,restrictions=restrictions)

        }

        Irfdraws[ii,,,1,jj]<-xx$irf1
        Irfdraws[ii,,,2,jj]<-xx$irf2

      }

    }

  }
  else{ # Parallel version

    # Register workers
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    irftmp <- array(0,dim=c(K,K,nhor,2))

    xtmp <- foreach(jj = 1:nLength) %dopar% {

      Alpha <- Alphadraws[,,,jj]
      Sigma <- Sigmadraws[,,,jj]
      tart <- tardraws[jj]
      thDelay <- deldraws[jj]

      irftmp <- tirf1(tvObj$mydata,Alpha,Sigma,tart,thVar,thDelay,NoLags,nhor,Intercept,bootrep,ident,restrictions,K)

    }

    # stop workers
    parallel::stopCluster(cl)

    for(jj in 1:nLength){

      Irfdraws[,,,,jj] <- xtmp[[jj]]

    }

  }

  irffinal <- array(0,dim=c(K,K,nhor,2,3))

  lowerquantile=min(irfquantiles)
  upperquantile=max(irfquantiles)
  for(ii in 1:K){
    for(jj in 1:K){
      for(kk in 1:nhor){
        irffinal[ii,jj,kk,1,1] <- quantile(Irfdraws[ii,jj,kk,1,],probs=0.5)
        irffinal[ii,jj,kk,2,1] <- quantile(Irfdraws[ii,jj,kk,2,],probs=0.5)

        irffinal[ii,jj,kk,1,2] <- quantile(Irfdraws[ii,jj,kk,1,],probs=lowerquantile)
        irffinal[ii,jj,kk,2,2] <- quantile(Irfdraws[ii,jj,kk,2,],probs=lowerquantile)

        irffinal[ii,jj,kk,1,3] <- quantile(Irfdraws[ii,jj,kk,1,],probs=upperquantile)
        irffinal[ii,jj,kk,2,3] <- quantile(Irfdraws[ii,jj,kk,2,],probs=upperquantile)
      }
    }
  }

  relist <- structure(list(irf=irffinal,irfhorizon=nhor,varnames=tvObj$varnames),class = "tvirf")

  return(relist)

}

irf.ftvar <- function(ftObj, nhor = 12, ncores = 1, irfquantiles = c(0.05,0.95), ident = 1, restrictions = NULL, bootrep = 50){

  Intercept  <- ftObj$Intercept
  betadraws  <- ftObj$Betadraws
  sigmadraws <- ftObj$Sigmadraws
  Ldraws     <- ftObj$Ldraws
  thVar      <- ftObj$thVar
  tardraws   <- ftObj$tardraws
  deldraws   <- ftObj$deldraws

  NoLags     <- ftObj$NoLags


  nreps <- dim(betadraws)[4]
  P     <- dim(sigmadraws)[1]
  dimXY <- dim(Ldraws)[1]

  IrfSmallDraws <- array(0,dim=c(P,P,nhor,2,nreps))
  IrfLargeDraws <- array(0,dim=c(P,dimXY,nhor,2,nreps))

  if(ncores == 1){ # No parallelization

    for(ii in 1:nreps){

      print(ii)

      Beta    <- betadraws[,,,ii]
      Sigma   <- sigmadraws[,,,ii]
      tart    <- tardraws[ii]
      thDelay <- deldraws[ii]
      L       <- Ldraws[,,,ii]

      xsplit <- splitVariables(y = ftObj$mydata, lags = NoLags, thDelay = thDelay, thresh = thVar,
                               tart = tart, intercept = Intercept)


      for(jj in 1:P){

        xx <- tirf(xsplit$ystar, xsplit$ytest, Beta[,,1], Beta[,,2],
                   Sigma[,,1], Sigma[,,2], tart, thVar, thDelay, NoLags,
                   nhor, Intercept = Intercept, shockvar = jj,
                   bootrep = bootrep)

        IrfSmallDraws[jj,,,1,ii] <- xx$irf1
        IrfSmallDraws[jj,,,2,ii] <- xx$irf2

        IrfLargeDraws[jj,,,1,ii] <- L[,,1] %*% IrfSmallDraws[jj,,,1,ii]
        IrfLargeDraws[jj,,,2,ii] <- L[,,2] %*% IrfSmallDraws[jj,,,2,ii]

      }

    }

  }
  else{ # Parallel Version

    # Register clusters
    # Warning: No check if number of cores are supported

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    irfDraws <- array(0,dim=c(P,P,nhor,2))

    xtmp <- foreach(ii = 1:nreps) %dopar% {

      Beta    <- betadraws[,,,ii]
      Sigma   <- sigmadraws[,,,ii]
      tart    <- tardraws[ii]
      thDelay <- deldraws[ii]

      irftmp <- tirf1(ftObj$mydata, Beta, Sigma, tart, thVar, thDelay, NoLags, nhor, Intercept,
                      bootrep, ident = 1, NULL, P)

    }

    # shut down workers
    parallel::stopCluster(cl)

    for(ii in  1:nreps){

      IrfSmallDraws[,,,,ii] <- xtmp[[ii]]
      L <- Ldraws[,,,ii]

      for(jj in 1:P){

        IrfLargeDraws[jj,,,1,ii] <- L[,,1] %*% IrfSmallDraws[jj,,,1,ii]
        IrfLargeDraws[jj,,,2,ii] <- L[,,2] %*% IrfSmallDraws[jj,,,2,ii]

      }

    }

  }

  # Prepare return values
  IrfSmallFinal <- array(0,dim=c(P,P,nhor,2,3))
  IrfLargeFinal <- array(0,dim=c(P,dimXY,nhor,2,3))
  irflower <- min(irfquantiles)
  irfupper <- max(irfquantiles)

  for(jj in 1:P){
    for(kk in 1:P){
      for(ll in 1:nhor){

        # Regime 1
        IrfSmallFinal[jj,kk,ll,1,1] <- median(IrfSmallDraws[jj,kk,ll,1,])
        IrfSmallFinal[jj,kk,ll,1,2] <- quantile(IrfSmallDraws[jj,kk,ll,1,],probs=irflower)
        IrfSmallFinal[jj,kk,ll,1,3] <- quantile(IrfSmallDraws[jj,kk,ll,1,],probs=irfupper)

        # Regime 2
        IrfSmallFinal[jj,kk,ll,2,1] <- median(IrfSmallDraws[jj,kk,ll,2,])
        IrfSmallFinal[jj,kk,ll,2,2] <- quantile(IrfSmallDraws[jj,kk,ll,2,],probs=irflower)
        IrfSmallFinal[jj,kk,ll,2,3] <- quantile(IrfSmallDraws[jj,kk,ll,2,],probs=irfupper)

      }
    }
  }

  for(jj in 1:P){
    for(kk in 1:dimXY){
      for(ll in 1:nhor){
        # Regime 1
        IrfLargeFinal[jj,kk,ll,1,1] <- mean(IrfLargeDraws[jj,kk,ll,1,])
        IrfLargeFinal[jj,kk,ll,1,2] <- quantile(IrfLargeDraws[jj,kk,ll,1,],probs=irflower)
        IrfLargeFinal[jj,kk,ll,1,3] <- quantile(IrfLargeDraws[jj,kk,ll,1,],probs=irfupper)

        # Regime 2
        IrfLargeFinal[jj,kk,ll,2,1] <- mean(IrfLargeDraws[jj,kk,ll,2,])
        IrfLargeFinal[jj,kk,ll,2,2] <- quantile(IrfLargeDraws[jj,kk,ll,2,],probs=irflower)
        IrfLargeFinal[jj,kk,ll,2,3] <- quantile(IrfLargeDraws[jj,kk,ll,2,],probs=irfupper)

      }
    }
  }


  relist <- structure(list(IrfLarge=IrfLargeFinal,IrfSmall=IrfSmallFinal,irfhorizon=nhor,varnames=ftObj$varnames),class = "tvirf")
  return(relist)
}

irf.favar <- function(fvObj, nhor = 12, ncores = 1, irfquantiles = c(0.05,0.95), ident = 1, restrictions = NULL){

  # Preliminaries

  Intercept  <- fvObj$Intercept
  Betadraws  <- fvObj$Betadraws
  Sigmadraws <- fvObj$Sigmadraws
  Ldraws     <- fvObj$Ldraws
  NoLags     <- fvObj$NoLags
  Intercept  <- fvObj$Intercept

  nreps      <- dim(Betadraws)[3]
  k          <- dim(Sigmadraws)[1]
  dimXY <- dim(Ldraws)[1]

  IrfSmallDraws   <- array(0,dim=c(k,k,nhor,nreps))
  IrfLargeDraws   <- array(0,dim=c(k,dimXY,nhor,nreps))

  if(ncores>1 && !require(parallel)){

    stop("The parallel package has to be installed")

  }

  if(ncores == 1){

    # No Parallelization
    for(ii in 1:nreps){

      Alpha <- Betadraws[,,ii]
      Sigma <- Sigmadraws[,,ii]
      L     <- Ldraws[,,ii]

      irf <- compirf(A=Alpha,Sigma=Sigma,NoLags=NoLags,intercept=Intercept,nhor = nhor)
      IrfSmallDraws[,,,ii] <- irf

      for(jj in 1:k){

        IrfLargeDraws[jj,,,ii] <- L[,]%*%irf[jj,,]

      } # End over expanding impulse-response functions
    }# end over loop
  }
  else{

    # Register workers
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    xtmp <- foreach(ii = 1:nreps) %dopar% {

      Alpha <- Betadraws[,,ii]
      Sigma <- Sigmadraws[,,ii]

      irf <- compirf(A=Alpha,Sigma=Sigma,NoLags=NoLags,intercept=Intercept,nhor = nhor)

    } # End over parallel loop

    for(ii in 1:nreps){

      L <- Ldraws[,,ii]
      IrfSmallDraws[,,,ii] <- xtmp[[ii]]

      for(jj in 1:k){

        IrfLargeDraws[jj,,,ii] <- L[,]%*%IrfSmallDraws[jj,,,ii]

      } # end over expanding impulse-response functions
    } # end over extracting irfs from the list
  }
  # Final computations
  IrfSmallFinal <- array(0,dim=c(k,k,nhor,3))
  IrfLargeFinal <- array(0,dim=c(k,dimXY,nhor,3))
  irflower <- min(irfquantiles)
  irfupper <- max(irfquantiles)

  for(jj in 1:k){
    for(kk in 1:k){
      for(ll in 1:nhor){

        IrfSmallFinal[jj,kk,ll,1] <- quantile(IrfSmallDraws[jj,kk,ll,],probs=0.5)
        IrfSmallFinal[jj,kk,ll,2] <- quantile(IrfSmallDraws[jj,kk,ll,],probs=irflower)
        IrfSmallFinal[jj,kk,ll,3] <- quantile(IrfSmallDraws[jj,kk,ll,],probs=irfupper)

      }
    }
  }
  for(jj in 1:k){
    for(kk in 1:dimXY){
      for(ll in 1:nhor){

        IrfLargeFinal[jj,kk,ll,1] <- quantile(IrfLargeDraws[jj,kk,ll,],probs=0.5)
        IrfLargeFinal[jj,kk,ll,2] <- quantile(IrfLargeDraws[jj,kk,ll,],probs=irflower)
        IrfLargeFinal[jj,kk,ll,3] <- quantile(IrfLargeDraws[jj,kk,ll,],probs=irfupper)

      }
    }
  }

  # Returning values
  relist <- structure(list(IrfLarge=IrfLargeFinal,irfhorizon=nhor,varnames=fvObj$varnames),class = "fvirf")
  return(relist)

}

tirf1 <- function(y,Alpha,Sigma,tart,thVar,thDelay,NoLags,nhor,Intercept,bootrep,ident,restrictions,K){



  xsplit <- splitVariables(y=y,lags=NoLags,thDelay=thDelay,thresh=thVar,tart=tart,intercept=Intercept)

  Irfdraws <- array(0,dim=c(K,K,nhor,2))

  for(ii in 1:K){
    if(ident==1){

      xx <- tirf(xsplit$ystar,xsplit$ytest,Alpha[,,1],Alpha[,,2],Sigma[,,1],Sigma[,,2],tart,thVar,thDelay,NoLags,nhor,Intercept,shockvar=ii,bootrep)

    }
    else if(ident==2){

      xx <- tirfsign(xsplit$ystar,xsplit$ytest,Alpha[,,1],Alpha[,,2],Sigma[,,1],Sigma[,,2],tart,thVar,thDelay,NoLags,nhor,Intercept,shockvar=ii,bootrep,restrictions=restrictions)

    }

    Irfdraws[ii,,,1]<-xx$irf1
    Irfdraws[ii,,,2]<-xx$irf2

  }

  return(Irfdraws)

}


