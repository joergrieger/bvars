irf.bvar <- function(bvObj, nhor = 12, ncores = 1,irfquantiles=c(0.05,0.95),ident=1,restrictions=NULL){

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

    # Register
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    xtmp <- foreach(ii = 1:nreps) %dopar% {
      print(ii)

      Alpha <- betadraws[,,ii]
      Sigma <- sigmadraws[,,ii]

      irf <- compirf(A=Alpha,Sigma=Sigma,NoLags=NoLags,intercept=intercept,nhor = nhor)

    }

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

irf.tvar <- function(tvObj, nhor=12, ncores=1,irfquantiles = c(0.05,0.95),ident=1,restrictions=NULL,bootrep=50){

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

      xsplit <- splitVariables(y=tvObj$mydata,lags=NoLags,thDelay=thDelay,thresh=thVar,tart=tart,intercept=Intercept)



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

    # Register
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    irftmp <- array(0,dim=c(K,K,nhor,2))

    xtmp <- foreach(jj = 1:nLength) %dopar% {




      Alpha <- Alphadraws[,,,jj]
      Sigma <- Sigmadraws[,,,jj]
      tart <- tardraws[jj]
      thDelay <- deldraws[jj]

      xsplit <- splitVariables(y=tvObj$mydata,lags=NoLags,thDelay=thDelay,thresh=thVar,tart=tart,intercept=Intercept)



      for(ii in 1:K){
        if(ident==1){

          xx <- tirf(xsplit$ystar,xsplit$ytest,Alpha[,,1],Alpha[,,2],Sigma[,,1],Sigma[,,2],tart,thVar,thDelay,NoLags,nhor,Intercept,shockvar=ii,bootrep)

        }
        else if(ident==2){

          xx <- tirfsign(xsplit$ystar,xsplit$ytest,Alpha[,,1],Alpha[,,2],Sigma[,,1],Sigma[,,2],tart,thVar,thDelay,NoLags,nhor,Intercept,shockvar=ii,bootrep,restrictions=restrictions)

        }

        irftmp[ii,,,1] <- xx$irf1
        irftmp[ii,,,2] <- xx$irf2
      }

    }
    print(xtmp[[1]])
    readline("Press [Enter] to continue")

  }




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

}


