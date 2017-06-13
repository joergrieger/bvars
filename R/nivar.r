#'@export
nivar <- function(mydata,lags=1,intercept=TRUE,coefprior=NULL,coefpriorvar=1,varprior=10,varpriordof=10,irfhorizon=16,irfquantiles=c(0.1,0.9),reps=500,burnin=100,stabletest=TRUE){
  y<- as.matrix(mydata)
  T   <- nrow(y)
  K   <- ncol(y)
  obs <- T-lags

  err <- .nierror(mydata,lags=lags,intercept=intercept,coefprior=coefprior,coefpriorvar=coefpriorvar,varprior=varprior,varpriordof=varpriordof,irfhorizon=irfhorizon,irfquantiles=irfquantiles,reps=reps,burnin=burnin,stabletest=stabletest)
  prior <- .niprior(mydata,lags,intercept,coefprior,coefpriorvar,varprior,varpriordof,irfhorizon,irfquantiles,reps,burnin,stabletest)
  lagdata <- .data(y,lags,intercept);
  results <- .nigibbs(lagdata$y,lagdata$x,lags,intercept,prior$coefprior,prior$coefpriorvar,prior$varprior,varpriordof,irfhorizon,irfquantiles,reps,burnin,stabletest)
  finalresults <- structure(list(prior=prior,betadist=results$Betadraws,Sigmadist=results$Sigmadraws,irf=results$irf),class="nivar")
  return(finalresults)
}

nivarSign <- function(mydata,lags=1,intercept=TRUE,coefprior=NULL,coefpriorvar=1,varprior=10,varpriordof=10,irfhorizon=16,irfquantiles=c(0.1,0.9),reps=500,burnin=100,Restrictions,stabletest=TRUE){
  y<- as.matrix(mydata)
  T   <- nrow(y)
  K   <- ncol(y)
  obs <- T-lags
  err <- .nierror(mydata,lags=lags,intercept=intercept,coefprior=coefprior,coefpriorvar=coefpriorvar,varprior=varprior,varpriordof=varpriordof,irfhorizon=irfhorizon,irfquantiles=irfquantiles,reps=reps,burnin=burnin,stabletest=stabletest)
  prior <- .niprior(mydata,lags,intercept,coefprior,coefpriorvar,varprior,varpriordof,irfhorizon,irfquantiles,reps,burnin,stabletest)
  lagdata <- .data(y,lags,intercept);
  results <- .nigibbsSign(lagdata$y,lagdata$x,lags,intercept,prior$coefprior,prior$coefpriorvar,prior$varprior,varpriordof,irfhorizon,irfquantiles,reps,burnin,Restrictions,stabletest)
  finalresults <- structure(list(prior=prior,betadist=results$Betadraws,Sigmadist=results$Sigmadraws,irf=results$irf),class="nivar")
  return(finalresults)
}

.nierror <- function(mydata,lags,intercept,coefprior,coefpriorvar,varprior,varpriordof,irfhorizon,irfquantiles,reps,burnin,stabletest){
  y<- as.matrix(mydata)
  T   <- nrow(y)
  K   <- ncol(y)
  obs <- T-lags
  constant=0
  if(intercept==TRUE) constant=1

  # some sanity checks
  if(anyNA(y)==TRUE){stop("Data contains N/As")}  # Test if data contains N/As
  if(K<2){stop("Only a univariate time series")}  # Test if data is univaraite
  if(T<lags){stop("Number of lags must be greater than the number of observations")} # Test if lags are smaller than the number of observations

  # test if provided priors have the correct form

  # Prior on coefficients
  print(coefprior)
  nrcoefpriorvar <- K*(K*lags+constant) # number of rows and columns for the prior variance of the coefficient
  if(!is.null(coefprior)){
    nrcoefprior <- nrow(coefprior)
    nccoefprior <- ncol(coefprior)
    nr <- K*lags+constant
    nc <- K
    if(nrcoefprior!=nr){stop("Incorrect number of rows in coefprior")}
    if(nccoefprior!=nc){stop("Incorrent number of columns in coefprior")}
  }
  if(!.isscalar(coefpriorvar)){
	nrtcoefpriorvar <- nrow(coefpriorvar)
    nctcoefpriorvar <- ncol(coefpriorvar)
    if(nrtcoefpriorvar!=nrcoefpriorvar){stop("Incorrect number of rows in coefpriorvar")}
    if(nctcoefpriorvar!=nrcoefpriorvar){stop("Incorrect number of columns in coefpriorvar")}
  }
  # Prior on Variance

  if(!.isscalar(varprior)){
    nrvarprior <- nrow(varprior)
    ncvarprior <- ncol(varprior)
    if(nrvarprior!=K){stop("Prior on variance must be of size KxK")}
    if(ncvarprior!=K){stop("Prior on variance must be of size KxK")}
  }
}
.niprior <- function(mydata,lags,intercept,coefprior,coefpriorvar,varprior,varpriordof,irfhorizon,irfquantiles,reps,burnin,stabletest){
  y<- as.matrix(mydata)
  T   <- nrow(y)
  K   <- ncol(y)
  obs <- T-lags
  # prior on mean
  constant=0
  if(intercept==TRUE){constant=1}
  if(is.null(coefprior)){
    coefprior <- array(0,dim=c(K*lags+constant,K))
  }
  nrcoefpriorvar <- K*(K*lags+constant) # number of rows and columns for a provided prior.
  if(.isscalar(coefpriorvar)){
    coefpriorvar <- coefpriorvar*.id(nrcoefpriorvar)
  }
  # prior on variances
  if(.isscalar(varprior)){
    varprior <- varprior*.id(K)
  }
  prior <- list(coefprior=coefprior,coefpriorvar=coefpriorvar,varprior=varprior)
  return(prior)
}
.nigibbs <- function(y,x,lags,intercept,coefprior,coefpriorvar,varprior,varpriordof,irfhorizon,irfquantiles,reps,burnin,stabletest){
  T   <- nrow(y)
  K   <- ncol(y)
  obs <- T-lags
  constant=0
  if(intercept==TRUE) constant=1
  # OLS estimates
  Aols <- solve(t(x)%*%x)%*%t(x)%*%y
  aols <- .vec(Aols)
  SSE  <- t(y-x%*%Aols)%*%(y-x%*%Aols)
  SIGMA_OLS <- SSE/(T-K+1)

  # Prepare Gibbs Sampling
  aprior    <- .vec(coefprior)
  Vprior    <- coefpriorvar
  vprior    <- varpriordof
  Sprior    <- varprior
  invSprior <- solve(Sprior)
  # Initialize parameters for Gibbs sampling
  alpha    <- aols
  Alpha    <- Aols
  SSEGibbs <- SSE
  SIGMA    <- SIGMA_OLS

  z <- .id(K)%x%x
  Betadraws <- array(0,dim=c(K*lags+constant,K,reps-burnin))
  Sigmadraws <- array(0,dim=c(K,K,reps-burnin))
  irfs <- array(0,dim=c(K,irfhorizon,K,reps-burnin))
  # Start Gibbs Sampling
  for(ii in 1:reps){
    print(ii)
    stable<-2 # only accept draws for VARS that are stable
    while(stable>1){
      variance <- solve(SIGMA)%x%.id(T)
      Vpost    <- solve(solve(Vprior)+t(z)%*%variance%*%z)
      apost    <- Vpost%*%(solve(Vprior)%*%aprior+t(z)%*%variance%*%.vec(y))
      alpha    <- mvrnorm(mu=apost,Sigma=Vpost) # Draw Alphas from a multivariate normal distribution
      Alpha1   <- matrix(alpha,ncol=K)
      # Draw variance
      vpost    <- obs+vprior
      Spost    <- Sprior+t(y-x%*%Alpha)%*%(y-x%*%Alpha)
      SIGMA1  <- solve(rWishart(1,vpost,Spost)[,,1]) # Draw Sigmas from a Wishart distribution

      if(intercept==TRUE){
        Alphatest <- Alpha1[2:(K*lags+1),]
      }
      else{
        Alphatest <- Alpha1
      }
      stable <- .stability(Alphatest,lags,K)
      if(stabletest==FALSE){stable=0} #
    }
    Alpha <- Alpha1
    SIGMA <- SIGMA1
    if(ii>burnin){
      Betadraws[,,ii-burnin]<-Alpha
      Sigmadraws[,,ii-burnin]<-SIGMA
      # Calculate Impulse-Response functions
      for(jj in 1:K){ # loop over the variables
        shock <- array(0,dim=c(K,irfhorizon+lags))
        shock[jj,lags+1]<-1
        cholsigma <- t(chol(SIGMA))
        yhat <- .irfsimu(beta=Alphatest,shocks=shock,sigma=cholsigma,lags=lags,horizon=irfhorizon,intercept=intercept,K=K)
        irfs[,,jj,ii-burnin]<-yhat
      }
    }
  }
  upperquantile <- max(irfquantiles)
  lowerquantile <- min(irfquantiles)
  irffinal <- array(0,dim=c(K,irfhorizon,K,3))
  for(ii in 1:K){
    for(jj in 1:K){
      for(kk in 1:irfhorizon){
        irffinal[jj,kk,ii,1]<-quantile(irfs[jj,kk,ii,],probs=0.5)
        irffinal[jj,kk,ii,2]<-quantile(irfs[jj,kk,ii,],probs=lowerquantile)
        irffinal[jj,kk,ii,3]<-quantile(irfs[jj,kk,ii,],probs=upperquantile)
      }
    }
  }
  return(list(Sigmadraws=Sigmadraws,Betadraws=Betadraws,irf=irffinal))
}

.nigibbsSign <- function(y,x,lags,intercept,coefprior,coefpriorvar,varprior,varpriordof,irfhorizon,irfquantiles,reps,burnin,Restrictions,stabletest){
  T   <- nrow(y)
  K   <- ncol(y)
  obs <- T-lags
  constant=0
  if(intercept==TRUE) constant=1
  # OLS estimates
  Aols <- solve(t(x)%*%x)%*%t(x)%*%y
  aols <- .vec(Aols)
  SSE  <- t(y-x%*%Aols)%*%(y-x%*%Aols)
  SIGMA_OLS <- SSE/(T-K+1)

  # Prepare Gibbs Sampling
  aprior    <- .vec(coefprior)
  Vprior    <- coefpriorvar
  vprior    <- varpriordof
  Sprior    <- varprior
  invSprior <- solve(Sprior)
  # Initialize parameters for Gibbs sampling
  alpha    <- aols
  Alpha    <- Aols
  SSEGibbs <- SSE
  SIGMA    <- SIGMA_OLS

  z <- .id(K)%x%x
  Betadraws <- array(0,dim=c(K*lags+constant,K,reps-burnin))
  Sigmadraws <- array(0,dim=c(K,K,reps-burnin))
  irfs <- array(0,dim=c(K,irfhorizon,K,reps-burnin))
  # Start Gibbs Sampling
  for(ii in 1:reps){
    print(ii)
    stable<-2 # only accept draws for VARS that are stable
    while(stable>1){
      variance <- solve(SIGMA)%x%.id(T)
      Vpost    <- solve(solve(Vprior)+t(z)%*%variance%*%z)
      apost    <- Vpost%*%(solve(Vprior)%*%aprior+t(z)%*%variance%*%.vec(y))
      alpha    <- mvrnorm(mu=apost,Sigma=Vpost) # Draw Alphas from a multivariate normal distribution
      Alpha1   <- matrix(alpha,ncol=K)
      # Draw variance
      vpost    <- obs+vprior
      Spost    <- Sprior+t(y-x%*%Alpha)%*%(y-x%*%Alpha)
      SIGMA1  <- solve(rWishart(1,vpost,Spost)[,,1]) # Draw Sigmas from a Wishart distribution

      if(intercept==TRUE){
        Alphatest <- Alpha1[2:(K*lags+1),]
      }
      else{
        Alphatest <- Alpha1
      }
      stable <- .stability(Alphatest,lags,K)
      if(stabletest==FALSE){stable=0} #
    }
    Alpha <- Alpha1
    SIGMA <- SIGMA1
    if(ii>burnin){
      Betadraws[,,ii-burnin]<-Alpha
      Sigmadraws[,,ii-burnin]<-SIGMA
      # Calculate Impulse-Response functions
      for(jj in 1:K){ # loop over the variables
        shock <- array(0,dim=c(K,irfhorizon+lags))
        shock[jj,lags+1]<-1
        cholsigma <- t(chol(SIGMA))
        # Implement Sign Restrictions
        SignRestriction <-FALSE
        while(!SignRestriction){
          qrmatrix <- matrix(rnorm(K*K),nrow=K)
          qrdecomp <- qr(qrmatrix)
          qrdecomp <- qr.Q(qrdecomp)
          testmatrix <- qrdecomp%*%cholsigma
          SignRestriction<-!.CheckSign(Restrictions,testmatrix)
        }
        cholsigma<-testmatrix
        yhat <- .irfsimu(beta=Alphatest,shocks=shock,sigma=cholsigma,lags=lags,horizon=irfhorizon,intercept=intercept,K=K)
        irfs[,,jj,ii-burnin]<-yhat
      }
    }
  }
  upperquantile <- max(irfquantiles)
  lowerquantile <- min(irfquantiles)
  irffinal <- array(0,dim=c(K,irfhorizon,K,3))
  for(ii in 1:K){
    for(jj in 1:K){
      for(kk in 1:irfhorizon){
        irffinal[jj,kk,ii,1]<-quantile(irfs[jj,kk,ii,],probs=0.5)
        irffinal[jj,kk,ii,2]<-quantile(irfs[jj,kk,ii,],probs=lowerquantile)
        irffinal[jj,kk,ii,3]<-quantile(irfs[jj,kk,ii,],probs=upperquantile)
      }
    }
  }
  return(list(Sigmadraws=Sigmadraws,Betadraws=Betadraws,irf=irffinal))
}
