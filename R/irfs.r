#
# internal function to calculate one path of a VAR model (without intercept) given beta,sigma and a shock
#
.irfsimu <- function(beta,sigma,shocks,horizon,lags,intercept,K){
  nl <- horizon+lags
  yhat <- as.matrix(array(0,dim=c(K,nl)))
  for(ii in 1:horizon){
    indxmax <- lags+ii-1
    indxmin <- ii
    ytmp <- yhat[,indxmax:indxmin]
    ytmp <- matrix(ytmp,nrow=1)
    tmp1 <- shocks[,ii+lags]%*%sigma
    tmp2 <- ytmp%*%beta+tmp1
    yhat[,ii+lags]<-tmp2

  }
  yhat <- yhat[,(lags+1):nl]
  return(yhat)
}

.tirf <- function(K,irfhorizon,lags,ytest,xstar,beta1,beta2,sigma1,sigma2,tart,thresh,thDelay,shockvar=1){
  csigma1 <- t(chol(sigma1))
  csigma2 <- t(chol(sigma2))

  irf1 <- 0  # The estimated impulse-response functions for regime 1
	irf2 <- 0  # The estimated impulse-response functions for regime 2
	T <- nrow(xstar)
	T <- T-thDelay

  for(jj in 1:T){
    # take estimates and reverse order
	  xvar0 <- .vec2matrix(lags,K,xstar[jj,])
		xvar1 <- xvar0[nrow(xvar0):1,]
		# additional running variable for the threshold variable
		yrun <- (matrix(0,nrow=irfhorizon))
		yrun[1:thDelay] <- ytest[jj:(jj+thDelay-1)]
		if(ytest[jj]<tart){
			 irf1x <- .gtirf(K,irfhorizon,lags,yrun,xvar1,beta1,beta2,csigma1,csigma2,tart,thresh,thDelay,shockvar)
			 irf1  <- irf1+irf1x
		 }
		 else{
			 irf2x <- .gtirf(K,irfhorizon,lags,yrun,xvar1,beta1,beta2,sigma1,sigma2,tart,thresh,thDelay,shockvar)
			 irf2  <- irf2+irf2x
		 }
		 #readline(prompt="Press [enter] to continue")

     }
	 return(list(irf1=irf1,irf2=irf2))
}
.tirfSign<-function(K,irfhorizon,lags,ytest,xstar,beta1,beta2,sigma1,sigma2,tart,thresh,thDelay,shockvar=1,Restrictions){

  csigma1 <- t(chol(sigma1))
  csigma2 <- t(chol(sigma2))

  irf1 <- 0  # The estimated impulse-response functions for regime 1
  irf2 <- 0  # The estimated impulse-response functions for regime 2

  T <- nrow(xstar)
  T <- T-thDelay

  # Sign restrictions using Haar-measure
  SignRestriction <-FALSE
  while(!SignRestriction){
    # Restrictions for first regime
    qrmatrix1     <- matrix(rnorm(K*K),nrow=K)
    qrdecomp1     <- qr(qrmatrix1)
    qrdecomp1     <- qr.Q(qrdecomp1)
    testmatrix1   <- qrdecomp1%*%csigma1
    SignRestriction <- !.CheckSign(Restrictions,testmatrix1)
  }
  SignRestriction <-FALSE
  while(!SignRestriction){
    # Restrictions for second regime
    qrmatrix2     <- matrix(rnorm(K*K),nrow=K)
    qrdecomp2     <- qr(qrmatrix2)
    qrdecomp2     <- qr.Q(qrdecomp2)
    testmatrix2   <- qrdecomp2%*%csigma2
    SignRestriction <- !.CheckSign(Restrictions,testmatrix2)
  }
  cholsigma1 <- testmatrix1
  cholsigma2 <- testmatrix2
  # Start simulating path
  for(jj in 1:T){
    # take estimates and reverse order
    xvar0 <- .vec2matrix(lags,K,xstar[jj,])
    xvar1 <- xvar0[nrow(xvar0):1,]
    # additional running variable for the threshold variable
    yrun <- (matrix(0,nrow=irfhorizon))
    yrun[1:thDelay] <- ytest[jj:(jj+thDelay-1)]
    if(ytest[jj]<tart){
      irf1x <- .gtirf(K,irfhorizon,lags,yrun,xvar1,beta1,beta2,cholsigma1,cholsigma2,tart,thresh,thDelay,shockvar)
      irf1  <- irf1+irf1x
    }
    else{
      irf2x <- .gtirf(K,irfhorizon,lags,yrun,xvar1,beta1,beta2,cholsigma1,cholsigma2,tart,thresh,thDelay,shockvar)
      irf2  <- irf2+irf2x
    }
    #readline(prompt="Press [enter] to continue")

  }
  return(list(irf1=irf1,irf2=irf2))



}
.gtirf <- function(K,irfhorizon,lags,yrun,xstar,beta1,beta2,csigma1,csigma2,tart,thresh,thDelay,shockvar){

     # Cholesky Decomposition of Variance-Covariance Matrix

	   xtmatrix <- array(0,dim=c(K,irfhorizon+lags))
	   ytmatrix <- array(0,dim=c(K,irfhorizon+lags))
	   ytmatrix[,1:lags]<-t(xstar)

	   # Shock matrices
	   uu1 <- array(0,dim=c(irfhorizon+lags,K))
	   uu2 <- array(0,dim=c(irfhorizon+lags,K))
	   uu1[lags+1,shockvar]<-1
	   uu2[lags+1,shockvar]<-1
	   # Start simulating path

	   for(ii in 1:irfhorizon){
		   yt <- ytmatrix[,ii:(ii+lags-1)]
		   xt <- xtmatrix[,ii:(ii+lags-1)]
		   xr <- array(0,dim=c(K*lags))
		   if(lags>1){
		     for(jj in 1:lags){
			     for(kk in 1:K){
				     #yr[(jj-1)*K+kk]<-yt[kk,lags-jj+1]
					 xr[(jj-1)*K+kk]<-xt[kk,lags-jj+1]
				 }
			 }
		 }
		 else{
			 xr <- xtmatrix[,ii]
		 }
		 # Get next periods values
		 if(yrun[ii]<tart){
		   temp <- t(xr%*%beta1)+t(uu1[ii+lags,]%*%csigma1)
			 ytmatrix[,ii+lags] <- ytmatrix[,ii+lags-1]+temp
			 xtmatrix[,ii+lags] <- temp
		 }
		 else{
			 temp <- t(xr%*%beta1)+t(uu1[ii+lags,]%*%csigma1)
			 ytmatrix[,ii+lags] <- ytmatrix[,ii+lags-1]+temp
			 xtmatrix[,ii+lags] <- temp
		 }
		 # Get next periods threshold variable
		 if(ii>=thDelay && ii<irfhorizon){
		     yrun[ii+1] <- ytmatrix[thresh,ii+lags+1-thDelay]
		 }

	 }
	 return(xtmatrix[,(lags+1):(lags+irfhorizon)])
 }

