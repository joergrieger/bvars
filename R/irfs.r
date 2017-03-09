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

.tirf <- function(K,T,irfhorizon,lags,Y,beta1,beta2,sigma1,sigma2,tart,thresh,thDelay,shockvar=1){
  irf <- 0
  for(jj in (lags+1):T){
    y1  <-matrix(Y[jj,],nrow=K)
    y0   <- y1[,ncol(y1):1]
    irfx <- .gtirf(K,irfhorizon,lags,y0,beta1,beta2,sigma1,sigma2,tart,thresh,thDelay,shockvar)
    irf <- irf+irfx
  }
  irfret <- irf/(T-(lags+1))
  return(irfret)
}
.gtirf <- function(K,irfhorizon,lags,Y,beta1,beta2,sigma1,sigma2,tart,thresh,thDelay,shockvar){
  csigma1 <- t(chol(sigma1))
  csigma2 <- t(chol(sigma2))
  ytmatrix <- array(0,dim=c(K,irfhorizon+lags))
  xtmatrix <- array(0,dim=c(K,irfhorizon+lags))
  ytmatrix[,1:lags]<-Y
  uu1 <- array(0,dim=c(irfhorizon+lags,K))
  uu1[lags+1,shockvar]<-1

  uu2 <- array(0,dim=c(irfhorizon+lags,K))
  uu2[lags+1,shockvar]<-1

  irf <-0

  for(ii in (lags+1):(irfhorizon+lags)){
    yt <- ytmatrix[,(ii-lags):(ii-1)]
    xt <- xtmatrix[,(ii-lags):(ii-1)]
    yr <- array(0,dim=c(K*lags))
    xr <- array(0,dim=c(K*lags))
    if(lags>1){
      for(jj in 1:lags){
        for(kk in 1:K){
          yr[(jj-1)*K+kk]<-yt[kk,lags-jj+1]
          xr[(jj-1)*K+kk]<-xt[kk,lags-jj+1]
        }
      }
    }
    else{
      yr <- ytmatrix[,ii-1]
      xr <- xtmatrix[,ii-1]
      #print(yr)
      #print(xr)
      #readline(prompt="Press [enter] to continue")
    }

    if(ytmatrix[thresh,ii-thDelay]<tart){
      # Regime 1
      ytmatrix[,ii]<-t(yr%*%beta1)+t(uu1[ii,]%*%csigma1)
      xtmatrix[,ii]<-t(xr%*%beta1)+t(uu1[ii,]%*%csigma1)
    }
    else
    {
      # Regime 2
      ytmatrix[,ii]<-t(yr%*%beta2)+t(uu2[ii,]%*%csigma2)
      xtmatrix[,ii]<-t(xr%*%beta2)+t(uu2[ii,]%*%csigma2)

    }
  }
  return(xtmatrix[,(lags+1):(irfhorizon+lags)])
}

