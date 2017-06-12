.tiprior<-function(y,lags,thMax,intercept,coefprior=NULL,coefpriorvar=10,varprior=1,RandomWalk=TRUE){
  T <- nrow(y)
  K <- ncol(y)
  obs <- T-max(lags,thMax)
  constant=0
  if(intercept==TRUE){constant=1}
  if(is.null(coefprior)){
    Aprior1 <- array(0,dim=c(K*lags+constant,K))
    Aprior2 <- array(0,dim=c(K*lags+constant,K))
    if(RandomWalk==TRUE){
      if(intercept==FALSE){
        for(ii in 1:K){
          Aprior1[ii,ii]=1
          Aprior2[ii,ii]=1
        }
      }
    }
    else if(intercept==TRUE){
      for(ii in 1:K){
        Aprior1[ii+1,ii]<-1
        Aprior2[ii+1,ii]<-1
      }
    }
    aprior1 <- matrix(Aprior1,ncol=1)
    aprior2 <- matrix(Aprior2,ncol=1)
  } else {
    Aprior1=coefprior
    Aprior2=coefprior
  }

  if(.isscalar(coefpriorvar)){
    n<-K*(K*lags+constant)
    Vprior1 <- coefpriorvar*.id(n)
    Vprior2 <- coefpriorvar*.id(n)
  } else {
    Vprior1 <- coefpriorvar
    Vprior2 <- coefpriorvar
  }
  if(.isscalar(varprior)){
    Sprior1 <- varprior*.id(K)
    Sprior2 <- varprior*.id(K)
  } else {
    Sprior1 <- varprior
    Sprior2 <- varprior
  }
  return(list(Aprior1=Aprior1,Aprior2=Aprior2,Vprior1=Vprior1,Vprior2=Vprior2,Sprior1=Sprior1,Sprior2=Sprior2))
}