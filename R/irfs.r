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
