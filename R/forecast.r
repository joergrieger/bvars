.forecast <- function(mydata,lags,beta,sigma,intercept,forecasthorizon){
  T <- nrow(mydata)
  K <- ncol(mydata)
  forecastdata <- matrix(0,ncol=(forecasthorizon+lags),nrow=K)
  if(intercept==TRUE){
    beta1 <- beta[2:(K*lags+1),]
    beta1_int <- beta[1,]
  }
  else{
    beta1 <- beta
    beta1_int <- numeric(K)
  }
  forecastdata[,1:lags] <- t(mydata[T:(T-lags+1),])
  mu=numeric(K)
  for(ii in 1:forecasthorizon){
	indxmax <- lags+ii-1
    indxmin <- ii
	ytmp <- forecastdata[,indxmax:indxmin]
	ytmp <- matrix(ytmp,nrow=1)
	tmp1 <- mvrnorm(mu=mu,Sigma=sigma)
	tmp2 <- beta1_int+ytmp%*%beta1+tmp1
	forecastdata[,ii+lags] <- tmp2
    
  }
  return(forecastdata[,(lags+1):(forecasthorizon+lags)])
}
