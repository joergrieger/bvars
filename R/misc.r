
#
#
# Function to create lagged data
#
#

lagdata <-function(y,lags,intercept=FALSE){
  T<-nrow(y)
  K<-ncol(y)
  obs <- T-lags
  x  <- embed(y,dimension=lags+1)[,-(1:K)]
  if(intercept==TRUE){
    x<-cbind(1,x)
  }
  yi <- y[(lags+1):T,]
  return(list(y=yi,x=x,obs=obs,T=T,K=K));
}


# function to test if a provided variable is a scalar
.isscalar <- function(x) is.atomic(x) && length(x) == 1L && !is.character(x) && Im(x)==0

vec2matrix <- function(lags,K,vec){
  matrix0 <- array(0,dim=c(lags,K))
  for(ii in 1:lags){
    for(jj in 1:K){
      matrix0[ii,jj] <- vec[(ii-1)*K+jj]
    }
  }
  return(matrix0)
}

invpd <- function(x){
  as.matrix(x)
  xncol <- ncol(x)
  temp <- diag(1,xncol)
  ipd <- mldivide(x,temp,pinv=TRUE)
  return(ipd)
}

demean <- function(x){
  x.nc <- ncol(x)
  for(i in 1:x.nc){
    x.mean <- mean(x[,i])
    x[,i]<-x[,i]-x.mean
  }
  return(x)
}


isscalar <- function(x) is.atomic(x) && length(x) == 1L && !is.character(x) && Im(x)==0

