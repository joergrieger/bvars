
# vectorizes a matrixs
.vec <- function(Mat){
  yret <- matrix(Mat,ncol=1)
  return(yret)
}
# creates an NxN identity matrix
.id <-function(N){
  return(diag(rep(1,N)))
}

# Checks if VAR is stable
.stability <- function(betadraw,lags,K){
  companionmatrix <- array(0,dim=c(K*lags,K*lags))
  Ai <- array(0,dim=c(K,K,lags))
  if(lags>1){
    for(jj in 1:lags){
      indxmin <- (jj-1)*K+1
      indxmax <- (jj-1)*K+K
      Ai[,,jj]<-betadraw[indxmin:indxmax,]
      companionmatrix[1:K,indxmin:indxmax]<-betadraw[indxmin:indxmax,]
    }
    indxmin <- K+1
    indxmax <- K*lags
    indxrep <- indxmax-indxmin
    for(jj in 0:indxrep){
      companionmatrix[indxmin+jj,jj+1]<-1
    }
    ev <- eigen(companionmatrix)$values
    ev <- Mod(ev)
    maxev <- max(abs(ev))
  }
  else{
    ev <- eigen(betadraw)$values
    ev <- Mod(ev)
    maxev <- max(abs(ev))
  }
  return(maxev)
}
# function to test if a provided variable is a scalar
.isscalar <- function(x) is.atomic(x) && length(x) == 1L && !is.character(x) && Im(x)==0
