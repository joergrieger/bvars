.hd <- function(y,lags,Beta,resid,Sigma){

  # Initialize Variables

  Fcomp <- companionmatrix(Beta)
  eps <- ginv(Sigma)%*%t(resid)
  obs <- nrow(Y)
  K <- ncol(Y)

  # Compute historical decomposition

  invABig <- matrix(0,K*lags,K)
  invABig[1:K,] <- Sigma

  Icomp <- cbind(diag(K),matrix(0,K,(lags-1)*K))
  HDshockBig <- array(0,dim=c(K*lags,(obs+1),K))
  HDshock <- array(0,dim=c(K,(obs+1),K))

  for(j in 1:K){
    epsBig <- matrix(0,K,(obs+1))
    epsBig[j,2:ncol(epsBig)] <- eps[j,]
    for(i in 2:(obs+1)){
      HDshock_big[,i,j] <- invABig %*% epsBig[,i]+Fcomp%*%HDshockBig[,(i-1),j]
      HDshock[,i,j] <- Icomp%*%HDshockBig[,i,j]
    }
  }
  for(i in 1:K){
    for(j in 1:K){
      HD.shock[,j,i] <- c(rep(NA,lags),HDshock[i,(2:dim(HDshock)[2]),j])
    }
  }
  return(HD.shock)
}
