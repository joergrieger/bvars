fevd.bvar <- function(bvObj,h=12){

  intercept <- bvObj$intercept
  betadraws <- bvObj$betadraws
  sigmadraws <- bvObj$sigmadraws
  NoLags <- bvObj$NoLags

  nreps <- dim(betadraws)[3]
  k <- dim(sigmadraws)[1]
  fe1 <- array(0,dim=c(k,k,nreps))

  for(ii in 1:nreps){
    beta  <- betadraws[,,ii]
    sigma <- sigmadraws[,,ii]
    if(intercept == TRUE){
      beta <- beta[-c(1),]
    }
    bigBeta <- companionmatrix(beta,NoLags)
    fe <- fevd1(bigBeta,sigma,h)
    fe1[,,ii] <- fe

  }
  decomposition <- array(NA,dim=c(k,k,4))
  for(ii in 1:k){
    for(jj in 1:k){

      decomposition[ii,jj,1] <- mean(fe1[ii,jj,])
      decomposition[ii,jj,2] <- median(fe1[ii,jj,])
      decomposition[ii,jj,3] <- quantile(fe1[ii,jj,],probs = c(0.05))
      decomposition[ii,jj,4] <- quantile(fe1[ii,jj,],probs = c(0.95))

    }
  }

  return(decomposition)
}

fevd1 <- function(beta,sigma,h){

  k  <- dim(sigma)[1]
  kp <- dim(beta)[1]
  J <- array(0,dim=c(k,kp))
  cholSigma <- t(chol(sigma))
  J[1:k,1:k] <- diag(1,k)
  TH <- array(0,dim=c(k,k,h))
  VC <- 0

  for(ii in 1:h){

    TH[,,ii] <- t(J %*% ( beta %^% (ii) )%*% t(J) %*% cholSigma)
    VC <- VC + TH[,,ii] * TH[,,ii]

  }

  xtmp <- apply(VC,2,sum)

  for(ii in 1:k){

    VC[ii,] <- VC[ii,]/xtmp

  }
  return(VC)

}

"%^%"<-function(A,n){
  if(n==1) A
  else {
    B<-A
    for(i in (2:n)){
      A <- A%*%B
      }
    }
  A
}
