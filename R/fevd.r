fevd.bvar <- function(bvarObj,h=10){

  nfevd       <- dim(bvarObj$betadraws)[3]
  fe <- 0

  for(ii in 1:nfevd){

    A     <- bvarObj$betadraws[,,ii]
    Sigma <- bvarObj$sigmadraws[,,ii]
    cholSigma <- t(chol(Sigma))
    dimSigma <- dim(Sigma)

    #
    # Create Companion matrix
    #

    # if intercept, remove coefficients for intercept
    if(bvarObj$intercept==TRUE){

      A <- A[-c(1),]

    }

    bigA <- companionmatrix(A,bvarObj$NoLags)
    bigJ <- array(0,dim=c(dim(bigA)[1],dimSigma[1]))
    bigJ[1:dimSigma[1],1:dimSigma[1]] <- diag(dimSigma[1])
    smallJ <- diag(dimSigma[1])

    # Calculate Mean-Squared prediction error

    mPsi <- array(0,dim=c(dimSigma[1],dimSigma[1],h))
    mPhi <- array(0,dim=c(dimSigma[1],dimSigma[1],h))
    for(jj in 1:h){

      mPsi[,,jj] <- t(bigJ) %*% mpower(bigA,jj) %*% bigJ
      mPhi[,,jj] <- mPsi[,,jj] %*% cholSigma

    }

    TH <- 0
    for(jj in 1:h){

      TH <- TH + (mPhi[,,jj] * mPhi[,,jj])

    }

    TH4 <- apply(TH,2,sum)

    # Calculate the forecast error variance

    fe1 <- array(0,dim=c(dimSigma[1],dimSigma[1]))
    for(jj in 1:dimSigma[1]){

      fe1[jj,] <- TH[jj,] / TH4

    }

    fe <- fe + fe1

  }

  return(fe/nfevd)


}

mpower <- function(x,n){

  mp <- x
  if(n > 1){
    for(ii in 1:(n-1)){

      mp <- mp %*% x

    }

  }
  return(mp)
}
