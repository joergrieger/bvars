hd.bvar <- function(bvObj){

  # Preliminary work
  nreps <- dim(bvObj$betadraws)[3]
  nT <- dim(bvObj$mydata)[1] - bvObj$NoLags

  # lagdata
  lg <- lagdata(bvObj$mydata,intercept = bvObj$intercept,lags=bvObj$NoLags)

  # Declare variables to store historical decomposition
  HD <- 0

  for(ii in 1:nreps){

    # Create residuals
    resid <- lg$y - lg$x%*%bv$betadraws[,,ii]

    # Create covariance matrix
    covar <- t(chol(bvObj$sigmadraws[,,ii])) # Identification using Cholesky decomposition
    nVar <- dim(covar)[1]
    seqCovar <- array(0,dim=c(nVar,nVar,nT))
    seqCovar[,,] <- covar

    # Create companion matrix
    if(bvObj$intercept){

      bta <- bvObj$betadraws[-c(1),,ii]

    }
    else{

      bta <- bvObj$betadraws[,,ii]

    }

    comp <- companionmatrix(bta,lags=bvObj$NoLags)
    nDim <- dim(comp)[1]
    seqCompanion <- array(0,dim=c(nDim,nDim,nT))
    seqCompanion[,,] <- comp

    # Create Historical decomposition

    HD <- HD + HDOnePath(seqCompanion,resid,seqCovar)


  }

  HD <- HD/nreps
  if(is.ts(bvObj$mydata)){

    timestamps <- time(bvObj$mydata)
    timestamps <- timestamps[-c(1:bvObj$NoLags)]

  }
  else{
    timestamps = NULL
  }

  retlist <- structure(list(hd=HD,varnames=colnames(bvObj$mydata),date=timestamps),class="histdecomp")
  return(retlist)

}

# function to calculate the historical decomposition for one path
# F <- sequence of companion matrices
# U <- sequence of reduced form residuals
# C <- sequence of identified covariance matrices
HDOnePath <- function(Fr,U,C){

  # Create 'structural structural shocks'
  nT   <- dim(U)[1]
  nVar <- dim(U)[2]

  eps  <- array(0,dim=c(nT,nVar^2))
  HDtemp <- array(NaN,dim=c(nVar^2,nT))
  bigeye <- array(0,dim=c(nVar,dim(Fr)[1]))
  bigeye[1:nVar,1:nVar] <- diag(1,nVar)


  for(ii in 1:nT){

    eps[ii,] <-   as.vector(repmat(mldivide(C[,,ii],U[ii,]),nVar,1))

  }

  HDtemp[,1] <- as.vector(C[,,1]) * eps[1,]

  for(jj in 2:nT){
    Ftrack <- diag(1,dim(Fr)[1])
    HDtemp[,jj] <- as.vector(C[,,jj]) * eps[jj,]

    for(kk in (jj-1):1){

      Ftrack <- Ftrack %*% Fr[,,kk+1]
      tmp1 <- as.vector(bigeye %*% Ftrack %*% t(bigeye) %*% C[,,kk]) * eps[kk,]
      HDtemp[,jj] <- HDtemp[,jj] + tmp1
    }


  }

  HD <- array(NaN,dim=c(nVar,nVar,nT))
  for(ii in 1:nVar){

    repma <- seq(ii,nVar^2,by=nVar)
    HD[ii,,] <- HDtemp[c(repma),]

  }

  return(HD)

}

