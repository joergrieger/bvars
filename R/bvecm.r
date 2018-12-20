detrend <- function(y,Intercept){
  if(Intercept == FALSE){
    resid <- y
  }
  else{
    nobs <- nrow(y)
    x     <- matrix(1,nrow=nobs,ncol=1)
    beta  <- solve(t(x)%*%x)%*%t(x)%*%y
    resid <- y - x%*%beta
  }

  return(resid)
}


johansen <- function(x,Intercept,NoLags=1){

  x   <- detrend(x,Intercept)
  dx  <- diff(x)
  z   <- lagdata(dx,NoLags)$x
  z   <- detrend(z,Intercept)
  dx  <- dx[-(1:NoLags),]
  dx  <- detrend(dx,Intercept)
  r0t <- dx - z %*% mldivide(z,dx)

  dx  <- lagdata(x,NoLags,intercept = FALSE)$x
  dx  <- detrend(dx[-c(1),],Intercept)
  rkt <- dx - z %*% mldivide(z,dx)

  skk <- t(rkt) %*% rkt / nrow(rkt)
  sk0 <- t(rkt) %*% r0t / nrow(r0t)
  s00 <- t(r0t) %*% r0t / nrow(r0t)
  sig <- sk0 %*% ginv(s00) %*% t(sk0)
  tmp <- ginv(skk)
  orig <- tmp %*% sig
  etmp <- eigen(orig)
  du   <- etmp$vectors
  au   <- etmp$values

  dt   <- du %*% solve(t(chol(t(du) %*% skk %*% du)))
  temp <-  solve(t(chol(t(du) %*% skk %*% du)))
  return(dt)
}

bvecm <- function(y,NoLags=1,Intercept=FALSE, r = 1,nreps = 200,burnin = 100,prior=4){

  # Check Input

  # Create cointegration vectors
  ecvec <- johansen(y,Intercept,NoLags)
  lagy  <- lagdata(y,NoLags,F)$x
  cointvec <- lagy %*% ecvec
  nr <- nrow(cointvec)
  nrend <- nr
  nrstart <- nr - NoLags+1
  cointvec <- cointvec[-c(nrstart:nrend),c(1:r)]
  dy <- diff(y)
  nr <- nrow(dy)
  nrend   <- nr
  nrstart <- nr - NoLags + 1
  dx <- dy[-c(nrstart:nrend),]
  dy <- dy[-c(1:NoLags),]
  dx <- cbind(dx,cointvec)

  # Declare variables for storage

  K <- ncol(dy)
  betadraws  <- array(0,dim=c(K*NoLags+r,K,nreps - burnin))
  sigmadraws <- array(0,dim=c(K,K,nreps - burnin))
  varnames   <- colnames(y)

  ###############################
  #
  # Create prior
  #
  ###############################

  if(prior == 1){
    # Independent Normal-Wishart Prior
    stop("Not implemented yet")
  }
  else if(prior == 2){
    stop("Not implemented yet")
  }
  else if(prior == 3){
    # Natural conjugate prior
    stop("Not implemented yet")
  }
  else if(prior == 4){
    # Uninformative prior
  }

  ##################################
  #
  # Initialize the MCMC algorithm
  #
  ##################################

  Aols <- solve(t(dx) %*% dx) %*% t(dx) %*% dy
  aols <- as.vector(Aols)
  resi <- dy - dx %*% Aols
  SSE  <- t(resi) %*% resi
  Sigma <- SSE/nrow(resi)

  ################################
  #
  # Start Gibbs Sampling
  #
  ################################
  for(irep in 1:nreps){
    if(prior == 1){
      stop("Not implemented")
    }
    else if(prior == 2){
      stop("Not implemented")
    }
    else if(prior == 3){
      stop("Not implemented")
    }
    else if(prior == 4){
      postdraw <- postun(y=dy,x=dx,Sigma=Sigma,stabletest=FALSE,Intercept=Intercept,NoLags=NoLags)
      Alpha <- postdraw$Alpha
      Sigma <- postdraw$Sigma
      print(Sigma)
    }

    ###############################
    #
    # Store Results
    #
    ###############################
    if(irep > burnin){
      betadraws[,,irep-burnin]  <- Alpha
      sigmadraws[,,irep-burnin] <- Sigma

    }
  }

}
