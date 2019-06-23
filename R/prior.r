# Set priors for Independent Normal-Wishart Prior
niprior <- function(K,NoLags,Intercept=TRUE,RandomWalk=TRUE,coefprior,coefpriorvar,varprior,varpriordof){

  #
  # Prior for coefficients
  #
  constant <- 0
  if(Intercept==TRUE) constant=1

  if(is.null(coefprior)){

    coefprior <- array(0,dim=c(K*NoLags+constant,K))

    if(RandomWalk==TRUE){

      coefprior[(1+constant):(K+constant),1:K] <- diag(1,K)

    }
  }
  else if(.isscalar(coefprior)){

    coefprior <- coefprior * array(1,dim=c(K*NoLags+constant,K))

  }

  if(.isscalar(coefpriorvar)){

    coefpriorvar <- coefpriorvar*diag(1,(K*(K*NoLags+constant)))

  }

  #
  # Prior on variance
  #
  if(.isscalar(varprior)){
    varprior <- diag(varprior,K)
  }
  return(list(coefprior=coefprior,coefpriorvar=coefpriorvar,varprior=varprior))
}

# Set parameters for Minnesota prior

# Input:
#
# y - Txm Matrix
# NoLags - Number of lags
# Intercept = T/F whether to include an intercept or not
# Random Walk - T/F or numeric
#
mbprior <- function(y,NoLags,Intercept=TRUE,RandomWalk=TRUE,lambda1=1,lambda2=1,lambda3=1,lambda4=2){
  y <- as.matrix(y)

  # Declare variables
  obs <- nrow(y)
  K   <- ncol(y)
  constant = 0

  if(Intercept==TRUE) constant=1

  #
  # Prior for coefficients
  #

  Aprior <- array(0,dim=c(K*NoLags+constant,K))

  if(is.numeric(RandomWalk)){

    rw = RandomWalk
    RandomWalk = TRUE

  }
  else{
    rw = 1
  }


  if(RandomWalk==TRUE){

    for(ii in 1:K){

      Aprior[(ii+constant),ii] <- rw

    }

  }
  aprior <- as.vector(Aprior)

  #
  # Prior for covariance matrix
  #
  sigmasq <- array(0,dim=c(K,1))

  for(ii in 1:K){

    Ylagi         <- embed(y[,ii],dimension=NoLags+1)[,-1]
    Yi            <- y[(NoLags+1):obs,ii]
    arest         <- lm(Yi~Ylagi-1)
    sigmasq[ii,1] <- summary(arest)$sigma

  }
  print(sigmasq)

  #
  # Covariance matrix for the prior
  #
  M <- K*NoLags+constant
  Vi <- array(0,dim=c(K*M,1))

  # without intercept
  for(ii in 1:K){ # loop over the ii-th equation
    for(jj in 1:NoLags){ #loop over the jj-th lag
      for(kk in 1:K){ #kk-th variable
        indx <- (ii-1)*(K*NoLags)+(jj-1)*K+kk
        if(ii==kk){
          Vi[indx,1] <- (lambda1)/(jj^lambda4)
        }
        else{
          Vi[indx,1] <- (lambda1*lambda2)/(jj^lambda4)*(sigmasq[ii,1]/sigmasq[kk,1])^2
        }
      }
    }
  }

  #
  # Add Covariance coefficients for intercepts
  #

  if(Intercept==TRUE){
    Vtmp <- array(0,dim=c(K*K*NoLags+K,1))
    for(ii in 1:K){
      coefinter <- lambda3*sigmasq[ii,1]^2
      indx <- (ii-1)*(K*NoLags)+ii
      Vtmp[indx,1] <- coefinter
      indxmin <- (ii-1)*(K*NoLags)+1
      indxmax <- (ii-1)*(K*NoLags)+K*NoLags
      Vtmp[(indx+1):(indx+(K*NoLags)),] <- Vi[indxmin:indxmax,]
    }
    Vi <- Vtmp
  }

  #
  # Create diagonal matrix
  #
  nr <- dim(Vi)[1]
  Vfinal <- array(0,dim=c(nr,nr))
  for(ii in 1:nr){
    Vfinal[ii,ii] <- Vi[ii,1]
  }
  return(list(aprior=aprior,Vmatrix=Vfinal))
}

ncprior <- function(K,NoLags,Intercept=TRUE,RandomWalk=TRUE,coefprior,coefpriorvar,varprior,varpriordof){
  #
  # Prior for coefficients
  #
  constant <- 0
  if(Intercept==TRUE) constant=1

  # prior for coefficients
  if(is.null(coefprior)){

    coefprior <- array(0,dim=c(K*NoLags+constant,K))

    if(RandomWalk==TRUE){

      coefprior[(1+constant):(K+constant),1:K] <- diag(1,K)
    }
  }
  if(.isscalar(coefprior)){

    coefprior <- coefprior * array(1,dim=c(K*NoLags + constant, K))

  }
  if(.isscalar(coefpriorvar)){
    coefpriorvar <- coefpriorvar*diag(1,K*NoLags+constant)
  }

  # prior for Variance-Covariance Matrix
  if(.isscalar(varprior)){
    varprior <- diag(varprior,K)
  }
  return(list(coefprior=coefprior,coefpriorvar=coefpriorvar,varprior=varprior))

}

dummyPrior <- function(y,lamda,tau,epsilon,p){

  N <- ncol(y)
  yy <- y[2:nrow(y),]
  xx <- y[1:(nrow(y)-1),]
  mu <- colMeans(y)
  sigma <- array(0,dim=c(ncol(yy)))
  delta <- array(0,dim=c(ncol(yy)))

  yy <- as.matrix(yy)
  xx <- as.matrix(xx)

  # Get AR(1) - coefficients
  for(ii in 1:N){

    x1 <- t(rbind(1,xx[,ii]))
    betaxy  <- solve(t(x1) %*% x1) %*% t(x1) %*% yy[,ii]
    resid   <- yy[,ii] - x1%*%betaxy
    if(abs(betaxy[2]) > 1) betaxy[2] = 1

    delta[ii] <- betaxy[2]
    sigma[ii] <- t(resid)%*%resid/nrow(yy)
  }
  if(lamda > 0){
    if(epsilon > 0){
      yd1 <- diag(sigma*delta)/lamda
      yd2 <- array(0,dim=c(N*(p-1),N))
      yd3 <- diag(sigma)
      yd <- rbind(yd1,yd2,yd3)

      jp <- diag(1:p)

      xxd1 <- jp%x%diag(sigma)/lamda
      xxd2 <- array(0,dim=c(N*p,1))
      xd1  <- cbind(xxd1,xxd2)

      xxd1 <- array(0,dim=c(N,(N * p) + 1))
      xxd2 <- cbind(array(0,dim=c(N,N*p)),epsilon)

      xd1  <- rbind(xd1,xxd1,xxd2)

    }
    else{
      yd1 <- diag(sigma*delta)/lamda

      yyd1 <- array(0,dim=c((N * (p - 1)), N))
      yyd2 <- diag(sigma)
      yd1 <- rbind(yd1,yyd1,yyd2)

      jp <- diag(1:p)
      xxd1 <- jp %x% diag(sigma)/lamda
      xxd2 <- array(0,dim=c(N, (N * p)))
      xd1 <- rbind(xxd1,xxd2)
    }
  }

  if(tau > 0){
    if(epsilon > 0){
      yd2  <- diag(delta * mu)/tau
      xd2  <- cbind(array(1,dim = c(1,p)) %x% yd2,0)
    }
    else{
      yd2 <- diag(delta * mu)/tau
      xd2 <- array(1,dim = c(1,p)) %x% yd2
    }
  }

  if(lamda > 0 && tau > 0){
    xd <- rbind(xd1,xd2)
    yd <- rbind(yd1,yd2)
  }
  else if(lamda > 0 && !(tau > 0 )){
    xd <- xd1
    yd <- yd1
  }
  else if( !(lamda > 0) && (tau > 0)){
    xd <- xd2
    yd <- yd2
  }
  else if(!(lamda > 0) && !(tau > 0)){
    xd <- NULL
    yd <- NULL
  }
  retlist <- list(xd=xd,yd=yd)
  return(retlist)

}
