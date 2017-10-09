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
mbprior <- function(y,NoLags,Intercept=TRUE,RandomWalk=TRUE,lambda1=1,lambda2=1,lambda3=1){
  y <- as.matrix(y)
  
  # Declare variables
  obs <- nrow(y)
  K <- ncol(y)
  constant=0
  
  if(Intercept==TRUE) constant=1
  
  #
  # Prior for coefficients
  #
  
  Aprior <- array(0,dim=c(K*NoLags+constant,K))
  
  if(RandomWalk==TRUE){
    for(ii in 1:K){
      Aprior[(ii+constant),ii] <- 1
    }
  }
  aprior <- as.vector(Aprior)
  print(aprior)
  
  #
  # Prior for covariance matrix
  #
  sigmasq <- array(0,dim=c(K,1))
  for(ii in 1:K){
    Ylagi <- embed(y[,ii],dimension=NoLags+1)[,-1]
    Yi <- y[(NoLags+1):obs,ii]
    arest <- lm(Yi~Ylagi-1)
    sigmasq[ii,1] <- summary(arest)$sigma
  }

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
          Vi[indx,1] <- lambda1/(jj^2)
        }
        else{
          Vi[indx,1] <- lambda2/(jj^2)*sigmasq[ii,1]/sigmasq[kk,1]
        }
      }
    }
  }
  print(Vi)
  #
  # Add Covariance coefficients for intercepts
  #
  
  if(intercept==TRUE){
    Vtmp <- array(0,dim=c(K*K*NoLags+K,1))
    for(ii in 1:K){
      coefinter <- lambda3*sigmasq[ii,1]
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
  if(.isscalar(coefpriorvar)){
    coefpriorvar <- coefpriorvar*diag(1,K*NoLags+constant)
  }

  # prior for Variance-Covariance Matrix
  if(.isscalar(varprior)){
    varprior <- diag(varprior,K)
  }
  return(list(coefprior=coefprior,coefpriorvar=coefpriorvar,varprior=varprior))
  
}