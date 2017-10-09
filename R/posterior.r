#################################################################################
#
# This R-file contains routines to draw posteriors for a given prior. Available
# priors are_
# 
# 1 - independent normal-wishart prior 
# 2 - Minnesota prior
# 3 - Natural Conjugate prior
# 4 - Uninformative prior 
# 
# Currently work in progress are
# 5 - dummy observations
# 6 - Stochastic search variable selection
#
#
# The general form of the function calls for priors 1-4 are 
#
# y - data
# x - lagged data
# aprior - prior on VAR-coefficients
# Vprior - prior on variance of VAR-coefficients
# vprior - prior on the shape of Variance
# Sprior - prior on scale of Variance
#
# stabletest - whether to accept only draws for the VAR-coefficients that lie outside the unit-root (TRUE/FALSE)
# Intercept  - whether the estimate includes an intercept or not (TRUE/FALSE)
# NoLags     - Number of lags of the VAR-model
#
#
#################################################################################


#
# Internal function to draw from a posterior distribution using an
# independent Normal-Wishart Prior. Parameters include:
#
postni <- function(y,x,aprior,Vprior,vprior,Sprior,Sigma,stabletest=FALSE,Intercept=TRUE,NoLags=2){
  # declare variables
  T <- nrow(y)
  K <- ncol(y)
  z <- diag(K)%x%x
  
  stable <- 2
  ii <- 1
  variance <- solve(Sigma)%x%diag(T)
  Vpost <- solve(solve(Vprior)+t(z)%*%variance%*%z)
  apost <- Vpost%*%(solve(Vprior)%*%aprior+t(z)%*%variance%*%as.vector(y))
  while(stable>1){
	
	#alpha <- mvrnorm(mu=apost,Sigma=Vpost)
	alpha <- apost+t(chol(Vpost))%*%rnorm(n=ncol(Vpost))
	#print(alpha)
	#readline(prompt="Press [enter] to continue")
	Alpha <- matrix(alpha,ncol=K)
	
	# Check if coefficients are good
	if(stabletest==TRUE){
		if(Intercept==TRUE){
			Alphatest <- Alpha[2:(K*NoLags+1),]
		}
		else{
			Alphatest <- Alpha
		}
		stable <- stability(betadraw = Alphatest,lags = NoLags)
	}
	else{
		stable <- 0
	}
  
  }
  # draw variance
  vpost <- vprior+T
  res <- y-x%*%Alpha
  Spost <- solve(Sprior+t(res)%*%res)
  Sigma <- solve(rWishart(1,vpost,Spost)[,,1])
  #print(Sigma)
  #readline(prompt="Press [enter] to continue")

  return(list(Alpha=Alpha,Sigma=Sigma))
}


#
# Function to draw from a posterior distribution using a Minnesota Prior 
#
# This function does not include prior parameters for the Variance-Covariance matrix.
#
postmb <- function(y,x,aprior,Vprior,Sigma,betaols,stabletest=FALSE,Intercept=TRUE,NoLags=2){
  
  # Declare variables
  K <- ncol(y)
  obs <- nrow(y)

  # Draw coefficients
  Vpost <- solve(solve(Vprior)+solve(Sigma)%x%(t(x)%*%x))
  apost <- Vpost%*%(solve(Vprior)%*%aprior+(solve(Sigma)%x%(t(x)%*%x))%*%betaols)
  
  stable <- 2
  while(stable>1){
    alpha <- mvrnorm(mu=apost,Sigma=Vpost)
    Alpha <- matrix(alpha,ncol=K)
    if(stabletest==TRUE){
	    if(Intercept==TRUE){
			Alphatest <- Alpha[2:(K*NoLags+1),]
		}
		else{
			Alphatest <- Alpha
		}
		stable <- stability(betadraw = Alphatest,lags = NoLags)
  
    }
	else{
	  stable <- 0
	}
  }
  
  
  # Draw Sigmas
  e <- y-x%*%Alpha
  scale <- t(e)%*%e/obs
  Sigma <- solve(rWishart(1,obs,scale)[,,1])
  
  # return values in a list
  return(list(Alpha=Alpha,Sigma=Sigma))

                                 
}

#
# Function to draw from a posterior distribution using an uninformative prior
#
postnc <- function(y,x,aprior,Vprior,vprior,Sprior,Sigma,stabletest=FALSE,Intercept=TRUE,NoLags=2){

  # Declare Variables
  K <- ncol(y)
  obs <- nrow(y)
  betaols <- solve(t(x)%*%x)%*%t(x)%*%y
  # Calculate posterior for coefficients
  Vpost <- solve(solve(Vprior)+(t(x)%*%x))
  Apost <- Vpost%*%(solve(Vprior)%*%aprior+t(x)%*%x%*%betaols)
  apost <- as.vector(Apost)
  
  # Calculate posterior for Variance-Covariance
  SSE <- t(y-x%*%betaols)%*%(y-x%*%betaols)
  Spost <- SSE+Sprior+t(betaols)%*%t(x)%*%x%*%betaols+t(aprior)%*%solve(Vprior)%*%aprior-t(Apost)%*%(solve(Vprior)+t(x)%*%x)%*%Apost
  vpost <- obs+vprior
  cova   <- Sigma%x%Vpost

  # Draw coefficients and variance-covariance
  stable <-2
  while(stable>1){
  
    alpha <- mvrnorm(mu=apost,Sigma=cova)
	Alpha <- matrix(alpha,ncol=K)
	
	# Check if coefficients are good
	if(stabletest==TRUE){
		if(Intercept==TRUE){
			Alphatest <- Alpha[2:(K*NoLags+1),]
		}
		else{
			Alphatest <- Alpha
		}
		stable <- stability(betadraw = Alphatest,lags = NoLags)
	}
	else{
		stable <- 0
	}
 }
 Sigma <- solve(rWishart(1,obs+vprior,solve(Spost))[,,1])
  
  
  
  return(list(Alpha=Alpha,Sigma=Sigma))
  
}

postun <- function(y,x,Sigma,stabletest=FALSE,Intercept=TRUE,NoLags=2){

  # uninformative prior 
  K <- ncol(y)
  obs <- nrow(y)
  
  # Calculate posterior
  Vpost   <- Sigma%x%solve(t(x)%*%x)
  betaols <- solve(t(x)%*%x)%*%t(x)%*%y
  sse     <- t(y-x%*%betaols)%*%(y-x%*%betaols)
  # Draw posterior
  stable <-2
  while(stable>1){

    alpha <- mvrnorm(mu=as.vector(betaols),Sigma=Vpost)
	Alpha <- matrix(alpha,ncol=K)
	
	# Check if coefficients are good
	if(stabletest==TRUE){
		if(Intercept==TRUE){
			Alphatest <- Alpha[2:(K*NoLags+1),]
		}
		else{
			Alphatest <- Alpha
		}
		stable <- stability(betadraw = Alphatest,lags = NoLags)
		
	}
	else{
		stable <- 0
	}
 }
 Sigma <- solve(rWishart(1,obs,solve(sse))[,,1])
 return(list(Alpha=Alpha,Sigma=Sigma))
  
}



stability <- function(betadraw,lags){
  K <- ncol(betadraw)
  companion <- companionmatrix(betadraw,lags)
  ev <- eigen(companion)$values
  maxev <- max(abs(ev))
  return(maxev)
}
# Function to create the Companion Matrix
companionmatrix <- function(Beta,lags){
  K <- ncol(Beta)
  #if(lags>1){
  #  xx <- cbind(diag(1,K*(lags-1)),matrix(0,K*(lags-1),K))
  #  PHI_Mat <- rbind(t(Beta),xx)
  #}
  #else{
  #  PHI_Mat <- Beta
  #}
  
  companion <- array(0,dim=c(K*lags,K*lags))
  Ai <- array(0,dim=c(K,K,lags))
  if(lags>1){
    for(jj in 1:lags){
      indxmin <- (jj-1)*K+1
      indxmax <- (jj-1)*K+K
      Ai[,,jj]<-Beta[indxmin:indxmax,]
      companion[1:K,indxmin:indxmax] <- Beta[indxmin:indxmax,]
    }
    indxmin <- K+1
    indxmax <- K*lags
    indxrep <- indxmax-indxmin
    for(jj in 0:indxrep){
      companion[indxmin+jj,jj+1]<-1
    }
  }
  else{
    companion <- Beta
  }
  
  return(companion)
}