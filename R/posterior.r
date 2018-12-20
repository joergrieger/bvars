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

postss <- function(y,x,SSEGibbs,SSE,omega,gammas,tau0,tau1,Sigma,aols,aprior,Intercept=TRUE,NoLags,stabletest=TRUE){
  # Preliminaries
  #print(omega)
  #readline(prompt="Press [enter] to continue")
  kappa0 <- 0.1
  kappa1 <- 10
  tau0 <- 0.1
  tau1 <- 10
  bi=0.01
  ai=0.01
  qij=0.5
  p_i=0.5

  K <- ncol(y)
  T <- nrow(y)
  NoRest <- ncol(y)*ncol(x)


  #
  # First: Draw Sigma
  #

  S <- array(list(),dim=c(K))
  for(kk2 in 1:K){
    S[[kk2]] <- SSEGibbs[1:kk2,1:kk2]
  }

  s <- array(list(),dim=(K-1))
  for(kk3 in 2:K){
    s[[kk3-1]] <- SSEGibbs[1:(kk3-1),kk3]
  }

  hh <- array(list(),dim=c(K-1))
  for(kk4 in 1:(K-1)){
    omeg <- omega[[kk4]]
	het <- hh[[kk4]]
	for(kkk in 1:dim(omeg)[1]){
	  if(omeg[kkk]==0){
	    het[kkk] <- kappa0
	  }
	  else{
	    het[kkk] <- kappa1
	  }
	}
	hh[[kk4]] <- het
  }

  Dj <- array(list(),dim=c(K-1))
  for(kk5 in 1:(K-1)){
    if(kk5==1){
	  Dj[[kk5]] <- diag(hh[[kk5]],1)
	}
	else{
	  Dj[[kk5]] <- diag(hh[[kk5]])
	}
  }
  DDj <- array(list(),dim=c(K-1))
  for(kk6 in 1:(K-1)){
    DD <- Dj[[kk6]]
	DDj[[kk6]] <- DD%*%DD
  }
  # Create B[i] matrix
  B <- array(list(),dim=c(K))
  for(rr in 1:K){
    if(rr==1){
	  B[[rr]] <- bi+0.5*SSEGibbs[rr,rr]
	}
	else if(rr>1){
	  si <- s[[rr-1]]
	  Si <- S[[rr-1]]
	  DiDi <- DDj[[rr-1]]
	  B[[rr]] <- bi+0.5*(SSEGibbs[rr,rr]-t(si)%*%solve(Si+solve(DiDi))%*%si)
	}
  }
  psi_ii_sq = array(0,dim=c(K))

  for(kk7 in 1:K){
    psi_ii_sq[kk7] <- rgamma(1,(ai+0.5*T),(B[[kk7]]))
  }

  #
  # Draw eta
  #

  eta <- array(list(),dim=c(K-1))
  for(kk8 in 1:(K-1)){
    si <- s[[kk8]]
	Si <- S[[kk8]]
	DiDi <- DDj[[kk8]]
	miuj <- -sqrt(psi_ii_sq[kk8+1])*(solve(Si+solve(DiDi))%*%si)
	Deltaj <- solve(Si+solve(DiDi))
	eta[[kk8]] <- miuj+t(chol(Deltaj))%*%rnorm(kk8)
  }

  #
  # Draw Omega
  #

  for(kk9 in 1:(K-1)){
    omegg <- omega[[kk9]]
	etag  <- eta[[kk9]]
    omegavec <- matrix(nrow=0,ncol=1)
	for(nn in 1:dim(omegg)[1]){
	  uij1 <- 1/(kappa0)*exp(-0.5*((etag[nn])^2)/(kappa0^2))*qij
	  uij2 <- 1/(kappa1)*exp(-0.5*((etag[nn])^2)/(kappa1^2))*(1-qij)
	  ost  <- max(0,uij1/(uij1+uij2))
	  ost  <- min(1,ost)
	  omegg <- rbinom(1,1,ost)
	  omegavec <- rbind(omegavec,omegg)
	}
	omega[[kk9]] <- omegavec
  }

  # Create PSI Matrix

  PSI_ALL <- array(0,dim=c(K,K))
  for(nn1 in 1:K){
    PSI_ALL[nn1,nn1] <- sqrt(psi_ii_sq[nn1])
  }
  for(nn2 in 1:(K-1)){
    etagg <- eta[[nn2]]
	for(nnn in 1:dim(etagg)[1]){
	  PSI_ALL[nnn,nn2+1] <- etagg[nnn]
	}
  }

  Sigma <- solve(PSI_ALL%*%t(PSI_ALL))

  #
  # End Drawing Sigma
  # Draw Alphas now
  #

  hi <- array(0,dim=c(NoRest))
  for(nn3 in 1:NoRest){

    if(gammas[nn3]==0){
	  hi[nn3] <- tau0#tau0[nn3]
	}
	else if(gammas[nn3]==1){
	  hi[nn3] <- tau1#tau1[nn3]
	}
  }

  D <- diag(as.vector(t(hi)%*%diag(1,NoRest)))
  #D <- diag(hi)

  DD <- D%*%D

  isig <- solve(Sigma)
  psixx <- isig%x%(t(x)%*%x)
  Vpost <- solve(psixx+solve(DD))
  apost <- Vpost%*%(psixx%*%aols+solve(DD)%*%aprior)
  stable <- 2
  while(stable>1){

    alpha <- mvrnorm(mu=apost,Sigma=Vpost)
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
  for(nn6 in 1:NoRest){
    ui1 <- 1/tau0*exp(-0.5*(alpha[nn6]/tau0)^2)*p_i
	ui2 <- 1/tau1*exp(-0.5*(alpha[nn6]/tau1)^2)*(1-p_i)

	gst <- min(1,ui1/(ui1+ui2))
	gst <- max(0,gst)
	if(ui1==Inf){gst=1}
	if(ui1+ui2<1e-7){gst=0.5}
	gammas[nn6] <- rbinom(1,1,gst)
  }

  SSEGibbs <- t(y-x%*%Alpha)%*%(y-x%*%Alpha)
  retlist <- list(Sigma=Sigma,Alpha=Alpha,gammas=gammas,omega=omega,SSEGibbs=SSEGibbs)
  return(retlist)


}


postdummy <- function(y,x,yd,xd,sigma,p,Intercept,stabletest){

  constant <- 0
  if(Intercept == TRUE){constant = 1}

  y0    <- rbind(y,yd)
  x0    <- rbind(x,xd)
  mstar <- as.vector(solve(t(x0) %*% x0) %*% t(x0) %*% y0)
  xx    <- t(x0) %*% x0
  ixx1  <- solve(xx,diag(1,ncol(xd)))

  vstar <- sigma %x% ixx1

  # Draw the beta
  stable <- 2
  while(stable > 1){
    Betadraw <- as.vector(mstar + rnorm(ncol(y0) * (ncol(y0) * p + constant)) %*% t(chol(vstar)))
    Alphatest <- matrix(Betadraw,ncol=ncol(y0))
    if(stabletest == TRUE){
      if(Intercept == TRUE){
        Alphatest <- Alphatest[-c(1),]
      }
      stable <- stability(betadraw = Alphatest,lags = p)
    }
    else{
      stable <- 0
    }
  }
  # Draw the Variance
  resid <- y0 - x0%*%matrix(Betadraw,ncol=ncol(y0))
  shape <- nrow(y0)
  scale <- t(resid)%*%resid/shape
  SigmaDraw <- solve(rWishart(1,shape,scale)[,,1])

  retlist <- list(Alpha = as.matrix(Betadraw,ncol=ncol(y0)),
                  Sigma = SigmaDraw)

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
