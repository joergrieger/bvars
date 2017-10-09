favar <- function(ydata,xdata,NoFactors=2,NoLags=2,RandomWalk=TRUE,prior=1,priorparam,alpha=0.01,beta=0.01,lipr=4,slowindex="",frotvar=NaN,nreps=20,burnin=10,irfhor=16,irfquantiles=c(0.05,0.95),stabletest=TRUE){
  
  #
  # Declare Variables and demean data
  #
  
  xdata <- as.matrix(xdata)
  ydata <- as.matrix(ydata)
  x.demeaned <- demean(xdata)
  y.demeaned <- demean(ydata)
  x <- x.demeaned
  y <- y.demeaned
  

  
  T <- nrow(y)
  N <- ncol(x)
  K <- ncol(y)
  P <- K+NoFactors
  
  irfSmallDraws <- array(0,dim=c(P,P,irfhor,NoRegimes,nreps-burnin))
  irfLargeDraws <- array(0,dim=c(P,(ncol(x)+ncol(y)),irfhor,NoRegimes,nreps-burnin))
  
  for(ii in 1:N){
    x[,ii] <- x[,ii]/sqrt(var(x[,ii]))
  }
  
  #
  #  Test for errors in input variables
  # 

  #
  # extract factors using pricing components
  #
  
  fac <- exfact(ydata=y,xdata=x,slowcode=slowindex,NoFactors=NoFactors)
  
  # put it in state-space form
  
  XY <- cbind(y,x)
  FY <- cbind(fac,y)
  L <- olssvd(XY,FY)
  
  e <- XY - FY%*%L
  Sigma <- t(e)%*%e/T
  
  # VAR equations
  
  FYlagged <- lagdata(FY,NoLags)
  FY.x <- FYlagged$x
  FY.y <- FYlagged$y
  Alpha <- solve(t(FY.x)%*%FY.x)%*%t(FY.x)%*%FY.y
  aols <- as.vector(Alpha)
  SSE <- t(FY.y-FY.x%*%Alpha)%*%(FY.y-FY.x%*%Alpha)
  SF <- SSE/(T-NoLags)
  
  
  #
  # Set up prior
  #
  
  # Prior on observation equation
  Liprvar <- lipr*diag(1,P)
  
  # Prior for VAR-Model
  if(prior==1){
    if(isempty(priorparam)){
	  stop("No prior parameters for Independent Normal-Wishart prior")
	}
	coefprior    <- priorparam[[1]]
	coefpriorvar <- priorparam[[2]]
	varprior     <- priorparam[[3]]
	varpriordof  <- priorparam[[4]]
    # Independent Normal-Wishart prior
    pr <- niprior(P,NoLags = NoLags,RandomWalk = RandomWalk,Intercept = FALSE,coefprior = coefprior,coefpriorvar = coefpriorvar,varprior = varprior,varpriordof = varpriordof)
    
    aprior <- as.vector(pr$coefprior)
    Vprior <- pr$coefpriorvar
    vprior <- varpriordof
    Sprior <- pr$varprior
  }
  else if(prior==2){
    if(isempty(priorparam)){
	  stop("No prior parameters for Independent Normal-Wishart prior")
	}
	else{
	  lambda1 <- priorparam[[1]]
	  lambda2 <- priorparam[[2]]
	  lambda3 <- priorparam[[3]]
	}
  
    # Minnesota prior
    pr <- mbprior(y=FY,NoLags = NoLags, intercept=FALSE,RandomWalk = RandomWalk)
    aprior <- pr$aprior
    Vprior <- pr$Vmatrix
  }
  else if(prior==3){
    if(isempty(priorparam)){
	  stop("No prior parameters for Natural conjugate prior")
	}
	coefprior    <- priorparam[[1]]
	coefpriorvar <- priorparam[[2]]
	varprior     <- priorparam[[3]]
	varpriordof  <- priorparam[[4]]
    # Natural conjugate prior 
	pr <- ncprior(P,NoLags = NoLags,RandomWalk = RandomWalk,Intercept = FALSE,coefprior = coefprior,coefpriorvar = coefpriorvar,varprior = varprior,varpriordof = varpriordof)
	aprior <- as.vector(pr$coefprior)
    Vprior <- pr$coefpriorvar
    vprior <- varpriordof
    Sprior <- pr$varprior
  }
  else if(prior==4){
    # uninformative prior, do nothing 
  }
  
  
  #
  # Start MCMC-Algorithm
  #
  L <- t(L)
  for(irep in 1:nreps){
    print(irep)
    
    # Step 1: Sample L and Sigma
    for(ii in 1:N){
      if(i>K){
        Li_postvar <- solve(solve(Liprvar)+Sigma[ii,ii]^(-1)*t(FY)%*%FY)
        Li_postmean <- Li_postvar%*%(Sigma[ii,ii]^(-1)*t(FY)%*%x[,ii])
        L[ii,1:P] <- t(Li_postmean)+rnorm(P)%*%chol(Li_postvar)
      }
      resi <- x[,ii]-FY%*%L[ii,]
      sh <- alpha/2+T/2
      sc <- beta/2+t(resi)%*%resi
      Sigma[ii,ii] <- rgamma(1,shape=sh,scale=sc)
    }
    
    # Step 2: Sample VAR-coefficients
    if(prior==1){
      postdraw <- postni(y=FY.y,x=FY.x,aprior,Vprior,vprior,Sprior,Sigma=SF,Intercept=FALSE,stabletest=stabletest)
      Alpha <- postdraw$Alpha
      SF <- postdraw$Sigma
    }
    else if(prior==2){
      postdraw <- postmb(y=FY.y,x=FY.x,Vprior=Vprior,aprior=aprior,Sigma=SF,betaols=aols)
      Alpha <- postdraw$Alpha
      SF <- postdraw$Sigma
    }
	else if(prior==3){
	  postdraw <- postnc(y=FY.y,x=FY.x,aprior,Vprior,vprior,Sprior,Sigma=SF,Intercept=FALSE,stabletest=stabletest)
	  Alpha <- postdraw$Alpha
      SF <- postdraw$Sigma
	}
	else if(prior==4){
	  postdraw <- postun(y=FY.y,x=FY.x,Sigma=SF,Intercept=FALSE,stabletest=stabletest)
	  Alpha <- postdraw$Alpha
      SF <- postdraw$Sigma
	}
    
    # compute impulse-response functions and save draws
    if(irep>burnin){
     irf <- compirf(A=Alpha,Sigma=SF,NoLags=NoLags,intercept=FALSE,nhor=irfhor)
     irfSmallDraws[,,,irep-burnin] <- irf
     for(ii in 1:P){
          irfLargeDraws[ii,,,irep-burnin] <- L[,,ireg]%*%irf[ii,,]
      }	 
    }
  }
  # Final computations
  irfSmallFinal <- array(0,dim=c(P,P,irfhor,3))
  irfLargeFinal <- array(0,dim=c(P,(ncol(x)+ncol(y)),irfhor,3))
  irflower <- min(irfquantiles)
  irfupper <- max(irfquantiles)
  for(jj in 1:P){
    for(kk in 1:P){
      for(ll in 1:irfhor){
        irfSmallfinal[jj,kk,ll,ii,1] <- mean(irfSmalldraws[jj,kk,ll,ii,])
        irfSmallfinal[jj,kk,ll,ii,2] <- quantile(irfSmalldraws[jj,kk,ll,ii,],probs=irflower)
        irfSmallfinal[jj,kk,ll,ii,3] <- quantile(irfSmalldraws[jj,kk,ll,ii,],probs=irfupper)
      }
    }
  }
  for(jj in 1:P){
    for(kk in 1:(ncol(x)+ncol(y))){
	  for(ll in 1:irfhor){
	    irfLargefinal[jj,kk,ll,ii,1] <- mean(irfLargedraws[jj,kk,ll,ii,])
		irfLargefinal[jj,kk,ll,ii,2] <- quantile(irfLargedraws[jj,kk,ll,ii,],probs=irflower)
		irfLargefinal[jj,kk,ll,ii,3] <- quantile(irfLargedraws[jj,kk,ll,ii,],probs=irfupper)
	  }
	}
  }
  
  
  
  
}

#
# Function to extract factors using principal components
# and rotate them
#

exfact <- function(ydata,xdata,slowcode="",NoFactors,frotvar=NaN){
  N <- ncol(xdata)
  
  print(dim(xdata))
  #
  # extract factors
  #
  
  factors <- extract(xdata,NoFactors)
  #
  # rotate factors as in Boivin et al (2005)
  #
  print("rotating data")
  if(!isempty(slowcode)){
    slow <- which(slowcode==1)
    xslow <- xdata[,slow]
    ftemp <- extract(xslow,NoFactors)
    fnr <- ncol(ftemp$fac)
    if(is.nan(frotvar)){
      frotvar=ncol(ydata)
    }
    Fr0 <- facrot(factors$fac,ydata[,frotvar],ftemp$fac)
  }
  else{
    Fr0 <- factors$fac
  }
  return(Fr0)
}

#
# extract factors using principal components
#

extract <- function(x,K){
  n <- ncol(x)
  x <- as.matrix(x)
  x.x <- t(x)%*%x
  evectors <- eigen(x.x)$vectors
  ret.evectors <- sqrt(n)*evectors[,1:K]
  fac <- x%*%ret.evectors/n
  return(list(lam=ret.evectors,fac=fac))
}

#
# Rotate factors
#

facrot <- function(F,Ffast,Fslow){
  
  Ffast <- as.matrix(Ffast)
  k1    <- ncol(Ffast)
  fm    <- cbind(matrix(1,nrow(Ffast),1), Ffast, Fslow)
  b     <- olssvd(F,fm)
  Fr <- F-Ffast%*%b[2:(k1+1),]
  return(Fr)
}

olssvd <- function(y,ly){
  duv          <- svd(t(ly)%*%ly)
  x.inv        <- duv$v%*%diag(1/duv$d)%*%t(duv$u)
  x.pseudo.inv <- x.inv %*% t(ly)
  b <- x.pseudo.inv%*%y
  return(b)
}
