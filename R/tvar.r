tvar <- function(mydata,NoLags=1,thMax=2,thVar=1,Intercept=TRUE,RandomWalk=TRUE,prior=1,coefprior=NULL,coefpriorvar=1,varprior=5,varpriordof=10,nreps=110,burnin=10,irfhorizon=16,irfquantiles=c(0.05,0.95),bootrep=10){
  
  #
  # Declare variables and some preliminary calculations
  #
  
  y <- as.matrix(mydata)
  T <- nrow(y)
  K <- ncol(y)
  constant <- 0
  if(Intercept==TRUE) constant <- 1
  
  thDelay <- thMax
  tard    <- seq(1:thMax)
  startest <- max(thMax,NoLags)
  ytest <- y[(startest+1-thDelay):(T-thDelay),thVar]
  tarmean <- mean(ytest)
  tarstandard <- sqrt(var(ytest))

  Alpha <- array(0,dim=c(K*NoLags+constant,K,2))
  Sigma <- array(0,dim=c(K,K,2))
  
  #
  # Declare variables to store results
  #
  Alphadraws <- array(0,dim=c(K*NoLags+constant,K,2,nreps-burnin))
  Sigmadraws <- array(0,dim=c(K,K,2,nreps-burnin))
  Irfdraws   <- array(0,dim=c(K,K,irfhorizon,2,nreps-burnin))
  tardraws   <- array(0,dim=c(nreps-burnin))
  deldraws   <- array(0,dim=c(nreps-burnin))
  NoRegimes  <- T-(startest+1)
  regimes    <- array(0,dim=c(NoRegimes,nreps-burnin))
  irffinal <- array(0,dim=c(K,K,irfhorizon,2,3))
  
  #
  # Create priors
  #
  
  if(prior==1){
    # Independent Normal-Wishart Prior
    pr <- niprior(K=K,NoLags=NoLags,RandomWalk=RandomWalk,Intercept=Intercept,coefprior=coefprior,coefpriorvar=coefpriorvar,varprior=varprior)
    
    aprior <- as.vector(pr$coefprior)
    Vprior <- pr$coefpriorvar
    vprior <- varpriordof
    Sprior <- pr$varprior
  }
  
  
  #
  # Initialization of the Gibbs-Sampler
  #
  tart <- tarmean
  xsplit <- splitVariables(y=y,lags=NoLags,thDelay=thDelay,thresh=thVar,tart=tarmean,intercept=Intercept)
  Alpha[,,1] <- solve(t(xsplit$x1)%*%xsplit$x1)%*%t(xsplit$x1)%*%xsplit$y1
  Alpha[,,2] <- solve(t(xsplit$x2)%*%xsplit$x2)%*%t(xsplit$x2)%*%xsplit$y2
  Sigma[,,1] <- t(xsplit$y1-xsplit$x1%*%Alpha[,,1])%*%(xsplit$y1-xsplit$x1%*%Alpha[,,1])
  Sigma[,,2] <- t(xsplit$y2-xsplit$x2%*%Alpha[,,2])%*%(xsplit$y2-xsplit$x2%*%Alpha[,,2])
  

  #
  # Start sampling
  #
  
  for(irep in 1:nreps){
    print(irep)
    
    #
    # Step 1: Split the data
    #
    xsplit <- splitVariables(y=y,lags=NoLags,thDelay=thDelay,thresh=thVar,tart=tart,intercept=Intercept)
    
    #
    # Step 2: Sample posteriors
    #
    if(prior==1){
      
      # First regime
      postdraw <- postni(xsplit$y1,xsplit$x1,aprior=aprior,Vprior=Vprior,vprior=vprior,Sprior=Sprior,Sigma[,,1])
      Alpha[,,1] <- postdraw$Alpha
      Sigma[,,1] <- postdraw$Sigma
      
      # Second regime
      postdraw <- postni(xsplit$y2,xsplit$x2,aprior=aprior,Vprior=Vprior,vprior=vprior,Sprior=Sprior,Sigma[,,2])
      Alpha[,,2] <- postdraw$Alpha
      Sigma[,,2] <- postdraw$Sigma
    }
    
    #
    # Step 3: Sample new threshold using a Random-Walk Metropolis-Hastings Algorithm
    #
    
    tarnew <- tart+rnorm(1,sd=tarstandard)
    l1post <- tarpost(xsplit$xstar,xsplit$ystar,Ytest=ytest,Alpha[,,1],Alpha[,,2],Sigma[,,1],Sigma[,,2],tarnew,NoLags,Intercept,tarmean,tarstandard)
    l2post <- tarpost(xsplit$xstar,xsplit$ystar,Ytest=ytest,Alpha[,,1],Alpha[,,2],Sigma[,,1],Sigma[,,2],tart,NoLags,Intercept,tarmean,tarstandard)
    acc <- min(1,exp(l1post$post-l2post$post))
    u <- runif(1)
    if(u < acc){
      tart=tarnew
    }
    tarmean=tart

    
    #
    # Step 4: Sample new delay parameter
    #
    prob <- matrix(0,nrow=thMax)
    for(jj in 1:thMax){
      split1 <- splitVariables(y=y,lags=NoLags,jj,thVar,tart,Intercept)
      x <- exptarpost(split1$xstar,split1$ystar,split1$ytest,Alpha[,,1],Alpha[,,2],Sigma[,,1],Sigma[,,2],tart,NoLags,Intercept,tarmean,tarstandard,ncrit=0.15)
      prob[jj,1] <- x$post
    }
    prob <- prob/sum(prob)
    thDelay <- sample(thMax,1,replace=FALSE,prob)
    
    #
    # Compute Impulse-Response functions and save data
    #
    if(irep>burnin){
      # Save draws of coefficients and variance
      Alphadraws[,,,irep-burnin] <- Alpha
      Sigmadraws[,,,irep-burnin] <- Sigma
      
      # Compute and save impulse response-functions
      if(Intercept==TRUE){
        beta1 <- Alpha[2:nrow(Alpha),,1]
        beta2 <- Alpha[2:nrow(Alpha),,2]
      }
      else{
        beta1 <- Alpha[,,1]
        beta2 <- Alpha[,,2]
      }
      for(ii in 1:K){
        xx <- tirf(xsplit$ystar,xsplit$ytest,Alpha[,,1],Alpha[,,2],Sigma[,,1],Sigma[,,2],tart,thVar,thDelay,NoLags,irfhorizon,Intercept,shockvar=ii,bootrep)
        #print(dim(xx$irf1))
        Irfdraws[ii,,,1,irep-burnin]<-xx$irf1
        Irfdraws[ii,,,2,irep-burnin]<-xx$irf2
      }
      
      # Regimes
      nT <- length(xsplit$e1)
      a  <- nT-NoRegimes
      regimes[,irep-burnin] <- xsplit$e1[(1+a):nT]
      
    }
  }
  
  #
  # Quantiles of Impulse-Response functions
  #
  lowerquantile=min(irfquantiles)
  upperquantile=max(irfquantiles)
  for(ii in 1:K){
    for(jj in 1:K){
      for(kk in 1:irfhorizon){
        irffinal[ii,jj,kk,1,1] <- mean(Irfdraws[ii,jj,kk,1,])
        irffinal[ii,jj,kk,2,1] <- mean(Irfdraws[ii,jj,kk,2,])
        
        irffinal[ii,jj,kk,1,2] <- quantile(Irfdraws[ii,jj,kk,1,],probs=lowerquantile)
        irffinal[ii,jj,kk,2,2] <- quantile(Irfdraws[ii,jj,kk,2,],probs=lowerquantile)
        
        irffinal[ii,jj,kk,1,3] <- quantile(Irfdraws[ii,jj,kk,1,],probs=upperquantile)
        irffinal[ii,jj,kk,2,3] <- quantile(Irfdraws[ii,jj,kk,2,],probs=upperquantile)
      }
    }
  }
  
  return(list(Alphadraws=Alphadraws,Sigmadraws=Sigmadraws,irf=irffinal))
  
}

####################################################################################
# Function to split a time series into two dependent on threshold and lag order
#
# Input variables:
# 
# y - time series
# lags - Order of lags
# thDelay - delay of threshold
# thresh - threshold variable
# tart - value of threshold
# intercept - icnlude intercept or not
###################################################################################

splitVariables <- function(y,lags,thDelay,thresh,tart,intercept){
  startest <- max(thDelay,lags)
  T <- nrow(y)
  K <- ncol(y)
  ytest <- y[(startest+1-thDelay):(T-thDelay),thresh]
  ystar <- y[(startest+1):T,]
  xstar <- embed(y,dimension=lags+1)[,-(1:K)]
  if(thDelay>lags){
    diff1=(thDelay-lags)+1
    xnr <- nrow(xstar)
    xstar <- xstar[diff1:xnr,]
  }
  
  e1 <- ytest < tart
  e2 <- ytest >=tart
  y1 <- ystar[e1,]
  y2 <- ystar[e2,]
  x1 <- xstar[e1,]
  x2 <- xstar[e2,]
  if(intercept==TRUE){
    x1 <- cbind(1,x1)
    x2 <- cbind(1,x2)
  }
  return(list(y1=y1,y2=y2,x1=x1,x2=x2,xstar=xstar,ytest=ytest,ystar=ystar,e1=e1))
}

################################################################
#
# Function to calculate the acceptance probability
# X - lagged time series
# Y - time series
# beta1, beta2 - VAR Coefficients
# sigma1,sigma2 - variance/covariance matrix
# tart - value of threshold variable
# lags - lag order
# intercept - does the series contain an intercept or not
# tarmean - mean
# tarstandard - standard deviation
# ncrit - percentage of observation each regime must have
#
################################################################

tarpost <- function(X,Ystar,Ytest,beta1,beta2,sigma1,sigma2,tart,lags,intercept=TRUE,tarmean,tarstandard,ncrit=0.15){
  e1 <- Ytest<tart
  e2 <- Ytest>=tart
  nc<-nrow(Ystar)
  # test if there are enough observations in each sample
  p1=0
  if(sum(e1)/nc < ncrit || sum(e2)/nc<ncrit){
    post=-Inf
    loglik1=-Inf
    loglik2=-Inf
    prior=-Inf
    
  }
  else{
    Y1=Ystar[e1,]
    Y2=Ystar[e2,]
    X1=X[e1,]
    X2=X[e2,]
    if(intercept==TRUE){
      X1=cbind(1,X1)
      X2=cbind(1,X2)
    }
    loglik1 <- loglike(beta1,sigma1,Y1,X1)$lik
    loglik2 <- loglike(beta2,sigma2,Y2,X2)$lik
    prior<-pnorm(tart,mean=tarmean,sd=tarstandard)
  }
  post <- (loglik1+loglik2+prior)
  return(list(post=post,lik1=loglik1,lik2=loglik2,prior=prior))
  
}

exptarpost <- function(X,Ystar,Ytest,beta1,beta2,sigma1,sigma2,tart,lags,intercept=TRUE,tarmean,tarstandard,ncrit=0.15){
  e1 <- Ytest<tart
  e2 <- Ytest>=tart
  nc<-nrow(Ystar)
  # test if there are enough observations in each sample
  if(sum(e1)/nc < ncrit || sum(e2)/nc<ncrit){
    post=0
    loglik1=0
    loglik2=0
    prior=0
    
  }
  else{
    Y1=Ystar[e1,]
    Y2=Ystar[e2,]
    X1=X[e1,]
    X2=X[e2,]
    if(intercept==TRUE){
      X1=cbind(1,X1)
      X2=cbind(1,X2)
    }
    loglik1 <- loglike(beta1,sigma1,Y1,X1)
    loglik2 <- loglike(beta2,sigma2,Y2,X2)
    prior<-pnorm(tart,mean=tarmean,sd=tarstandard)
    ds1 <- loglik1$ds
    ds2 <- loglik2$ds
    sterm1 <- loglik1$sterm
    sterm2 <- loglik2$sterm
    post <- ds1*ds2*exp(-0.5*(sterm1+sterm2))
    
  }
  
  
  return(list(post=post,lik1=loglik1,lik2=loglik2,prior=prior))
  
}

loglike <- function(beta1,sigma,Y,X){
  T=nrow(Y)
  N=ncol(Y)
  v=Y-X%*%beta1
  sterm=0
  isigma <- solve(sigma)
  for(ii in 1:T){
    sterm<-sterm+t(v[ii,])%*%isigma%*%v[ii,]
  }
  dsigma <- log(det(isigma))
  lik <- (T/2)*dsigma-0.5*sterm
  return(list(lik=lik,ds=(T/2)*dsigma,sterm=sterm))
  
}