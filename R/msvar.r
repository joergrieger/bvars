msvar <- function(mydata,NoLags=1,NoRegimes=2,RandomWalk=TRUE,coefprior=NULL,coefpriorvar=10,prior=1,varprior=10,varpriordof=100,alphaprior=2,Intercept=TRUE,nreps=500,burnin=100,irfhor=40,irfquantiles=c(0.05,0.95),stabletest=TRUE){
  
  #
  # Preliminaries
  #
  
  y <- as.matrix(mydata)
  T <- nrow(y)
  K <- ncol(y)
  obs <- T-NoLags
  h <- NoRegimes
  
  constant=0
  if(Intercept==TRUE) constant=1
  
  #
  # Declare variables for storing results
  #
  
  Betadraws <- array(0,dim=c(K*NoLags+constant,K,h,nreps-burnin))
  Sigmadraws <- array(0,dim=c(K,K,h,nreps-burnin))
  irfdraws   <- array(0,dim=c(K,K,irfhor,h,nreps-burnin))
  
  ###############################################################
  #
  # Define Priors
  # We use an independent Normal-Wishart prior for the VAR-
  # coefficient in all regimes. The prior for the transition-
  # probability matrix is a Dirichlet prior.
  # All regimes have the same prior.
  #
  #############################################################
  
  if(prior==1){
    pr <- niprior(K=K,NoLags=NoLags,RandomWalk=RandomWalk,Intercept=Intercept,coefprior=coefprior,coefpriorvar=coefpriorvar,varprior=varprior)
    
    aprior <- as.vector(pr$coefprior)
    Vprior <- pr$coefpriorvar
    vprior <- varpriordof
    Sprior <- pr$varprior
    
  }
  
  
  # Prior on probabilities of transition matrix
  alphaprior <- alphaprior*diag(1,h)+matrix(2,h,h)
  stransmat <- matrix(0,h,h)
  
  ##################################################
  #
  # Getting starting values for the Gibbs sampler.
  # Here we use the OLS-Values from a VAR with a 
  # single regime.
  #
  ##################################################
  
  laggedData <- lagdata(y,lags=NoLags,intercept=Intercept)
  ylagged    <- laggedData$y
  xlagged    <- laggedData$x
  betaest    <- solve(t(xlagged)%*%xlagged)%*%t(xlagged)%*%ylagged
  err        <- (ylagged-xlagged%*%betaest)
  sig        <- t(err)%*%err/T
  Sigma      <- array(0,dim=c(nrow(sig),nrow(sig),h))
  Alpha      <- array(0,dim=c(K*NoLags+constant,K,h))
  
  for(ii in 1:h){
    Alpha[,,ii]  <- betaest
    Sigma[,,ii]  <- sig
  }
  
  ###########################################
  #
  # Start of the MCMC algorithm
  #
  ##########################################
  
  for(irep in 1:nreps){
    print(irep)
    
    #
    # Step 1: Draw S_t
    # Run Hamilton Filter first, then use results from hamilton filter to draw states
    #
    
    # Get posterior for Dirichlet prior and draw transition matrix
    alpha <- alphaprior + stransmat - 1
    transmat <- t(apply(alpha,1,rdirichlet,n=1))
    
    # Run Hamilton Filter and draw states
    filteredprob <- hamiltonfilter(Alpha,Sigma,transmat,ylagged,xlagged,h=NoRegimes)
    stt <- getst(filteredprob$fprob,transmat,h=NoRegimes)
    
    
    # Count transitions
    nseq <- countseq(stt,h=NoRegimes)
    stransmat <- nseq
    
    # Draw VAR coefficients and Variance for all states
    # Here we use the Independent Normal-Wishart prior
    for(ii in 1:h){
      xlaggedfilt <- xlagged[stt==ii,]
      ylaggedfilt <- ylagged[stt==ii,]
      if(prior==1){
        postdraw <- postni(y=ylaggedfilt,x=xlaggedfilt,aprior=aprior,Vprior=Vprior,vprior=vprior,Sprior=Sprior,Sigma=Sigma[,,ii])
      }
      
      Alpha[,,ii] <- postdraw$Alpha
      Sigma[,,ii] <- postdraw$Sigma

    }
    #readline(prompt="Press [enter] to continue")
    if(irep>burnin){
      # Save draws
      Betadraws[,,,irep-burnin] <- Alpha
      Sigmadraws[,,,irep-burnin] <- Sigma
      
      # Compute impulse-response functions for each regime
      for(ii in 1:h){
        irf <- compirf(Alpha[,,ii],Sigma[,,ii],NoLags=NoLags,intercept=Intercept,nhor=irfhor)
        irfdraws[,,,ii,irep-burnin] <- irf
      }
    }
  }
  
  #
  # Compute mean and quantiles of irfs
  #
  
  irffinal <- array(0,dim=c(K,K,irfhor,h,3))
  irflower <- min(irfquantiles)
  irfupper <- max(irfquantiles)
  for(ii in 1:h){
    for(jj in 1:K){
      for(kk in 1:K){
        for(ll in 1:irfhor){
          irffinal[jj,kk,ll,ii,1] <- mean(irfdraws[jj,kk,ll,ii,])
          irffinal[jj,kk,ll,ii,2] <- quantile(irfdraws[jj,kk,ll,ii,],probs=irflower)
          irffinal[jj,kk,ll,ii,3] <- quantile(irfdraws[jj,kk,ll,ii,],probs=irfupper)
        }
      }
    }
  }
  
  return(list(Betadraws=Betadraws,Sigmadraws=Sigmadraws,irf=irffinal))
}