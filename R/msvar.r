msvar <- function(mydata,NoLags=1,NoRegimes=2,RandomWalk=TRUE,Intercept=TRUE,coefprior=NULL,prior=1,priorparam,alphaprior=2,nreps=500,burnin=100,irfhor=40,irfquantiles=c(0.05,0.95),ident=1,restrictions=NULL,stabletest=TRUE){

  #
  # Preliminaries
  #

  y <- as.matrix(mydata)
  T <- nrow(y)
  K <- ncol(y)
  obs <- T-NoLags
  h <- NoRegimes
  transmat <- array(0,dim=c(h,h))

  constant=0
  if(Intercept==TRUE) constant=1

  #
  # Declare variables for storing results
  #

  Betadraws <- array(0,dim=c(K*NoLags+constant,K,h,nreps-burnin))
  Sigmadraws <- array(0,dim=c(K,K,h,nreps-burnin))
  irfdraws   <- array(0,dim=c(K,K,irfhor,h,nreps-burnin))
  sttsave    <- array(0,dim=c(obs,nreps-burnin))
  trmsave    <- array(0,dim=c(h,h,nreps-burnin))

  ###############################################################
  #
  # Define Priors
  # The user has the choice between several priors on the VAR-
  # coefficients:
  #
  # 1 - independent Normal-Wishart prior
  # 2 - Minnesota Prior (not implemented)
  # 3 - Natural Conjugate prior
  # 4 - Uninformative Prior
  # The prior for the transition-
  # probability matrix is a Dirichlet prior.
  # All regimes have the same prior.
  #
  #############################################################

  if(prior==1){
    if(isempty(priorparam)){
	  stop("No prior parameters for Independent Normal-Wishart prior")
	}
	coefprior    <- priorparam[[1]]
	coefpriorvar <- priorparam[[2]]
	varprior     <- priorparam[[3]]
	varpriordof  <- priorparam[[4]]

    # Independent Normal-Wishart prior
    pr <- niprior(K=K,NoLags=NoLags,RandomWalk=RandomWalk,Intercept=Intercept,coefprior=coefprior,coefpriorvar=coefpriorvar,varprior=varprior)

    aprior <- as.vector(pr$coefprior)
    Vprior <- pr$coefpriorvar
    vprior <- varpriordof
    Sprior <- pr$varprior

  }
  else if(prior==3){
    if(isempty(priorparam)){
	  stop("No prior parameters for Natural conjugate prior")
	}
	coefprior    <- priorparam[[1]]
	coefpriorvar <- priorparam[[2]]
	varprior     <- priorparam[[3]]
	varpriordof  <- priorparam[[4]]

    # Natural Conjugate prior
    pr <- ncprior(K=K,NoLags=NoLags,RandomWalk=RandomWalk,Intercept=Intercept,
	              coefprior=coefprior,coefpriorvar=coefpriorvar,varprior=varprior)

    aprior <- pr$coefprior
    Vprior <- pr$coefpriorvar
    vprior <- varpriordof
    Sprior <- pr$varprior
  }
  else if(prior==4){
    # Uninformative prior, do nothing
  }
  else if(prior == 6){
    #SSVS prior
    xx <- lagdata(mydata,lags=NoLags,intercept=Intercept)
    xlagged <- xx$x
    ylagged <- xx$y
    betaest    <- solve(t(xlagged) %*% xlagged) %*% t(xlagged) %*% ylagged
    err        <- (ylagged - xlagged %*% betaest)
    sig        <- t(err) %*% err / (T - NoLags)
    omega    <- array(list(),dim=c(K-1,NoRegimes))
    for(kk1 in 1:(K-1)){
      for(ii in 1:NoRegimes){
        omega[[kk1,ii]] <- array(1,dim=c(kk1))
      }
    }
    if(Intercept==TRUE){
      constant=1
    }
    else{
      constant=0
    }

    NoRest <- K*(K*NoLags+constant)
    gammas <- array(1,dim=c(NoRest,NoRegimes))
    si <- sig%x%(t(xlagged)%*%xlagged)
    tau0 <- 0.1*si
    tau1 <-  10*si
    aprior <- array(0,dim=c(NoRest))
    SSEGibbs <- array(0,dim=c(K,K,2))

    for(ii in 1:2){
      SSEGibbs[,,ii] <- t(err)%*%err
    }

  }


  # Prior on probabilities of transition matrix
  if(!is.null(alphaprior)){
	# prior for alpha is not empty
	if(isscalar(alphaprior)){
		# if value for alpha prior is a scalar
		alphaprior <- alphaprior*matrix(1,h,h)
	}
	else{
	  if(nrow(alphaprior)!=h || ncol(alphaprior)!=h){
		stop("wrong number of elements in prior for transition probability matrix")
	  }
	}
  }
  else{
	# Set alphaprior to some standardized values
	Tmax <- T/h
	Toffdiag <- as.integer(Tmax/(h))
	Tdiag    <- Tmax-Toffdiag
	Toffdiag <- max(as.integer(Toffdiag/2),2)
	Tdiag    <- max(as.integer(Tdiag/2),2)
	alphaprior <- matrix(Toffdiag,h,h) + diag(Tdiag,h)
  }

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
  sig        <- t(err)%*%err/(T-NoLags)
  Sigma      <- array(0,dim=c(nrow(sig),nrow(sig),h))
  Alpha      <- array(0,dim=c(K*NoLags+constant,K,h))

  for(ii in 1:h){
    Alpha[,,ii]  <- betaest/(1.5^(ii-1))
    Sigma[,,ii]  <- sig/(1.5^(ii-1))
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
    alpha <- alphaprior + stransmat - matrix(1,h,h)
    transmat <- t(apply(alpha,1,rdirichlet,n=1))


    # Run Hamilton Filter and draw states
    filteredprob <- hamiltonfilter(Alpha,Sigma,transmat,ylagged,xlagged,h=NoRegimes)
	  stt <- getst(filteredprob$fprob,transmat,h=NoRegimes)

    # Count transitions
    nseq <- countseq(stt,h=NoRegimes)
    stransmat <- nseq

    # Draw VAR coefficients and Variance for all states
    for(ii in 1:h){
      xlaggedfilt <- xlagged[stt==ii,]
      ylaggedfilt <- ylagged[stt==ii,]

      if(prior==1){

        # Independent Normal-Wishart
        postdraw <- postni(y=ylaggedfilt,x=xlaggedfilt,aprior=aprior,Vprior=Vprior,vprior=vprior,Sprior=Sprior,Sigma=Sigma[,,ii],stabletest=stabletest,NoLags=NoLags,Intercept=Intercept)

      }
      else if(prior==3){

        # Natural Conjugate Prior
        postdraw <- postnc(y=ylaggedfilt,x=xlaggedfilt,aprior=aprior,Vprior=Vprior,vprior=vprior,Sprior=Sprior,Sigma=Sigma[,,ii],stabletest=stabletest,NoLags=NoLags,Intercept=Intercept)

      }
      else if(prior==4){

        # Uninformative prior
        postdraw <- postun(y=ylaggedfilt,x=xlaggedfilt,Sigma=Sigma[,,ii],stabletest=stabletest,NoLags=NoLags,Intercept=Intercept)

      }


      Alpha[,,ii] <- postdraw$Alpha
      Sigma[,,ii] <- postdraw$Sigma


    }
	# Identification of regimes
	sorder <- order(Sigma[1,1,])
	Sigma <- Sigma[,,rev(sorder)]
	Alpha <- Alpha[,,rev(sorder)]

	# Compute Impulse-Response Functions and store all important results

    if(irep>burnin){
      # Save draws for VAR-coefficients and Variance-Covariance matrix
      Betadraws[,,,irep-burnin] <- Alpha
      Sigmadraws[,,,irep-burnin] <- Sigma

	  # save draws for states and transition matrix
	  sttsave[,irep-burnin]  <- stt
	  trmsave[,,irep-burnin] <- transmat

      # Compute impulse-response functions for each regime
      for(ii in 1:h){
	    if(ident==1){
	      irf <- compirf(Alpha[,,ii],Sigma[,,ii],NoLags=NoLags,intercept=Intercept,nhor=irfhor)
	    }
	    else if(ident==2){
		  irf <- compirfsign(Alpha[,,ii],Sigma[,,ii],NoLags=NoLags,intercept=Intercept,nhor=irfhor,restrictions=restrictions)
	    }
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
  varnames = colnames(mydata)
  retlist <- structure(list(Betadraws=Betadraws,Sigmadraws=Sigmadraws,irf=irffinal,stt=sttsave,trm=trmsave,
                            NoRegimes = NoRegimes, mydata=mydata, varnames = varnames),class="msvar")
  return(retlist)
}
