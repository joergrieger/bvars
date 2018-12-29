# bvar - estimates a vector autoregressive model using bayesian inference
#
# mydata - multivariate time series
# NoLags - Number of lags
# Intercept - whether the model has an intercept (True/False)
# prior

#' @export
bvar <- function(mydata,NoLags=1,Intercept=TRUE,RandomWalk=TRUE,prior=1,priorparam,irfhorizon=16,irfquantiles=c(0.05,0.95),ident=1,Restrictions=NULL,nreps=110,burnin=10,stabletest=TRUE){
  ###############################
  #
  # Declare Variables
  #
  ###############################

  Y <- as.matrix(mydata)
  T <- nrow(Y)
  K <- ncol(Y)
  obs <- T-NoLags
  constant <- 0
  if(Intercept == TRUE) constant=1

  # Variables for storage
  betadraws <- array(0,dim=c(K*NoLags+constant,K,nreps-burnin))
  sigmadraws <- array(0,dim=c(K,K,nreps-burnin))
  irfdraws <- array(0,dim=c(K,K,irfhorizon,nreps-burnin))
  irffinal <- array(0,dim=c(K,K,irfhorizon,3))
  varnames <- colnames(mydata)

  if(is.ts(mydata)){

    dates = time(mydata)

  }
  else{

    dates=NULL

  }

  ##############################
  #
  # Check if input is correct
  #
  ##############################
  if(prior>6){
    stop("Invalid choice for prior")
  }
  if(ident>2){

    stop("Invalid choice for identification of structural shocks")

  }

  ##############################
  #
  # Create prior
  #
  ##############################
  if(prior==1){
    # Independent Normal-Wishart Prior
	if(isempty(priorparam)){
	  stop("No prior parameters for Independent Normal-Wishart prior")
	}
	coefprior    <- priorparam[[1]]
	coefpriorvar <- priorparam[[2]]
	varprior     <- priorparam[[3]]
	varpriordof  <- priorparam[[4]]

    pr <- niprior(K = K,NoLags = NoLags,RandomWalk = RandomWalk,Intercept = Intercept,coefprior = coefprior,
                  coefpriorvar = coefpriorvar, varprior = varprior)

    aprior <- as.vector(pr$coefprior)
    Vprior <- pr$coefpriorvar
    vprior <- varpriordof
    Sprior <- pr$varprior
  }
  else if(prior==2){

    # Minnesota Prior
    if(isempty(priorparam)){
	  stop("No prior parameters for Minnesota Prior")

    }


    else{

      lambda1 <- priorparam[[1]]
      lambda2 <- priorparam[[2]]
      lambda3 <- priorparam[[3]]

    }

    pr <- mbprior(y = mydata, NoLags = NoLags, Intercept = Intercept, RandomWalk = RandomWalk, lambda1 = lambda1,
                  lambda2 = lambda2, lambda3 = lambda3)
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
    # Natural Conjugate Prior
    pr <- ncprior(K = K, NoLags = NoLags, RandomWalk = RandomWalk, Intercept = Intercept, coefprior = coefprior,
                  coefpriorvar = coefpriorvar, varprior = varprior)

    aprior <- matrix(pr$coefprior,ncol=1)
    Vprior <- pr$coefpriorvar
    vprior <- varpriordof
    Sprior <- pr$varprior
  }
  else if(prior == 4){
    # Uninformative prior, do nothing
  }
  else if(prior == 5){
    # Dummy prior
    # to implement
  }
  else if(prior == 6){

    # SSVS prior

    # Prior on Coefficients

    xx <- lagdata(mydata,lags=NoLags,intercept=Intercept)
    xlagged <- xx$x
    ylagged <- xx$y
    betaest    <- solve(t(xlagged)%*%xlagged)%*%t(xlagged)%*%ylagged
    err        <- (ylagged-xlagged%*%betaest)
    sig        <- t(err)%*%err/(T-NoLags)
    omega    <- array(list(),dim=c(K-1,1))

    for(kk1 in 1:(K-1)){

      for(ii in 1:1){

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
    gammas <- array(1,dim=c(NoRest,1))
    si <- sig%x%(t(xlagged)%*%xlagged)
    tau0 <- 0.1*si
    tau1 <-  10*si
    aprior <- array(0,dim=c(NoRest))
    SSEGibbs <- array(0,dim=c(K,K))
    SSEGibbs <- t(err) %*% err

  }


  ######################################
  #
  # Initialize the MCMC algorithm
  #
  ######################################

  # lag data
  dat <- lagdata(Y,lags=NoLags,intercept=Intercept)
  y.lagged <- dat$y
  x.lagged <- dat$x

  # OLS estimates
  Aols <- solve(t(x.lagged)%*%x.lagged)%*%t(x.lagged)%*%y.lagged
  aols <- as.vector(Aols)
  resi <- y.lagged-x.lagged%*%Aols
  SSE  <- t(resi)%*%resi
  Sigma <- SSE/T

  ######################################
  #
  # Start Gibbs Sampling
  #
  ######################################
  for(irep in 1:nreps){
    print(irep)
    #readline(prompt="Press [enter] to continue")
    if(prior==1){
      postdraw <- postni(y=y.lagged,x=x.lagged,aprior=aprior,Vprior=Vprior,vprior=vprior,Sprior=Sprior,Sigma=Sigma,stabletest=TRUE,Intercept=Intercept,NoLags=NoLags)
      Alpha <- postdraw$Alpha
      Sigma <- postdraw$Sigma
    }
    else if(prior==2){
      postdraw <- postmb(y=y.lagged,x=x.lagged,Vprior=Vprior,aprior=aprior,Sigma=Sigma,betaols=aols)
      Alpha <- postdraw$Alpha
      Sigma <- postdraw$Sigma
    }
	  else if(prior==3){

	    postdraw <- postnc(y=y.lagged,x=x.lagged,aprior=aprior,Vprior=Vprior,vprior=vprior,Sprior=Sprior,Sigma=Sigma,stabletest=TRUE,Intercept=Intercept,NoLags=NoLags)
	    Alpha <- postdraw$Alpha
	    Sigma <- postdraw$Sigma

	  }
	  else if(prior==4){
	    postdraw <- postun(y=y.lagged,x=x.lagged,Sigma=Sigma,stabletest=TRUE,Intercept=Intercept,NoLags=NoLags)
	    Alpha <- postdraw$Alpha
	    Sigma <- postdraw$Sigma
	  }
    else if(prior == 6){

      #aols <- as.vector(solve(t(xsplit$x1)%*%xsplit$x1)%*%t(xsplit$x1)%*%xsplit$y1)
      aols <- as.vector(solve(t(x.lagged)%*%x.lagged)%*%t(x.lagged)%*%y.lagged)
      postdraw <- postss(y=y.lagged,x=x.lagged,SSEGibbs=SSEGibbs[,],omega=omega,gammas=gammas,tau0=tau0,tau1=tau1,Sigma=Sigma,aprior=aprior,aols=aols,Intercept=Intercept,NoLags=NoLags,stabletest=stabletest)

      Alpha <- postdraw$Alpha
      Sigma <- postdraw$Sigma

      SSEGibbs <- postdraw$SSEGibbs
      gammas   <- postdraw$gammas
      omega    <- postdraw$omega

    }


    if(irep>burnin){
      # compute and save impulse-response functions
      if(ident==1){
        # Recursive identification
        irf <- compirf(A=Alpha,Sigma=Sigma,NoLags=NoLags,intercept=Intercept,nhor=irfhorizon)
      }
      else if(ident==2){
        # Identification using Sign restrictions
		if(is.null(Restrictions)){
		  stop("No Restrictions provided")
		}
		irf <- compirfsign(A=Alpha,Sigma=Sigma,NoLags=NoLags,intercept=Intercept,nhor=irfhorizon,restrictions=Restrictions)
      }
      irfdraws[,,,irep-burnin] <- irf
      #plot(irf[1,1,],type="l")

      # save draws
      betadraws[,,irep-burnin] <- Alpha
      sigmadraws[,,irep-burnin] <- Sigma

    }
  }
  # Final computations
  irffinal <- array(0,dim=c(K,K,irfhorizon,3))
  irflower <- min(irfquantiles)
  irfupper <- max(irfquantiles)
  for(jj in 1:K){
    for(kk in 1:K){
      for(ll in 1:irfhorizon){

        irffinal[jj,kk,ll,1] <- median(irfdraws[jj,kk,ll,])
        irffinal[jj,kk,ll,2] <- quantile(irfdraws[jj,kk,ll,],probs=irflower)
        irffinal[jj,kk,ll,3] <- quantile(irfdraws[jj,kk,ll,],probs=irfupper)

      }
    }
  }

  relist <- structure(list(type=prior,intercept=Intercept,betadraws=betadraws,sigmadraws=sigmadraws,irfdraws=irffinal,varnames=varnames,
                           NoLags=NoLags,mydata=mydata),class="bvar")
  return(relist)

}
