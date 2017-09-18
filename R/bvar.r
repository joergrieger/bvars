#' @export
bvar <- function(mydata,NoLags=1,Intercept=TRUE,prior=1,RandomWalk=TRUE,coefprior=NULL,coefpriorvar=1,varprior=1,varpriordof=100,lambda1=1,lambda2=1,lambda3=100,irfhorizon=16,irfquantiles=c(0.1,0.9),ident=1,Restrictions=NULL,nreps=110,burnin=10,forcasthorizon=10,stabletest=TRUE){
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
  if(Intercept==TRUE) constant=1

  # Variables for storage
  betadraws <- array(0,dim=c(K*NoLags+constant,K,nreps-burnin))
  sigmadraws <- array(0,dim=c(K,K,nreps-burnin))
  irfdraws <- array(0,dim=c(K,K,irfhorizon,nreps-burnin))
  irffinal <- array(0,dim=c(K,K,irfhorizon,3))

  ##############################
  #
  # Check if input is correct
  #
  ##############################
  if(prior>2){
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
    pr <- niprior(K=K,NoLags=NoLags,RandomWalk=RandomWalk,Intercept=Intercept,coefprior=coefprior,coefpriorvar=coefpriorvar,varprior=varprior)

    aprior <- as.vector(pr$coefprior)
    Vprior <- pr$coefpriorvar
    vprior <- varpriordof
    Sprior <- pr$varprior
  }
  else if(prior==2){
    # Minnesota Prior
    pr <- mbprior(y=mydata,NoLags=NoLags,intercept=Intercept,RandomWalk=RandomWalk,lambda1=lambda1,lambda2=lambda2,lambda3=lambda3)
    aprior <- pr$aprior
    Vprior <- pr$Vmatrix

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
      postdraw <- postni(y=y.lagged,x=x.lagged,aprior=aprior,Vprior=Vprior,vprior=vprior,Sprior=Sprior,Sigma=Sigma)
      Alpha <- postdraw$Alpha
      Sigma <- postdraw$Sigma
    }
    else if(prior==2){
      postdraw <- postmb(y=y.lagged,x=x.lagged,Vprior=Vprior,aprior=aprior,Sigma=Sigma,betaols=aols)
      Alpha <- postdraw$Alpha
      Sigma <- postdraw$Sigma
    }

    if(irep>burnin){
      # compute and save impulse-response functions
      if(ident==1){
        # Recursive identification
        irf <- compirf(A=Alpha,Sigma=Sigma,NoLags=NoLags,intercept=Intercept,nhor=irfhorizon)
      }
      else if(ident==2){
        # Identification using Sign restrictions

      }
      irfdraws[,,,irep-burnin] <- irf
      #plot(irf[1,1,],type="l")

      # save draws
      betadraws[,,irep-burnin] <- Alpha
      sigmadraws[,,irep-burnin] <- Sigma

    }
  }
  relist <- list(type=prior,betadraws=betadraws,sigmadraws=sigmadraws,irfdraws=irfdraws)
  return(relist)

}
