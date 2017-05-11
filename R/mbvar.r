#' @export
mbvar <- function(mydata,lags=1,intercept=TRUE,RandomWalk=TRUE,lambda1=1,lambda2=1,lambda3=10,irfhorizon=16,reps=200,burnin=100,irfquantiles=c(0.05,0.9),verbose=TRUE,stabletest=FALSE){
  y<-as.matrix(mydata)
  # Create prior

  lagdata <- .data(y,lags,intercept)
  if(verbose==TRUE){
    print("Creating Prior")
  }
  prior <- .prior(y,lags=lags,intercept=intercept,RandomWalk=RandomWalk,lambda1=lambda1,lambda2=lambda2,lambda3=lambda3)
  if(verbose==TRUE){
    print("Starting the Gibbs Sampling")
  }
  results <- .gibbsSign(lagdata$y,lagdata$x,lagdata$obs,lagdata$K,intercept,prior$Vmatrix,prior$aprior,lags=lags,reps,burnin,irfhorizon,irfquantiles,verbose,stabletest)
  finalresults <- structure(list(aprior=prior$aprior,Vprior=prior$Vmatrix,betadist=results$betad,sigmadist=results$sigma,irf=results$irf),class="mbvar")
  return(finalresults)
}

mbvarSign <- function(mydata,lags=1,intercept=TRUE,RandomWalk=TRUE,lambda1=1,lambda2=1,lambda3=10,irfhorizon=16,reps=200,burnin=100,irfquantiles=c(0.05,0.9),verbose=TRUE,Restrictions,stabletest=FALSE){
  y<-as.matrix(mydata)
  # Create prior

  lagdata <- .data(y,lags,intercept)
  if(verbose==TRUE){
    print("Creating Prior")
  }
  prior <- .prior(y,lags=lags,intercept=intercept,RandomWalk=RandomWalk,lambda1=lambda1,lambda2=lambda2,lambda3=lambda3)
  if(verbose==TRUE){
    print("Starting the Gibbs Sampling")
  }
  results <- .gibbsSign(lagdata$y,lagdata$x,lagdata$obs,lagdata$K,intercept,prior$Vmatrix,prior$aprior,lags=lags,reps,burnin,irfhorizon,irfquantiles,verbose,Restrictions,stabletest)
  finalresults <- structure(list(aprior=prior$aprior,Vprior=prior$Vmatrix,betadist=results$betad,sigmadist=results$sigma,irf=results$irf),class="mbvar")
  return(finalresults)
}

.data <-function(y,lags,intercept=intercept){
  T<-nrow(y)
  K<-ncol(y)
  obs <- T-lags
  x  <- embed(y,dimension=lags+1)[,-(1:K)]
  if(intercept==TRUE){
    x<-cbind(1,x)
  }
  yi <- y[(lags+1):T,]
  return(list(y=yi,x=x,obs=obs,T=T,K=K));
}
 # Creates Minnesota-Prior
 #
 # Creates prior for means and Variance
 #
 # For a random walk prior, set coefficient for 1-period lag of a variable to 1
 # and every other coefficient to 0
 #
 # For a non-random walk prior set all coefficients to 0
 #
.prior <- function(y,lags=1,intercept=TRUE,RandomWalk=FALSE,lambda1=1,lambda2=1,lambda3=1){
  obs <-nrow(y)
  K <- ncol(y) # number of variabkes
  M <- K*lags  # number of rows for prior = K*lags or K*lags+1 (if there is an intercept)
  if(intercept==TRUE){M<-M+1}

  Aprior <- diag(0,M,K) # create an MxK matrix for the prior on the means

  if(RandomWalk==FALSE){
    aprior <- .vec(Aprior)
  }
  else if(RandomWalk==TRUE){

    if(intercept==FALSE){
      for(ii in 1:K){
        Aprior[ii,ii]=1
      }
    }
    else if(intercept==TRUE){
      for(ii in 1:K){
        Aprior[ii+1,ii]<-1
      }
    }
    aprior <- .vec(Aprior)
  }

  #
  # Create prior for Covariance matrix
  #

  # Run AR(p) model on each equation,ignoring constants
  sigmasq <- matrix(0,K,1)
  print(dim(y))
  for(ii in 1:K){
    Ylagi  <- embed(y[,ii],dimension=lags+1)[,-1]
    print(ii)
    Yi     <- y[(lags+1):obs,ii]
    # OLS estimation of AR(p) model without intercept
    arest <- lm(Yi~Ylagi-1)
    lmsummary <- summary(arest)
    sigmasq[ii,1] <- lmsummary$sigma
    #print(arest$residuals)
  }

  # Create Covariance Matrix for the prior
  # First: Without the intercept
  Vi <- array(0,dim=c(M*M*lags,1))
  for(ii in 1:K){ # Loop over the ii-th equation
    for(jj in 1:lags){ ##j-th lag
      for(kk in 1:K){ ## kk-th Variable
        indx <- (ii-1)*(K*lags)+(jj-1)*K+kk
        if(ii==kk){
          Vi[indx,1]<-lambda1/(jj^2)
        }
        else{
          Vi[indx,1]<-lambda2/(jj^2)*sigmasq[ii,1]/sigmasq[kk,1]
        }
      }
    }
  }
  #
  # Now: Add Covariance coefficients for intercepts
  #
  if(intercept==TRUE){
    Vtmp <- array(0,dim=c(K*K*lags+K,1))
    for(ii in 1:K){
      coefinter <- lambda3*sigmasq[ii,1]
      indx <- (ii-1)*(K*lags)+ii
      Vtmp[indx,1]<-coefinter
      indxmin <- (ii-1)*(K*lags)+1
      indxmax <- (ii-1)*(K*lags)+K*lags
      Vtmp[(indx+1):(indx+(K*lags)),]<-Vi[indxmin:indxmax,]
    }
    Vi <- Vtmp
  }
  # Create diagonal Matrix
  nr <- dim(Vi)[1]
  Vfinal <- array(0,dim=c(nr,nr))
  for(ii in 1:nr){
    Vfinal[ii,ii]<-Vi[ii,1]
  }
  return(list(aprior=aprior,Vmatrix=Vfinal))
}

.gibbs <- function(y,x,obs,K,intercept,Vprior,aprior,lags,reps,burnin,irfhorizon,irfquantiles,verbose,stabletest){
  # OLS estimates
  Betaols  <- solve(t(x)%*%x)%*%t(x)%*%y
  betaols  <- matrix(Betaols,ncol=1)

  epsi     <- y-x%*%Betaols
  sigmaols <- t(epsi)%*%epsi

  # Get posterior
  tmp1       <- solve(sigmaols)%x%(t(x)%*%x)
  Vposterior <- solve(solve(Vprior)+tmp1)
  tmp1       <- (solve(sigmaols)%x%(t(x)%*%x))%*%betaols
  aposterior <- Vposterior%*%(solve(Vprior)%*%aprior+tmp1)
  Aposterior <- matrix(aposterior,ncol=K)
  cholsigma  <- t(chol(sigmaols))

  #
  # Calculate Impulse-Response function using Gibbs Sampling
  #

  irfs <- array(0,dim=c(K,irfhorizon,K,reps))
  if(intercept==TRUE){
    betadist <- array(0,dim=c(K*lags+1,K,reps-burnin))
  }
  else{
    betadist <- array(0,dim=c(K*lags,K,reps-burnin))
  }
  sigmadist <- array(0,dim=c(K,K,reps-burnin))

  Sigma=solve(sigmaols)
  if(verbose==TRUE){print("Calculating impulse-response function using Gibbs sampling")}
  # Start Gibbs Sampling
  for(ii in 1:reps){
    if(verbose==TRUE){
      print(ii)
    }
    stable=2

    Vpost <- solve((solve(Vprior)+solve(Sigma)%x%(t(x)%*%x)))
    apost <- Vpost%*%(solve(Vprior)%*%aprior+(solve(Sigma)%x%(t(x)%*%x))%*%betaols)
    while(stable>1){
      # Draw betas
      betadraw1 <- mvrnorm(mu=apost,Sigma=Vpost)
      betadraw1 <- matrix(betadraw1,ncol=K)
      # Draw Sigmas
      e <- y-x%*%betadraw1
      scale <- (t(e)%*%e)/(obs+1)
      Sigma1 <- solve(rWishart(1,df=obs,scale)[,,1])
      if(intercept==TRUE){
        nrb <- nrow(betadraw1)
        betatest <- betadraw1[2:nrb,]
      }
      stable<-.stability(betatest,lags,K)
    }
    betadraw <- betadraw1
    Sigma    <- Sigma1

    betadist[,,ii-burnin]<-betadraw1
    sigmadist[,,ii-burnin]<-Sigma1
    #print(stable)
    #readline(prompt="Press [enter] to continue")
    if(ii>burnin){
      cholsigma <- t(chol(Sigma))
      for(jj in 1:K){ # loop over the variables
        shock <- array(0,dim=c(K,irfhorizon+lags))
        shock[jj,lags+1]<-1
        yhat <- .irfsimu(beta=betatest,sigma=cholsigma,shocks=shock,horizon=irfhorizon,lags=lags,intercept=intercept,K=K)
        irfs[,,jj,ii-burnin]<-yhat
      }
    }
  }
  upperquantile <- max(irfquantiles)
  lowerquantile <- min(irfquantiles)
  irffinal <- array(0,dim=c(K,irfhorizon,K,3))
  for(ii in 1:K){
    for(jj in 1:K){
      for(kk in 1:irfhorizon){
        irffinal[jj,kk,ii,1]<-quantile(irfs[jj,kk,ii,],probs=0.5)
        irffinal[jj,kk,ii,2]<-quantile(irfs[jj,kk,ii,],probs=lowerquantile)
        irffinal[jj,kk,ii,3]<-quantile(irfs[jj,kk,ii,],probs=upperquantile)
      }
    }
  }
  retlist <- list(irf=irffinal,betad=betadist,sigma=sigmadist)
  return(retlist)
}

.gibbsSign <- function(y,x,obs,K,intercept,Vprior,aprior,lags,reps,burnin,irfhorizon,irfquantiles,verbose,Restrictions,stabletest){
  # OLS estimates
  print("Testing Minnesota Prior and Sign Restrictions")
  Betaols  <- solve(t(x)%*%x)%*%t(x)%*%y
  betaols  <- matrix(Betaols,ncol=1)

  epsi     <- y-x%*%Betaols
  sigmaols <- t(epsi)%*%epsi

  # Get posterior
  tmp1       <- solve(sigmaols)%x%(t(x)%*%x)
  Vposterior <- solve(solve(Vprior)+tmp1)
  tmp1       <- (solve(sigmaols)%x%(t(x)%*%x))%*%betaols
  aposterior <- Vposterior%*%(solve(Vprior)%*%aprior+tmp1)
  Aposterior <- matrix(aposterior,ncol=K)
  cholsigma  <- t(chol(sigmaols))

  #
  # Calculate Impulse-Response function using Gibbs Sampling
  #

  irfs <- array(0,dim=c(K,irfhorizon,K,reps))
  if(intercept==TRUE){
    betadist <- array(0,dim=c(K*lags+1,K,reps-burnin))
  }
  else{
    betadist <- array(0,dim=c(K*lags,K,reps-burnin))
  }
  sigmadist <- array(0,dim=c(K,K,reps-burnin))

  Sigma=solve(sigmaols)
  if(verbose==TRUE){print("Calculating impulse-response function using Gibbs sampling")}
  # Start Gibbs Sampling
  for(ii in 1:reps){
    if(verbose==TRUE){
      print(ii)
    }
    stable=2

    Vpost <- solve((solve(Vprior)+solve(Sigma)%x%(t(x)%*%x)))
    apost <- Vpost%*%(solve(Vprior)%*%aprior+(solve(Sigma)%x%(t(x)%*%x))%*%betaols)
    while(stable>1){
      # Draw betas
      betadraw1 <- mvrnorm(mu=apost,Sigma=Vpost)
      betadraw1 <- matrix(betadraw1,ncol=K)
      # Draw Sigmas
      e <- y-x%*%betadraw1
      scale <- (t(e)%*%e)/(obs+1)
      Sigma1 <- solve(rWishart(1,df=obs,scale)[,,1])
      if(intercept==TRUE){
        nrb <- nrow(betadraw1)
        betatest <- betadraw1[2:nrb,]
      }
      stable<-.stability(betatest,lags,K)
    }
    betadraw <- betadraw1
    Sigma    <- Sigma1

    betadist[,,ii-burnin]<-betadraw1
    sigmadist[,,ii-burnin]<-Sigma1
    #print(stable)
    #readline(prompt="Press [enter] to continue")
    if(ii>burnin){
      cholsigma <- t(chol(Sigma))
      SignRestriction <-FALSE
      while(!SignRestriction){
        qrmatrix <- matrix(rnorm(K*K),nrow=K)
        qrdecomp <- qr(qrmatrix)
        testmatrix <- qrdecomp$qr%*%cholsigma
        SignRestriction<-!.CheckSign(Restrictions,testmatrix)
      }

      for(jj in 1:K){ # loop over the variables
        shock <- array(0,dim=c(K,irfhorizon+lags))
        shock[jj,lags+1]<-1
        yhat <- .irfsimu(beta=betatest,sigma=cholsigma,shocks=shock,horizon=irfhorizon,lags=lags,intercept=intercept,K=K)
        irfs[,,jj,ii-burnin]<-yhat
      }
    }
  }
  upperquantile <- max(irfquantiles)
  lowerquantile <- min(irfquantiles)
  irffinal <- array(0,dim=c(K,irfhorizon,K,3))
  for(ii in 1:K){
    for(jj in 1:K){
      for(kk in 1:irfhorizon){
        irffinal[jj,kk,ii,1]<-quantile(irfs[jj,kk,ii,],probs=0.5)
        irffinal[jj,kk,ii,2]<-quantile(irfs[jj,kk,ii,],probs=lowerquantile)
        irffinal[jj,kk,ii,3]<-quantile(irfs[jj,kk,ii,],probs=upperquantile)
      }
    }
  }
  retlist <- list(irf=irffinal,betad=betadist,sigma=sigmadist)
  return(retlist)
}

