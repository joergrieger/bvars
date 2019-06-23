favar <- function(ydata,xdata,NoFactors = 2,NoLags = 2,RandomWalk = TRUE,prior = 1,priorm = 1,priorparam,alpha = 0.01,
                  beta = 0.01,lipr = 4,slowindex = "",frotvar = NaN,nreps = 20,burnin = 10,irfhor = 16,
                  irfquantiles = c(0.05,0.95),stabletest = TRUE, Intercept = FALSE){



  #
  #  Test for errors in input variables
  #

  #
  # Declare Variables and normalize data
  #

  varnames <- c(colnames(ydata),colnames(xdata))



  xdata <- scale(as.matrix(xdata))
  ydata <- scale(as.matrix(ydata))
  #x.demeaned <- demean(xdata)
  #y.demeaned <- demean(ydata)
  #x <- x.demeaned
  #y <- y.demeaned

  T <- nrow(ydata)
  N <- ncol(xdata)
  K <- ncol(ydata)
  P <- K+NoFactors

  #Intercept <- FALSE
  constant  <- 0
  if(Intercept == TRUE) constant  <- 1

  irfSmallDraws <- array(0,dim=c(P, P, irfhor, nreps - burnin))
  irfLargeDraws <- array(0,dim=c(P,(ncol(xdata) + ncol(ydata)), irfhor, nreps - burnin))

  Alphadraws <- array(0,dim=c(P * NoLags + constant, P, nreps - burnin))
  Sigmadraws <- array(0,dim=c(P, P ,nreps - burnin))
  Ldraws     <- array(0,dim=c(ncol(ydata) + ncol(xdata), P, nreps - burnin))

  #
  # extract factors using pricing components
  #

  fac <- exfact(ydata=ydata,xdata=xdata,slowcode=slowindex,NoFactors=NoFactors)

  # put it in state-space form

  XY <- cbind(ydata,xdata)
  FY <- cbind(ydata,fac)
  L  <- olssvd(XY,FY)

  e <- XY - FY%*%L
  Sigma <- t(e)%*%e/T

  # VAR equations

  FYlagged <- lagdata(FY,NoLags,intercept=Intercept)
  FY.x <- FYlagged$x
  FY.y <- FYlagged$y

  Alpha <- solve(t(FY.x) %*% FY.x) %*% t(FY.x) %*% FY.y
  aols <- as.vector(Alpha)
  SSE <- t(FY.y - FY.x %*% Alpha) %*% (FY.y - FY.x %*% Alpha)
  SF <- SSE/(T - NoLags)


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
    pr <- niprior(P,NoLags = NoLags,RandomWalk = RandomWalk,Intercept = Intercept,coefprior = coefprior,coefpriorvar = coefpriorvar,
                  varprior = varprior,varpriordof = varpriordof)

    aprior <- as.vector(pr$coefprior)
    Vprior <- pr$coefpriorvar
    vprior <- varpriordof
    Sprior <- pr$varprior

  }
  else if(prior==2){

    if(isempty(priorparam)){
	  stop("No prior parameters for the Minnesota Prior")

    }
    else{

      lambda1 <- priorparam[[1]]
      lambda2 <- priorparam[[2]]
      lambda3 <- priorparam[[3]]

    }
    # Minnesota prior
    pr <- mbprior(y=FY,NoLags = NoLags, Intercept=Intercept,RandomWalk = RandomWalk)
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

    pr <- ncprior(P,NoLags = NoLags,RandomWalk = RandomWalk,Intercept = FALSE,coefprior = coefprior,
                  coefpriorvar = coefpriorvar,varprior = varprior,varpriordof = varpriordof)
    aprior <- as.vector(pr$coefprior)
    Vprior <- pr$coefpriorvar
    vprior <- varpriordof
    Sprior <- pr$varprior

  }
  else if(prior==4){
    # uninformative prior, do nothing
  }
  else if(prior == 5){
    # Dummy - observation prior
  }
  else if(prior == 6){
    # SSVS - Prior
    xx <- lagdata(FY,lags=NoLags,intercept=Intercept)
    xlagged <- xx$x
    ylagged <- xx$y

    betaest <- solve(t(xlagged)%*%xlagged)%*%t(xlagged)%*%ylagged
    err     <- (ylagged - xlagged%*%betaest)
    sig     <- t(err)%*%err/(T-NoLags)
    omega   <- array(list(),dim=c(P-1))
    for(kk1 in 1:(P-1)){
      omega[[kk1]] <- array(1,dim=c(kk1))
    }
    if(Intercept == TRUE){
      constant = 1
    }
    else{
      constant = 0
    }

    NoRest <- P * (P * NoLags + constant)
    gammas <- array(0.5,dim=c(NoRest))
    si     <- sig%x%(t(xlagged)%*%xlagged)
    tau0   <- 0.01*si
    tau1   <- 1*si
    aprior <- array(0,dim=c(NoRest))
    SSEGibbs <- (array(0,dim=c(P,P)))
    SSEGibbs <- t(err)%*%err

  }


  if(priorm == 2){

    gammam <- array(0.5,dim=c(P,N))
    tau2 <- 0.0000001
    c2   <- 9/tau2

  }


  #
  # Start MCMC-Algorithm
  #
  L <- t(L)
  for(irep in 1:nreps){
    print(irep)

    # Step 1: Sample L and Sigma
    if(priorm == 1){
      # Normal priors for measurement equation
      for(ii in 1:N){
        if(ii > K){

          Li_postvar  <- solve(solve(Liprvar) + Sigma[ii,ii]^(-1) * t(FY) %*% FY)
          Li_postmean <- Li_postvar %*% (Sigma[ii,ii]^(-1) * t(FY) %*% xdata[,ii])
          L[ii,1:P]   <- t(Li_postmean) + rnorm(P) %*% chol(Li_postvar)

        }

        resi <- xdata[,ii] - FY %*% L[ii,]
        sh <- alpha/2 + T/2
        sc <- beta/2 + t(resi) %*% resi
        Sigma[ii,ii] <- rgamma(1,shape=sh,scale=sc)

      }
    }
    else if(priorm == 2){
      # SSVS-Prior for measurement equation
      for(ii in 1:N){
        if(ii > K){

          # Sample the betas
          VBeta <- diag(gammam[,ii] * c2 * tau2 + ( 1 - gammam[,ii] ) * tau2)
          DBeta <- solve(t(FY)%*%FY) * Sigma[ii,ii]^(-1)
          dBeta <- t(FY)%*%x[,ii] * Sigma[ii,ii]^(-1)
          HBeta <- t(chol(DBeta))
          xTemp <- rnorm(P)%*%HBeta
          L[ii,1:P] <- t(DBeta %*% dBeta) + (rnorm(P)%*%HBeta)

          # Sample the gammas
          for(jj in 1:P){
            numerator <- pnorm(L[ii,jj], mean = 0, sd = sqrt(c2*tau2))
            denominator <- numerator + pnorm(L[ii,jj],mean=0,sd=sqrt(tau2))
            prob <- numerator/denominator
            gammam[jj,ii] <- 0.5 * sign(runif(1) - prob) + 0.5
          }

          # Sample the variance
          resid <- xdata[,ii] - FY %*% L[ii,]
          sh <- alpha/2 + T/2
          sc <- beta/2 + t(resid) %*% resid
          Sigma[ii,ii] <- rgamma(1,shape=sh,scale=sc)

        }
      }
    }

    # Step 2: Sample VAR-coefficients for the state equation
    if(prior==1){
      postdraw <- postni(y=FY.y,x=FY.x,aprior,Vprior,vprior,Sprior,Sigma=SF,Intercept=Intercept,stabletest=stabletest)
      Alpha <- postdraw$Alpha
      SF <- postdraw$Sigma
    }
    else if(prior==2){

      postdraw <- postmb(y=FY.y,x=FY.x,Vprior=Vprior,aprior=aprior,Sigma=SF,betaols=aols,Intercept = Intercept)
      Alpha <- postdraw$Alpha
      SF <- postdraw$Sigma

    }
    else if(prior==3){

      postdraw <- postnc(y=FY.y,x=FY.x,aprior,Vprior,vprior,Sprior,Sigma=SF,Intercept = Intercept,stabletest=stabletest)
      Alpha <- postdraw$Alpha
      SF <- postdraw$Sigma

    }
    else if(prior==4){

      postdraw <- postun(y=FY.y,x=FY.x,Sigma=SF,Intercept=Intercept,stabletest=stabletest,Intercept=Intercept)
      Alpha <- postdraw$Alpha
      SF <- postdraw$Sigma

    }
    else if(prior == 6){

      aols     <- as.vector(solve(t(FY.x) %*% FY.x) %*% t(FY.x) %*% FY.y)
      postdraw <- postss(y=FY.y,x=FY.y,SSEGibbs=SSEGibbs, omega = omega,
                         gammas = gammas, tau0 = tau0, tau1 = tau1,
                         Sigma = SF, aols = aols, Intercept = Intercept,
                         NoLags = NoLags, stabletest = stabletest)

      Alpha <- postdraw$Alpha
      SF    <- postdraw$Sigma

      SSEGibbs <- postdraw$SSEGibbs
      gammas   <- postdraw$gammas
      omega    <- postdraw$omega

    }


    # compute impulse-response functions and save draws
    if(irep>burnin){

      # Save draws
      Alphadraws[,,irep-burnin] <- Alpha
      Sigmadraws[,,irep-burnin] <- SF
      Ldraws[,,irep-burnin]     <- L

      irf <- compirf(A=Alpha,Sigma=SF,NoLags=NoLags,intercept=FALSE,nhor=irfhor)
      irfSmallDraws[,,,irep-burnin] <- irf

      for(ii in 1:P){
        irfLargeDraws[ii,,,irep-burnin] <- L[,]%*%irf[ii,,]
      }
    }
  }
  # Final computations
  irfSmallFinal <- array(0,dim=c(P,P,irfhor,3))
  irfLargeFinal <- array(0,dim=c(P,(ncol(xdata)+ncol(ydata)),irfhor,3))
  irflower <- min(irfquantiles)
  irfupper <- max(irfquantiles)
  for(jj in 1:P){
    for(kk in 1:P){
      for(ll in 1:irfhor){

        irfSmallFinal[jj,kk,ll,1] <- mean(irfSmallDraws[jj,kk,ll,])
        irfSmallFinal[jj,kk,ll,2] <- quantile(irfSmallDraws[jj,kk,ll,],probs=irflower)
        irfSmallFinal[jj,kk,ll,3] <- quantile(irfSmallDraws[jj,kk,ll,],probs=irfupper)

      }
    }
  }
  for(jj in 1:P){
    for(kk in 1:(ncol(xdata)+ncol(ydata))){
      for(ll in 1:irfhor){

        irfLargeFinal[jj,kk,ll,1] <- mean(irfLargeDraws[jj,kk,ll,])
        irfLargeFinal[jj,kk,ll,2] <- quantile(irfLargeDraws[jj,kk,ll,],probs=irflower)
        irfLargeFinal[jj,kk,ll,3] <- quantile(irfLargeDraws[jj,kk,ll,],probs=irfupper)

      }
    }
  }

  # Return results in a favar object
  retlist <- structure(list(Betadraws = Alphadraws,Sigmadraws = Sigmadraws,Ldraws = Ldraws,
                            NoLags = NoLags,mydata = FY,
                            Intercept = Intercept,varnames = varnames),class="favar")

  return(retlist)

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
  #if(!isempty(slowcode)){
  #  slow <- which(slowcode==1)
  #  xslow <- xdata[,slow]
  #  ftemp <- extract(xslow,NoFactors)
  #  fnr <- ncol(ftemp$fac)
  #  if(is.nan(frotvar)){
  #    frotvar=ncol(ydata)
  #  }
  #  Fr0 <- facrot(factors$fac,ydata[,frotvar],ftemp$fac)
  #}
  #else{
    Fr0 <- factors$fac
  #}
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
