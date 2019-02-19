ftvar <- function(mydata,factors, NoFactors = 1, NoLags = 1, slowindex = "", frotvar, thMax = 4, thVar = 1,
                  RandomWalk = TRUE, Intercept = TRUE, prior = 1, priorm = 1, priorparam = NULL, Lipr = 4,
                  nreps = 200,burnin = 100,irfhorizon = 20, irfquantiles = c(0.05,0.95), bootrep = 10,
                  ncrit = 0.2,stabletest = TRUE, thin = 10){
  #
  # Preliminaries
  #

  varnames <- c(colnames(mydata),colnames(factors))

  constant <- 0
  if(Intercept==TRUE){
    constant <- 1
  }
  xdata <- as.matrix(factors)
  ydata <- as.matrix(mydata)
  x <- demean(xdata)
  y <- demean(ydata)

  T <- nrow(y)
  N <- ncol(x)
  K <- ncol(y)
  P <- K+NoFactors

  # normalize data
  x <- scale(x)

  #
  # Check if input is correct
  #

  #
  # Declare Variables for storage
  #
  startest <- max(thMax,NoLags)
  xsave <- (nreps - burnin)/thin
  Alphadraws <- array(0,dim=c(P*NoLags+constant,P,2,xsave))
  Sigmadraws <- array(0,dim=c(P,P,2,nreps-burnin))
  irfdraws   <- array(0,dim=c(P,P,irfhorizon,2,xsave))
  Ldraws     <- array(0,dim=c(ncol(y)+ncol(x),P,2,xsave))
  irfSmalldraws <- array(0,dim=c(P,P,irfhorizon,2,xsave))
  irfLargedraws  <- array(0,dim=c(P,ncol(y)+ncol(x),irfhorizon,2,xsave))
  tardraws   <- array(0,dim=c(xsave))
  deldraws   <- array(0,dim=c(xsave))
  NoRegimes  <- T-(startest+1)
  regimes    <- array(0,dim=c(NoRegimes,xsave))
  gammamdraws <- array(NA,dim=c(P,N+K,2,xsave))

  #
  # Extract factors and put it in state-space form
  #
  print("extracting factors")
  fac <- exfact(ydata=y,xdata=x,slowcode=slowindex,NoFactors=NoFactors)
  print("putting it into state-space form")
  XY  <- cbind(y,x)
  FY  <- cbind(y,fac)
  Li  <- olssvd(XY,FY)

  res <- XY-FY%*%Li
  Sig <- t(res)%*%res/T

  L <- array(0,dim=c(ncol(XY),P,2))
  Sigma <- array(0,dim=c(ncol(XY),ncol(XY),2))
  for(ii in 1:2){
    Sigma[,,ii] <- Sig
    L[,,ii] <- Li
  }
  if(priorm == 2){
    if(priorm == 2){
      gammam <- array(0.5,dim=c(P,ncol(XY),2))
      tau2 <- 0.001
      c2   <- 1/tau2
    }

  }

  #
  # Set priors
  #

  # Var coefficients
  if(prior==1){
    # Independent Normal-Wishart prior
    if(isempty(priorparam)){
	  stop("No prior parameters for Independent Normal-Wishart prior")
	}
	coefprior    <- priorparam[[1]]
	coefpriorvar <- priorparam[[2]]
	varprior     <- priorparam[[3]]
	varpriordof  <- priorparam[[4]]

	pr <- niprior(P,NoLags,RandomWalk=RandomWalk,Intercept=Intercept,coefprior=coefprior,coefpriorvar=coefpriorvar,varprior=varprior)

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

    pr <- mbprior(y=FY,NoLags=NoLags,Intercept=Intercept,RandomWalk=RandomWalk,lambda1=lambda1,lambda2=lambda2,lambda3=lambda3)

    aprior <- pr$aprior
    Vprior <- pr$Vmatrix

  }
  else if(prior == 3){

    # Natural conjugate prior
    if(isempty(priorparam)){

	  stop("No prior parameters for Natural conjugate prior")

    }

    coefprior    <- priorparam[[1]]
    coefpriorvar <- priorparam[[2]]
    varprior     <- priorparam[[3]]
    varpriordof  <- priorparam[[4]]
    print(.isscalar(coefprior))


    pr <- ncprior(P,NoLags,RandomWalk=RandomWalk,Intercept=Intercept,coefprior=coefprior,coefpriorvar=coefpriorvar,varprior=varprior)

    aprior <- pr$coefprior
    Vprior <- pr$coefpriorvar
    vprior <- varpriordof
    Sprior <- pr$varprior

  }
  else if(prior==4){
    # Uninformative prior, do nothing
  }
  else if(prior == 5){
    # Dummy observation prior, to be implemented
    if(isempty(priorparam)){
      stop("No prior parameters for Dummy observation prior")
    }

    tau     <- priorparam[[1]]
    lamda   <- priorparam[[2]]
    epsilon <- priorparam[[3]]

    if(abs(epsilon) > 0 && Intercept == FALSE){
      stop("If intercept is not included, epsilon can not be greater than zero")
    }

    dp <- dummyPrior(y = FY, p = NoLags, tau = tau, epsilon = epsilon, lamda = lamda)
    xd <- dp$xd
    yd <- dp$yd

  }
  else if(prior == 6){
    # SSVS - Prior
    xx <- lagdata(FY,lags=NoLags,intercept=Intercept)
    xlagged <- xx$x
    ylagged <- xx$y

    betaest <- solve(t(xlagged)%*%xlagged)%*%t(xlagged)%*%ylagged
    err     <- (ylagged - xlagged%*%betaest)
    sig     <- t(err)%*%err/(T-NoLags)
    omega   <- array(list(),dim=c(P-1,2))

    for(kk1 in 1:(P-1)){
      omega[[kk1,ii,1]] <- array(1,dim=c(kk1))
      omega[[kk1,ii,2]] <- array(1,dim=c(kk1))
    }

    if(Intercept == TRUE){
      constant = 1
    }
    else{
      constant = 0
    }

    NoRest <- P * (P * NoLags + constant)
    gammas <- array(0.5,dim=c(NoRest))
    si     <- sig%X%(t(xlagged)%*%xlagged)
    tau0   <- 0.01*si
    tau1   <- 1*si
    aprior <- array(0,dim=c(NoRest))
    SSEGibbs <- (array(0,dim=c(P,P,2)))
    SSEGibbs[,,1] <- t(err)%*%err
    SSEGibbs[,,2] <- t(err)%*%err

  }


  # observation equation
  Liprvar <- Lipr*diag(1,P)
  alpha   <- 0.01
  beta    <- 0.01

  #
  # Initialize the Gibbs sampler
  #

  thDelay     <- thMax
  tard        <- seq(1:thMax)
  startest    <- max(thMax,NoLags)
  ytest       <- y[(startest+1-thDelay):(T-thDelay),thVar]
  tarmean     <- mean(ytest)
  tarstandard <- sqrt(var(ytest))
  tart        <- tarmean
  thx         <- thVar+NoFactors

  xsplit <- splitVariables(y=FY,lags=NoLags,thDelay=thDelay,thresh=thx,tart=tart,intercept=Intercept)

  Beta <- array(0,dim=c(P*NoLags+constant,P,2))
  SF   <- array(0,dim=c(P,P,2))

  Beta[,,1] <- solve(t(xsplit$x1) %*% xsplit$x1) %*% t(xsplit$x1) %*% xsplit$y1
  Beta[,,2] <- solve(t(xsplit$x2) %*% xsplit$x2) %*% t(xsplit$x2) %*% xsplit$y2

  SF[,,1]   <- t(xsplit$y1 - xsplit$x1 %*% Beta[,,1]) %*% (xsplit$y1 - xsplit$x1 %*% Beta[,,1])
  SF[,,2]   <- t(xsplit$y2 - xsplit$x2 %*% Beta[,,2]) %*% (xsplit$y2 - xsplit$x2 %*% Beta[,,2])

  #
  # Start MCMC algorithm
  #

  isave <- 0
  for(irep in 1:nreps){
    print(irep)

    # Step 1: split states

    xsplit <- splitVariables(y = FY,lags=NoLags,thDelay=thDelay,thresh=(thVar+NoFactors),tart=tart,intercept=FALSE)

    # Step 2: Sample L and Sigma for both regimes (measurementt equation)

    # regime 1
    nr <- nrow(XY)
    nr2 <- nrow(as.matrix(xsplit$ytest,ncol=1))


    rdiff <- nr-nr2+1
    XYred <- XY[rdiff:nr,]
    XYsplit <- XYred[xsplit$e1,]
    if(priorm == 1){

      # Regime 1
      for(ii in 1:ncol(XYsplit)){
        if(ii>K){

          Lipostvar <- solve(solve(Liprvar)+Sigma[ii,ii,1]^(-1)*t(xsplit$y1)%*%xsplit$y1)
          Lipostmean <- Lipostvar%*%(Sigma[ii,ii,1]^(-1)*t(xsplit$y1)%*%XYsplit[,ii])
          L[ii,1:P,1] <- t(Lipostmean)+rnorm(P)%*%chol(Lipostvar)

        }
        resi <- XYsplit[,ii]-xsplit$y1%*%L[ii,,1]
        sh   <- alpha/2+nrow(xsplit$y1)/2
	      sc   <- beta/2+t(resi)%*%resi/2
	      Sigma[ii,ii,1] <- rgamma(1,shape=sh,scale=sc)
      }

      # Regime 2
      XYsplit <- XYred[!xsplit$e1,]
      for(ii in 1:ncol(XYsplit)){
        if(ii>K){
          Lipostvar <- solve(solve(Liprvar)+Sigma[ii,ii,2]^(-1)*t(xsplit$y2)%*%xsplit$y2)
          Lipostmean <- Lipostvar%*%(Sigma[ii,ii,2]^(-1)*t(xsplit$y2)%*%XYsplit[,ii])
          L[ii,1:P,2] <- t(Lipostmean)+rnorm(P)%*%chol(Lipostvar)
        }

        resi <- XYsplit[,ii]-xsplit$y2%*%L[ii,,2]
        sh   <- alpha/2+nrow(xsplit$y2)/2
        sc   <- beta/2+t(resi)%*%resi/2
        Sigma[ii,ii,2] <- rgamma(1,shape=sh,scale=sc)

      }
    }
    else if(priorm == 2){

      # Regime 1
      XYsplit <- XYred[xsplit$e1,]
      FYred   <- FY[rdiff:nr,]
      for(ii in 1:ncol(XYsplit)){
        if(ii > K){


         # Sample the betas
          VBeta <- diag(gammam[,ii,1] * c2 * tau2 + ( 1 - gammam[,ii,1] ) * tau2)
          DBeta <- solve(t(xsplit$y1)%*%xsplit$y1) * Sigma[ii,ii,1]^(-1)
          dBeta <- t(xsplit$y1)%*%XYsplit[,ii] * Sigma[ii,ii,1]^(-1)
          HBeta <- t(chol(DBeta))
          xTemp <- rnorm(P)%*%HBeta
          L[ii,1:P,1] <- t(DBeta %*% dBeta) + (rnorm(P)%*%HBeta)


          # Sample the gammas
          for(jj in 1:P){
            numerator <- pnorm(L[ii,jj,1], mean = 0, sd = sqrt(c2*tau2))
            denominator <- numerator + pnorm(L[ii,jj,1],mean=0,sd=sqrt(tau2))
            prob <- numerator/denominator
            gammam[jj,ii,1] <- 0.5 * sign(runif(1) - prob) + 0.5
          }

          # Sample the variance
          resid <- XYsplit[,ii] - xsplit$y1 %*% L[ii,,1]
          sh <- alpha/2 + T/2
          sc <- beta/2 + t(resid) %*% resid
          Sigma[ii,ii,1] <- rgamma(1,shape=sh,scale=sc)

        }

      }

      # Regime 2
      XYsplit <- XYred[!xsplit$e1,]
      for(ii in 1:ncol(XYsplit)){
        if(ii > K){

          # Sample the betas
          VBeta <- diag(gammam[,ii,2] * c2 * tau2 + ( 1 - gammam[,ii,2] ) * tau2)
          DBeta <- solve(t(xsplit$y2) %*% xsplit$y2) * Sigma[ii,ii,2]^(-1)
          dBeta <- t(xsplit$y2) %*% XYsplit[,ii] * Sigma[ii,ii,2]^(-1)
          HBeta <- t(chol(DBeta))
          xTemp <- rnorm(P)%*%HBeta
          L[ii,1:P,2] <- t(DBeta %*% dBeta) + (rnorm(P)%*%HBeta)


          # Sample the gammas
          for(jj in 1:P){
            numerator <- pnorm(L[ii,jj,2], mean = 0, sd = sqrt(c2*tau2))
            denominator <- numerator + pnorm(L[ii,jj,2],mean=0,sd=sqrt(tau2))
            prob <- numerator/denominator
            gammam[jj,ii,2] <- 0.5 * sign(runif(1) - prob) + 0.5
          }


          # Sample the variance
          resid <- XYsplit[,ii] - xsplit$y2 %*% L[ii,,2]
          sh <- alpha/2 + T/2
          sc <- beta/2 + t(resid) %*% resid
          Sigma[ii,ii,2] <- rgamma(1,shape=sh,scale=sc)

        }
      }
    }

    # Step 3: sample var coefficients of the state equation


    xsplit <- splitVariables(y=FY,lags=NoLags,thDelay=thDelay,thresh=(thVar+NoFactors),tart=tart,intercept=Intercept)


    if(prior==1){

      # Independent Normal-Wishart prior
      # First regime

      postdraw <- postni(xsplit$y1,xsplit$x1,aprior=aprior,Vprior=Vprior,vprior=vprior,Sprior=Sprior,Sigma=SF[,,1],Intercept=Intercept,stabletest=stabletest,NoLags=NoLags)
      Beta[,,1] <- postdraw$Alpha
      SF[,,1]   <- postdraw$Sigma


      # Second regime
      postdraw <- postni(xsplit$y2,xsplit$x2,aprior=aprior,Vprior=Vprior,vprior=vprior,Sprior=Sprior,Sigma=SF[,,2],Intercept=Intercept,stabletest=stabletest,NoLags=NoLags)
      Beta[,,2] <- postdraw$Alpha
      SF[,,2]   <- postdraw$Sigma

    }
    else if(prior == 2){
      # Minnesota Prior

      # First Regime
      Aols <- solve(t(xsplit$x1)%*%xsplit$x1)%*%t(xsplit$x1)%*%xsplit$y1
      aols <- as.vector(Aols)
      postdraw  <- postmb(y=xsplit$y1,x=xsplit$x1,Vprior=Vprior,aprior=aprior,Sigma=SF[,,1],betaols=aols,stabletest=stabletest,Intercept=Intercept)
      Beta[,,1] <- postdraw$Alpha
      SF[,,1]   <- postdraw$Sigma

      # Second Regime
      Aols <- solve(t(xsplit$x2)%*%xsplit$x2)%*%t(xsplit$x2)%*%xsplit$y2
      aols <- as.vector(Aols)
      postdraw  <- postmb(y=xsplit$y2,x=xsplit$x2,Vprior=Vprior,aprior=aprior,Sigma=SF[,,2],betaols=aols,stabletest=stabletest,Intercept=Intercept)
      Beta[,,2] <- postdraw$Alpha
      SF[,,2]   <- postdraw$Sigma

    }
    else if(prior==3){
      # Natural Conjugate prior
      # First regime
      postdraw <- postnc(xsplit$y1,xsplit$x1,aprior=aprior,Vprior=Vprior
                         ,vprior=vprior,Sprior=Sprior,Sigma=SF[,,1]
                         ,Intercept=Intercept,stabletest=stabletest
                         ,NoLags=NoLags)


      Beta[,,1] <- postdraw$Alpha
      SF[,,1]   <- postdraw$Sigma


      # Second regime

      postdraw <- postnc(xsplit$y2,xsplit$x2,aprior=aprior,Vprior=Vprior
                         ,vprior=vprior,Sprior=Sprior,Sigma=SF[,,2]
                         ,Intercept=Intercept,stabletest=stabletest
                         ,NoLags=NoLags)


      Beta[,,2] <- postdraw$Alpha
      SF[,,2]   <- postdraw$Sigma
    }
    else if(prior==4){

      # uninformative prior

      # First regime

      postdraw <- postun(xsplit$y1,xsplit$x1,Sigma=SF[,,1],Intercept=Intercept
                         ,stabletest=stabletest,NoLags=NoLags)

      Beta[,,1] <- postdraw$Alpha
      SF[,,1]   <- postdraw$Sigma


      # Second regime

      postdraw <- postun(xsplit$y2,xsplit$x2,Sigma=SF[,,2],Intercept=Intercept
                         ,stabletest=stabletest,NoLags=NoLags)

      Beta[,,2] <- postdraw$Alpha
      SF[,,2]   <- postdraw$Sigma
    }
    else if(prior == 5){
      # First regime
      postdraw <- postdummy(y = xsplit$y1, x = xsplit$x1, yd = yd, xd = xd, sigma = SF[,,1], p = NoLags, Intercept = Intercept
                            , stabletest = stabletest)

      Beta[,,1] <- postdraw$Alpha
      SF[,,1]   <- postdraw$Sigma

      # Second regime
      postdraw <- postdummy(y = xsplit$y2, x = xsplit$x2, yd = yd, xd = xd, sigma = SF[,,2], p = NoLags, Intercept = Intercept
                            ,stabletest = stabletest)

      Beta[,,2] <- postdraw$Alpha
      SF[,,2]   <- postdraw$Sigma

    }
    else if(prior == 6){


      # First regime

      aols <- as.vector(solve(t(xsplit$x1)%*%xsplit$x1)%*%t(xsplit$x1)%*%xsplit$y1)
      postdraw <- postss(y=xsplit$y1,x=xsplit$x1,SSEGibbs=SSEGibbs[,,1],omega=omega[,1],
                         gammas=gammas[,1],tau0=tau0,tau1=tau1,Sigma=Sigma[,,1],
                         aprior=aprior,aols=aols,Intercept=Intercept,NoLags=NoLags,
                         stabletest=stabletest)


      Beta[,,1] <- postdraw$Alpha
      SF[,,1] <- postdraw$Sigma
      SSEGibbs[,,1] <- postdraw$SSEGibbs
      gammas[,1] <- postdraw$gammas
      omega[,1] <- postdraw$omega



      # Second regime

      aols <- as.vector(solve(t(xsplit$x2)%*%xsplit$x2)%*%t(xsplit$x2)%*%xsplit$y2)

      postdraw <- postss(y=xsplit$y2,x=xsplit$x2,SSEGibbs=SSEGibbs[,,2],omega=omega[,2],
                         gammas=gammas[,2],tau0=tau0,tau1=tau1,Sigma=Sigma[,,2],
                         aprior=aprior,aols=aols,Intercept=Intercept,NoLags=NoLags,
                         stabletest=stabletest)


      Beta[,,2] <- postdraw$Alpha
      SF[,,2] <- postdraw$Sigma


      SSEGibbs[,,2] <- postdraw$SSEGibbs
      gammas[,2] <- postdraw$gammas
      omega[,2] <- postdraw$omega
    }


    # Step 4: sample new threshold
    tarnew <- tart+rnorm(1,sd=tarstandard)
    l1post <- tarpost(xsplit$xstar,xsplit$ystar,Ytest=ytest,Beta[,,1],Beta[,,2]
                      ,SF[,,1],SF[,,2],tarnew,NoLags,intercept=Intercept
                      ,tarmean,tarstandard,ncrit=ncrit)


    l2post <- tarpost(xsplit$xstar,xsplit$ystar,Ytest=ytest,Beta[,,1],Beta[,,2]
                      ,SF[,,1],SF[,,2],tart,NoLags,intercept=Intercept
                      ,tarmean,tarstandard,ncrit=ncrit)


    acc <- min(1,exp(l1post$post-l2post$post))
    u <- runif(1)
    if(is.na(acc)){acc = 0}


    if(u<acc){
      tart=tarnew
    }
    #tarmean=tart
    # Step 5: Sample new delay parameter
    prob <- matrix(0,nrow=thMax)
    for(jj in 1:thMax){

      split1 <- splitVariables(y=FY,lags=NoLags,jj,thVar+NoFactors,tart,intercept=Intercept)
      x <- exptarpost(split1$xstar,split1$ystar,split1$ytest,Beta[,,1],Beta[,,2],SF[,,1],SF[,,2],
                      tart,NoLags,intercept=Intercept,tarmean,tarstandard,ncrit=ncrit)
      prob[jj,1] <- x$post

    }
    #print(prob)
    mprob <- max(prob)
    prob <- exp(prob-mprob)
    prob <- prob/sum(prob)


    if(anyNA(prob)){
      prob <- matrix(1/thMax,nrow=thMax)
    }
    thDelay <- sample(thMax,1,replace=FALSE,prob)


    # Store results after burnin-period

    xmod <- mod(irep-burnin,thin)

    if(irep>burnin && xmod == 0){
      isave <- isave + 1

      # Store draws for Beta and Sigma and L

      Sigmadraws[,,,isave] <- SF
      Alphadraws[,,,isave] <- Beta
      Ldraws[,,,isave] <- L
      tardraws[isave] <- tart
      deldraws[isave] <- thDelay

      if(priorm==2){

        gammamdraws[,,1,isave] <- gammam[,,1]
        gammamdraws[,,2,isave] <- gammam[,,2]

      }


      # Regimes

      nT <- length(xsplit$e1)
      a  <- nT-NoRegimes
      regimes[,isave] <- xsplit$e1[(1+a):nT]

	  }	# end storing results
  } # End loop over mcmc algorithm


  # Prepare return values

	retlist <- structure(list(Betadraws = Alphadraws,Sigmadraws = Sigmadraws,Ldraws = Ldraws, deldraws = deldraws,
	                          regimes = regimes,tardraws = tardraws, gammam = gammamdraws, NoLags = NoLags,mydata = FY,
	                          Intercept = Intercept, thVar = thVar,varnames = varnames),class="ftvar")

	return(retlist)


}
