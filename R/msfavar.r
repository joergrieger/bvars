msfavar <- function(y,x,NoFactors=1,NoLags=2,slowindex="",frotvar,RandomWalk=TRUE,
                    coefprior=NULL,prior=1,priorparam=NULL,Lipr=4,NoRegimes=2,
                    alphaprior=5,nreps=200,burnin=100,irfhor=40,ident=1,
                    restrictions=NULL,Intercept=FALSE,irfQuantiles = c(0.16,0.84),stabletest=TRUE){
  #
  # Preliminaries
  #
  xdata <- scale(as.matrix(x))
  ydata <- scale(as.matrix(y))
  #x.demeaned <- demean(xdata)
  #y.demeaned <- demean(ydata)
  #x <- x.demeaned
  #y <- y.demeaned



  T <- nrow(y)
  N <- ncol(x)
  K <- ncol(y)
  P <- K+NoFactors

  #
  # Check if input is correct
  #

  #
  # Declare Variables for Storage
  #

  Alphadraws <- array(0,dim=c(P*NoLags,P,NoRegimes,nreps-burnin))
  Sigmadraws <- array(0,dim=c(P,P,NoRegimes,nreps-burnin))
  irfdraws   <- array(0,dim=c(P,P,irfhor,NoRegimes,nreps-burnin))
  Ldraws     <- array(0,dim=c(ncol(y)+ncol(x),P,NoRegimes,nreps-burnin))
  sttdraws   <- array(0,dim=c(T-NoLags,nreps-burnin))
  irfSmallDraws <- array(0,dim=c(P,P,irfhor,NoRegimes,nreps-burnin))
  irfLargeDraws <- array(0,dim=c(P,(ncol(x)+ncol(y)),irfhor,NoRegimes,nreps-burnin))

  # Declare Variables for final storage

  irfDraw  <- array(0,dim=c(P,(ncol(x)+ncol(y)),irfhor,NoRegimes,3))

  #
  # Extract principal components
  #

  fac <- exfact(ydata=ydata,xdata=xdata,slowcode=slowindex,NoFactors=NoFactors)

  # put it in state-space form

  XY <- cbind(ydata,xdata)
  FY <- cbind(ydata,fac)
  Li <- olssvd(XY,FY)

  e   <- XY - FY%*%Li
  Sig <- t(e)%*%e/T

  L <- array(0,dim=c(ncol(XY),P,NoRegimes))
  Sigma <- array(0,dim=c(ncol(XY),ncol(XY),NoRegimes))

  for(ii in 1:NoRegimes){
    L[,,ii] <- Li
    Sigma[,,ii] <- Sig
  }

  #
  # Set priors
  #

  # VAR-coefficients
  if(prior==1){
    # Independent Normal-Wishart prior

    if(isempty(priorparam)){

      stop("No prior parameters for Independent Normal-Wishart prior")
    }


    coefprior    <- priorparam[[1]]
    coefpriorvar <- priorparam[[2]]
    varprior     <- priorparam[[3]]
    varpriordof  <- priorparam[[4]]
    pr <- niprior(P,NoLags = NoLags,RandomWalk = RandomWalk,Intercept = FALSE,coefprior = coefprior,coefpriorvar = coefpriorvar,varprior = varprior,varpriordof = varpriordof)

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
  else if(prior==3){
    # Natural conjugate prior
	if(isempty(priorparam)){
	  stop("No prior parameters for Natural-conjugate prior")
	}
	coefprior    <- priorparam[[1]]
	coefpriorvar <- priorparam[[2]]
	varprior     <- priorparam[[3]]
	varpriordof  <- priorparam[[4]]

    pr <- ncprior(P,NoLags = NoLags,RandomWalk = RandomWalk,Intercept = FALSE,coefprior = coefprior,coefpriorvar = coefpriorvar,varprior = varprior,varpriordof = varpriordof)

    aprior <- as.vector(pr$coefprior)
    Vprior <- pr$coefpriorvar
    vprior <- varpriordof
    Sprior <- pr$varprior
  }
  else if(prior==4){
    # Uninformative prior
  }


  # observation equation

  Liprvar <- Lipr*diag(1,P)
  alpha <- 0.01
  beta  <- 0.01

  # Prior on transition-probability matrix
  Alphaprior <- alphaprior*matrix(1,NoRegimes,NoRegimes)
  Stransmat <- matrix(0,NoRegimes,NoRegimes)

  #
  # Initialize MCMC algorithm
  #

  Beta <- array(0,dim=c(P*NoLags,P,NoRegimes))
  SF   <- array(0,dim=c(P,P,NoRegimes))

  FYlagged <- lagdata(FY,NoLags,intercept=FALSE)
  XYred <- XY[(NoLags+1):nrow(XY),]
  FYred <- FY[(NoLags+1):nrow(XY),]
  FYx <- FYlagged$x
  FYy <- FYlagged$y

  Phi  <- solve(t(FYx)%*%FYx)%*%t(FYx)%*%FYy
  resi <- FYy-FYx%*%Phi
  sse  <- t(resi)%*%resi/T
  for(ii in 1:NoRegimes){
    Beta[,,ii] <- Phi
    SF[,,ii] <- sse
  }

  #
  # Start MCMC Algorithm
  #
  for(irep in 1:nreps){
    print(irep)

    # Step 1: Sample states
    transmat      <- t(apply(Alphaprior+Stransmat-1,1,rdirichlet,n=1))
    filteredprob  <- hamiltonfilter(Beta,SF,transmat,FYy,FYx,h=NoRegimes)
    stt           <- getst(filteredprob$fprob,transmat,h=NoRegimes)
    Stransmat     <- countseq(stt,h=NoRegimes)

    # Next step: Inference in the FAVAR state-by-state
    for(ireg in 1:NoRegimes){
      FYreg <- FYred[stt==ireg,]
      XYreg <- XYred[stt==ireg,]
      FYxreg <- FYx[stt==ireg,]
      FYyreg <- FYy[stt==ireg,]

      # Step 2: Sample L and  Sigma
      for(ii in 1:ncol(XY)){
        if(ii>K){

          Lipostvar <- solve(solve(Liprvar)+Sigma[ii,ii,ireg]^(-1)*t(FYreg)%*%FYreg)
          Lipostmean <-Lipostvar%*%( Sigma[ii,ii,ireg]^(-1)*t(FYreg)%*%XYreg[,ii])
          L[ii,1:P,ireg] <- t(Lipostmean)+rnorm(P)%*%chol(Lipostvar)

        }
        resi <- XYreg[,ii]-FYreg%*%L[ii,,ireg]
        sh <- alpha/2+nrow(FYreg)/2
        sc <- beta/2+t(resi)%*%resi/2
        Sigma[ii,ii,ireg] <- rgamma(1,shape=sh,scale=sc)
      }

      # Step 3: Sample VAR coefficients
      if(prior==1){

        postdraw <- postni(y=FYyreg,x=FYxreg,aprior,Vprior,vprior,Sprior,Sigma=SF[,,ireg],stabletest=TRUE,Intercept=Intercept,NoLags=NoLags)

      }
      else if(prior == 2){

        Aols <- solve(t(FYxreg)%*%FYxreg)%*%t(FYxreg)%*%FYyreg
        aols <- as.vector(Aols)
        postdraw  <- postmb(y=FYyreg,x=FYxreg,Vprior=Vprior,aprior=aprior,Sigma=SF[,,ireg],betaols=aols,stabletest=stabletest,Intercept=Intercept)


      }
      else if(prior == 3){

        postdraw <- postnc(y=FYyreg,x=FYxreg,aprior,Vprior,vprior,Sprior,Sigma=SF[,,ireg],stabletest=TRUE,Intercept=Intercept,NoLags=NoLags)

      }
      else if(prior==4){

        postdraw <- postun(y=FYyreg,x=FYxreg,Sigma=SF[,,ireg],stabletest=TRUE,Intercept=Intercept,NoLags=NoLags)

      }
      # Order by means
      #sorder <- order(Sigma[1,1,])
      #Sigma <- Sigma[,,rev(sorder)]
      #Alpha <- Alpha[,,rev(sorder)]

      Beta[,,ireg] <- postdraw$Alpha
      SF[,,ireg]   <- postdraw$Sigma
    } # End Loop over regimes

    # Identify Regimes
    idVar = 3
    idRegimes <- array(0,dim=c(NoRegimes))

    for(ij in 1:NoRegimes){

      idRegimes[ij] <- mean(FYred[stt==ij,idVar])

    }
    sorder <- order(idRegimes)
    Beta <- Beta[,,rev(sorder)]
    SF   <- SF[,,rev(sorder)]


    # compute Impulse-Response functions and save variables
    if(irep>burnin){
      for(ireg in 1:NoRegimes){
        # Save values from draws
        Alphadraws[,,ireg,irep-burnin] <- Beta[,,ireg]
        Sigmadraws[,,ireg,irep-burnin] <- SF[,,ireg]

        Ldraws[,,ireg,irep-burnin] <- L[,,ireg]

        # compute and save impulse-response functions
        sttdraws[,irep-burnin] <- stt

        # Compute Impulse-response functions

        if(ident==1){

          # identification of structural shocks using cholesky decomposition
          irf <- compirf(A=Beta[,,ireg],Sigma=SF[,,ireg],NoLags=NoLags,intercept=Intercept,nhor = irfhor)

        }
        else if(ident==2){

          # Identification of structural shocks using sign restrictions
          irf <- compirf(A=Beta[,,ireg],Sigma=SF[,,ireg],NoLags=NoLags,intercept=FALSE,nhor = irfhor,restrictions=restrictions)

        }
        irfSmallDraws[,,,ireg,irep-burnin] <- irf
        # Compute IRFs for all variables approximated by factors
        for(ii in 1:P){

          irfLargeDraws[ii,,,ireg,irep-burnin] <- L[,,ireg]%*%irf[ii,,]

        }


      }

    }  # End Saving

  } # End loop over Gibbs Sampler

  # Final Calculations for impulse-response functions

  for(ii in 1:NoRegimes){
    for(jj in 1:irfhor){
      for(kk in 1:P){
        for(ll in 1:(N+K)){
          irfDraw[kk,ll,jj,ii,1] <- median(irfLargeDraws[kk,ll,jj,ii,])
          irfDraw[kk,ll,jj,ii,2] <- quantile(irfLargeDraws[kk,ll,jj,ii,],probs=min(irfQuantiles))
          irfDraw[kk,ll,jj,ii,3] <- quantile(irfLargeDraws[kk,ll,jj,ii,],probs=max(irfQuantiles))
        }
      }
    }

  }
  retList <- list(Beta=Alphadraws,Sigma=Sigmadraws,stt=sttdraws,irfLarge=irfDraw)
  return(retList)
}
