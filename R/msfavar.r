msfavar <- function(y,x,NoFactors=1,NoLags=2,slowindex,frotvar,RandomWalk=TRUE,coefprior=NULL,coefpriorvar=1,varprior=1,varpriordof=10,Lipr=4,NoRegimes=2,alphaprior=5,nreps=200,burnin=100,irfhor=40){
  #
  # Preliminaries
  #
  xdata <- as.matrix(x)
  ydata <- as.matrix(y)
  x.demeaned <- demean(xdata)
  y.demeaned <- demean(ydata)
  x <- x.demeaned
  y <- y.demeaned
  
  
  
  T <- nrow(y)
  N <- ncol(x)
  K <- ncol(y)
  P <- K+NoFactors
  
  for(ii in 1:N){
    x[,ii] <- x.demeaned[,ii]/sqrt(var(x[,ii]))
  }
  
  
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
  #
  # Extract principal components
  #
  
  fac <- exfact(ydata=y,xdata=x,slowcode=slowindex,NoFactors=NoFactors)
  
  # put it in state-space form
  
  XY <- cbind(x,y)
  FY <- cbind(fac,y)
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
  pr <- niprior(P,NoLags = NoLags,RandomWalk = RandomWalk,Intercept = FALSE,coefprior = coefprior,coefpriorvar = coefpriorvar,varprior = varprior,varpriordof = varpriordof)
  
  aprior <- as.vector(pr$coefprior)
  Vprior <- pr$coefpriorvar
  vprior <- varpriordof
  Sprior <- pr$varprior
  
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
    transmat <-  t(apply(Alphaprior+Stransmat-1,1,rdirichlet,n=1))
    filteredprob <- hamiltonfilter(Beta,SF,transmat,FYy,FYx,h=NoRegimes)
    stt <- getst(filteredprob$fprob,transmat,h=NoRegimes)
    Stransmat <- countseq(stt,h=NoRegimes)
    #print(Stransmat)
    #readline(prompt="Press [enter] to continue")
    
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
      
      postdraw <- postni(y=FYyreg,x=FYxreg,aprior,Vprior,vprior,Sprior,Sigma=SF[,,ireg])
      
      Beta[,,ireg] <- postdraw$Alpha
      SF[,,ireg] <- postdraw$Sigma
      
      # compute Impulse-Response functions and save variables
      if(irep>burnin){
        # Save values from draws
        Alphadraws[,,ireg,irep-burnin] <- Beta[,,ireg]
        Sigmadraws[,,ireg,irep-burnin] <- SF[,,ireg]
        
        Ldraws[,,ireg,irep-burnin] <- L[,,ireg]
        
        # compute and save impulse-response functions
        irf <- compirf(A=Beta[,,ireg],Sigma=SF[,,ireg],NoLags=NoLags,intercept=FALSE,nhor = irfhor)
        irfSmallDraws[,,,ireg,irep-burnin] <- irf
        for(ii in 1:P){
          irfLargeDraws[ii,,,ireg,irep-burnin] <- L[,,ireg]%*%irf[ii,,]
        }
        

      }
    }
    
    #
    # Save draws
    #
    if(irep>burnin){
      sttdraws[,irep-burnin] <- stt
    }
    
  }
  
}