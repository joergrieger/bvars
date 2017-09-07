nifavar <- function(ydata,xdata,NoFactors=2,NoLags=1,RandomWalk=TRUE,coefpriorvar=10,varprior=10,varpriordof=20,alpha=0.01,beta=0.01,lipr=4,method="PC",slowindex="",frotvar=NaN,nreps=11000,burnin=1000,irfhor=40,irfquantiles=c(0.05,0.95)){
  
  xdata <- as.matrix(xdata)
  ydata <- as.matrix(ydata)
  
  y.demeaned <- demean(ydata)
  x.demeaned <- demean(xdata)

  X <- x.demeaned
  Y <- y.demeaned
  
  T <- nrow(Y)
  N <- ncol(X)
  M <- ncol(Y)
  P <- M+NoFactors
  #
  # Sanity checks for input variables
  #
  # Check if both time series have the same length
  if(nrow(Y) != nrow(X)){
    stop("Both time series must have the same length")
  }
  # check for factor rotation
  if(!is.nan(frotvar) && !isempty(slowcode) && frotvar>ncol(Y)){
    stop("invalid variables for factor rotation")
  }
  
  # Extract factors
  # The user can choose between 2 methods:
  # 1. Extraction of Factors using Principal Components
  # 2. Extraction of Factors using MCMC methods
  
  # Extract factors using principal components
  if(method=="PC"){
    # Standardize
    for(i in 1:N){
      X[,i]<-X[,i]/sqrt(var(X[,i]))
    }
    factors <- extract(X,NoFactors)
    # Add Factor rotations
    if(!isempty(slowcode)){
      slow=which(slowindex==1)
      xslow <- X[,slow]
      ftemp <- extract(xslow,NoFactors)
      fnr <- ncol(ftemp$fac)
	  if(is.nan(frotvar)){
	    frotvar=ncol(Y)
	  }
      Fr0 <- facrot(factors$fac,Y[,frotvar],ftemp$fac)

    }
    else {
      Fr0 <- factors$fac
    }
  }
  else if(method=="MCMC"){
    print("Not implemented")
  }
  #
  # Now, put it into State-Space representation
  #
  
  XY <- cbind(X,Y)
  FY <- cbind(Fr0,Y)
  L <- olssvd(XY,FY)
  

  e <- XY - FY%*%L
  SIGMA <- t(e)%*%e/T
  # VAR equations
  FY_lagged <- .data(FY,NoLags)
  FY.X <- FY_lagged$x
  FY.Y <- FY_lagged$y
  
  PHI <- solve(t(FY.X)%*%FY.X)%*%(t(FY.X)%*%FY.Y)
  SSE <- t(FY.Y-FY.X%*%PHI)%*%(FY.Y-FY.X%*%PHI)
  S_F <- SSE/(T-NoLags)
  # Put everything into a VAR(1)-form
  if(NoLags>1){
    S_F_Mat <- cbind(S_F,matrix(0,P,P*(NoLags-1)))
    S_F_Mat <- rbind(S_F_Mat,P*(NoLags-1),P*NoLags)
    xx <- cbind(diag(1,P*(NoLags-1)),matrix(0,P*(NoLags-1),P))
    PHI_Mat <- rbind(t(PHI),xx)
  }
  else{
    S_F_Mat <- S_F
    PHI_mat <- t(PHI)
  }
  # Setting Priors for Factor Loadings
  Li_prvar <- lipr*diag(1,P)
  
  
  # Setting Priors for the Independent Normal-Wishart prior
  coefprior <- array(0,dim=c(P*NoLags,P))
  if(RandomWalk==TRUE){
    for(ii in 1:P){
      coefprior[ii,ii]<-1
    }
  }

  nrcoefpriorvar <- P*(P*NoLags)
  coefpriorvar1 <- coefpriorvar*.id(nrcoefpriorvar)
  varprior1 <- varprior*.id(P)
  varpriordof1 <- T+varpriordof
  
  sampledResults <- .fnimcmc(FY=FY,X=X,Xlag=FY.X,Ylag=FY.Y,NoFactors=NoFactors,L=t(L),P=P,NoLags=NoLags,nreps=nreps,burnin=burnin,Li_prvar=Li_prvar,SIGMA=SIGMA,alpha=alpha,beta=beta,
                              coefprior=coefprior,coefpriorvar=coefpriorvar1,varprior=varprior1,varpriordof=varpriordof1,nhor=irfhor)
  irffinal <- array(0,dim=c((ncol(Y)+NoFactors),(ncol(Y)+ncol(X)),irfhor,3))
  lower_quantile <- min(irfquantiles)
  upper_quantile <- max(irfquantiles)
  for(ii in 1:irfhor){
    for(jj in 1:(ncol(X)+ncol(Y))){
      for(kk in 1:(ncol(Y)+NoFactors)){
        irffinal[kk,jj,ii,2] <- quantile(sampledResults$irf[,kk,jj,ii],probs=lower_quantile)
        irffinal[kk,jj,ii,3] <- quantile(sampledResults$irf[,kk,jj,ii],probs=upper_quantile)
        irffinal[kk,jj,ii,1] <- mean(sampledResults$irf[,kk,jj,ii])
      }
    }
  }
  return(betadraws=sampledResults$Betadraws,sigmadraws=sampledResults$Sigmadraws,irf=irffinal)
}

.fnimcmc <- function(FY,X,Xlag,Ylag,NoFactors=1,L,P,NoLags=1,nreps,burnin,Li_prvar,SIGMA,alpha,beta,coefprior,coefpriorvar,varprior,varpriordof,nhor=40){
  # Preliminaries
  N <- ncol(X)
  T <- nrow(FY)
  K <- ncol(Ylag)
  bigj <- matrix(0,P,P*NoLags)
  bigj[1:P,1:P] <- diag(1,P)
  imp <- array(0,dim=c(nreps-burnin,K,ncol(X)+ncol(Ylag)-NoFactors,nhor))
  Betadraws  <- array(0,dim=c(K*NoLags,K,nreps-burnin))
  Sigmadraws <- array(0,dim=c(K,K,nreps-burnin))
  # OLS estimates
  PHI <- solve(t(Xlag)%*%Xlag)%*%(t(Xlag)%*%Ylag)
  SSE <- t(Ylag-Xlag%*%PHI)%*%(Ylag-Xlag%*%PHI)
  S_F <- SSE/(T-NoLags)
  if(NoLags>1){
    S_F_Mat <- cbind(S_F,matrix(0,P,P*(NoLags-1)))
    S_F_Mat <- rbind(S_F_Mat,P*(NoLags-1),P*NoLags)
    xx <- cbind(diag(1,P*(NoLags-1)),matrix(0,P*(NoLags-1),P))
    PHI_Mat <- rbind(t(PHI),xx)
  }
  else{
    S_F_Mat <- S_F
    PHI_Mat <- t(PHI)
  }
  # Preliminaries for the Gibbs Sampler
  aprior <- .vec(coefprior)
  Vprior <- coefpriorvar
  vpior  <- varpriordof
  Sprior <- varprior
  invSprior <- solve(Sprior)
  Alpha <- PHI
  z <- .id(P)%x%Xlag
  for(irep in 1:nreps){
    print(irep)
    # Step 1: Sample L and SIGMA
    for(i in 1:N){
      if(i>K){
        Li_postvar   <- solve(solve(Li_prvar)+(SIGMA[i,i])^(-1)*t(FY)%*%FY)
        Li_postmean  <- Li_postvar%*%((SIGMA[i,i])^(-1)*t(FY)%*%X[,i])
        Lidraw <- t(Li_postmean)+rnorm(P)%*%chol(Li_postvar)
        L[i,1:P]<-Lidraw
      }
      ed <- X[,i]-FY%*%L[i,]
      S_1 <- alpha/2+T/2
      S_2 <- beta/2+t(ed)%*%ed/2
      Sidraw <- rgamma(1,shape=S_1,scale=1/S_2)
      SIGMA[i,i] <- Sidraw
    }
    # Step 2: Draw VAR coefficients using an Normal-Independent Wishart Prior
    variance <- solve(S_F)%x%.id(nrow(Ylag))
    Vpost    <- solve(solve(Vprior)+t(z)%*%variance%*%z)
    apost    <- Vpost%*%(solve(Vprior)%*%aprior+t(z)%*%variance%*%.vec(Ylag))
    alpha    <- mvrnorm(mu=apost,Sigma=Vpost)
    Alpha1   <- matrix(alpha,ncol=K)
    # Draw Variance
    vpost    <- T+Vprior
    Spost    <- solve(invSprior+t(Ylag-Xlag%*%Alpha)%*%(Ylag-Xlag%*%Alpha))
    S_F1     <- solve(rWishart(1,vpost,Spost)[,,1]) # Draw Variance from a Wishart distribution
    Alpha    <- Alpha1
    S_F      <- S_F1
    # Save draws
    if(irep>burnin){
      Betadraws[,,irep-burnin] <- Alpha
      Sigmadraws[,,irep-burnin] <- S_F
      # Calculate Impulse-Response Function
      if(NoLags>1){
        S_F_Mat <- cbind(S_F,matrix(0,P,P*(NoLags-1)))
        S_F_Mat <- rbind(S_F_Mat,P*(NoLags-1),P*NoLags)
        xx <- cbind(diag(1,P*(NoLags-1)),matrix(0,P*(NoLags-1),P))
        PHI_Mat <- rbind(t(PHI),xx)
      }
      else{
        S_F_Mat <- S_F
        PHI_mat <- t(PHI)
      }
      biga <- PHI_Mat
      bigai <- biga
      shock <- t(chol(S_F))
      impresp <- matrix(0,P,P*nhor)
      impresp[1:P,1:P]<-shock
      for(j in 1:(nhor-1)){
        impresp[,(j*P+1):((j+1)*P)] <- (bigj%*%bigai%*%t(bigj)%*%shock)
        bigai <- bigai%*%biga
      }
      imp_m <- array(0,dim=c(P,P,nhor))
      jj <- 0
      for(ii in 1:P){
        for(ij in 1:nhor){
          #jj <- jj+P
          jj <- ii + P*(ij-1)
          imp_m[ii,,ij] <- impresp[,jj]
        }
      imp[irep-burnin,ii,,]<-L%*%imp_m[ii,,]
      }
    } # Save IRFs
    
  } # End loop over Gibbs sampler
  return(list(Betadraws=Betadraws,Sigmadraws=Sigmadraws,irf=imp))
}