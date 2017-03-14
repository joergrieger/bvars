#' @importFrom stats embed lm pnorm quantile rWishart rnorm runif var
#' @export

tvar <- function(mydata,lags=1,thDelay=1,thresh=1,tarscale=0.5,tarstandard=NULL,intercept=TRUE,coefprior=NULL,coefpriorvar=1,varprior=1,varpriordof=1,irfhorizon=16,irfquantiles=c(0.05,0.95),reps=300,burnin=100,stabletest=FALSE){
  y <- as.matrix(mydata)
  T <- nrow(y)
  K <- ncol(y)
  obs <- T-max(thDelay,lags)
  .tierror(mydata,lags=1,thDelay,thresh,tarscale,tarstandard,intercept,coefprior,coefpriorvar,varprior,varpriordof,irfhorizon,irfquantiles,reps,burnin,stabletest)
  prior <- .tiprior(y,lags,thDelay,thresh,tarscale,tarstandard,intercept,coefprior,coefpriorvar,varprior,varpriordof,irfhorizon,irfquantiles,reps,burnin,stabletest)
  results <- .tigibbs(y,lags,thDelay,thresh,tarscale,tarstandard,intercept,prior$Aprior1,prior$Aprior2,prior$Vprior1,prior$Vprior2,prior$Sprior1,prior$Sprior2,varpriordof,irfhorizon,irfquantiles,reps,burnin,stabletest,prior)
  retresults <-  structure(list(prior=prior,results=results),class="tvar")
}
.tierror <-function(mydata,lags=1,thDelay=1,thresh=1,tarscale=0.5,tarstandard=NULL,intercept=TRUE,coefprior=NULL,coefpriorvar=1,varprior=1,varpriordof=1,irfhorizon=16,irfquantiles=c(0.05,0.95),reps=300,burnin=100,stabletest=FALSE){
  if(thDelay > lags){stop("Number of lags has to be greater than the delay")}

}

.tiprior<-function(y,lags,thDelay,thresh,tarscale,tarstandard,intercept,coefprior,coefpriorvar,varprior,varpriordof,irfhorizon,irfquantiles,reps,burnin,stabletest){
  T <- nrow(y)
  K <- ncol(y)
  obs <- T-max(lags,thDelay)
  constant=0
  if(intercept==TRUE){constant=1}
  if(is.null(coefprior)){
    Aprior1 <- array(0,dim=c(K*lags+constant,K))
    Aprior2 <- array(0,dim=c(K*lags+constant,K))
    aprior1 <- matrix(Aprior1,ncol=1)
    aprior2 <- matrix(Aprior2,ncol=1)
  } else {
    Aprior1=coefprior
    Aprior2=coefprior
  }

  if(.isscalar(coefpriorvar)){
    n<-K*(K*lags+constant)
    Vprior1 <- coefpriorvar*.id(n)
    Vprior2 <- coefpriorvar*.id(n)
  } else {
    Vprior1 <- coefpriorvar
    Vprior2 <- coefpriorvar
  }
  if(.isscalar(varprior)){
    Sprior1 <- varprior*.id(K)
    Sprior2 <- varprior*.id(K)
  } else {
    Sprior1 <- varprior
    Sprior2 <- varprior
  }
  return(list(Aprior1=Aprior1,Aprior2=Aprior2,Vprior1=Vprior1,Vprior2=Vprior2,Sprior1=Sprior1,Sprior2=Sprior2))
}

.tigibbs <- function(y,lags,thDelay,thresh,tarscale,tarstandard,intercept,Aprior1,Aprior2,Vprior1,Vprior2,Sprior1,Sprior2,varpriordof,irfhorizon,irfquantiles,reps,burnin,stabletest,prior){
  # vectorize prior
  aprior1 <- matrix(Aprior1,ncol=1)
  aprior2 <- matrix(Aprior2,ncol=1)
  # Some preliminary Calculations
  startest <- max(thDelay,lags)
  T <- nrow(y)
  K <- ncol(y)
  constant=0;
  if(intercept==TRUE) constant=1

  ytest <- y[(startest+1-thDelay):(T-thDelay),thresh]
  ystar <- y[(startest+1):T,]
  tarmean <- mean(ytest)
  tart <- tarmean
  xstar <- embed(y,dimension=lags+1)[,-(1:K)]
  # Initialize the Gibbs Sampler
  e1 <- ytest < tart
  e2 <- ytest >=tart
  y1 <- ystar[e1,]
  y2 <- ystar[e2,]
  x1 <- xstar[e1,]
  x2 <- xstar[e2,]
  if(intercept==TRUE){
    x1 <- cbind(1,x1)
    x2 <- cbind(1,x2)
  }
  # OLS estimation
  beta01 <- solve(t(x1)%*%x1)%*%t(x1)%*%y1
  beta02 <- solve(t(x2)%*%x2)%*%t(x2)%*%y2

  SSE1 <- t(y1-x1%*%beta01)%*%(y1-x1%*%beta01)
  SSE2 <- t(y2-x2%*%beta02)%*%(y2-x2%*%beta02)
  obs1 <- nrow(y1)
  obs2 <- nrow(y2)
  SIGMA1Ols <- SSE1/obs1
  SIGMA2Ols <- SSE2/obs2
  SIGMA1 <- SIGMA1Ols
  SIGMA2 <- SIGMA2Ols
  ALPHAdraws1 <- array(0,dim=c(K*lags+constant,K,reps-burnin))
  ALPHAdraws2 <- array(0,dim=c(K*lags+constant,K,reps-burnin))
  SIGMAdraws1 <- array(0,dim=c(K,K,reps-burnin))
  SIGMAdraws2 <- array(0,dim=c(K,K,reps-burnin))
  tardraws    <- array(0,dim=c(reps-burnin))
  irf1draws   <- array(0,dim=c(K,irfhorizon,K,reps-burnin))
  irf2draws   <- array(0,dim=c(K,irfhorizon,K,reps-burnin))
  if(is.null(tarstandard)){
    tarstandard <- tarscale*sqrt(var(ytest))
  }
  # Do the Gibbs Sampling
  for(ii in 1:reps){
    cat(ii)
    cat("\n")
    # Step 1: Split the Data
    e1 <- ytest<tart
    e2 <- ytest>=tart
    y1 <- ystar[e1,]
    y2 <- ystar[e2,]
    x1 <- xstar[e1,]
    x2 <- xstar[e2,]
    if(intercept==TRUE){
      x1 <- cbind(1,x1)
      x2 <- cbind(1,x2)
    }
    obs1 <- nrow(y1)
    obs2 <- nrow(y2)
    # Step 2: get coefficients
    # Step 2a: Coefficients for the first regime
    z1 <- .id(K)%x%x1
    variance1 <- solve(SIGMA1)%x%.id(obs1)
    Vpost01   <- solve(solve(Vprior1)+t(z1)%*%variance1%*%z1)
    apost01    <- Vpost01%*%(solve(Vprior1)%*%aprior1+t(z1)%*%variance1%*%matrix(y1,ncol=1))
    stable <- 2

    while(stable>1){
      alpha01    <- mvrnorm(mu=apost01,Sigma=Vpost01)
      Alpha01    <- matrix(alpha01,ncol=K)
      if(intercept==TRUE){
        Alphatest <- Alpha01[2:(K*lags+1),]
      }
      else{
        Alphatest <- Alpha01
      }
      stable <- .stability(Alphatest,lags,K)
      if(stabletest==FALSE){stable=0}
    }
    vpost01    <- obs1+varpriordof
    Spost01    <- Sprior1+t(y1-x1%*%Alpha01)%*%(y1-x1%*%Alpha01)
    SIGMA1     <- solve(rWishart(1,vpost01,Spost01)[,,1])
    ALPHA1 <- Alpha01

    # Step 2b: Coefficients for the second regime
    z2 <- .id(K)%x%x2
    variance2 <- solve(SIGMA2)%x%.id(obs2)
    Vpost02   <- solve(solve(Vprior2)+t(z2)%*%variance2%*%z2)
    apost02   <- Vpost02%*%(solve(Vprior2)%*%aprior2+t(z2)%*%variance2%*%matrix(y2,ncol=1))
    stable <- 2
    while(stable>1){
      alpha02 <-mvrnorm(mu=apost02,Sigma=Vpost02)
      Alpha02 <-matrix(alpha02,ncol=K)
      if(intercept==TRUE){
        Alphatest <- Alpha02[2:(K*lags+1),]
      }
      else{
        Alphatest <- Alpha02
      }
      stable <- .stability(Alphatest,lags,K)
      if(stabletest==FALSE){stable=0}
    }
    vpost02 <- obs2+varpriordof
    Spost02 <- Sprior2+t(y2-x2%*%Alpha02)%*%(y2-x2%*%Alpha02)
    SIGMA2  <- solve(rWishart(1,vpost02,Spost02)[,,1])
    ALPHA2  <- Alpha02

    # Step 3: Draw new
    tarnew <- tart+rnorm(1,sd=tarstandard)
    l1post <- .tarpost(xstar,Ystar=ystar,Ytest=ytest,ALPHA1,ALPHA2,SIGMA1,SIGMA2,tarnew,lags,intercept=intercept,tarmean,tarstandard,ncrit=0.15)
    l2post <- .tarpost(xstar,Ystar=ystar,Ytest=ytest,ALPHA1,ALPHA2,SIGMA1,SIGMA2,tart,lags,intercept=intercept,tarmean,tarstandard,ncrit=0.15)
    acc=min(1,exp(l1post$post-l2post$post))
    u <- runif(1)
    if(u<acc){
      tart=tarnew
    }
    tarmean=tart
    if(ii>burnin){
      ALPHAdraws1[,,ii-burnin]<-ALPHA1
      ALPHAdraws2[,,ii-burnin]<-ALPHA2
      SIGMAdraws1[,,ii-burnin]<-SIGMA1
      SIGMAdraws2[,,ii-burnin]<-SIGMA2
      tardraws[ii-burnin]<-tart
      # Calculating Impulse-response function for each regime
      if(intercept==TRUE){
        beta1 <- ALPHA1[2:(K*lags+1),]
        beta2 <- ALPHA2[2:(K*lags+1),]
        x11   <- x1[,2:(K*lags+1)]
        x22   <- x2[,2:(K*lags+1)]
      }
      for(ll in 1:K){
        tirf1 <- .tirf(K,obs1,irfhorizon,lags,x11,beta1,beta2,SIGMA1,SIGMA2,tart,thresh,thDelay,shockvar=ll)
        tirf2 <- .tirf(K,obs2,irfhorizon,lags,x22,beta1,beta2,SIGMA1,SIGMA2,tart,thresh,thDelay,shockvar=ll)
        irf1draws[,,ll,ii-burnin]<-tirf1
        irf2draws[,,ll,ii-burnin]<-tirf2
      }

    }
    if(ii==burnin){
      cat("Burn-in finished\n")
    }
  }
  for(ii in 1:K){
    for(jj in 1:K){
      for(kk in 1:irfhorizon){
        irffinal1[jj,kk,ii,1]<-quantile(irf1draws[jj,kk,ii,],probs=0.5)
        irffinal1[jj,kk,ii,2]<-quantile(irf1draws[jj,kk,ii,],probs=lowerquantile)
        irffinal1[jj,kk,ii,3]<-quantile(irf1draws[jj,kk,ii,],probs=upperquantile)
        irffinal2[jj,kk,ii,1]<-quantile(irf2draws[jj,kk,ii,],probs=0.5)
        irffinal2[jj,kk,ii,2]<-quantile(irf2draws[jj,kk,ii,],probs=lowerquantile)
        irffinal2[jj,kk,ii,3]<-quantile(irf2draws[jj,kk,ii,],probs=upperquantile)
      }
    }
  }
  return(list(Sigma1=SIGMAdraws1,Sigma2=SIGMAdraws2,beta1=ALPHAdraws1,beta2=ALPHAdraws2,thDraws=tardraws,irf1=irffinal1,irf2=irffinal2))
}

#'
#' Helper function for the creation.
.tarpost <- function(X,Ystar,Ytest,beta1,beta2,sigma1,sigma2,tart,lags,intercept=TRUE,tarmean,tarstandard,ncrit=0.15){
  e1 <- Ytest<tart
  e2 <- Ytest>=tart
  nc<-nrow(Ystar)
  # test if there are enough observations in each sample
  if(sum(e1)/nc < ncrit || sum(e2)/nc<ncrit){
    post=-Inf
    loglik1=-Inf
    loglik2=-Inf
    prior=-Inf
  }
  else{
    Y1=Ystar[e1,]
    Y2=Ystar[e2,]
    X1=X[e1,]
    X2=X[e2,]
    if(intercept==TRUE){
      X1=cbind(1,X1)
      X2=cbind(1,X2)
    }
    loglik1 <- .loglike(beta1,sigma1,Y1,X1)
    loglik2 <- .loglike(beta2,sigma2,Y2,X2)
    prior<-pnorm(tart,mean=tarmean,sd=tarstandard)
  }
  post=loglik1+loglik2+prior
  return(list(post=post,lik1=loglik1,lik2=loglik2,prior=prior))

}
.loglike <- function(beta1,sigma,Y,X){
  T=nrow(Y)
  N=ncol(Y)
  v=Y-X%*%beta1
  sterm=0
  isigma <- solve(sigma)
  for(ii in 1:T){
    sterm<-sterm+t(v[ii,])%*%isigma%*%v[ii,]
  }
  dsigma <- det(isigma)
  lik <- (T/2)*dsigma-0.5*sterm

}
