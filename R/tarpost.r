.tarpost <- function(X,Ystar,Ytest,beta1,beta2,sigma1,sigma2,tart,lags,intercept=TRUE,tarmean,tarstandard,ncrit=0.15){
  e1 <- Ytest<tart
  e2 <- Ytest>=tart
  nc<-nrow(Ystar)
  # test if there are enough observations in each sample
  p1=0
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
    loglik1 <- .loglike(beta1,sigma1,Y1,X1)$lik
    loglik2 <- .loglike(beta2,sigma2,Y2,X2)$lik
    prior<-pnorm(tart,mean=tarmean,sd=tarstandard)
  }
  post <- (loglik1+loglik2+prior)
  return(list(post=post,lik1=loglik1,lik2=loglik2,prior=prior))

}

.exptarpost <- function(X,Ystar,Ytest,beta1,beta2,sigma1,sigma2,tart,lags,intercept=TRUE,tarmean,tarstandard,ncrit=0.15){
  e1 <- Ytest<tart
  e2 <- Ytest>=tart
  nc<-nrow(Ystar)
  # test if there are enough observations in each sample
  if(sum(e1)/nc < ncrit || sum(e2)/nc<ncrit){
    post=0
    loglik1=0
    loglik2=0
    prior=0
	
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
	ds1 <- loglik1$ds
	ds2 <- loglik2$ds
	sterm1 <- loglik1$sterm
	sterm2 <- loglik2$sterm
	post <- ds1*ds2*exp(-0.5*(sterm1+sterm2))
	
  }
  
  
  return(list(post=post,lik1=loglik1,lik2=loglik2,prior=prior))

}

.loglike <- function(beta1,sigma,Y,X){
  T=nrow(Y)
  N=ncol(Y)
  v=Y-X%*%beta1
  sterm=0
  isigma <- (sigma)
  for(ii in 1:T){
    sterm<-sterm+t(v[ii,])%*%isigma%*%v[ii,]
  }
  dsigma <- det(isigma)
  lik <- (T/2)*dsigma-0.5*sterm
  return(list(lik=lik,ds=(T/2)*dsigma,sterm=sterm))

}