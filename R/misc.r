#' @export
#' @title lagdata - function to lag data
#' @param mydata time series
#' @param nolags integer, number of lags
#' @param intercept logical, whether the lagged series has an intercept or not
#' @return list with the lagged time series and the shortened original series

lagdata <-function(mydata, nolags, intercept = FALSE){

  Time <- nrow(mydata)
  K    <- ncol(mydata)
  obs  <- Time - nolags
  x    <- embed(mydata,dimension = nolags + 1)[,-(1:K)]
  if(intercept == TRUE){

    x <- cbind(1,x)

  }
  yi <- mydata[(nolags + 1):Time,]

  return(list(y = yi,x = x,obs = obs,Time = Time, K = K));
}

#' @export
#' @title stability - test stability of estimate
#' @param betadraw K*nolags x K - matrix with the draw to test for stability, with K being the number of variables. Draw should not include the intercept
#' @param nolags number of lags
#' @return double with the maximum eigenvalue of the companion matrix


stability <- function(betadraw,nolags){


  companion <- companionmatrix(betadraw,nolags)

  K <- ncol(betadraw)
  ev <- eigen(companion)$values
  maxev <- max(abs(ev))

  return(maxev)

}

#' @export
#' @title stability - Function to rewrite the VAR-estimate in companion form.
#' @param betadraw K*nolags x K - matrix with the draw to rewrite in companion form, with K being the number of variables. Draw should not include the intercept
#' @param nolags number of lags
#' @return a K * nolags x K * nolags-matrix of the VAR-estimate in companion form

companionmatrix <- function(betadraw,nolags){

  K <- ncol(betadraw)

  companion <- array(0, dim=c( K * nolags, K * nolags ))

  Ai <- array(0,dim=c( K, K, nolags ))

  if(nolags > 1){

    for(jj in 1:nolags){

      indxmin  <- ( jj - 1 ) * K + 1
      indxmax  <- ( jj - 1 ) * K + K

      Ai[,,jj] <- betadraw[indxmin:indxmax,]

      companion[1:K,indxmin:indxmax] <- betadraw[indxmin:indxmax,]
    }

    indxmin <- K + 1
    indxmax <- K * nolags
    indxrep <- indxmax - indxmin
    for(jj in 0:indxrep){

      companion[indxmin + jj, jj + 1] <- 1

    }
  }
  else{

    companion <- betadraw

  }

  return(companion)

}

# Function to check whether a method for a class exists exists or not
check_exist_method <- function(func,object){

  meth <- methods(class = class(object))
  tmp  <- gsub(paste0(".",class(object)),"",meth)
  return(sum(tmp == func) > 0)

}


invpd <- function(x){
  as.matrix(x)
  xncol <- ncol(x)
  temp <- diag(1,xncol)
  ipd <- pracma::mldivide(x,temp,pinv=TRUE)
  return(ipd)
}

#' function to split a series into two
#' @param y data
#' @param lags number of lags
#' @param thDelay delay of threshold
#' @param thresh threshold variable
#' @param tart value of threshold
#' @param intercept whether both series should include an intercept

splitVariables <- function(y,lags,thDelay,thresh,tart,intercept){

  startest <- max(thDelay,lags)

  T <- nrow(y)
  K <- ncol(y)

  ytest <- y[(startest+1-thDelay):(T-thDelay),thresh]
  ystar <- y[(startest+1):T,]
  xstar <- embed(y,dimension=lags+1)[,-(1:K)]

  if(thDelay>lags){

    diff1=(thDelay-lags)+1
    xnr <- nrow(xstar)
    xstar <- xstar[diff1:xnr,]

  }

  e1 <- ytest <  tart
  e2 <- ytest >= tart
  y1 <- ystar[e1,]
  y2 <- ystar[e2,]
  x1 <- xstar[e1,]
  x2 <- xstar[e2,]
  if(intercept==TRUE){
    x1 <- cbind(1,x1)
    x2 <- cbind(1,x2)
  }

  return(list(y1=y1,y2=y2,x1=x1,x2=x2,xstar=xstar,ytest=ytest,ystar=ystar,e1=e1))

}

################################################################
#
# Function to calculate the acceptance probability
# X - lagged time series
# Y - time series
# beta1, beta2 - VAR Coefficients
# sigma1,sigma2 - variance/covariance matrix
# tart - value of threshold variable
# lags - lag order
# intercept - does the series contain an intercept or not
# tarmean - mean
# tarstandard - standard deviation
# ncrit - percentage of observation each regime must have
#
################################################################

tarpost <- function(X,Ystar,Ytest,beta1,beta2,sigma1,sigma2,tart,lags,intercept=TRUE,tarmean,tarstandard,ncrit=0.15){

  e1 <- Ytest < tart
  e2 <- Ytest >= tart
  nc <- nrow(Ystar)

  # test if there are enough observations in each sample

  p1 = 0

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

    loglik1 <- loglike(beta1,sigma1,Y1,X1)$lik
    loglik2 <- loglike(beta2,sigma2,Y2,X2)$lik

    prior<-pnorm(tart,mean=tarmean,sd=tarstandard)

  }
  post <- (loglik1+loglik2+prior)
  return(list(post=post,lik1=loglik1,lik2=loglik2,prior=prior))

}

exptarpost <- function(X,Ystar,Ytest,beta1,beta2,sigma1,sigma2,tart,lags,intercept=TRUE,tarmean,tarstandard,ncrit=0.15){
  e1 <- Ytest < tart
  e2 <- Ytest >=tart
  nc <- nrow(Ystar)
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

    loglik1 <- loglike(beta1,sigma1,Y1,X1)
    loglik2 <- loglike(beta2,sigma2,Y2,X2)

    #prior <- pnorm(tart,mean=tarmean,sd=tarstandard)

    isig <- 1/tarstandard
    res <- tart-tarmean
    expterm <- -0.5*res*isig*res

    constant <- log(1/(2*sqrt(pi)))-0.5*log(tarstandard)

    prior <- constant + expterm

    ds1 <- loglik1$ds
    ds2 <- loglik2$ds

    sterm1 <- loglik1$sterm
    sterm2 <- loglik2$sterm

    post   <- loglik1$lik + loglik2$lik + prior
    #post <- ds1*ds2*exp(-0.5*(sterm1+sterm2))
  }

  return(list(post=post,lik1=loglik1,lik2=loglik2,prior=prior))

}


#' @title calculate the log-likelihood function for a VAR
#' @param beta VAR coefficients
#' @param sigma variance-covariance matrix
#' @param Y lhs
#' @param X rhs

loglike <- function(beta,sigma,Y,X){

  obs <- nrow(Y)
  K   <- ncol(Y)

  residuals <- Y - X %*% beta
  sterm = 0
  isigma <- solve(sigma)
  dsigma <- log(det(isigma))

  for(ii in 1:obs){

    resid <- matrix(residuals[ii,],nrow = 1)
    sterm <- sterm + resid %*% isigma %*% t(resid)

  }


  lik <- (obs/2) * dsigma - 0.5 * sterm

  return(list(lik = lik,ds = (obs / 2) * dsigma, sterm = sterm))

}



tirf <- function(y,ytest,beta1,beta2,sigma1,sigma2,tar,thVar,thDelay,nolags,irfhor,intercept=TRUE,shockvar,bootrep=50){
  T <- nrow(y)
  endest <- max(thDelay,nolags)
  startdel <- endest - thDelay
  startlag <- endest - nolags
  Tmax <- T-endest
  e1 <- 0
  e2 <- 0
  irf1 <- 0
  irf2 <- 0
  for(ii in 1:Tmax){

    y0 <- y[(startlag+ii):(endest+ii-1),]
    y1 <- y[(startdel+ii):(endest+ii-1),]

    if(ytest[ii]<=tar){

      e1 <- e1+1
      xx <- tirfsimu(y0,y1,beta1,beta2,sigma1,sigma2,tar,thVar,thDelay,nolags,irfhor,intercept,shockvar=shockvar,bootrep)
      irf1 <- irf1+xx

    }
    else{

      e2 <- e2+1
      xx <- tirfsimu(y0,y1,beta1,beta2,sigma1,sigma2,tar,thVar,thDelay,nolags,irfhor,intercept,shockvar=shockvar,bootrep)
      irf2 <- irf2+xx

    }
  }
  irf1 <- irf1/e1
  irf2 <- irf2/e2
  return(list(irf1=irf1,irf2=irf2))
}

tirfsimu <- function(y0,y1,beta1,beta2,sigma1,sigma2,tar,thVar,thDelay,NoLags,irfhor,Intercept=TRUE,shockvar=1,bootrep=50){

  K <- ncol(y0)
  constant <- 0
  if(Intercept==TRUE) constant=1

  csigma1 <- t(chol(sigma1))
  csigma2 <- t(chol(sigma2))
  d1 <- diag(csigma1)
  d2 <- diag(csigma2)
  csx1 <- csigma1/d1
  csx2 <- csigma2/d2

  saveshock   <- 0
  savenoshock <- 0

  for(irep in 1:bootrep){
    yhatnoshock <- array(0,dim=c(irfhor+NoLags,K))
    yhatnoshock[1:NoLags,] <- y0

    yhatshock <- array(0,dim=c(irfhor+NoLags,K))
    yhatshock[1:NoLags,] <- y0

    ystarnoshock <- array(0,dim=c(irfhor+thDelay,K))
    ystarnoshock[1:thDelay,] <- y1

    ystarshock <- array(0,dim=c(irfhor+thDelay,K))
    ystarshock[1:thDelay,] <- y1
    for(ii in 1:irfhor){
      fi <- NoLags+ii
      xhatnoshock <- matrix(0,nrow=1,ncol=(K*NoLags+constant))
      xhatshock   <- matrix(0,nrow=1,ncol=(K*NoLags+constant))

      for(ji in 1:NoLags){
        xhatnoshock[1,((ji-1)*K+1+constant):(ji*K+constant)]<-yhatnoshock[(fi-ji),]
        xhatshock[1,((ji-1)*K+1+constant):(ji*K+constant)]<-yhatshock[(fi-ji),]
      }
      if(Intercept==TRUE){
        xhatnoshock[1,1] <- 1
        xhatshock[1,1] <- 1
      }

      # Next step for shocked variable
      e1 <- ystarshock[ii,thVar]<=tar
      e2 <- ystarshock[ii,thVar]>tar

      if(ii==1){

        uu1 <- matrix(0,nrow=1,ncol=K)
        uu2 <- matrix(0,nrow=1,ncol=K)

        uu1[1,shockvar] <- 0.1
        uu2[1,shockvar] <- 0.1

        yyshock <- (xhatshock%*%beta1+uu1%*%csx1)*e1+(xhatshock%*%beta2+uu2%*%csx2)*e2

      }
      else{

        uu1 <- rnorm(K,0,0.1)
        uu2 <- rnorm(K,0,0.1)

        yyshock <- (xhatshock%*%beta1+uu1%*%csigma1)*e1+(xhatshock%*%beta2+uu2%*%csigma2)*e2

      }

      # Next step for unshocked variable
      e1 <- ystarnoshock[ii,thVar]<=tar
      e2 <- ystarnoshock[ii,thVar]>tar

      if(ii==1){

        uu1 <- matrix(0,nrow=1,ncol=K)
        uu2 <- matrix(0,nrow=1,ncol=K)

        yynoshock <- (xhatnoshock%*%beta1+uu1%*%csx1)*e1+(xhatnoshock%*%beta2+uu2%*%csx2)*e2

      }
      else{

        uu1 <- rnorm(K,0,0.1)
        uu2 <- rnorm(K,0,0.1)

        yynoshock <- (xhatnoshock%*%beta1+uu1%*%csigma1)*e1+(xhatnoshock%*%beta2+uu2%*%csigma2)*e2

      }

      ystarshock[thDelay+ii,] <- yyshock
      yhatshock[NoLags+ii,] <- yyshock

      ystarnoshock[thDelay+ii,] <- yynoshock
      yhatnoshock[NoLags+ii,] <- yynoshock

      #readline(prompt="Press [enter] to continue")
    } # finish loop for one path
    saveshock <- saveshock+yhatshock[(NoLags+1):(irfhor+NoLags),]
    savenoshock <- savenoshock+yhatnoshock[(NoLags+1):(irfhor+NoLags),]

  } # loop over bootstraps

  # save variables

  saveshock <- saveshock/bootrep
  savenoshock <- savenoshock/bootrep
  irf <- saveshock-savenoshock
  return(t(irf))
}
