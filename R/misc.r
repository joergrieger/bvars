#' @importFrom Rdpack reprompt

#' @title lagdata - function to lag data
#' @param mydata time series
#' @param nolags integer, number of lags
#' @param intercept logical, whether the lagged series has an intercept or not
#' @return list with the lagged time series and the shortened original series
#' @importFrom stats embed
#' @noRd

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

#' @title stability - test stability of estimate
#' @param betadraw K*nolags x K - matrix with the draw to test for stability, with K being the number of variables. Draw should not include the intercept
#' @param nolags number of lags
#' @return double with the maximum eigenvalue of the companion matrix
#' @noRd
stability <- function(betadraw,nolags){


  companion <- companionmatrix(betadraw,nolags)

  K <- ncol(betadraw)
  ev <- eigen(companion)$values
  maxev <- max(abs(ev))

  return(maxev)

}


#' @title stability - Function to rewrite the VAR-estimate in companion form.
#' @param betadraw K*nolags x K - matrix with the draw to rewrite in companion form, with K being the number of variables. Draw should not include the intercept
#' @param nolags number of lags
#' @return a K * nolags x K * nolags-matrix of the VAR-estimate in companion form
#' @noRd
#'
companionmatrix <- function(betadraw,nolags){

  K <- ncol(betadraw)

  companion <- array(0, dim=c( K * nolags, K * nolags ))

  Ai <- array(0,dim=c( K, K, nolags ))

  if(nolags > 1){

    for(jj in 1:nolags){

      indxmin  <- ( jj - 1 ) * K + 1
      indxmax  <- ( jj - 1 ) * K + K

      Ai[,,jj] <- betadraw[indxmin:indxmax,]

      companion[1:K,indxmin:indxmax] <- t(betadraw[indxmin:indxmax,])
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
#' @importFrom utils methods
#' @noRd
check_exist_method <- function(func,object){

  meth <- utils::methods(class = class(object))
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
#' @noRd
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
#' @noRd
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



tirf <- function(y,ytest,beta1,beta2,sigma1,sigma2,tar,thVar,thDelay,nolags,irfhor,intercept=TRUE,shockvar,bootrep=50,id_obj){
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

    if(ytest[ii] <= tar){

      e1 <- e1+1
      xx <- tirfsimu(y0,y1,beta1,beta2,sigma1,sigma2,tar,thVar,thDelay,nolags,irfhor,intercept,shockvar=shockvar,bootrep, id_obj)
      irf1 <- irf1+xx

    }
    else{

      e2 <- e2+1
      xx <- tirfsimu(y0,y1,beta1,beta2,sigma1,sigma2,tar,thVar,thDelay,nolags,irfhor,intercept,shockvar=shockvar,bootrep,id_obj)
      irf2 <- irf2+xx

    }
  }
  irf1 <- irf1/e1
  irf2 <- irf2/e2
  return(list(irf1=irf1,irf2=irf2))
}



tirfsimu <- function(y0,y1,beta1,beta2,sigma1,sigma2,tar,thVar,thDelay,NoLags,irfhor,Intercept=TRUE,shockvar=1,bootrep=50,id_obj){

  K <- ncol(matrix(y0,nrow=NoLags))

  constant <- 0
  if(Intercept==TRUE) constant=1

  csigma1 <- structural(id_obj,beta1,sigma1) # t(chol(sigma1))
  csigma2 <- structural(id_obj,beta2,sigma2) #t(chol(sigma2))
  d1 <- diag(csigma1)
  d2 <- diag(csigma2)
  csx1 <- csigma1/d1
  csx2 <- csigma2/d2

  saveshock   <- 0
  savenoshock <- 0

  for(irep in 1:bootrep){

    yhatnoshock <- array(0,dim=c(irfhor+NoLags,K))
    yhatnoshock[1:NoLags,] <- matrix(y0,nrow = NoLags)

    yhatshock <- array(0,dim=c(irfhor+NoLags,K))
    yhatshock[1:NoLags,] <- matrix(y0, nrow = NoLags)

    ystarnoshock <- array(0,dim=c(irfhor+thDelay,K))
    ystarnoshock[1:thDelay,] <- y1

    ystarshock <- array(0,dim=c(irfhor+thDelay,K))
    ystarshock[1:thDelay,] <- y1

    for(ii in 1:irfhor){
      fi <- NoLags+ii
      xhatnoshock <- matrix(0,nrow=1,ncol=(K*NoLags+constant))
      xhatshock   <- matrix(0,nrow=1,ncol=(K*NoLags+constant))

      for(ji in 1:NoLags){

        xhatnoshock[1,((ji-1)*K+1+constant):(ji * K + constant)] <- yhatnoshock[(fi-ji),]
        xhatshock[1,((ji-1)*K+1+constant):(ji*K + constant)] <- yhatshock[(fi-ji),]

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

        yyshock <- (xhatshock %*% beta1 + uu1 %*% csx1) * e1 + (xhatshock %*% beta2 + uu2 %*% csx2) * e2

      }
      else{

        uu1 <- rnorm(K,0,0.1)
        uu2 <- rnorm(K,0,0.1)

        yyshock <- (xhatshock %*% beta1 + uu1 %*% csigma1 ) * e1 + (xhatshock %*% beta2 + uu2 %*% csigma2) * e2

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

#' @title VMA representation of VAR estimate
#' @param Alpha an nxnxp array with the VAR coefficients
#' @param vma_order order of the vector moving average process
#' @param ar_order the order of the vector autoregressive process
#' @param k the number of variables in the VAR model
#' @return an nxnxvma_order array of the moving average coefficients
#' @noRd
var2vma <- function(Alpha,vma_order,ar_order,k){


  vma <- array(0,dim=c(k,k,vma_order))

  vma[,,1] <- Alpha[,,1]

  if(vma_order > 1){

    for(ii in 2:vma_order){

      if(ii <= ar_order){

        vma[,,ii] <- Alpha[,,ii] %*% vma[,,ii - 1]

      }
    }
  }
  return(vma)
}

#' @export
#' @title parallel compuation of irfs of tvars
#' @param y data
#' @param Alpha coefficients
#' @param Sigma variance-covariance matrices
#' @param tart Value of threshold
#' @param thVar threshold-Variable
#' @param thDelay delay of threshold
#' @param nolags number of lags in the model
#' @param nhor length of the impulse-response function
#' @param intercept intercept or no intercept
#' @param bootrep number of bootstrap replications
#' @param id_obj S3 object for identification of the model
#' @param K number of variables

tirf1 <- function(y,Alpha,Sigma,tart,thVar,thDelay,nolags,nhor,intercept,bootrep,id_obj,K){


  print(thDelay)
  xsplit <- splitVariables(y = y,lags = nolags, thDelay = thDelay, thresh = thVar, tart = tart, intercept = intercept)

  Irfdraws <- array(0,dim=c(K,K,nhor,2))

  for(ii in 1:K){

    xx <- tirf(y = xsplit$ystar, ytest = xsplit$ytest,
               beta1 = Alpha[,,1], beta2 =Alpha[,,2],
               sigma1 = Sigma[,,1], sigma2 = Sigma[,,2],
               tar = tart, thVar = thVar, thDelay = thDelay,
               nolags = nolags, irfhor = nhor, intercept = intercept,
               shockvar=ii,bootrep = bootrep,
               id_obj = id_obj)
    Irfdraws[ii,,,1]<-xx$irf1
    Irfdraws[ii,,,2]<-xx$irf2

  }

  return(Irfdraws)

}

getmode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
