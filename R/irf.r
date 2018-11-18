compirf <- function(A,Sigma,NoLags,intercept=TRUE,nhor){
  K <- nrow(Sigma)
  bigj <- matrix(0,K,K*NoLags)
  bigj[1:K,1:K] <- diag(K)
  if(intercept==TRUE){
    B <- A[2:(nrow(A)),]
  }
  else{
    B <- A
  }
  #if(NoLags>1){
  #  S_F_Mat <- cbind(Sigma,matrix(0,K,K*(NoLags-1)))
  #  S_F_Mat <- rbind(S_F_Mat,matrix(0,K*(NoLags-1),K*NoLags))
  #  xx <- cbind(diag(1,K*(NoLags-1)),matrix(0,K*(NoLags-1),K))
  #  PHI_Mat <- rbind(t(B),xx)
  #}
  #else{
  #  S_F_Mat <- Sigma
  #  PHI_Mat <- B
  #}
  PHI_Mat <- companionmatrix(B,NoLags)
  biga <- PHI_Mat
  bigai <- biga
  shock <- t(chol(Sigma))
  impresp <- matrix(0,K,K*nhor)
  impresp[1:K,1:K] <- shock
  for(ii in 1:(nhor-1)){
	#print(bigai)
	#readline(prompt="Press [enter] to continue")
    impresp[,(ii*K+1):((ii+1)*K)] <- (bigj%*%bigai%*%t(bigj)%*%shock)
    bigai <- bigai%*%biga
  }
  imp <- array(0,dim=c(K,K,nhor))
  jj <- 0
  for(ii in 1:K){
    for(ij in 1:nhor){
      jj <- ii+K*(ij-1)
      imp[ii,,ij] <- impresp[,jj]
    }
  }
  return(imp)
}

compirfsign <- function(A,Sigma,NoLags,intercept=TRUE,nhor,restrictions){
  K <- ncol(Sigma)
  SignRestriction <- FALSE
  cholsigma <- t(chol(Sigma))
  while(!SignRestriction){
    qrmatrix <- matrix(rnorm(K*K),nrow=K)
    qrdecomp <- qr(qrmatrix)
    qrdecomp <- qr.Q(qrdecomp)
    testmatrix <- qrdecomp%*%cholsigma
    SignRestriction <- !CheckSign(Restrictions,testmatrix)
  }
  Sigma <- testmatrix
  irf <- compirf(A,Sigma,NoLags,intercept,nhor)
  return(irf)
}

CheckSign <- function(RestrictionMatrix,TestMatrix){
  # Check if RestrictionMatrix and TestMatrix are of the same sign
  Test1 <- dim(as.matrix(RestrictionMatrix))
  Test2 <- dim(as.matrix(TestMatrix))
  Test <- Test1==Test2
  if(!Test[1]){
    stop("Matrix with sign restrictions and test matrix do not have the same size")
  }
  if(!Test[2]){
    stop("Matrix with sign restrictions and test matrix do not have the same size")
  }
  n1 <- Test1[1]
  n2 <- Test1[2]
  TestFail=FALSE
  for(ii in 1:n1){
    for(jj in 1:n2){
      if(!is.na(RestrictionMatrix[ii,jj])){
        if(RestrictionMatrix[ii,jj]<0){
          if(TestMatrix[ii,jj]>0){
            TestFail=TRUE
          }
        }
        if(RestrictionMatrix[ii,jj]>0){
          if(TestMatrix[ii,jj]<0){
            TestFail=TRUE
          }
        }
      }
    }
  }
  return(TestFail)
}


tirf <- function(y,ytest,beta1,beta2,sigma1,sigma2,tar,thVar,thDelay,NoLags,irfhor,Intercept=TRUE,shockvar,bootrep=50){
  T <- nrow(y)
  endest <- max(thDelay,NoLags)
  startdel <- endest - thDelay
  startlag <- endest - NoLags
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
      xx <- tirfsimu(y0,y1,beta1,beta2,sigma1,sigma2,tar,thVar,thDelay,NoLags,irfhor,Intercept,shockvar=shockvar,bootrep)
      irf1 <- irf1+xx
    }
    else{
      e2 <- e2+1
      xx <- tirfsimu(y0,y1,beta1,beta2,sigma1,sigma2,tar,thVar,thDelay,NoLags,irfhor,Intercept,shockvar=shockvar,bootrep)
      irf2 <- irf2+xx
    }
  }
  irf1 <- irf1/e1
  irf2 <- irf2/e2
  return(list(irf1=irf1,irf2=irf2))
}
tirfsign <- function(y,ytest,beta1,beta2,sigma1,sigma2,tar,thVar,thDelay,NoLags,irfhor,Intercept=TRUE,shockvar,bootrep=50,restrictions){
  K <- ncol(Sigma)
  SignRestriction <- FALSE
  cholsigma <- t(chol(Sigma))
  while(!SignRestriction){
    qrmatrix <- matrix(rnorm(K*K),nrow=K)
    qrdecomp <- qr(qrmatrix)
    qrdecomp <- qr.Q(qrdecomp)
    testmatrix <- qrdecomp%*%cholsigma
    SignRestriction <- !CheckSign(Restrictions,testmatrix)
  }
  Sigma <- testmatrix
  irf <- tirf(y,ytest,beta1,beta2,sigma1,sigma2,tar,thVar,thDelay,NoLags,irfhor,Intercept=TRUE,shockvar,bootrep=50)
  return(irf)

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

  saveshock <- 0
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

