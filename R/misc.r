
# vectorizes a matrixs
.vec <- function(Mat){
  yret <- matrix(Mat,ncol=1)
  return(yret)
}
# creates an NxN identity matrix
.id <-function(N){
  return(diag(rep(1,N)))
}

# Checks if VAR is stable
.stability <- function(betadraw,lags,K){
  companion <- companionmatrix(betadraw,lags)
  ev <- eigen(companion)$values
  maxev <- max(abs(ev))
  #companion <- array(0,dim=c(K*lags,K*lags))
  #Ai <- array(0,dim=c(K,K,lags))
  #if(lags>1){
  #  for(jj in 1:lags){
  #    indxmin <- (jj-1)*K+1
  #    indxmax <- (jj-1)*K+K
  #    Ai[,,jj]<-betadraw[indxmin:indxmax,]
  #    companionmatrix[1:K,indxmin:indxmax]<-betadraw[indxmin:indxmax,]
  #  }
  #  indxmin <- K+1
  #  indxmax <- K*lags
  #  indxrep <- indxmax-indxmin
  #  for(jj in 0:indxrep){
  #    companionmatrix[indxmin+jj,jj+1]<-1
  #  }
  #  ev <- eigen(companionmatrix)$values
  #  ev <- Mod(ev)
  #  maxev <- max(abs(ev))
  #}
  #else{
  #  ev <- eigen(betadraw)$values
  #  ev <- Mod(ev)
  #  maxev <- max(abs(ev))
  #}
  return(maxev)
}
# Function to create the Companion Matrix
companionmatrix <- function(Beta,lags){
  companion <- array(0,dim=c(K*lags,K*lags))
  Ai <- array(0,dim=c(K,K,lags))
  if(lags>1){
    for(jj in 1:lags){
      indxmin <- (jj-1)*K+1
      indxmax <- (jj-1)*K+K
      Ai[,,jj]<-betadraw[indxmin:indxmax,]
      companion[1:K,indxmin:indxmax] <- betadraw[indxmin:indxmax,]
    }
    indxmin <- K+1
    indxmax <- K*lags
    indxrep <- indxmax-indxmin
    for(jj in 0:indxrep){
      companion[indxmin+jj,jj+1]<-1
    }
  }
  else{
    companion <- Beta
  }
  
  return(companion)
}
# function to test if a provided variable is a scalar
.isscalar <- function(x) is.atomic(x) && length(x) == 1L && !is.character(x) && Im(x)==0

# Function to check sign restrictions
# Returns False if sign-restrictions are satisfied and TRUE if sign restrictions are not satisfied
.CheckSign <- function(RestrictionMatrix,TestMatrix){
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

.vec2matrix <- function(lags,K,vec){
     matrix0 <- array(0,dim=c(lags,K))
     for(ii in 1:lags){
	     for(jj in 1:K){
		     matrix0[ii,jj] <- vec[(ii-1)*K+jj]
		 }
	 }
    return(matrix0)
}
detrend1 <- function(y,q){
  vague <- 1e+4
  f <- array(0,dim=c(2,2))
  f[1,1] <- 2
  f[1,2] <- -1
  f[2,1] <- 1
  x <- array(0,dim=c(2,2))
  p <- vague*array(1,dim=c(2,2))
  xf <- array(0,dim=c(nrow(y,1)))
  for(i in 1:nrow(y)){
    x <- f*x
    p <- f*p*t(f)
    p[1,1] <- p[1,1]+q
    if(is.nan(y[i])==FALSE){
      h <- p[1,1]+1
      e <- y[i]-x[1]
      k <- p[,1]/h
      x <- x+k%*%e
      p <- p-(k%*%p[1,])
    }
    xf[i]<-x[1]
  }
  xc <- y-xf
  return(list(xc=xc,xf=xf))
}

transform <- function(x,tcode){
  # tcodes:
  # 1. Level
  # 2. First Difference
  # 3. Second Difference
  # 4. Log-Level
  # 5. Log-First-Difference
  # 6. Log-Second-Difference
  # 7. Detrend Log Using 1-sided HP detrending for Monthly data
  # 8. Detrend Log Using 1-sided HP detrending for Quarterly data
  
  T <- nrow(x)
  print(T)
  small <- 1e-6
  relvarm <- 0.00000075
  relvarq <- 0.000625

  n <- dim(x)[1]
  y <- array(0,dim=c(T,1))
  if(tcode==1){
    y=x
  }
  else if(tcode==3){
    y[3:n] <- x[3:n]-2*x[2:(n-1)]+x[1:(n-2)]
  }
  else if(tcode==2){
    y[2:n] <- x[2:n]-x[1:(n-1)]
  }
  else if(tcode==4){
    if(min(x)<small){
      y <- NaN
    }
    else{
      y <- log(x)
    }
  }
  else if(tcode==5){
    if(min(x)<small){
      y <- NaN
    }
    else{
      x <- log(x)
      y[2:n,1] <- x[2:n,1]-x[1:(n-1),1]
    }
  }
  else if(tcode==6){
    if(min(x)<small){
      y <- NaN
    }
    else{
      y[3:n] <- x[3:n]-2*x[2:(n-1)]+x[1:(n-2)]
    }
  }
  else if(tcode==7){
    if(min(x)<small){
      y <- NaN
    }
    else{
      x <- log(x)
      fr <- detrend1(x,relvarq)
      y <- fr$xc
      tl <- fr$xf
    }
  }
  else if(tcode==8){
    if(min(x)<small){
      y <- NaN
    }
    else{
      x <- log(x)
      fr <- detrend1(x,relvarm)
      y <- fr$xc
      tl <- fr$xf
    }
  }
  else if(tcode==16){
    if(min(x)<small){
      y <- NaN
    }
    else{
      x <- log(x)
      y[3:n]<-x[3:n]-2*x[2:(n-1)]+x[1:(n-2)]
    }
  }
  else if(tcode==17){
    if(min(x)<small){
      y <- NaN
    }
    else{
      y[14:n] <- x[14:n]-x[13:(n-1)]-x[2:(n-12)]+x[1:(n-13)]
    }
  }
  return(y)
}
