 #
 # Calculate ergodic probabilities of the transition matrix
 # See eq. 4.44-4.49 in Kim/Nelson
 #

ergodicprob <- function(transmat,h){
	temp1 <- diag(1,h)-transmat
	temp2 <- matrix(1,nrow=1,ncol=h)
	A <- rbind(temp1,temp2)
	EN <- matrix(0,ncol=1,nrow=(h+1))
	EN[(h+1),1] <- 1
	ett <- solve(t(A)%*%A)%*%t(A)[,(h+1)]
	return(ett)
}

hamiltonfilter = function(BETA,SIGMA,transmat,y,x,h=2){
  # Preliminaries
  T <- nrow(y)
  transmat <- t(transmat)
  #ett <- ergodicprob(transmat,h)
  ett <- 1/h*array(1,dim=c(h,1))

  #
  # Calculate the inverse of the Variance-Covariance Matrix for each regime
  #

  invSigma <- array(0,dim=c(nrow(SIGMA[,,1]),nrow(SIGMA[,,1]),h))
  detSigma <- array(0,dim=c(h,1))
  for(ii in 1:h){

    invSigma[,,ii] <- invpd(SIGMA[,,ii])
    detSigma[ii,1] <- det(SIGMA[,,ii])

  }

  # Filter Forward

  lik <- 0
  fprob <- array(0,dim=c(T,h))
  resi <- array(0,dim=c(T,ncol(y),h))
  neta <- array(0,dim=c(h,1))

  for(it in 1:T){

    for(ii in 1:h){


      em <- y[it,]-x[it,]%*%BETA[,,ii]
      neta[ii,] <- log(1/sqrt(detSigma[ii,1])) + (-0.5*(em%*%invSigma[,,ii]%*%t(em)))

    }

    # Hamilton Filter
    ett1 <- ett*exp(neta)
    fit <- sum(ett1,na.rm=TRUE)
    ett <- (transmat%*%ett1)/fit
	  fprob[it,] <- t(ett)

	  if(anyNA(fprob[it,])){ fprob[it,] <- 1/h}

	  if(fit > 0){

      lik <- lik+log(fit)

    }
    else{
      lik <- lik - 10
    }
  }
  return(list(fprob=fprob,lik=lik))
}

getst <- function(fprob,transmat,ncrit=0.15,h=2){
  T <- nrow(fprob)
  ST <- array(1,dim=c(T,1))
  pr_tr <- t(transmat)
  # Time T
  check <- 2

  while(check>1){

    ST[T,1] <- sample(x=seq(1:h),size=1,replace=FALSE,prob=fprob[T,])#bingen(fprob[T,],h)
    p1 <- array(0,dim=c(h,1))

    for(ii in (T-1):1){
      for(jj in 1:h){

        p1[jj,1] <- pr_tr[ST[ii+1],jj]%*%fprob[ii,jj]

      }

      p <- p1/sum(p1)
      ST[ii] <- sample(x=seq(1:h),size=1,replace=FALSE,prob=p) #bingen(p,h)

    }
    checkncrit <- TRUE
    for(ii in 1:h){
      checkncrit <- checkncrit && (sum(ST==ii)/T>ncrit)
    }
    if(checkncrit){
      check=0
    }
  }
  return(ST)
}

countseq <- function(x,h){
  x <- as.matrix(x)
  nt <- nrow(x)
  stt <- matrix(0,h,h)
  for(ii in 2:nt){
    stt[x[ii-1],x[ii]] <- stt[x[ii-1],x[ii]]+1
  }
  return(stt)
}
# Draw next step in markov chain
bingen <- function(transprob,h){
  r <- runif(1)
  ret <- 1
  tpcdf <- array(0,dim=c(h,1))
  for(ii in 1:h){
    tpcdf[ii,1]<-sum(transprob[1:ii])
  }
  return(min(which(tpcdf>r)))
}
