hamiltonfilter = function(BETA,SIGMA,transmat,y,x,h=2){
  # Preliminaries
  T <- nrow(y)
  # Initialise the filter (ergodic probabilities)
  temp1 <- diag(1,h)-transmat
  temp2 <- matrix(1,nrow=1,ncol=h)
  A <- rbind(temp1,temp2)
  EN <- matrix(0,ncol=1,nrow=(h+1))
  EN[(h+1),1] <- 1
  ett <- solve(t(A)%*%A)%*%t(A)%*%EN
  
  # Filter Forward
  lik <- 0
  fprob <- array(0,dim=c(T,h))
  #readline(prompt="Press [enter] to continue")
  for(it in 1:T){
    neta <- array(0,dim=c(h,1))
    for(ii in 1:h){
      em <- y[it,]-x[it,]%*%BETA[,,ii]
      neta[ii,] <- 1/sqrt(det(SIGMA[,,ii]))*exp(-0.5*(em%*%invpd(SIGMA[,,ii])%*%t(em)))
    }
    # Hamilton Filter
    ett1 <- ett*neta
    #print(ett1)
    #readline(prompt="Press [enter] to continue")
    fit <- sum(ett1,na.rm=TRUE)
    ett <- (transmat%*%ett1)/fit
    fprob[it,] <- t(ett1/fit)
    if(fit > 0){
      lik <- lik+log(fit)
    }
    else{
      lik <- lik - 10
    }
  }
  #print(fprob)
  #readline(prompt="Press [enter] to continue")
  return(list(fprob=fprob,lik=lik))
}

getst <- function(fprob,transmat,ncrit=0.05,h=2){
  T <- nrow(fprob)
  ST <- array(1,dim=c(T,1))
  pr_tr <- transmat
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
      #print(sum(ST==ii))
      checkncrit <- checkncrit && (sum(ST==ii)/T>ncrit)
    }
    #print(checkncrit)
    #print(transmat)
    #readline(prompt="Press [enter] to continue")
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