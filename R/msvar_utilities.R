#' @title Hamiltonfilter
#' @param BETA 3D-array of VAR-coefficients
#' @param SIGMA 3D-array of variance-covariance matrices
#' @param transmat an hxh matrix with the transition probabilities
#' @param y lhs of data
#' @param x rhs of data
#' @param h number of regimes
#' @noRd
hamiltonfilter = function(BETA,SIGMA,transmat,y,x,h=2){

  # Preliminaries
  T <- nrow(y)
  transmat <- t(transmat)

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


      em <- y[it,] - x[it,] %*% BETA[,,ii]

      neta[ii,] <- log(1 / sqrt(detSigma[ii,1])) + (-0.5 * (em %*% invSigma[,,ii] %*% t(em)))

    }

    # Hamilton Filter
    ett1 <- ett * exp(neta)
    fit <- sum(ett1,na.rm=TRUE)
    ett <- (transmat %*% ett1)/fit

    fprob[it,] <- t(ett)

    if(anyNA(fprob[it,])){ fprob[it,] <- 1 / h}

    if(fit > 0){

      lik <- lik + log(fit)

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

    ST[T,1] <- sample(x=seq(1:h),size = 1,replace = FALSE, prob = fprob[T,])

    p1 <- array(0,dim=c(h,1))

    for(ii in (T-1):1){

      for(jj in 1:h){

        p1[jj,1] <- pr_tr[ST[ii+1],jj]%*%fprob[ii,jj]

      }

      p <- p1 / sum(p1)
      ST[ii] <- sample(x = seq(1:h),size = 1,replace = FALSE,prob = p) #bingen(p,h)

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

#' @title counts the sequence of transitions
#' @param h the number of regimes
#' @param x a vector of length T with the regimes
#' @return an hxh matrix with the counts
#' @noRd
countseq <- function(x,h){

  x <- as.matrix(x)
  nt <- nrow(x)
  stt <- matrix(0,h,h)

  for(ii in 2:nt){

    stt[x[ii-1],x[ii]] <- stt[x[ii-1],x[ii]]+1

  }

  return(stt)
}
