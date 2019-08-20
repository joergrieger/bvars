#' @export
#' @title extract factors from a time series
#' @param factordata data from which the factors should be extracted
#' @param no_factors number of factors that are going to be extracted
#' @return matrix with the extracted factors
#'
get_factors <- function(factordata,no_factors){

  nv <- ncol(factordata) # number of variables in factordata
  nObs <- nrow(factordata)
  x.x <- t(factordata) %*% factordata
  evectors <- eigen(x.x)$vectors
  ret.evectors <- sqrt(nObs)*evectors[,1:no_factors]

  fac <- factordata %*% ret.evectors/nObs

  return(fac)

}
#' @export
#' @title draw posterior for measurement equation using a normal-gamma prior
#' @param li_prvar prior on variance for coefficients
#' @param fy 'independent' variables
#' @param xy 'dependent' variables
#' @param K number of variables in ts model
#' @param N total number of variables in factor data
#' @param P Number variables in ts model plus number of factors
#' @param Sigma previous draw of variance-covariance matrix
#' @param L previous draw of coefficients
#' @param alpha,beta prior on variances
#'
draw_posterior_normal <- function(li_prvar,fy,xy,K,P,N,Sigma,L,alpha,beta){

  for(ii in 1:(N + K)){
    if(ii > K){

      Li_postvar   <- solve(solve(li_prvar) + Sigma[ii,ii]^(-1) * t(fy) %*% fy)
      Li_postmean  <- Li_postvar %*% (Sigma[ii,ii]^(-1) * t(fy) %*% xy[,ii])
      L[ii,1:P]    <- t(Li_postmean) + rnorm(P) %*% t(chol(Li_postvar))

    }
    resi <- xy[,ii] - fy %*% L[ii,]
    sh   <- alpha/2 + T/2
    sc   <- beta/2  + t(resi) %*% resi
    Sigma[ii,ii] <- rgamma(1,shape=sh,scale=sc)

  }

  return( list(L=L,Sigma=Sigma) )

}
#' @export
#' @title Draw posterior for measurement equation using an SSVS-prior
#' @param fy
#' @param fy 'independent' variables
#' @param xy 'dependent' variables
#' @param K number of variables in ts model
#' @param N total number of variables in factor data
#' @param P Number variables in ts model plus number of factors
#' @param L previous draw of coefficients
#' @param Sigma previous draw of variances
#' @param tau2 variance of coefficients
#' @param c2 factor for tau2
#' @param gammam previous draw of gammas
#' @param alpha,beta priors for variances
#'
draw_posterior_ssvs <- function(fy,xy,K,P,N,Sigma,tau2,c2,gammam,alpha,beta,L){

  for(ii in 1:(N + K)){
    if(ii > K){

      # Sample betas
      VBeta <- diag(gammam[,ii] * c2 * tau2 + (1-gammam[,ii]) * tau2)
      DBeta <- solve(t(fy) %*% fy * Sigma[ii,ii]^(-1))
      dbeta <- t(fy) %*% xy[,ii] * Sigma[ii,ii]^(-1)
      HBeta <- t(chol(DBeta))

      L[ii,1:P] <- t(DBeta %*% dbeta) + (rnorm(P) %*% HBeta)

      # Sample the gammas
      for(jj in 1:P){

        numerator <- pnorm(L[ii,jj],mean=0,sd=sqrt(c2 * tau2))
        denominator <- numerator + pnorm(L[ii,jj],mean=0,sd=sqrt(tau2))
        prob <- numerator / denominator
        gammam[jj,ii] <- 0.5*sign(runif(1)-prob)+0.5

      }

      # Sample the variance
      resid <- xy[,ii] - fy %*% L[ii,]
      sh <- alpha/2 + T/2
      sc <- beta/2 + t(resid) %*% resid
      Sigma[ii,ii] <- rgamma(1, shape = sh, scale = sc)

    }

  }

  return(list(Sigma = Sigma, L=L,gammam=gammam))
}

#' @title linear regression using single value decomposition
#' @param y dependent variable
#' @param ly independent variable

olssvd <- function(y,ly){

  duv   <- svd(t(ly) %*% ly)
  x_inv <-duv$v %*% diag(1/duv$d) %*% t(duv$u)
  x_pseudo_inv <- x_inv %*% t(ly)
  b <- x_pseudo_inv %*% y
  return(b)
}
