#' @title Initialize the mcmc sampler for conjugate Normal-Wishart-Prior
#'
#' @param priorObj an object containing the prior
#' @param yLagged lagged data lhs
#' @param xLagged lagged data rhs
#' @param ... currently not used
#'
initialize_mcmc.cnw <- function(priorObj,yLagged,xLagged,...){

  cat("Initialize MCMC sampler for Conjugate Normal-Wishart prior\n")

  AOLS      <- solve(t(xLagged) %*% xLagged) %*% t(xLagged) %*% yLagged
  residuals <- yLagged - xLagged %*% AOLS
  SSE       <- t(residuals) %*% residuals
  Sigma     <- SSE / nrow(yLagged)

  return(list(Alpha = AOLS, Sigma = Sigma, addinfo = NULL))

}

#' @title Initialize the mcmc sampler for uninformative prior
#'
#' @param priorObj an object containing the prior
#' @param yLagged lagged data lhs
#' @param xLagged lagged data rhs
#' @param ... currently not used
#'

initialize_mcmc.unf <- function(priorObj,yLagged,xLagged,...){

  cat("Initialize MCMC sampler for uninformative prior\n")

  AOLS      <- solve(t(xLagged) %*% xLagged) %*% t(xLagged) %*% yLagged
  residuals <- yLagged - xLagged %*% AOLS
  SSE       <- t(residuals) %*% residuals
  Sigma     <- SSE / nrow(yLagged)

  return(list(Alpha = AOLS, Sigma = Sigma, addinfo = NULL))

}

#' @title Initialize the mcmc sampler for Minnesota Prior
#'
#' @param priorObj an object containing the prior
#' @param yLagged lagged data lhs
#' @param xLagged lagged data rhs
#' @param ... currently not used
#'
#'
initialize_mcmc.minnesota <- function(priorObj,yLagged,xLagged,...){

  cat("Initialize MCMC sampler for Minnesota prior\n")

  AOLS      <- solve(t(xLagged) %*% xLagged) %*% t(xLagged) %*% yLagged
  AOLS      <- solve(t(xLagged) %*% xLagged) %*% t(xLagged) %*% yLagged
  residuals <- yLagged - xLagged %*% AOLS
  SSE       <- t(residuals) %*% residuals
  Sigma     <- SSE / nrow(yLagged)

  return(list(Alpha = AOLS, Sigma = Sigma, addinfo = NULL))

}

#' @title Initialize the mcmc sampler for SSVS-Prior
#'
#' @param priorObj an object containing the prior
#' @param yLagged lagged data lhs
#' @param xLagged lagged data rhs
#' @param ... currently not used
initialize_mcmc.ssvs <- function(priorObj,yLagged,xLagged,...){

  cat("Initialize MCMC sampler for SSVS prior\n")

  K         <- ncol(yLagged)
  obs       <- nrow(yLagged)
  intercept <- priorObj$intercept
  norest    <- K * (K * priorObj$nolags + intercept)

  # Initial values on parameters for tightness on coefficients and variance-covariance

  gammas    <- array(1,dim=c(norest,1)) # coefficients
  omega     <- array(list(),dim=c(K-1,1)) # variance-covariance

  for(kk1 in 1:(K - 1)){

    for( ii in 1:1){

      omega[[kk1,ii]] <- array(1,dim=c(kk1))

    }
  }

  # Initial values for coefficients and variance-covariance matrix

  Alpha     <- solve( t(xLagged) %*% xLagged ) %*% t(xLagged) %*% yLagged

  residuals <- yLagged - xLagged %*% Alpha
  Sigma       <- t(residuals) %*% residuals / obs

  addInfo = list(gammas = gammas,
                 omega  = omega,
                 SSEGibbs = Sigma)

  retlist <- list(Alpha = Alpha,
                  Sigma = Sigma,
                  addInfo = addInfo)

  return(retlist)


}

