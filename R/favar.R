#' @export
#' @title Factor-Augmented Vector Autoregression
#' @param data data that is not going to be reduced to factors
#' @param factordata data that is going to be reduced to its factors
#' @param nreps total number of draws
#' @param burnin number of burn-in draws.
#' @param nthin thinning parameter
#' @param priorObj An S3 object containing information about the prior.
#' @param priorm Selects the prior on the measurement equation, 1=Normal-Gamma Prior and 2=SSVS prior.
#' @param alpha,beta prior on the variance of the measurement equation
#' @param tau2 variance of the coefficients in the measurement equation (only used if priorm=2)
#' @param c2 factor for the variance of the coefficients (only used if priorm=2)
#' @param li_prvar prior on variance of coefficients (only used if priorm = 1)
#' @param stabletest boolean, check if a draw is stationary or not
favar <- function(data,priorObj,factordata,nreps,burnin,alpha,beta,tau2,c2,li_prvar,priorm,stabletest = TRUE,nthin=1){

  # normalize data
  scaled_data <- scale(data)
  scaled_factordata <- scale(factordata)

  # Variables
  no_lags    <- priorObj$nolags
  no_factors <- priorObj$nofactors
  intercept  <- priorObj$intercept

  nObs <- nrow(scaled_data)
  N    <- ncol(factordata)
  K    <- ncol(data)
  P    <- K + no_factors



  # Declare Variables for storage

  # extract factors and join series

  fac <- get_factors(factordata,no_factors)
  xy  <- cbind(data,factordata)
  fy  <- cbind(data,fac)

  L <- olssvd(xy,fy)

  resids <- xy - fy %*% L
  Sigma <- t(resids) %*% resids

  L <- t(L)

  # Prior on the measurement equation
  gammam  <- array(0.5,dim=c(P,N))
  Liprvar <- li_prvar * diag(1,P)

  # Initialize the MCMC algorithm
  fy_lagged <- lagdata(fy,nolags=no_lags,intercept=intercept)
  draw <- initialize_mcmc(priorObj,fy_lagged$y,fy_lagged$x)

  # Start the MCMC sampler

  for(ireps in 1:nreps){
    print(ireps)

    # Draw posterior on measurement equation
    if(priorm == 1){

      draw_measurement <- draw_posterior_normal(Liprvar,fy,xy,K,P,N,Sigma,L,alpha,beta)
      L <- draw_measurement$L
      Sigma <- draw_measurement$Sigma

    }
    else if(priorm == 2){

      draw_measurement <- draw_posterior_ssvs(fy,xy,K,P,N,Sigma,tau2,c2,gammam,alpha,beta,L)
      L <- draw_measurement$L
      Sigma <- draw_measurement$Sigma
      gammam <- draw_measurement$gammam

    }

    # Draw posterior for state equation
    draw <- draw_posterior(priorObj, fy_lagged$y, fy_lagged$x, previous = draw, stabletest = stabletest)

  } # End loop over MCMC sampler

}
