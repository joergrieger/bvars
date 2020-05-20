#' @export
#' @title Factor-Augmented Vector Autoregression
#' @description Function to estimate factor-augmented vector autoregressions using a 2-step procedure.
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
#' @return A S3-object of the class favar.
#' @details Estimates a favar-model using a 2-step procedure. In the first step the factors are extracted from the series using principal components. In the second step, a VAR-model of order p of both the factor series and other variables is estimated. To uniquely identify the favar all series are normalized with mean 0 and standard deviation of 1. Furthermore, the variance-covariance matrix is assumed to be diagonal. The VAR-model is of the form
#' \deqn{\left[\begin{array}{c}Y_t\\ F_t\end{array}\right]=\Phi(L)\left[\begin{array}{c}Y_{t-1}\\ F_{t-1}\end{array}\right]+w_t}
#' and the observation equation takes the form
#' \deqn{X_t=\Lambda^fF_t+\Lambda^yY_t+e_t}
#' Since a model with \eqn{\tilde{\Lambda}^f=\Lambda^fH} and \eqn{\tilde{F}_t=H^{-1}F_t} are observationally equivalent to eqn{\Lambda,F} we impose the standard normalization restriction implicit in the principal components. One interpretation of the factors is that they are a diffusion index as in Stock and Watson (1998)
#' @references Bernanke, Ben S., Jean Boivin and Piotr Eliasz, Measuring the effects of monetary policy: a factor-augmented vector autoregressive (favar) approach
#' @references Stock, James and Mark Watson, Diffusion Indexes, NBER Working Paper No. 6702, 1998
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
  if(intercept == TRUE){

    constant = 1

  }
  else{

    constant = 0

  }
  storevar   <- floor( ( nreps - burnin ) / nthin )
  addInfo    <- array(list(),dim=c(storevar))

  Alphadraws <- array(NA,dim=c(P * priorObj$nolags + constant,P,storevar))
  Sigmadraws <- array(NA,dim=c(P,P,storevar))
  Ldraws     <- array(NA,dim=c(N+K,P,storevar))
  Sigma_measure <- array(NA,dim=c(N+K,N+K,storevar))
  gammam_draws <- array(NA,dim=c(P,N+K,storevar))

  # extract factors and join series

  fac <- get_factors(factordata,no_factors)
  xy  <- cbind(data,factordata)
  fy  <- cbind(data,fac)

  L <- olssvd(xy,fy)

  resids <- xy - fy %*% L
  Sigma <- t(resids) %*% resids

  L <- t(L)
  print(dim(L))

  # Prior on the measurement equation
  gammam  <- array(0.5,dim=c(P,N+K))
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

    # Store results
    if(ireps > burnin && (ireps - burnin) %% nthin == 0){

      Alphadraws[,,( ireps - burnin ) / nthin] <- draw$Alpha
      Sigmadraws[,,( ireps - burnin ) / nthin] <- draw$Sigma
      addInfo[[(ireps - burnin) / nthin]] <- draw$addInfo
      Ldraws[,,(ireps - burnin) / nthin] <- L
      Sigma_measure[,,(ireps - burnin) / nthin] <- Sigma

      if(priorm == 2){

        gammam_draws[,,( ireps - burnin)/nthin ] <- gammam

      }

    }

  } # End loop over MCMC sampler

  # Store results

  # general information
  general_information <- list(intercept = priorObj$intercept,
                              nolags    = priorObj$nolags,
                              nofactors = no_factors,
                              nreps     = nreps,
                              burnin    = burnin,
                              nthin     = nthin)

  # Information about the data

  if(sum(class(data) == "xts") > 0){

    tstype = "xts"
    var_names = colnames(data)

  }
  else if(sum(class(data) == "ts") > 0){

    tstype = "ts"
    var_names = colnames(data)

  }
  else if(is.matrix(data)){

    tstype = "matrix"
    var_names = colnames(data)

  }

  data_info <- list(type         = tstype,
                    var_names    = var_names,
                    data         = data,
                    no_variables = K)

  # Information about the data used to extract the factors

  if(sum(class(factordata) == "xts") > 0){

    tstype = "xts"
    var_names = colnames(factordata)

  }
  else if(sum(class(factordata) == "ts") > 0){

    tstype = "ts"
    var_names = colnames(factordata)

  }
  else if(is.matrix(factordata)){

    tstype = "matrix"
    var_names = colnames(factordata)

  }

  factordata_info <- list(type         = tstype,
                          var_names    = var_names,
                          data         = factordata,
                          no_variables = K)


  # The results of the mcmc draws

  draw_info <- list(Alpha = Alphadraws,
                    Sigma = Sigmadraws,
                    Sigma_measure = Sigma_measure,
                    L = Ldraws,
                    additional_info = addInfo )

  # Return information
  ret_object <- structure(list(general_info = general_information,
                               data_info    = data_info,
                               factordata_info = factordata_info,
                               mcmc_draws   = draw_info ),
                          class = "favar")

  return(ret_object)

}
