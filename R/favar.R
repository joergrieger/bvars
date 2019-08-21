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

#' @export
#' @title Function to calculate irfs
#' @param obj an S3 object of class favar
#' @param id_obj an S3 object with information about identifiaction of the model
#' @param nhor horizon of the impulse-response function
#' @param irfquantiles quantiles for the impulse-response functions
#' @param ncores number of cores used
#' @param ... currently not used
#'
#' @return returns an S3-object of the class fvirf
irf.favar <- function(obj,id_obj,nhor=12,ncores=1,irfquantiles = c(0.05,0.95),...){

  # Preliminaries
  intercept <- obj$general_info$intercept
  Betadraws <- obj$mcmc_draws$Alpha
  Sigmadraws <- obj$mcmc_draws$Sigma
  Ldraws <- obj$mcmc_draws$L
  nolags <- obj$general_info$nolags

  nreps <- dim(Betadraws)[3]
  k     <- dim(Sigmadraws)[1]
  dimXY <- dim(Ldraws)[1]

  irf_small_draws <- array(0,dim=c(k,k,nhor,nreps))
  irf_large_draws <- array(0,dim=c(k,dimXY,nhor,nreps))

  if(ncores>1 && !requireNamespace("foreach",quietly=TRUE)){

    stop("The foreach package cannot be loaded.")

  }

  if(ncores == 1){

    for(ii in 1:nreps){

      Alpha <- Betadraws[,,ii]
      Sigma <- Sigmadraws[,,ii]
      L     <- Ldraws[,,ii]

      irf <- compirf(Alpha = Alpha, Sigma = Sigma, id_obj = id_obj, nolags = nolags, intercept = intercept, nhor = nhor)
      irf_small_draws[,,,ii] <- irf

      for(jj in 1:k){

        irf_large_draws[jj,,,ii] <- L[,] %*% irf[jj,,]

      }


    }

  }
  else{

    # Get impulse-response for VAR-system
    # Register workers
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    `%dopar%` <- foreach::`%dopar%`

    xtmp <- foreach::foreach(ii = 1:nreps) %dopar% {

      Alpha <- Betadraws[,,ii]
      Sigma <- Sigmadraws[,,ii]

      irf <- compirf(Alpha = Alpha, Sigma = Sigma, id_obj = id_obj, nolags = nolags, intercept = intercept, nhor = nhor)

    } # End getting IRFs

    # Transform it to larger system
    for(ii in 1:nreps){

      L <- Ldraws[,,ii]
      irf_small_draws <- xtmp[[ii]]

      for(jj in 1:k){

        irf_large_draws[jj,,,ii] <- L[,] %*% irf_small_draws

      }

    }# end transformation

  }# End loop over parallel version

  # Store values
  IrfSmallFinal <- array(0,dim=c(k,k,nhor,3))
  IrfLargeFinal <- array(0,dim=c(k,dimXY,nhor,3))
  irflower <- min(irfquantiles)
  irfupper <- max(irfquantiles)

  for(jj in 1:k){
    for(kk in 1:k){
      for(ll in 1:nhor){

        IrfSmallFinal[jj,kk,ll,1] <- quantile(irf_small_draws[jj,kk,ll,],probs=0.5)
        IrfSmallFinal[jj,kk,ll,2] <- quantile(irf_small_draws[jj,kk,ll,],probs=irflower)
        IrfSmallFinal[jj,kk,ll,3] <- quantile(irf_small_draws[jj,kk,ll,],probs=irfupper)

      }
    }
  }
  for(jj in 1:k){
    for(kk in 1:dimXY){
      for(ll in 1:nhor){

        IrfLargeFinal[jj,kk,ll,1] <- quantile(irf_large_draws[jj,kk,ll,],probs=0.5)
        IrfLargeFinal[jj,kk,ll,2] <- quantile(irf_large_draws[jj,kk,ll,],probs=irflower)
        IrfLargeFinal[jj,kk,ll,3] <- quantile(irf_large_draws[jj,kk,ll,],probs=irfupper)

      }
    }
  }

  # Returning values
  relist <- structure(list(irf = IrfLargeFinal,irfhorizon = nhor,varnames = obj$data_info$varnames,
                           factor_varnames = obj$factordata_info$varnames),class = "fvirf")
  return(relist)

}
