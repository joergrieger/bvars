#' @export
#' @title Function to draw one single path for the impulse-response functions.
#' @param Alpha the estimate of the VAR-coefficients
#' @param Sigma the estimate of the Variance-covariance matrix
#' @param id_obj S3-object for identifying structural shocks
#' @param nolags Number of lags in the model
#' @param intercept Whether the model has an intercept or not
#' @param nhor the horizon of the impulse-response functions

compirf <- function(Alpha, Sigma,id_obj, nolags, intercept = TRUE, nhor){

  # Preliminaries

  K <- nrow(Sigma)
  bigj <- matrix(0, K, K * nolags)
  bigj[1:K,1:K] <- diag(K)

  if(intercept == TRUE){

    B <- Alpha[-c(1),]

  }
  else{

    B <- Alpha

  }

  PHI_Mat <- companionmatrix(B,nolags)
  biga <- PHI_Mat
  bigai <- biga

  # identify the model

  shock <- structural(id_obj, Alpha, Sigma)

  impresp <- matrix(0,K,K * nhor)
  impresp[1:K,1:K] <- shock

  for(ii in 1:(nhor-1)){

    impresp[,(ii * K + 1):(( ii + 1) * K)] <- (bigj %*% bigai %*% t(bigj) %*% shock)
    bigai <- bigai %*% biga

  }

  imp <- array(0,dim=c(K,K,nhor))
  jj <- 0

  for(ii in 1:K){
    for(ij in 1:nhor){

      jj <- ii + K * (ij - 1)
      imp[ii,,ij] <- impresp[,jj]

    }
  }

  return(imp)

}

#' @export
#' @title Compute Impulse-Response Functions
#' @param obj an S3 object of class bvar
#' @param id_obj an S3 object with information about identifiaction of the model
#' @param nhor horizon of the impulse-response function
#' @param irfquantiles quantiles for the impulse-response functions
#' @param ncores number of cores used
#' @param bootrep bootstrap replication for generalized impulse-response functions
#' @param ... currently not used
#' @return returns an S3-object of the class bvirf, fvirf, msirf or tirf
#' @details The function irf computes impulse-response functions for estimatet linear and non-linear VAR models. For Threshold-VARs and Factor-augmented Threshold-VARs generalized impulse-response functions are computed. The parameter bootrep determines how many bootstrap replications are going to be computed.
#' @rdname irf
irf <- function(obj,id_obj,nhor=12,ncores=1,irfquantiles = c(0.05,0.95),...) UseMethod("irf")

#' @importFrom stats quantile
#' @importFrom stats median
#' @export
#' @rdname irf
irf.bvar  <- function(obj, id_obj, nhor = 12, ncores = 1, irfquantiles = c(0.05,0.95),...){

  # Declare variables

  intercept  <- obj$general_info$intercept
  Alphadraws <- obj$mcmc_draws$Alpha
  Sigmadraws <- obj$mcmc_draws$Sigma

  nolags <- obj$general_info$nolags
  nreps  <- ( obj$general_info$nreps - obj$general_info$burnin ) / obj$general_info$nthin
  K      <- dim(Sigmadraws[,,1])[1]

  irfdraws <- array(0,dim=c(K,K,nhor,nreps))
  irffinal <- array(0,dim=c(K,K,nhor,3))

  #chck <- requireNamespace("doParallel") && requireNamespace("foreach") && requireNamespace("parallel")

  if(ncores>1 && !requireNamespace("foreach",quietly=TRUE)){

    stop("The foreach package cannot be loaded.")

  }

  if(ncores == 1){ # No parallelization

    for(ii in 1:nreps){

      Alpha <- Alphadraws[,,ii]
      Sigma <- Sigmadraws[,,ii]

      irf <- compirf(Alpha = Alpha, Sigma = Sigma, id_obj = id_obj, nolags = nolags, intercept = intercept, nhor = nhor)
      irfdraws[,,,ii] <- irf

    }

  }
  else{ # parallel version

    # Register workers
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    `%dopar%` <- foreach::`%dopar%`
    xtmp <- foreach::foreach(ii = 1:nreps) %dopar% {

      Alpha <- Alphadraws[,,ii]
      Sigma <- Sigmadraws[,,ii]

      irf <- compirf(Alpha = Alpha, Sigma = Sigma, id_obj = id_obj, nolags = nolags, intercept = intercept, nhor = nhor)

    }

    # stop workers
    parallel::stopCluster(cl)

    for(ii in 1:nreps){

      irfdraws[,,,ii] <- xtmp[[ii]]

    }

  }

  irflower <- min(irfquantiles)
  irfupper <- max(irfquantiles)
  for(jj in 1:K){
    for(kk in 1:K){
      for(ll in 1:nhor){

        irffinal[jj,kk,ll,1] <- stats::median(irfdraws[jj,kk,ll,])
        irffinal[jj,kk,ll,2] <- stats::quantile(irfdraws[jj,kk,ll,],probs=irflower)
        irffinal[jj,kk,ll,3] <- stats::quantile(irfdraws[jj,kk,ll,],probs=irfupper)

      }
    }
  }

  relist <- structure(list(irfdraws     = irffinal,
                           irfhorizon   = nhor,
                           varnames     = obj$data_info$var_names,
                           no_variables = obj$data_info$no_variables),
                      class = "bvirf")

  return(relist)

}


#' @export
#' @rdname irf
irf.favar <- function(obj, id_obj, nhor=12, ncores=1, irfquantiles = c(0.05,0.95),...){

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

#' @export
#' @rdname irf
irf.ftvar <- function(obj, id_obj, nhor=12, ncores = 1, irfquantiles = c(0.05,0.95), bootrep = 50,...){
  nLength    <- dim(obj$mcmc_draws$Alpha)[4]
  Alphadraws <- obj$mcmc_draws$Alpha
  Sigmadraws <- obj$mcmc_draws$Sigma
  Ldraws     <- obj$mcmc_draws$L
  thVar      <- obj$general_info$thVar
  tardraws   <- obj$mcmc_draws$tardraws
  deldraws   <- obj$mcmc_draws$deldraws
  nolags     <- obj$general_info$nolags
  intercept  <- obj$general_info$intercept

  P <- dim(Sigmadraws[,,,1])[1]
  dimxy <- dim(Ldraws)[1]

  irf_small_draws <- array(0,dim=c(P,P,nhor,2,nLength))
  irf_large_draws <- array(0,dim=c(P,dimxy,nhor,2,nLength))

  # get data
  factors <- obj$factordata_info$factors
  data    <- obj$data_info$data

  fy <- cbind(data,factors)

  if(ncores>1 && !requireNamespace("foreach",quietly = TRUE)){

    stop("The foreach package cannot be loaded.")

  }


  if(ncores == 1){ # No parallelization
    for(ii in 1:nLength){
      print(ii)
      Beta  <- Alphadraws[,,,ii]
      Sigma <- Sigmadraws[,,,ii]
      tart  <- tardraws[ii]
      L     <- Ldraws[,,,ii]
      thDelay <- deldraws[ii]

      xsplit <- splitVariables(y = fy, lags = nolags, thDelay = thDelay, thresh = thVar,
                               tart = tart, intercept = intercept)

      for(jj in 1:P){

        xx <- tirf(y = xsplit$ystar,ytest = xsplit$ytest,
                   beta1 = Beta[,,1],beta2 = Beta[,,2],
                   sigma1 = Sigma[,,1], sigma2 = Sigma[,,2],
                   tar = tart, thVar = thVar,thDelay = thDelay,
                   nolags = nolags, irfhor = nhor, intercept = intercept,
                   shockvar = jj,
                   bootrep = bootrep,
                   id_obj = id_obj)

        irf_small_draws[jj,,,1,ii] <- xx$irf1
        irf_small_draws[jj,,,2,ii] <- xx$irf2

        irf_large_draws[jj,,,1,ii] <- L[,,1] %*% irf_small_draws[jj,,,1,ii]
        irf_large_draws[jj,,,2,ii] <- L[,,2] %*% irf_small_draws[jj,,,2,ii]

      }


    }
  }
  else{

    # Register clusters
    # Warning: No check if number of cores are supported

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    irf_draws <- array(0,dim=c(P,P,nhor,2))
    `%dopar%` <- foreach::`%dopar%`

    xtmp <- foreach::foreach(ii = 1:nLength) %dopar% {

      Beta  <- Alphadraws[,,,ii]
      Sigma <- Sigmadraws[,,,ii]
      tart  <- tardraws[ii]
      L     <- Ldraws[,,,ii]
      thDelay <- deldraws[ii]

      irftmp <- tirf1(y = fy,
                      Alpha = Beta, Sigma = Sigma,
                      tart = tart, thVar = thVar, thDelay = thDelay,
                      nolags = nolags, nhor = nhor, intercept = intercept,
                      bootrep, id_obj = id_obj, P)

    } # End of parallel loop
    # shut down workers
    parallel::stopCluster(cl)

    for(ii in 1:nLength){

      irf_small_draws[,,,,ii] <- xtmp[[ii]]
      L <- Ldraws[,,,ii]

      for(jj in 1:P){
        irf_large_draws[jj,,,1,ii] <- L[,,1] %*% irf_small_draws[jj,,,1,ii]
        irf_large_draws[jj,,,2,ii] <- L[,,2] %*% irf_small_draws[jj,,,2,ii]
      }
    }# End projecting small system onto larger system
  }# End getting impulse-response-functions

  # Prepare return values
  IrfSmallFinal <- array(0,dim=c(P,P,nhor,2,3))
  IrfLargeFinal <- array(0,dim=c(P,dimxy,nhor,2,3))
  irflower <- min(irfquantiles)
  irfupper <- max(irfquantiles)

  for(jj in 1:P){
    for(kk in 1:P){
      for(ll in 1:nhor){

        # Regime 1
        IrfSmallFinal[jj,kk,ll,1,1] <- median(irf_small_draws[jj,kk,ll,1,])
        IrfSmallFinal[jj,kk,ll,1,2] <- quantile(irf_small_draws[jj,kk,ll,1,],probs=irflower)
        IrfSmallFinal[jj,kk,ll,1,3] <- quantile(irf_small_draws[jj,kk,ll,1,],probs=irfupper)

        # Regime 2
        IrfSmallFinal[jj,kk,ll,2,1] <- median(irf_small_draws[jj,kk,ll,2,])
        IrfSmallFinal[jj,kk,ll,2,2] <- quantile(irf_small_draws[jj,kk,ll,2,],probs=irflower)
        IrfSmallFinal[jj,kk,ll,2,3] <- quantile(irf_small_draws[jj,kk,ll,2,],probs=irfupper)

      }
    }
  }

  for(jj in 1:P){
    for(kk in 1:dimxy){
      for(ll in 1:nhor){
        # Regime 1
        IrfLargeFinal[jj,kk,ll,1,1] <- mean(irf_large_draws[jj,kk,ll,1,])
        IrfLargeFinal[jj,kk,ll,1,2] <- quantile(irf_large_draws[jj,kk,ll,1,],probs=irflower)
        IrfLargeFinal[jj,kk,ll,1,3] <- quantile(irf_large_draws[jj,kk,ll,1,],probs=irfupper)

        # Regime 2
        IrfLargeFinal[jj,kk,ll,2,1] <- mean(irf_large_draws[jj,kk,ll,2,])
        IrfLargeFinal[jj,kk,ll,2,2] <- quantile(irf_large_draws[jj,kk,ll,2,],probs=irflower)
        IrfLargeFinal[jj,kk,ll,2,3] <- quantile(irf_large_draws[jj,kk,ll,2,],probs=irfupper)

      }
    }
  }

  relist <- structure(list(irf=IrfLargeFinal,irfhorizon=nhor,varnames=obj$data_info$varnames,factor_varnames = obj$factordata_info$var_names),class = "tvirf")
  return(relist)
}

#' @export
#' @rdname irf
irf.msvar <- function(obj, id_obj, nhor = 12, ncores = 1, irfquantiles = c(0.05,0.95),...){

  #
  # Extract information Variables
  #

  intercept <- obj$general_info$intercept
  nolags    <- obj$general_info$nolags
  nreps     <- obj$general_info$nreps
  burnin    <- obj$general_info$burnin
  nthin     <- obj$general_info$nthin
  noregimes <- obj$general_info$noregimes

  Alphadraws <- obj$mcmc_draws$Alpha
  Sigmadraws <- obj$mcmc_draws$Sigma

  K <- dim(obj$mcmc_draws$Sigma)[2]

  totdraws <- floor((nreps - burnin) / nthin)

  irfdraws <- array(NA,dim=c(K,K,nhor,noregimes,totdraws))
  irffinal <- array(NA, dim =c(K,K,nhor,3,noregimes))

  # Write a function check_parallel

  if(ncores>1 && !requireNamespace("foreach",quietly = TRUE)){

    stop("The foreach package cannot be loaded.")

  }

  if(ncores == 1){ # No Parallelization

    for(jj in 1:noregimes){ # Loop over regimes
      for(ii in 1:totdraws){ # Loop over all draws

        Alpha <- Alphadraws[,,jj,ii]
        Sigma <- Sigmadraws[,,jj,ii]

        irf <- compirf(Alpha = Alpha, Sigma = Sigma, id_obj = id_obj,nolags = nolags, intercept = intercept, nhor = nhor)
        irfdraws[,,,jj,ii] <- irf


      } # End loop over all draws

    } # End loop over regimes

  }
  else{

    for(jj in 1:noregimes){ # Loop over regimes, not parallelized

      cl <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cl)

      `%dopar%` <- foreach::`%dopar%`
      xtmp <- foreach::foreach(ii = 1:totdraws) %dopar% {

        Alpha <- Alphadraws[,,jj,ii]
        Sigma <- Sigmadraws[,,jj,ii]

        irf <- compirf(Alpha = Alpha, Sigma = Sigma, id_obj = id_obj, nolags = nolags, intercept = intercept, nhor = nhor)

      } # End parallel loop

      # stop workers
      parallel::stopCluster(cl)

      for(ii in 1:totdraws){

        irfdraws[,,,jj,ii] <- xtmp[[ii]]

      }


    }

  }

  irflower <- min(irfquantiles)
  irfupper <- max(irfquantiles)
  for(ii in 1:noregimes){

    for(jj in 1:K){
      for(kk in 1:K){
        for(ll in 1:nhor){

          irffinal[jj,kk,ll,1,ii] <- median(irfdraws[jj,kk,ll,ii,])
          irffinal[jj,kk,ll,2,ii] <- quantile(irfdraws[jj,kk,ll,ii,],probs=irflower)
          irffinal[jj,kk,ll,3,ii] <- quantile(irfdraws[jj,kk,ll,ii,],probs=irfupper)

        }
      }
    }

  }

  relist <- structure(list(irfdraws     = irffinal,
                           irfhorizon   = nhor,
                           varnames     = obj$data_info$var_names,
                           no_variables = obj$data_info$no_variables,
                           noregimes    = noregimes),
                      class = "msirf")
  return(relist)

}

#' @export
#' @rdname irf
irf.tvar  <- function(obj, id_obj, nhor = 12, ncores = 1, irfquantiles = c(0.05,0.95), bootrep=50,...){

  nLength    <- (obj$general_info$nreps - obj$general_info$burnin) / obj$general_info$nthin

  Alphadraws <- obj$mcmc_draws$Alpha
  Sigmadraws <- obj$mcmc_draws$Sigma
  thVar      <- obj$general_info$thVar

  tardraws   <- obj$mcmc_draws$tardraws
  deldraws   <- obj$mcmc_draws$deldraws

  nolags     <- obj$general_info$nolags
  intercept  <- obj$general_info$intercept

  if(ncores>1 && !requireNamespace("foreach",quietly = TRUE)){

    stop("The foreach package cannot be loaded.")

  }


  K <- dim(Sigmadraws[,,,1])[1]

  irfdraws <- array(0,dim=c(K,K,nhor,2,nLength))

  if(ncores == 1){

    for(jj in 1:nLength){
      print(jj)


      Alpha   <- Alphadraws[,,,jj]
      Sigma   <- Sigmadraws[,,,jj]
      tart    <- tardraws[jj]
      thDelay <- deldraws[jj]

      xsplit <- splitVariables(y = obj$data_info$data,lags = nolags, thDelay = thDelay, thresh = thVar, tart = tart, intercept = intercept)


      for(ii in 1:K){

        xx <- tirf(xsplit$ystar,xsplit$ytest,Alpha[,,1],Alpha[,,2],Sigma[,,1],Sigma[,,2],tart,thVar,thDelay,nolags,nhor,intercept,shockvar=ii,bootrep,id_obj)

        irfdraws[ii,,,1,jj] <- xx$irf1
        irfdraws[ii,,,2,jj] <- xx$irf2

      }

    }
  }
  else{ # Parallel version

    # Register workers
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    irftmp <- array(0,dim=c(K,K,nhor,2))

    `%dopar%` <- foreach::`%dopar%`
    xtmp <- foreach::foreach(jj = 1:nLength) %dopar% {

      Alpha <- Alphadraws[,,,jj]
      Sigma <- Sigmadraws[,,,jj]
      tart <- tardraws[jj]
      thDelay <- deldraws[jj]

      irftmp <- tirf1(y = obj$data_info$data,
                      Alpha = Alpha,Sigma = Sigma,
                      tart = tart,thVar = thVar, thDelay = thDelay,
                      nolags = nolags, nhor = nhor, intercept = intercept,
                      bootrep = bootrep, K = K, id_obj = id_obj)

    }

    # stop workers
    parallel::stopCluster(cl)

    for(jj in 1:nLength){

      irfdraws[,,,,jj] <- xtmp[[jj]]

    }

  }

  irffinal <- array(0,dim=c(K,K,nhor,2,3))

  lowerquantile=min(irfquantiles)
  upperquantile=max(irfquantiles)
  for(ii in 1:K){
    for(jj in 1:K){
      for(kk in 1:nhor){
        irffinal[ii,jj,kk,1,1] <- quantile(irfdraws[ii,jj,kk,1,],probs=0.5)
        irffinal[ii,jj,kk,2,1] <- quantile(irfdraws[ii,jj,kk,2,],probs=0.5)

        irffinal[ii,jj,kk,1,2] <- quantile(irfdraws[ii,jj,kk,1,],probs=lowerquantile)
        irffinal[ii,jj,kk,2,2] <- quantile(irfdraws[ii,jj,kk,2,],probs=lowerquantile)

        irffinal[ii,jj,kk,1,3] <- quantile(irfdraws[ii,jj,kk,1,],probs=upperquantile)
        irffinal[ii,jj,kk,2,3] <- quantile(irfdraws[ii,jj,kk,2,],probs=upperquantile)
      }
    }
  }

  relist <- structure(list(irf=irffinal,irfhorizon=nhor,varnames=obj$general_info$var_names),class = "tvirf")

  return(relist)

}

