#' @export
#' @title bvar - estimate a bayesian vector autoregressive model
#' @param mydata the time series used for estimating the VAR model
#' @param priorObj a S3-object containing information about the prior
#' @param stabletest logical, flag to test whether a draw is stable or not
#' @param nreps number of draws for the mcmc sampler
#' @param burnin number of burnin-draws
#' @param nthin thinning parameter
#' @return returns an S3 object of the class "bvar" with the following fields
#'
#' `general_info` list with general information about the model
#'
#' `intercept` whether the model has an intercept or not
#'
#'  `nolags` number of lags in the model
#'
#'   `nreps` total number of draws
#'
#'   `burnin` number of burn-in draws
#'
#'   `nthin` the thinning parameter
#'
#'   `data_info` information about the data
#'
#'   `type` type of the data object (can be ts, xts or matrix)
#'
#'   `var_names` variable names
#'
#'   `mydata` the data itself
#'
#'   `mcmc_draws` the draws from the mcmc algortithm
#'
#'   `Alpha` an (K * p + Intercept) x K x (nreps - burnin) / nthin matrix with the draws for the VAR-coefficients. With K being the number of variables, p the number of lags and Intercept is 1 if the model has an intercept and 0 otherwise.
#'
#'   `Sigma` an K x K x (nreps - burnin) / nthin - matrix with the draws of the Variance-Covariance matrix
#'
#'   `additional_info` an array of length (nreps - burnin) / nthin of lists with any additional information returned by the posterior.

bvar <- function(mydata,priorObj,stabletest = FALSE, nreps = 15000, burnin = 5000, nthin = 1){

  # Declare Variables

  K    <- ncol(mydata)
  Time <- nrow(mydata)
  obs  <- Time - priorObj$nolags
  constant <- 0
  if(priorObj$intercept) constant <- 1

  # Variables for storage
  storevar   <- floor( ( nreps - burnin ) / nthin )
  addInfo    <- array(list(),dim=c(storevar))
  Alphadraws <- array(NA,dim=c(K * priorObj$nolags + constant,K,storevar))
  Sigmadraws <- array(NA,dim=c(K,K,storevar))

  # Initialize MCMC algorithm
  tmp <- lagdata(mydata, nolags = priorObj$nolags, intercept = priorObj$intercept)
  draw <- initialize_mcmc(priorObj,tmp$y,tmp$x)

  # start the MCMC sampler

  for(ireps in 1:nreps){

    if(ireps %% 1000 == 0){

      cat("draw no.",ireps,"\n")

    }

    draw <- draw_posterior(priorObj, tmp$y, tmp$x, previous = draw, stabletest = stabletest)

    if(ireps > burnin){

      if( (ireps - burnin) %% nthin == 0){

        Alphadraws[,,( ireps - burnin ) / nthin] <- draw$Alpha
        Sigmadraws[,,( ireps - burnin ) / nthin] <- draw$Sigma
        addInfo[[(ireps - burnin) / nthin]] <- draw$addInfo

      }

    }


  }

  # Store values

  general_information <- list(intercept = priorObj$intercept,
                              nolags    = priorObj$nolags,
                              nreps     = nreps,
                              burnin    = burnin,
                              nthin     = nthin)


  # Information about the data

  if(sum(class(mydata) == "xts") > 0){

    tstype = "xts"
    var_names = colnames(mydata)

  }
  else if(sum(class(mydata) == "ts") > 0){

    tstype = "ts"
    var_names = colnames(mydata)

  }
  else if(is.matrix(mydata)){

    tstype = "matrix"
    var_names = colnames(mydata)

  }

  data_info <- list(type         = tstype,
                    var_names    = var_names,
                    data         = mydata,
                    no_variables = K)

  draw_info <- list(Alpha = Alphadraws,
                    Sigma = Sigmadraws,
                    additional_info = addInfo )

  # Return information
  ret_object <- structure(list(general_info = general_information,
                               data_info    = data_info,
                               mcmc_draws   = draw_info ),
                          class = "bvar")

  return(ret_object)

}


#' @export
#' @title Function to calculate irfs
#' @param obj an S3 object of class bvar
#' @param id_obj an S3 object with information about identifiaction of the model
#' @param nhor horizon of the impulse-response function
#' @param irfquantiles quantiles for the impulse-response functions
#' @param ncores number of cores used
#' @param ... currently not used
#'
#' @return returns an S3-object of the class
irf.bvar  <- function(obj,id_obj, nhor = 12, ncores = 1, irfquantiles = c(0.05,0.95),...){

  # Declare variables

  intercept  <- obj$general_info$intercept
  Alphadraws <- obj$mcmc_draws$Alpha
  Sigmadraws <- obj$mcmc_draws$Sigma

  nolags <- obj$general_info$nolags
  nreps  <- ( obj$general_info$nreps - obj$general_info$burnin ) / obj$general_info$nthin
  K      <- dim(Sigmadraws[,,1])[1]

  irfdraws <- array(0,dim=c(K,K,nhor,nreps))
  irffinal <- array(0,dim=c(K,K,nhor,3))

  if(ncores>1 && !require(doParallel)){

    stop("The parallel package has to be installed")

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

    xtmp <- foreach(ii = 1:nreps) %dopar% {

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

        irffinal[jj,kk,ll,1] <- median(irfdraws[jj,kk,ll,])
        irffinal[jj,kk,ll,2] <- quantile(irfdraws[jj,kk,ll,],probs=irflower)
        irffinal[jj,kk,ll,3] <- quantile(irfdraws[jj,kk,ll,],probs=irfupper)

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

