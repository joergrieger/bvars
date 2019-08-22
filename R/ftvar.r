#' @export
#' @title bayesian estimation of threshold VAR
#' @param mydata data
#' @param factordata data to extract the factors from
#' @param priorObj S3 object of the prior
#' @param thMax maximum delay of threshold variable
#' @param thVar the threshold variable
#' @param nreps total number of mcmc draws
#' @param burnin number of burn-in draws
#' @param nthin thinning parameter
#' @param stabletest test for stability
#' @param priorm Selects the prior on the measurement equation, 1=Normal-Gamma Prior and 2=SSVS prior.
#' @param alpha,beta prior on the variance of the measurement equation
#' @param tau2 variance of the coefficients in the measurement equation (only used if priorm=2)
#' @param c2 factor for the variance of the coefficients (only used if priorm=2)
#' @param li_prvar prior on variance of coefficients (only used if priorm = 1)
#'
ftvar <- function(mydata,factordata,priorObj,thMax,thVar,nreps,burnin,nthin,stabletest,alpha,beta,tau2,c2,li_prvar,priorm){

  if(priorObj$intercept){

    constant = 1

  }
  else{

    constant = 0

  }

  mydata     <- scale(mydata)
  factordata <- scale(factordata)

  intercept = priorObj$intercept

  nobs <- nrow(mydata)
  N <- ncol(factordata)
  K <- ncol(mydata)
  P <- K + priorObj$nofactors
  #
  #  Declare variables for the thresholding
  #

  thDelay  <- thMax # the delay
  tard     <- seq(1:thMax)
  startest <- max(thMax,priorObj$nolags)
  ytest    <- mydata[(startest + 1 - thDelay):(nobs - thDelay),thVar]
  tarmean  <- mean(ytest)
  tarstandard <- sqrt(stats::var(ytest))

  #
  # Declar Variables for storage
  #
  startest <- max(thMax,priorObj$nolags)


  # Draws for VAR-coefficient
  Alphadraws <- array(NA,dim = c(P * priorObj$nolags + intercept, P, 2, (nreps - burnin) / nthin))
  Sigmadraws <- array(NA, dim = c(P, P, 2, (nreps - burnin) / nthin))
  addInfo    <- array(list(),dim=c(2, (nreps - burnin) / nthin))

  # Measurement equation
  Ldraws     <- array(0,dim=c(ncol(mydata)+ncol(factordata),P,2,(nreps - burnin)/nthin))
  Sigma_m_draws <- array(0,dim=c(ncol(mydata)+ncol(factordata),ncol(mydata)+ncol(factordata),2,(nreps - burnin)/nthin))
  tardraws   <- array(NA,dim = c((nreps-burnin) / nthin))
  deldraws   <- array(NA,dim = c((nreps-burnin) / nthin))
  totRegimes <- nobs - (startest + 1)
  regimes    <- array(NA,dim = c(nobs - (startest + 1), (nreps - burnin) / nthin))
  gammamdraws <- array(NA,dim=c(P,N+K,2,(nreps - burnin) / nthin))

  #
  # get factors and put the model in appropriate form
  # Also get initial values for mcmc algor
  #
  factors <- get_factors(factordata,priorObj$nofactors)
  xy      <- cbind(mydata,factordata)
  fy      <- cbind(mydata,factors)
  Li      <- olssvd(xy,fy)

  res <- xy - fy %*% Li
  sse <- t(res) %*% res / nobs

  L             <- array(0,dim=c(ncol(xy),P,2))
  Sigma_measure <- array(0,dim=c(ncol(xy),ncol(xy),2))

  for(ii in 1:2){

    Sigma_measure[,,ii] <- sse
    L[,,ii] <- Li

  }

  # Initialize mcmc algorithm
  tart <- tarmean

  xsplit <- splitVariables(y = fy, lags = priorObj$nolags, thDelay = thDelay, thresh = thVar, tart = tart, intercept = intercept)
  prev <- array(list(),dim=c(2))

  prev[[1]] <- initialize_mcmc(priorObj,xsplit$y1,xsplit$x1)
  prev[[2]] <- initialize_mcmc(priorObj,xsplit$y2,xsplit$x2)

  # Initial values for the measurement equation
  # Prior on the measurement equation
  gammam  <- array(0.5,dim=c(P,N+K))
  Liprvar <- li_prvar * diag(1,P)

  for(ireps in 1:nreps){
    print(ireps)


    # split variables and make fy and xy the same length for both regimes
    fy_split <- splitVariables(y = fy, lags = priorObj$nolags, thDelay = thDelay, thresh = thVar, tart = tart, intercept = intercept)
    nr  <- nrow(xy)
    nr2 <- nrow(as.matrix(fy_split$ytest,ncol = 1))
    rdiff <- nr - nr2 + 1
    xyred <- xy[rdiff:nr,]


    # Draw posterior for measurement equation
    # first regime
    xy_split <- xyred[fy_split$e1,]
    draw_measurement <- draw_posterior_normal(li_prvar = Liprvar, fy = fy_split$y1, xy = xy_split,
                                              K = K,P = P,N = N,Sigma = Sigma_measure[,,1],L = L[,,1], alpha, beta)
    L[,,1] <- draw_measurement$L
    Sigma_measure[,,1] <- draw_measurement$Sigma



    # second regime
    xy_split <- xyred[!fy_split$e1,]
    draw_measurement <- draw_posterior_normal(Liprvar,fy_split$y2,xy_split,K,P,N,Sigma_measure[,,2],L[,,2],alpha,beta)
    L[,,2] <- draw_measurement$L
    Sigma_measure[,,2] <- draw_measurement$Sigma

    # Sample VAR-parameters

    prev[[1]] <- draw_posterior(priorObj,fy_split$y1,fy_split$x1,previous = prev[[1]], stabletest = stabletest)
    prev[[2]] <- draw_posterior(priorObj,fy_split$y2,fy_split$x2,previous = prev[[2]], stabletest = stabletest)

    # Sample new threshold
    tarnew <- tart + rnorm(1, sd = tarstandard) # new suggestion for threshold


    l1post <- tarpost(fy_split$xstar, fy_split$ystar, Ytest = fy_split$ytest,
                      prev[[1]]$Alpha, prev[[2]]$Alpha,
                      prev[[1]]$Sigma, prev[[2]]$Sigma,
                      tarnew, priorObj$nolags, intercept, tarmean, tarstandard)

    l2post <- tarpost(fy_split$xstar, fy_split$ystar, Ytest = fy_split$ytest,
                      prev[[1]]$Alpha, prev[[2]]$Alpha,
                      prev[[1]]$Sigma, prev[[2]]$Sigma,
                      tart, priorObj$nolags, intercept, tarmean, tarstandard)
    acc <- min(1, exp(l1post$post - l2post$post))

    if(is.na(acc)){acc = 0}

    u <- runif(1)
    if(u < acc) tart = tarnew

    # Sample new delay parameters

    prob <- matrix(0,nrow = thMax)
    # Loop over all potential threshold delays
    for(jj in 1:thMax){

      split1 <- splitVariables(y = fy, lags = priorObj$nolags, jj, thVar, tart, intercept)

      x <- exptarpost(split1$xstar,split1$ystar, split1$ytest,
                      prev[[1]]$Alpha, prev[[2]]$Alpha,
                      prev[[1]]$Sigma, prev[[2]]$Sigma,
                      tart, priorObj$nolags, intercept, tarmean,
                      tarstandard, ncrit = 0.05)

      prob[jj,1] <- exp(x$post)


    }

    prob <- prob/sum(prob)

    if(anyNA(prob)){

      prob <- matrix(1/thMax,nrow=thMax)

    }

    thDelay <- sample(thMax,1,replace=FALSE,prob)

    # Store values

    if(ireps > burnin && (ireps - burnin) %% nthin == 0){

      for(ii in 1:2){

        Alphadraws[,,ii,( ireps - burnin ) / nthin] <- prev[[ii]]$Alpha
        Sigmadraws[,,ii,( ireps - burnin ) / nthin] <- prev[[ii]]$Sigma

        if(!is.null((prev[[ii]]$addInfo))){

          addInfo[[ii, (ireps - burnin) / nthin]]    <- prev[[ii]]$addInfo

        }

        # Draws for measurement equation
        Ldraws[,,ii,(ireps - burnin) / nthin] <- L[,,ii]
        Sigma_m_draws[,,ii,(ireps - burnin) / nthin] <- Sigma_measure[ii]

      }

      # Regimes
      nT <- length(xsplit$e1)
      a  <- nT-totRegimes
      regimes[ ,(ireps- burnin) / nthin] <- xsplit$e1[(1+a):nT]

      deldraws[(ireps - burnin) / nthin] <- thDelay
      tardraws[(ireps - burnin) / nthin] <- tart



    } # End Storing values
  } # End loop over mcmc

  general_information <- list(intercept = priorObj$intercept,
                              nolags    = priorObj$nolags,
                              nofactors = priorObj$nofactors,
                              nreps     = nreps,
                              burnin    = burnin,
                              nthin     = nthin,
                              thVar     = thVar,
                              thMax     = thMax)
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
                          factors      = factors,
                          no_variables = K)

  draw_info <- list(Alpha = Alphadraws,
                    Sigma = Sigmadraws,
                    L     = Ldraws,
                    Sigmam = Sigma_measure,
                    additional_info = addInfo,
                    regimes = regimes,
                    deldraws = deldraws,
                    tardraws = tardraws)
  # Return information
  ret_object <- structure(list(general_info = general_information,
                               data_info    = data_info,
                               factordata_info = factordata_info,
                               mcmc_draws   = draw_info ),
                          class = "ftvar")

  return(ret_object)




} # End function

#' @export
#' @title impulse-response functions for threshold VAR-Models
#' @param obj an S3 object of class bvar
#' @param id_obj an S3 object with information about identifiaction of the model
#' @param nhor horizon of the impulse-response function
#' @param irfquantiles quantiles for the impulse-response functions
#' @param ncores number of cores used
#' @param bootrep the number of bootstrap replications
#' @param ... currently not used
#'
#' @return returns an S3-object of the class tvirf
irf.ftvar <- function(obj, id_obj, nhor=12, irfquantiles = c(0.05,0.95), ncores = 1,bootrep = 50,...){
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

