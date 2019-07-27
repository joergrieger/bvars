#' @export
#' @title bayesian estimation of threshold VAR
#' @param mydata data
#' @param priorObj S3 object of the prior
#' @param thMax maximum delay of threshold variable
#' @param thVar the threshold variable
#' @param nreps total number of mcmc draws
#' @param burnin number of burn-in draws
#' @param nthin thinning parameter
#' @param stabletest test for stability
tvar <- function(mydata,priorObj,thMax=2,thVar=1,nreps = 1100,burnin=100,nthin=1,stabletest = TRUE){

  #
  # Declare general variables
  #

  obs <- nrow(mydata)
  K   <- ncol(mydata)
  nolags    <- priorObj$nolags
  intercept <- priorObj$intercept

  #
  #  Declare variables for the thresholding
  #

  thDelay  <- thMax # the delay
  tard     <- seq(1:thMax)
  startest <- max(thMax,nolags)
  ytest    <- mydata[(startest + 1 - thDelay):(obs - thDelay),thVar]
  tarmean  <- mean(ytest)
  tarstandard <- sqrt(var(ytest))


  #
  # Declare variables to store results
  #

  Alphadraws <- array(NA,dim = c(K * nolags + intercept, K, 2, (nreps - burnin) / nthin))
  Sigmadraws <- array(NA, dim = c(K, K, 2, (nreps - burnin) / nthin))
  addInfo    <- array(list(),dim=c(2, (nreps - burnin) / nthin))
  tardraws   <- array(NA,dim = c((nreps-burnin) / nthin))
  deldraws   <- array(NA,dim = c((nreps-burnin) / nthin))
  totRegimes <- obs - (startest + 1)
  regimes    <- array(NA,dim = c(obs - (startest + 1), (nreps - burnin) / nthin))

  # Initialize MCMC-sampler

  tart <- tarmean

  xsplit <- splitVariables(y = mydata, lags = nolags, thDelay = thDelay, thresh = thVar, tart = tart, intercept = intercept)
  prev <- array(list(),dim=c(2))

  prev[[1]] <- initialize_mcmc(priorObj,xsplit$y1,xsplit$x1)
  prev[[2]] <- initialize_mcmc(priorObj,xsplit$y2,xsplit$x2)

  #
  # The loop over the gibbs sampler
  #

  for(ireps in 1:nreps){

    if(ireps %% 1 == 0){

      cat("draw no.",ireps,"\n")

    }

    #
    # Step 1: split variables
    #

    xsplit <- splitVariables(y = mydata, lags = nolags, thDelay = thDelay, thresh = thVar, tart = tart, intercept = intercept)

    #
    # Step 2: Sample posteriors
    #

    prev[[1]] <- draw_posterior(priorObj,xsplit$y1,xsplit$x1,previous = prev[[1]], stabletest = stabletest)
    prev[[2]] <- draw_posterior(priorObj,xsplit$y2,xsplit$x2,previous = prev[[2]], stabletest = stabletest)



    #
    # Step 3: sample new threshold using Random-Walk Metropolis-Hastings Algorithm
    #

    tarnew <- tart + rnorm(1,sd = tarstandard) # proposal for new threshold value


    l1post <- tarpost(xsplit$xstar, xsplit$ystar, Ytest = xsplit$ytest,
                      prev[[1]]$Alpha, prev[[2]]$Alpha,
                      prev[[1]]$Sigma, prev[[2]]$Sigma,
                      tarnew, nolags, intercept, tarmean, tarstandard)

    l2post <- tarpost(xsplit$xstar, xsplit$ystar, Ytest = xsplit$ytest,
                      prev[[1]]$Alpha, prev[[2]]$Alpha,
                      prev[[1]]$Sigma, prev[[2]]$Sigma,
                      tart, nolags, intercept, tarmean, tarstandard)



    acc <- min(1,exp(l1post$post - l2post$post))
    u   <- runif(1)

    if(u < acc){

      # accept proposal

      tart = tarnew

    }

    tarmean = tart

    #
    # Step 4: Sample new delay parameter
    #

    prob <- matrix(0,nrow = thMax)


    # Loop over all potential threshold delays
    for(jj in 1:thMax){

      split1 <- splitVariables(y = mydata, lags = nolags, jj, thVar, tart, intercept)

      x <- exptarpost(split1$xstar,split1$ystar, split1$ytest,
                      prev[[1]]$Alpha, prev[[2]]$Alpha,
                      prev[[1]]$Sigma, prev[[2]]$Sigma,
                      tart, nolags, intercept, tarmean,
                      tarstandard, ncrit = 0.05)

      prob[jj,1] <- exp(x$post)


    }

    prob <- prob/sum(prob)

    if(anyNA(prob)){

     prob <- matrix(1/thMax,nrow=thMax)

    }
    if(ireps > burnin){

      if( (ireps - burnin) %% nthin == 0){

        for(ii in 1:2){

          Alphadraws[,,ii, (ireps - burnin) / nthin] <- prev[[ii]]$Alpha
          Sigmadraws[,,ii, (ireps - burnin) / nthin] <- prev[[ii]]$Sigma

          if(!is.null((prev[[ii]]$addInfo))){

            addInfo[[ii, (ireps - burnin) / nthin]]    <- prev[[ii]]$addInfo

          }


        }

        # Regimes
        nT <- length(xsplit$e1)
        a  <- nT-totRegimes
        regimes[ ,(ireps- burnin) / nthin] <- xsplit$e1[(1+a):nT]

        deldraws[(ireps - burnin) / nthin] <- thDelay
        tardraws[(ireps - burnin) / nthin] <- tart

      }

    }

  }

  #
  # Final storage of data
  #

  general_information <- list(intercept = priorObj$intercept,
                              nolags    = priorObj$nolags,
                              nreps     = nreps,
                              burnin    = burnin,
                              nthin     = nthin,
                              thVar     = thVar)


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
                    additional_info = addInfo,
                    regimes = regimes,
                    deldraws = deldraws,
                    tardraws = tardraws)

  # Return information
  ret_object <- structure(list(general_info = general_information,
                               data_info    = data_info,
                               mcmc_draws   = draw_info ),
                          class = "tvar")

  return(ret_object)


}

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

irf.tvar  <- function(obj,id_obj, nhor = 12, ncores = 1, irfquantiles = c(0.05,0.95), bootrep=50,...){

  nLength    <- (obj$general_info$nreps - obj$general_info$burnin) / obj$general_info$nthin

  Alphadraws <- obj$mcmc_draws$Alpha
  Sigmadraws <- obj$mcmc_draws$Sigma
  thVar      <- obj$general_info$thVar

  tardraws   <- obj$mcmc_draws$tardraws
  deldraws   <- obj$mcmc_draws$deldraws

  nolags     <- obj$general_info$nolags
  intercept  <- obj$general_info$intercept


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

        xx <- tirf(xsplit$ystar,xsplit$ytest,Alpha[,,1],Alpha[,,2],Sigma[,,1],Sigma[,,2],tart,thVar,thDelay,nolags,nhor,intercept,shockvar=ii,bootrep)

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

    xtmp <- foreach(jj = 1:nLength) %dopar% {

      Alpha <- Alphadraws[,,,jj]
      Sigma <- Sigmadraws[,,,jj]
      tart <- tardraws[jj]
      thDelay <- deldraws[jj]

      irftmp <- tirf(obj$mydata,Alpha,Sigma,tart,thVar,thDelay,nolags,nhor,intercept,bootrep,K)

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
