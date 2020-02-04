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
#'
#' @importFrom stats runif
#' @importFrom stats var
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
  tarstandard <- sqrt(stats::var(ytest))


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

    tarnew <- tart + stats::rnorm(1,sd = tarstandard) # proposal for new threshold value


    l1post <- tarpost(xsplit$xstar, xsplit$ystar, Ytest = xsplit$ytest,
                      prev[[1]]$Alpha, prev[[2]]$Alpha,
                      prev[[1]]$Sigma, prev[[2]]$Sigma,
                      tarnew, nolags, intercept, tarmean, tarstandard)

    l2post <- tarpost(xsplit$xstar, xsplit$ystar, Ytest = xsplit$ytest,
                      prev[[1]]$Alpha, prev[[2]]$Alpha,
                      prev[[1]]$Sigma, prev[[2]]$Sigma,
                      tart, nolags, intercept, tarmean, tarstandard)



    acc <- min(1,exp(l1post$post - l2post$post))
    u   <- stats::runif(1)

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

    thDelay <- sample(thMax,1,replace=FALSE,prob)

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


