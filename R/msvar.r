#' @export
#' @title estimate regime-switching models with fixed transition probabilities
#' @param mydata data
#' @param priorObj  S3 object containing information about the prior used.
#' @param stabletest test for stability of each draw of the VAR-coefficients
#' @param noregimes Number of regimes
#' @param nreps Total number of draws
#' @param burnin number of burn-in draws
#' @param nthin thinning parameter

msvar <- function(mydata,priorObj,stabletest = FALSE, noregimes = 2, nreps = 15000, burnin = 10000, nthin = 1){

  # Preliminaries

  stransmat <- matrix(0,noregimes,noregimes)
  alphaprior <- matrix(10,noregimes,noregimes)

  K    <- ncol(mydata)
  Time <- nrow(mydata)
  obs  <- Time - priorObj$nolags

  constant <- 0
  if(priorObj$intercept) constant <- 1

  # Variables for storage
  storevar   <- floor( ( nreps - burnin ) / nthin )
  addInfo    <- array(list(),dim=c(noregimes, storevar))
  Alphadraws <- array(NA,dim=c(K * priorObj$nolags + constant, K, noregimes, storevar))
  Sigmadraws <- array(NA,dim=c(K,K,noregimes, storevar))

  Alpha2 <- array(NA,dim=c(K * priorObj$nolags + constant, K ,noregimes))
  Sigma2 <- array(NA,dim=c(K,K,noregimes))

  # Initialize the Gibbs sampler for the VAR coefficients
  tmp <- lagdata(mydata, nolags = priorObj$nolags, intercept = priorObj$intercept)
  yLagged <- tmp$y
  xLagged <- tmp$x

  prev <- array(list(),dim=c(noregimes))

  for(ii in 1:noregimes){

    prev[[ii]] <- initialize_mcmc(priorObj,tmp$y,tmp$x)

  }

  for(ireps in 1:nreps){

    if(ireps %% 10 == 0){

      cat("draw no.",ireps,"\n")

    }

    # Get posterior for Dirichlet prior on transition probabilities
    alpha    <- alphaprior + stransmat - matrix(1,noregimes,noregimes)
    transmat <- t(apply(alpha,1,extraDistr::rdirichlet,n=1))

    # Run Hamilton Filter and draw States
    for(ii in 1:noregimes){

      Alpha2[,,ii] <- prev[[ii]]$Alpha
      Sigma2[,,ii] <- prev[[ii]]$Sigma

    }

    filteredprob <- hamiltonfilter(Alpha2,Sigma2,transmat,yLagged, xLagged, h=noregimes)
    stt <- getst(filteredprob$fprob,transmat,h=noregimes)

    # count transitions
    nseq <- countseq(stt,noregimes)
    stransmat <- nseq

    # Draw VAR coefficients and Variance-Covariance for all states

    for(ii in 1:noregimes){

      yLagged_filt <- yLagged[stt==ii, ]
      xLagged_filt <- xLagged[stt==ii, ]

      prev[[ii]] <- draw_posterior(priorObj, yLagged,xLagged, previous = prev[[ii]], stabletest = stabletest)

    }

    if(ireps > burnin){

      if( (ireps - burnin) %% nthin == 0){

        for(ii in 1:noregimes){

          Alphadraws[,,ii, (ireps - burnin) / nthin] <- prev[[ii]]$Alpha
          Sigmadraws[,,ii, (ireps - burnin) / nthin] <- prev[[ii]]$Sigma

          if(!is.null((prev[[ii]]$addInfo))){

            addInfo[[ii, (ireps - burnin) / nthin]]    <- prev[[ii]]$addInfo

          }


        }

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
                              noregimes = noregimes)


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
                          class = "msvar")

  return(ret_object)


}

#' @export
#' @title impulse-response functions for regime-switching models
#' @param obj an S3 object of class bvar
#' @param id_obj an S3 object with information about identifiaction of the model
#' @param nhor horizon of the impulse-response function
#' @param irfquantiles quantiles for the impulse-response functions
#' @param ncores number of cores used
#' @param ... currently not used
#'
#' @return returns an S3-object of the class msirf

irf.msvar <- function(obj,id_obj,nhor,irfquantiles = c(0.05,0.95),ncores = 1,...){

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

  if(ncores>1 && !require(doParallel)){

    stop("The parallel package has to be installed")

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



      xtmp <- foreach(ii = 1:totdraws) %dopar% {

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
