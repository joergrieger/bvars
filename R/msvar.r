
#' Estimate regime-switching models with fixed transition probabilities.
#'
#' The function msvar estimates a regime-switching models with fixed transition probabilities. To estimate the msvar-model the user has to provide the data in mydata, which can be a simple TxK matrix or a ts or xts object. Furthermore the user provides a prior via priorObj. The user can also specify the number of regimes with the parameter noregimes. However, in order to provide a good estimate of the msvar model the number of regimes shouldn't be too high. The logical parameter stabletest tells the function whether to check the eigenvalue of the associated companion matrices of the model. If TRUE, the model will draw coefficients until the largest eigenvalue of the companion matrix is smaller than one. The total number of draws is given by the parameter nreps and the number of retained draws is (nreps-burnin)/nthin.
#'
#' @param mydata data
#' @param priorObj  S3 object containing information about the prior used.
#' @param stabletest logical, test for stability of each draw of the VAR-coefficients
#' @param noregimes Number of regimes
#' @param nreps Total number of draws
#' @param burnin number of burn-in draws
#' @param nthin thinning parameter
#' @seealso \code{\link{bvar}} for BVAR-Models and \code{\link{tvar}} for threshold models.
#' @export

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
  Regimedraws <- array(NA,dim=c(obs,storevar))
  trans_mat_draws  <- array(dim=c(noregimes,noregimes,storevar))

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

    if(ireps %% 1000 == 0){

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

        # Store draw of coefficients

        for(ii in 1:noregimes){

          Alphadraws[,,ii, (ireps - burnin) / nthin] <- prev[[ii]]$Alpha
          Sigmadraws[,,ii, (ireps - burnin) / nthin] <- prev[[ii]]$Sigma

          if(!is.null((prev[[ii]]$addInfo))){

            addInfo[[ii, (ireps - burnin) / nthin]]    <- prev[[ii]]$addInfo

          }


        }
        # Store draws of regimes, smoothed probabilities and transition probabilities
        trans_mat_draws[,,(ireps-burnin)/nthin] <- transmat
        Regimedraws[,(ireps-burnin)/nthin] <- stt

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
                    additional_info = addInfo,
                    transmat = trans_mat_draws,
                    regimes = Regimedraws)

  # Return information
  ret_object <- structure(list(general_info = general_information,
                               data_info    = data_info,
                               mcmc_draws   = draw_info ),
                          class = "msvar")

  return(ret_object)


}





