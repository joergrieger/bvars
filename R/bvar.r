
#' Estimate a bayesian vector autoregressive model.
#'
#' This function is the main function to estimate a bayesian VAR model for the TxK series mydata. To estimate a bayesian VAR model to user first has to choose a prior and select the parmameters for it and submit it to the function via priorObj. It should be noted that the data submitted to the bvar function and the prior have to be the same. The logical parameter stabletest tells the function whether to check if a draw of coefficients is stable, i.e. if the largest eigenvalue of the companion matrix smaller than 1. Furthermore, the parameters nreps, burnin and nthin determine the number of mcmc-draws and the how many of the draws are retained. The number of retained draws is (nreps - burnin)/nthin.
#'
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
#'
#' @seealso \code{\link{msvar}} for regime switching models and \code{\link{tvar}} for threshold models.
#' @export
#' @importFrom stats frequency
#' @importFrom stats is.ts
#' @importFrom stats ts
#' @details  This is the main function for estimating a Bayesian Vectorautoregressive model. The user has to supply the data by mydata and a previously defined prior via priorObj. Several standard priors such as the Minnesota prior or the Independent Normal-Wishart are provided and the user has to parameterize them.


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



