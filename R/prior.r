#
#' @export
#' @title sets up conjugate Normal-Wishart prior
#' @param mydata a TxK xts-object needed for setting up the prior
#' @param nolags integer Number of lags in the VAR model
#' @param factordata for factor models additional time series for the factors are needed (not yet implemented)
#' @param no_factors number of factors in a factor model (not yet implemented)
#' @param coefprior double or () matrix with the prior for the VAR-coefficients. If only a scalar variable is provided the prior will be set to
#' @param coefpriorvar double or () matrix with the prior on the variance of the VAR-coefficients. If only a scalar is provided the prior will be set to diag(1,K)
#' @param varprior double or () matrix with the prior on Variance-Covariance matrix.
#' @param varpriordof integer. The degree of freedom for prior on the Variance-Covariance matrix.
#' @param intercept logical whether the VAR model has an intercept (TRUE) or not (FALSE)
#' @return returns an S3-object of class "cnw"
#' @details
#' This sets up the Normal-Wishart prior to use for bvar-estimation. The conjugate Normal-Wishart prior has the following form
#' \deqn{\alpha|\Sigma\sim N(\underline{\alpha},\Sigma\otimes \underline{V})}
#' and
#' \deqn{\Sigma^{-1}\sim W(\underline{S}^{-1},\underline{v})}
#' with \eqn{\underline{\alpha},\underline{V},\underline{v}}, and \eqn{\underline{S}} the prior hyperparameters chosen by the user.
#' @family priors
#' @author Joerg Rieger
#' @references K. Rao Kadiyala and Sune Karlsson, Numerical Methods for Estimation and Inference in Bayesian VAR-Models, Journal of Applied Econometrics 12(2), 99-132
#' @references Gary Koop and Dimitris Korobilis (2010), Bayesian Multivariate Time Series Methods for Empirical Macroeconomics, Foundations and Trends in Econometrics 3(4), 267-358

set_prior_cnw <- function(mydata = NULL, factordata = NULL, no_factors = 0, coefprior = NULL,
                          coefpriorvar = NULL, varprior = NULL, varpriordof = NULL,
                          nolags = 1, intercept = TRUE){

  # Input checks

  # Are their any NA elements in the data
  if(anyNA(mydata)){

    stop("There are NA elements in the data")

  }

  # Check for inconsistencies between factor data and no. of factors
  if(is.null(factordata) && no_factors > 0){

    stop("No data for factor model but number of factors > 0")

  }

  # Preliminary calculations

  K <- dim(mydata)[2] + no_factors # get dimension of time series

  constant = 0
  if(intercept) constant = 1

  # Prior on coefficients

  coefprior    <- array( 1,dim=c(K * nolags + constant,K )) * coefprior # coefficients
  coefpriorvar <- diag( 1,K * nolags + constant ) * coefpriorvar # variance on coefficients

  # prior for variance-covariance matrix
  varprior    <- diag(1,K) * varprior
  varpriordof <- varpriordof

  pr <- list(type         = "cnw",
             nolags       = nolags,
             nofactors    = no_factors,
             intercept    = intercept,
             coefprior    = coefprior,
             coefpriorvar = coefpriorvar,
             varprior     = varprior,
             varpriordf   = varpriordof)

  pr <- structure(pr, class = "cnw")

  return(pr)

}

#' @export
#' @title set up uninformative prior
#' @param mydata data
#' @param factordata data for factor models  (not yet implemented)
#' @param no_factors number of factors (not yet implemented)
#' @param nolags number of lags
#' @param intercept whether the model has an intercept
#' @return returns an S3 object of the class "unf"
#' @family priors
#' @author Joerg Rieger
#' @references K. Rao Kadiyala and Sune Karlsson, Numerical Methods for Estimation and Inference in Bayesian VAR-Models, Journal of Applied Econometrics 12(2), 99-132
#' @references Gary Koop and Dimitris Korobilis (2010), Bayesian Multivariate Time Series Methods for Empirical Macroeconomics, Foundations and Trends in Econometrics 3(4), 267-358


set_prior_uninformative <- function(mydata=NULL,factordata=NULL,no_factors=0,nolags=1,intercept=TRUE){

  pr <- list(type         = "unf",
             nolags       = nolags,
             nofactors    = no_factors,
             intercept    = intercept)

  pr <- structure(pr, class = "unf")

  return(pr)

}

#' @export
#' @title set up Minnesota Prior
#' @param mydata data
#' @param factordata data for factor models
#' @param no_factors number of factors (not yet implemented)
#' @param nolags number of lags (not yet implemented)
#' @param intercept whether the model has an intercept
#' @param lambda1 hyperparameter 1
#' @param lambda2 hyperparameter 2
#' @param lambda3 hyperparameter 3
#' @param lambda4 hyperparameter 4
#' @return returns an S3 object of the class ""minnesota"
#'
#' @family priors
#' @author Joerg Rieger
#' @importFrom stats lm
#' @references K. Rao Kadiyala and Sune Karlsson, Numerical Methods for Estimation and Inference in Bayesian VAR-Models, Journal of Applied Econometrics 12(2), 99-132
#' @references Gary Koop and Dimitris Korobilis (2010), Bayesian Multivariate Time Series Methods for Empirical Macroeconomics, Foundations and Trends in Econometrics 3(4), 267-358
#' @references Thomas Doan, Robert Litterman, Christopher A. Sims (1984), Forecasting and conditional projection using realistic prior distributions, Econometric Reviews 3(1), 1-100

set_prior_minnesota <- function(mydata,factordata=NULL,no_factors=0,nolags,intercept=TRUE,lambda1=1,lambda2=1,lambda3=1,lambda4=2){
  mydata <- as.matrix(mydata)

  if(is.null(factordata) && no_factors > 0){
    stop("Please provide data from which the factors can be extracted from.")
  }

  # get factor data
  if(no_factors > 0){

    factors <- get_factors(factordata,no_factors)
    data = cbind(mydata,factors)

  }
  else{

    data = mydata

  }

  # Declare variables
  obs <- nrow(data)
  K   <- ncol(data)

  constant = 0
  if(intercept==TRUE) constant=1

  #
  # Prior for coefficients
  #

  Aprior <- array(0,dim=c(K * nolags + constant,K))
  aprior <- as.vector(Aprior)

  #
  # Prior for covariance matrix
  #

  sigmasq <- array(0,dim=c(K,1))

  for(ii in 1:K){

    Ylagi         <- stats::embed(data[,ii],dimension = nolags + 1)[,-1]
    Yi            <- data[(nolags + 1):obs,ii]
    arest         <- stats::lm(Yi~Ylagi-1)
    sigmasq[ii,1] <- summary(arest)$sigma

  }

  #print(sigmasq)

  #
  # Covariance matrix for the prior
  #

  M <- K * nolags + constant
  Vi <- array(0,dim=c( K * M, 1))

  # without intercept
  for(ii in 1:K){ # loop over the ii-th equation

    for(jj in 1:nolags){ #loop over the jj-th lag

      for(kk in 1:K){ #kk-th variable

        indx <- (ii - 1) * (K * nolags) + (jj - 1) * K + kk

        if(ii==kk){

          Vi[indx,1] <- lambda1/(jj^lambda4)

        }
        else{
          Vi[indx,1] <- (lambda1 * lambda2)/(jj ^ lambda4) *(sigmasq[ii,1] / sigmasq[kk,1])^2
        }
      }
    }
  }

  #
  # Add Covariance coefficients for intercepts
  #

  if(intercept==TRUE){

    Vtmp <- array(0,dim=c(K * K * nolags + K, 1))

    for(ii in 1:K){

      coefinter    <- lambda3 * sigmasq[ii,1]^2
      indx         <- (ii-1)*(K*nolags)+ii
      Vtmp[indx,1] <- coefinter
      indxmin      <- (ii - 1) * (K * nolags) + 1
      indxmax      <- (ii - 1) * (K * nolags) + K * nolags

      Vtmp[(indx + 1):(indx + (K * nolags)),] <- Vi[indxmin:indxmax,]

    }
    Vi <- Vtmp
  }

  #
  # Create diagonal matrix
  #
  nr <- dim(Vi)[1]

  Vfinal <- array(0,dim=c(nr,nr))

  for(ii in 1:nr){

    Vfinal[ii,ii] <- Vi[ii,1]

  }

  pr <- list(type         = "Minnesota",
             nolags       = nolags,
             nofactors   = no_factors,
             intercept    = intercept,
             Aprior       = Aprior,
             Vprior       = Vfinal)

  pr <- structure(pr, class = "minnesota")

  return(pr)
}


#' @export
#' @title set up Stochastic Search Variable Selection Prior
#' @param mydata data
#' @param factordata data for factor models
#' @param no_factors number of factors (not yet implemented)
#' @param nolags number of lags (not yet implemented)
#' @param intercept whether the model has an intercept
#' @param tau parameter for prior on coefficients
#' @param kappa parameter for prior on variance-covariance matrix
#' @return returns an S3 object of the class ""minnesota"
#' @family priors
#' @author Joerg Rieger
#' @references K. Rao Kadiyala and Sune Karlsson, Numerical Methods for Estimation and Inference in Bayesian VAR-Models, Journal of Applied Econometrics 12(2), 99-132
#' @references Gary Koop and Dimitris Korobilis (2010), Bayesian Multivariate Time Series Methods for Empirical Macroeconomics, Foundations and Trends in Econometrics 3(4), 267-358
#' @references Dimitris Korobilis (2013), VAR Forecasting using Bayesian Variable Selection, Journal of Applied Econometrics 28(2),204-230
#'
set_prior_ssvs <- function(mydata, factordata = NULL, no_factors = 0, nolags, intercept=TRUE, tau, kappa){

  #tmp <- lagdata(mydata = mydata, nolags = nolags, intercept = intercept)
  #xLagged <- tmp$x
  #yLagged <- tmp$y
  if(is.null(factordata) && no_factors >0){

    stop("Please provide additional data for factor model")

  }


  K <- ncol(mydata) + no_factors

  norest <- K * (K * nolags + intercept)
  aprior <- array(0, dim =c(norest))

  tau0 <- tau
  tau1 <- 1/tau

  kappa0 <- kappa
  kappa1 <- 1/kappa

  pr <- list(type         = "SSVS",
             nolags       = nolags,
             nofactors    = no_factors,
             intercept    = intercept,
             tau0         = tau0,
             tau1         = tau1,
             kappa0       = kappa0,
             kappa1       = kappa1,
             aprior       = aprior)

  pr <- structure(pr, class = "ssvs")

}
