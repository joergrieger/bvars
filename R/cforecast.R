#' @export
#' @title Conditional forecasts
#' @param obj Estimated model
#' @param forecastHorizon forecast horizon
#' @param cfconds A \eqn{v\times s}-matrix of linear restrictions with v as the number of restrtictions and \eqn{s=h\times n}.
#' @param id_obj Identification
#' @param interval forecast bands
#' @param ... not used
#' @return returns an S3 object of the class fcbvar
#' @details
#' Conditional forecasts are forecasts conditional on given values for a subset of variables and are obtained by pre-determining the path of certain variables. Waggoner and Zha (1999) show that the distribution of future shocks is normal with
#' \deqn{\eta\sim N(\bar{\eta},\bar\Gamma})}
#' where
#' \deqn{\bar{\eta}=R'(RR')^{-1}r}
#' and
#' \deqn{\Gamma=I-R'(RR')^{-1}R}
#' with \eqn{\eta} the \eqn{s\times1}-vector of structural shocks and \eqn{r} is the vector of differences between predicred and conditional values.
#'
#'
#' @rdname cforecast
#' @references Waggoner, Daniel F. and Tao Zha, Conditional Forecasts in Dynamic Multivariate Models, The Review of Economics and Statistics, Vol. 81, No. 4 (Nov 1999), pp. 639-651

cforecast <- function(obj,forecastHorizon,id_obj,cfconds,interval = c(0.05,0.95),...) UseMethod("cforecast")



#' @export
#' @rdname cforecast
cforecast.bvar <- function(obj,forecastHorizon,id_obj,cfconds,interval = c(0.05,0.95),...){

  # Step one: unconditional forecasts
  # Preliminary Calculations
  nVariables       <- dim(obj$data_info$data)[2] # Number of variables
  nLength          <- dim(obj$data_info$data)[1] # Length of time series
  nLags            <- obj$general_info$nolags  # Number of lags
  nForecasts       <- dim(obj$mcmc_draws$Alpha)[3] # Number of forecasts, depends on the number sampled posteriors
  mForecastStorage <- array(NA,dim=c(forecastHorizon,nVariables,nForecasts)) # Matrix to storage forecasts
  cdForecast       <- array(NA,dim=c(forecastHorizon,nVariables))
  fc_tmp  <- array(NA,dim=c(forecastHorizon+nLags,nVariables)) # Matrix to temporary store forecasts

  # If the data is of class ts, get start date and frequency of data

  if(stats::is.ts(obj$data_info$data)){

    tsStart          <- stats::start(obj$data_info$data)
    tsFrequency      <- stats::frequency(obj$data_info$data)

  }

  tsColnames <- colnames(obj$data_info$data)

  # Check if mydata is a time series object
  # If it is a time series object get frequency, etc.
  if(stats::is.ts(obj$data_info$data)){

    nFreq <- frequency(obj$data_info$data)
    nFirstDate <- min(stats::time(obj$data_info$data))
    nFirstYear <- floor(nFirstDate)
    nFirstMonthQuarter <- (nFirstDate-nFirstYear) * nFreq
    nLastDate  <- max(stats::time(obj$data_info$data))

  }

  # Do the forecasts

  for(ii in 1:nForecasts){
    fc_tmp  <- array(NA,dim=c(forecastHorizon+nLags,nVariables)) # Matrix to temporary store forecasts
    fc_tmp[1:nLags,] <- (obj$data_info$data[nLength:(nLength-nLags+1),])

    # Step 1: Unconditional forecast for one path
    for(jj in 1:forecastHorizon){
      nstart <- jj
      nEnd   <- jj + nLags - 1
      y <- fc_tmp[nstart:nEnd,]
      randDraw <- rnorm(nVariables) %*% t(chol(obj$mcmc_draws$Sigma[,,ii]))

      if(obj$general_info$intercept == TRUE){

        y <- c(1,t(y))
        tempForecast <- y %*% obj$mcmc_draws$Alpha[,,ii] + randDraw

      }
      else{

        y <- t(y)
        tempForecast <- y %*% obj$mcmc_draws$Alpha[,,ii] + randDraw

      }

      # Store the unconditional
      fc_tmp[jj + nLags,] <- tempForecast
    }

    # Step 2: Calculate Impulse-Response Function
    ortirf <- compirf(Alpha = obj$mcmc_draws$Alpha[,,ii], Sigma = obj$mcmc_draws$Sigma[,,ii], id_obj = id_obj,
                      nolags = nLags, intercept = obj$general_info$intercept, nhor = forecastHorizon)

    # Step 3: create R matrix and r vector
    R <- vector()
    r <- vector()
    for(jj in 1:forecastHorizon){
      for(kk in 1:nVariables){
        if(!is.na(cfconds[jj,kk])){

          # r-vector condition - unconditional forecast
          r_tmp <- matrix(cfconds[jj,kk] - fc_tmp[jj + nLags,kk],nrow=1)
          r <- rbind(r,r_tmp)
          # R-Matrix
          R_tmp <- matrix(0,nrow=1,ncol = nVariables * forecastHorizon)
          for(ll in 1:jj){
            R_tmp[1,((ll - 1) * nVariables + 1):(ll * nVariables)] <- ortirf[kk,,(jj - ll + 1)]
          }
          R <- rbind(R,R_tmp)
        }
      }
    }

    # Step 4: Draw shocks from the Waggoner-Zha distribution
    Rdim <-  dim(R)
    Q <- Rdim[1]
    K <- Rdim[2]

    Rsvd <- svd(R,nu=nrow(R),nv=ncol(R))
    U <- Rsvd$u
    V <- Rsvd$v
    S <- Rsvd$d

    P  <- diag(S)
    V1 <- V[,1:Q]
    V2 <- V[,-c(1:Q)]

    eta <- matrix(pracma::mrdivide(V1,P) %*% t (U) %*% r + V2 %*% rnorm(n = K - Q),ncol=forecastHorizon)

    # Step 5: Obtain the conditional forecasts
    for(jj in 1:forecastHorizon){

      temp <- 0

      for(kk in 1:jj){

        temp <- temp + ortirf[,,jj - kk + 1] %*% eta[,kk]

      }

      cdForecast[jj,] <- fc_tmp[nLags + jj,] + t(temp)

    }

    # Step 6: Save forecasts
    mForecastStorage[,,ii] <- cdForecast

  }
  # Summarize results
  forecastMean  <- array(0,dim=c(forecastHorizon,nVariables))
  forecastUpper <- array(0,dim=c(forecastHorizon,nVariables))
  forecastLower <- array(0,dim=c(forecastHorizon,nVariables))

  for(ii in 1:nVariables){
    for(jj in 1:forecastHorizon){

      forecastMean[jj,ii] <- mean(mForecastStorage[jj,ii,])
      forecastLower[jj,ii] <- quantile(mForecastStorage[jj,ii,],probs=max(interval))
      forecastUpper[jj,ii] <- quantile(mForecastStorage[jj,ii,],probs=min(interval))

    }
  }

  OriginalPath <- array(NA,dim=c(nLength + forecastHorizon, nVariables))
  forecastFinalMean <- OriginalPath
  forecastFinalUpper <- OriginalPath
  forecastFinalLower <- OriginalPath

  OriginalPath[1:nLength,] <-obj$data_info$data
  forecastFinalMean[(nLength + 1):(nLength + forecastHorizon),] <- forecastMean
  forecastFinalUpper[(nLength + 1):(nLength + forecastHorizon),] <- forecastUpper
  forecastFinalLower[(nLength + 1):(nLength + forecastHorizon),] <- forecastLower

  if(is.ts(obj$data_info$data)){

    forecastFinalMean  <- ts(forecastFinalMean,start=tsStart,frequency=tsFrequency)
    forecastFinalUpper <- ts(forecastFinalUpper,start=tsStart,frequency=tsFrequency)
    forecastFinalLower <- ts(forecastFinalLower,start=tsStart,frequency=tsFrequency)
    OriginalPath       <- ts(OriginalPath,start=tsStart,frequency=tsFrequency)

  }

  colnames(forecastFinalMean)  <- colnames(obj$data_info$data)
  colnames(forecastFinalUpper) <- colnames(obj$data_info$data)
  colnames(forecastFinalLower) <- colnames(obj$data_info$data)
  colnames(OriginalPath)       <- colnames(obj$data_info$data)

  retList <- structure(list(forecast = forecastFinalMean, Upper = forecastFinalUpper, Lower = forecastFinalLower,
                            Original = OriginalPath),class="fcbvar")
  return(retList)

}

