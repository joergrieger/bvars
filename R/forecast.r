#' @export
#' @title forecasts for a bayesian VAR model
#' @param obj S3 object of the class bvar
#' @param forecastHorizon Forecast Horizon
#' @param interval forecast bands
#' @param ... currently not used
#' @return returns an S3 object of the class fcbvar
#' @rdname forecast
forecast <- function (obj,forecastHorizon = 16,interval = c(0.95,0.05),...) UseMethod("forecast")

#' @export
#' @importFrom stats start
#' @importFrom stats time
#' @rdname forecast
forecast.bvar <- function(obj,forecastHorizon = 16,interval = c(0.95,0.05),...){

  # Preliminary Calculations
  nVariables       <- dim(obj$data_info$data)[2] # Number of variables
  nLength          <- dim(obj$data_info$data)[1] # Length of time series
  nLags            <- obj$general_info$nolags  # Number of lags
  nForecasts       <- dim(obj$mcmc_draws$Alpha)[3] # Number of forecasts, depends on the number sampled posteriors
  mForecastStorage <- array(NA,dim=c(forecastHorizon+nLags,nVariables,nForecasts)) # Matrix to storage forecasts

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
    mForecastStorage[1:nLags,,ii] <- (obj$data_info$data[nLength:(nLength-nLags+1),])

    for(jj in 1:forecastHorizon){

      nStart <- jj
      nEnd   <- jj+nLags-1
      y <- mForecastStorage[nStart:nEnd,,ii]

      randDraw <- rnorm(nVariables) %*% t(chol(obj$mcmc_draws$Sigma[,,ii]))

      if(obj$general_info$intercept == TRUE){

        y <- c(1,t(y))

        tempForecast <- y %*% obj$mcmc_draws$Alpha[,,ii] + randDraw

      }
      else{

        y <- t(y)
        tempForecast <- y %*% obj$mcmc_draws$Alpha[,,ii] + randDraw

      }
      # Storing the forecast
      mForecastStorage[jj+nLags,,ii] <- tempForecast
    }

  }

  # Remove initial values
  mForecast <- mForecastStorage[-c(1:nLags),,]
  forecastMean  <- array(0,dim=c(forecastHorizon,nVariables))
  forecastUpper <- array(0,dim=c(forecastHorizon,nVariables))
  forecastLower <- array(0,dim=c(forecastHorizon,nVariables))


  for(ii in 1:nVariables){
    for(jj in 1:forecastHorizon){

      forecastMean[jj,ii]  <- mean(mForecast[jj,ii,])
      forecastUpper[jj,ii] <- quantile(mForecast[jj,ii,],probs=min(interval))
      forecastLower[jj,ii] <- quantile(mForecast[jj,ii,],probs=max(interval))

    }
  }

  #forecastFinalMean  <- rbind(as.matrix(bvarObj$mydata),forecastMean)
  #forecastFinalUpper <- rbind(as.matrix(bvarObj$mydata),forecastUpper)
  #forecastFinalLower <- rbind(as.matrix(bvarObj$mydata),forecastLower)

  OriginalPath       <- array(NA,dim=c(nLength + forecastHorizon, nVariables))
  forecastFinalMean  <- array(NA,dim=c(nLength + forecastHorizon, nVariables))
  forecastFinalUpper <- array(NA,dim=c(nLength + forecastHorizon, nVariables))
  forecastFinalLower <- array(NA,dim=c(nLength + forecastHorizon, nVariables))

  OriginalPath[1:nLength,] <- obj$data_info$data
  forecastFinalMean[(nLength + 1):(nLength + forecastHorizon),] <- forecastMean
  forecastFinalUpper[(nLength + 1):(nLength + forecastHorizon),] <- forecastUpper
  forecastFinalLower[(nLength + 1):(nLength + forecastHorizon),] <- forecastLower


  if(stats::is.ts(obj$data_info$data)){

    forecastFinalMean  <- stats::ts(forecastFinalMean,start=tsStart,frequency=tsFrequency)
    forecastFinalUpper <- stats::ts(forecastFinalUpper,start=tsStart,frequency=tsFrequency)
    forecastFinalLower <- stats::ts(forecastFinalLower,start=tsStart,frequency=tsFrequency)
    OriginalPath       <- stats::ts(OriginalPath, start = tsStart, frequency = tsFrequency)

  }

  colnames(forecastFinalMean)  <- colnames(obj$data_info$data)
  colnames(forecastFinalUpper) <- colnames(obj$data_info$data)
  colnames(forecastFinalLower) <- colnames(obj$data_info$data)
  colnames(OriginalPath)       <- colnames(obj$data_info$data)

  retList <- structure(list(forecast = forecastFinalMean, Upper = forecastFinalUpper, Lower = forecastFinalLower,
                            Original = OriginalPath),class="fcbvar")
  return(retList)

}

#'@export
#'@rdname forecast
forecast.msvar <- function(obj,forecastHorizon = 16, interval = c(0.05,0.95),...){

  # Preliminary Calculations
  nVariables <- dim(obj$data_info$data)[2]  # Number of variables
  nLength    <- dim(obj$data_info$data)[1]  # Length of time series
  nLags      <- obj$general_info$nolags     # Number of lags of the model
  nForecasts <- dim(obj$mcmc_draws$Alpha)[4] # total number of forecasts
  noregimes  <- obj$general_info$noregimes

  mForecastStorage <- array(NA,dim=c(forecastHorizon+nLags,nVariables,nForecasts))


  if(is.ts(obj$data_info$data)){
    tsStart     <- start(obj$data_info$data)
    tsFrequency <- frequency(obj$data_info$data)
  }
  tsColnames <- colnames(obj$data_info$data)


  # Weights for the forecasts
  fcWeights <- array(1/noregimes,dim=c(noregimes))
  # Do the forecasts

  for( ii in 1:nForecasts){
    mForecastStorage[1:nLags,,ii] <- obj$data_info$data[nLength:(nLength - nLags +1),]
    for(jj in 1:forecastHorizon){

      tmpForecast <- 0
      nStart <- jj
      nEnd   <- jj + nLags - 1

      y <- mForecastStorage[nStart:nEnd,,ii]

      # Prepare data
      if(obj$general_info$intercept == TRUE){

        y <- c(1,t(y))

      }
      else{

        y <- t(y)

      }
      #Forecast over all regimes than average it
      for(kk in 1:noregimes){

        randDraw <- stats::rnorm(nVariables) %*% t(chol(obj$mcmc_draws$Sigma[,,kk,ii]))

        tmp <- y %*% obj$mcmc_draws$Alpha[,,kk,ii] + randDraw
        tmpForecast <- tmpForecast + fcWeights[kk]  * tmp



      }
      # Update the forecast weights

      # Storing the forecast
      mForecastStorage[jj+nLags,,ii] <- tmpForecast

    }

  }
  # Preparing the data for returning
  forecastMean   <- array(0,dim=c(forecastHorizon,nVariables))
  forecastUpper  <- array(0,dim=c(forecastHorizon,nVariables))
  forecastLower  <- array(0,dim=c(forecastHorizon,nVariables))

  # Remove initial values
  mForecast <- mForecastStorage[-c(1:nLags),,]

  for(ii in 1:nVariables){
    for(jj in 1:forecastHorizon){

      forecastMean[jj,ii]   <- mean(mForecast[jj,ii,])
      forecastUpper[jj,ii] <- quantile(mForecast[jj,ii,],probs = max(interval))
      forecastLower[jj,ii]  <- quantile(mForecast[jj,ii,],probs = min(interval))

    }

  }

  OriginalPath        <- array(NA,dim=c(nLength + forecastHorizon,nVariables))
  forecastFinalMean   <- array(NA,dim=c(nLength + forecastHorizon,nVariables))
  forecastFinalUpper  <- array(NA,dim=c(nLength + forecastHorizon,nVariables))
  forecastFinalLower  <- array(NA,dim=c(nLength + forecastHorizon,nVariables))

  OriginalPath[1:nLength,] <- obj$data_info$data
  forecastFinalMean[(nLength + 1):(nLength + forecastHorizon),]  <- forecastMean
  forecastFinalUpper[(nLength + 1):(nLength + forecastHorizon),] <- forecastUpper
  forecastFinalLower[(nLength + 1):(nLength + forecastHorizon),] <- forecastLower

  if(is.ts(obj$data_info$data)){

    forecastFinalMean  <- ts(forecastFinalMean,start=tsStart,frequency=tsFrequency)
    forecastFinalUpper <- ts(forecastFinalUpper,start=tsStart,frequency=tsFrequency)
    forecastFinalLower <- ts(forecastFinalLower,start=tsStart,frequency=tsFrequency)
    OriginalPath       <- ts(OriginalPath, start = tsStart, frequency = tsFrequency)

  }

  colnames(forecastFinalMean)  <- colnames(obj$data_info$data)
  colnames(forecastFinalUpper) <- colnames(obj$data_info$data)
  colnames(forecastFinalLower) <- colnames(obj$data_info$data)
  colnames(OriginalPath)       <- colnames(obj$data_info$data)

  retList <- structure(list(forecast = forecastFinalMean, Upper = forecastFinalUpper, Lower = forecastFinalLower,
                            Original = OriginalPath),class="fcbvar")
  return(retList)

}

#' @export
#' @rdname forecast
forecast.tvar <- function(obj,forecastHorizon,interval= c(0.05,0.95),...){


  # Preliminary Calculations
  nVariables       <- dim(obj$data_info$data)[2] # Number of variables
  nLength          <- dim(obj$data_info$data)[1] # Length of time series
  nLags            <- obj$general_info$nolags  # Number of lags
  nForecasts       <- dim(obj$mcmc_draws$Alpha)[4] # Number of forecasts, depends on the number sampled posteriors

  nl1 <- max(obj$general_info$thMax,obj$gneral_info$nolags)
  mForecastStorage <- array(NA,dim=c(forecastHorizon+nl1,nVariables,nForecasts)) # Matrix to storage forecasts

  print(dim(mForecastStorage))


  if(is.ts(obj$data_info$data)){

    tsStart          <- start(obj$data_info$data)
    tsFrequency      <- frequency(obj$data_info$data)

  }

  for(ii in 1:nForecasts){

    mForecastStorage[1:nl1,,ii] <- (obj$data_info$data[nLength:(nLength-nl1+1),])

    for(jj in 1:forecastHorizon){

      nStart <- jj
      nEnd   <- jj+nLags-1
      y <- mForecastStorage[nStart:nEnd,,ii]
      thCheck <- mForecastStorage[jj + nl1 - obj$mcmc_draws$deldraws[ii] ,obj$general_info$thVar, ii]

      if(thCheck < obj$mcmc_draws$tardraws[ii]){

        # First Regime
        Alpha <- obj$mcmc_draws$Alpha[,,1,ii]
        Sigma <- obj$mcmc_draws$Sigma[,,1,ii]

      }
      else{

        # Second Regime
        Alpha <- obj$mcmc_draws$Alpha[,,2,ii]
        Sigma <- obj$mcmc_draws$Sigma[,,2,ii]

      }

      randDraw <- stats::rnorm(nVariables) %*% t(chol(Sigma))

      if(obj$general_info$intercept){

        y <- c(1,t(y))
        tempForecast <- y %*% Alpha + randDraw

      }
      else{
        y <- t(y)
        tempForecast <- y %*% Alpha + randDraw
      }

      # Storing the forecast
      mForecastStorage[jj+nl1,,ii] <- tempForecast

    }

  }
  # Remove initial values
  mForecast <- mForecastStorage[-c(1:nl1),,]
  forecastMean  <- array(0,dim=c(forecastHorizon,nVariables))
  forecastUpper <- array(0,dim=c(forecastHorizon,nVariables))
  forecastLower <- array(0,dim=c(forecastHorizon,nVariables))
  for(ii in 1:nVariables){
    for(jj in 1:forecastHorizon){

      forecastMean[jj,ii]  <- mean(mForecast[jj,ii,])
      forecastUpper[jj,ii] <- quantile(mForecast[jj,ii,],probs=max(interval))
      forecastLower[jj,ii] <- quantile(mForecast[jj,ii,],probs=min(interval))

    }
  }

  #forecastFinalMean  <- rbind(as.matrix(tvarObj$mydata),forecastMean)
  #forecastFinalUpper <- rbind(as.matrix(tvarObj$mydata),forecastUpper)
  #forecastFinalLower <- rbind(as.matrix(tvarObj$mydata),forecastLower)
  OriginalPath       <- array(NA,dim=c(nLength + forecastHorizon, nVariables))
  forecastFinalMean  <- array(NA,dim=c(nLength + forecastHorizon, nVariables))
  forecastFinalUpper <- array(NA,dim=c(nLength + forecastHorizon, nVariables))
  forecastFinalLower <- array(NA,dim=c(nLength + forecastHorizon, nVariables))

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
